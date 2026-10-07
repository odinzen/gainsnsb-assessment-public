"""Post-process the corrected-solid run: diagnostics, posterior summary, database update,
posterior-predictive invariants and the three-way transferability table.

Steps (each writes into results/ and bayesian/):
  diag    convergence, burn-in, flattened chain, medians, credible intervals, correlations
  tdb     write the posterior medians into databases/GaInSnSb.tdb
  ppc     ternary and Ga-Sn eutectics over posterior draws, full-solid model
  three   Table 2: liquid fitted to thermochemistry only, phase equilibria only, jointly
Usage: python finalize_full.py diag tdb ppc three
"""
import json
import os
import re
import sys
from multiprocessing import Pool

import numpy as np
from scipy.optimize import least_squares, minimize, minimize_scalar

import assess_data as E
import emcee_full as M
import forward_full as F

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
TDB = os.path.join(HERE, "..", "databases", "GaInSnSb.tdb")
LABELS = ["a0", "b0", "a1", "b1", "a2", "b2", "Lt"]
K = 273.15


def tau_sokal(x):
    x = x - x.mean(); n = len(x)
    f = np.fft.fft(x, n=2 * n)
    acf = np.fft.ifft(f * np.conj(f))[:n].real
    acf /= acf[0]
    tau = 1.0
    for m in range(1, n):
        tau = 1.0 + 2.0 * np.sum(acf[1:m])
        if m >= 5 * tau:
            break
    return max(tau, 1.0)


def diag():
    raw = np.load(os.path.join(HERE, "emcee_full_chain_raw.npy"))      # walkers, steps, dim
    nw, ns, nd = raw.shape
    # walker-averaged chain for tau (emcee 3 convention), per-walker mean as a cross-check
    tau_avg = np.array([tau_sokal(raw[:, :, k].mean(axis=0)) for k in range(nd)])
    burn = min(int(max(10 * tau_avg.max(), 0.2 * ns)), ns // 2)
    post = raw[:, burn:, :]
    tau_w = np.array([np.mean([tau_sokal(post[w, :, k]) for w in range(nw)]) for k in range(nd)])
    flat = post.reshape(-1, nd)
    np.save(os.path.join(HERE, "emcee_full_chain.npy"), flat)
    med = np.median(flat, axis=0)
    np.save(os.path.join(HERE, "posterior_median.npy"), med)
    q = np.percentile(flat, [16, 50, 84], axis=0)
    lnp = np.load(os.path.join(HERE, "emcee_full_lnp_raw.npy"))
    out = {
        "walkers": nw, "steps": ns, "burn_in": burn, "n_samples": int(flat.shape[0]),
        "steps_over_tau": round(float(ns / tau_w.max()), 1),
        "autocorr_time": dict(zip(LABELS, np.round(tau_w, 1).tolist())),
        "ess": dict(zip(LABELS, [int(flat.shape[0] / t) for t in tau_w])),
        "median": dict(zip(LABELS, med.tolist())),
        "ci68_plus": dict(zip(LABELS, (q[2] - q[1]).tolist())),
        "ci68_minus": dict(zip(LABELS, (q[1] - q[0]).tolist())),
        "corr": np.round(np.corrcoef(flat.T), 2).tolist(),
        "frac_Lt_ge_0": float(np.mean(flat[:, 6] >= 0)),
        "max_lnp": float(lnp.max()),
    }
    json.dump(out, open(os.path.join(RES, "posterior_summary.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "corr"}, indent=1))


def tdb():
    med = np.load(os.path.join(HERE, "posterior_median.npy"))
    txt = open(TDB).read()
    for o, (a, b) in enumerate([(med[0], med[1]), (med[2], med[3]), (med[4], med[5])]):
        txt, n = re.subn(rf"(PARAMETER G\(LIQUID,GA,SN;{o}\) 298\.15 )[^;]*;", rf"\g<1>{a:+.1f}{b:+.4f}*T;", txt)
        assert n == 1
    txt, n = re.subn(r"(PARAMETER G\(LIQUID,GA,IN,SN;0\) 298\.15 )[^;]*;", rf"\g<1>{med[6]:+.1f};", txt)
    assert n == 1
    open(TDB, "w").write(txt)
    print("wrote medians into", TDB)


def ternary_eutectic(theta):
    f = lambda z: F.liquidus((z[0], z[1], 1 - z[0] - z[1]), theta, 270.0, 300.0) if 0.05 < z[0] < 0.95 and 0.02 < z[1] < 0.9 else 400.0
    f2 = lambda z: f(z) if np.isfinite(f(z)) else 400.0
    r = minimize(f2, np.array([0.752, 0.154]), method="Nelder-Mead", options=dict(xatol=5e-4, fatol=5e-3, maxiter=300))
    return r.fun - K, r.x[0], r.x[1]


def gasn_eutectic(theta):
    f = lambda xs: F.liquidus((1 - xs, 0.0, xs), theta, 280.0, 320.0)
    r = minimize_scalar(lambda xs: f(xs) if np.isfinite(f(xs)) else 400.0, bounds=(0.04, 0.12), method="bounded",
                        options=dict(xatol=2e-4))
    return r.fun - K, r.x


def _ppc_one(th):
    return list(ternary_eutectic(th)) + list(gasn_eutectic(th))


def ppc(n=240):
    flat = np.load(os.path.join(HERE, "emcee_full_chain.npy"))
    idx = np.random.default_rng(3).choice(len(flat), n, replace=False)
    with Pool(6) as pool:
        rows = np.array(pool.map(_ppc_one, [flat[i] for i in idx]))
    med = np.load(os.path.join(HERE, "posterior_median.npy"))
    m = _ppc_one(med)
    pct = lambda a: [float(np.percentile(a, p)) for p in (16, 50, 84)]
    out = {"n_draws": n, "median_model": {"tern_T_C": m[0], "tern_xGa": m[1], "tern_xIn": m[2], "gasn_T_C": m[3], "gasn_xSn": m[4]},
           "tern_T_C_16_50_84": pct(rows[:, 0]), "gasn_T_C_16_50_84": pct(rows[:, 3]),
           "tern_xGa_16_50_84": pct(rows[:, 1]), "tern_xIn_16_50_84": pct(rows[:, 2]), "gasn_xSn_16_50_84": pct(rows[:, 4])}
    json.dump(out, open(os.path.join(RES, "posterior_predictive.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))


def thermo_resid(th):
    a0, b0, a1, b1, a2, b2 = th[:6]
    r = []
    for j, T in enumerate(E.KAT_T):
        aT = (a0 + b0 * T, a1 + b1 * T, a2 + b2 * T)
        for i, xg in enumerate(E.KAT_XGA):
            r.append((E.gasn_activity(xg, T, aT)[0] - E.A_GA[i, j]) / E.SIG_A)
    for xs, h in zip(E.ZIV_XSN, E.ZIV_HMIX):
        r.append((E.hmix(xs, (a0, a1, a2)) - h) / E.SIG_H)
    r += list(np.array(th[1:6:2]) / 4.0)       # the same weak prior on the b terms
    return np.array(r)


def phase_resid(th):
    xs = M.GASN_EUT[0]
    Te = F.liquidus((1 - xs, 0.0, xs), th, 283.0, 400.0)
    r = [((Te if np.isfinite(Te) else 400.0) - M.GASN_EUT[1]) / E.SIG_EUT]
    for (x, TK), sig in zip(M.TERN, M.SIG_TERN):
        Tl = F.liquidus(x, th, max(276.0, TK - 40), TK + 60)
        r.append(((Tl if np.isfinite(Tl) else TK + 60) - TK) / sig)
    r += list(np.array(th[1:6:2]) / 4.0)
    return np.array(r)


def scores(th):
    a0, b0, a1, b1, a2, b2, Lt = th
    ares = [E.gasn_activity(xg, T, (a0 + b0 * T, a1 + b1 * T, a2 + b2 * T))[0] - E.A_GA[i, j]
            for j, T in enumerate(E.KAT_T) for i, xg in enumerate(E.KAT_XGA)]
    hres = [E.hmix(xs, (a0, a1, a2)) - h for xs, h in zip(E.ZIV_XSN, E.ZIV_HMIX)]
    tres = [F.liquidus(x, th, max(276.0, TK - 40), TK + 60) - TK for x, TK in M.TERN]
    Te, xe = gasn_eutectic(th)
    rms = lambda a: float(np.sqrt(np.nanmean(np.square(a))))
    return {"aGa_rms": rms(ares), "H_rms": rms(hres), "gasn_eut_C": Te, "gasn_eut_xSn": xe,
            "tern_liq_rms_C": rms(tres), "n_tern_finite": int(np.isfinite(tres).sum())}


def three():
    med = np.load(os.path.join(HERE, "posterior_median.npy"))
    c1 = least_squares(thermo_resid, med[:6], x_scale=np.array([300, 0.3, 300, 0.3, 300, 0.3])).x
    c1 = np.concatenate([c1, [med[6]]])
    c2 = least_squares(phase_resid, med, x_scale=np.array([300, 0.3, 300, 0.3, 300, 0.3, 800]),
                       diff_step=np.array([0.02, 0.3, 0.05, 0.3, 0.05, 0.3, 0.02]), max_nfev=200).x
    out = {}
    for name, th in [("thermochemical_only", c1), ("phase_equilibrium_only", c2), ("joint", med)]:
        out[name] = {"params": dict(zip(LABELS, np.round(th, 4).tolist())), **scores(th)}
    json.dump(out, open(os.path.join(RES, "three_way.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))


def eutectic():
    """Why the calculated eutectic sits below the measured one: driving forces of the solids at the
    measured eutectic, the liquidus there, and what a smaller ternary term would do."""
    med = np.load(os.path.join(HERE, "posterior_median.npy"))
    xE, TE = (0.764, 0.144, 0.092), 10.7 + K
    mu = F.mu_liquid(xE, TE, med)
    df = {}
    for p in F.SOLIDS:
        ph = F.SOL[p]
        mus = np.array([mu[c] for c in ph.consts])
        df[p] = float((F.GRIDS[p] @ mus - ph.G(TE, F.GRIDS[p], med)).max())
    scan = []
    for Lt in [float(med[6]), -3000.0, -2000.0]:
        th = med.copy(); th[6] = Lt
        r = [F.liquidus(x, th, max(276.0, TK - 40), TK + 60) - TK for x, TK in M.TERN]
        scan.append({"Lt": Lt, "eutectic_C": ternary_eutectic(th)[0], "section_rms_C": float(np.sqrt(np.nanmean(np.square(r[:14])))),
                     "liquidus_at_measured_eutectic_C": float(r[14] + 10.7)})
    out = {"driving_force_at_measured_eutectic_J": df, "Lt_scan": scan}
    json.dump(out, open(os.path.join(RES, "eutectic_diagnostic.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    os.makedirs(RES, exist_ok=True)
    for step in sys.argv[1:]:
        {"diag": diag, "tdb": tdb, "ppc": ppc, "three": three, "eutectic": eutectic}[step]()
