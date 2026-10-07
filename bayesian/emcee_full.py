"""Bayesian joint assessment of the Ga-In-Sn liquid against the shipped solid descriptions.

Same data, likelihood and priors as emcee_assess.py, with two changes: the forward model is
forward_full (every phase read from databases/GaInSnSb.tdb, solids with their Ga solubility), and
all 17 ternary liquidus points enter the likelihood instead of a subset, and only the
measured Ga activities are fitted (Katayama's Sn activities are integrated from them, not independent). The chain runs long enough
for about 50 autocorrelation times.
"""
import os
import sys
import time

import numpy as np
import emcee

import assess_data as E
import forward_full as F

HERE = os.path.dirname(os.path.abspath(__file__))
TERN = E.TERN_LIQUIDUS                        # 14 Evans-Prince section points, eutectic, 2 further points
# Each liquidus datum carries its own uncertainty: the section points are read off Evans-Prince Fig. 2
# (+/-2.5 C), the eutectic is their chemically analysed invariant, 10.7 +/- 0.3 C (sigma 0.5 C allows
# for the composition uncertainty). Set SIG_EUT_TERN = 2.5 to recover the earlier uniform weighting.
SIG_EUT_TERN = float(os.environ.get("SIG_EUT_TERN", 0.5))
SIG_TERN = [E.SIG_TL] * 14 + [SIG_EUT_TERN] + [E.SIG_TL] * 2
assert len(SIG_TERN) == len(TERN)
GASN_EUT = (0.085, 20.5 + 273.15)             # Ga-Sn eutectic, 8.5 at.% Sn (Evans-Prince Table 1)


def log_like(theta):
    a0, b0, a1, b1, a2, b2, Lt = theta
    ll = 0.0
    for j, T in enumerate(E.KAT_T):
        aT = (a0 + b0 * T, a1 + b1 * T, a2 + b2 * T)
        for i, xg in enumerate(E.KAT_XGA):
            aGa, _ = E.gasn_activity(xg, T, aT)   # a_Sn in Katayama Table 2 is Gibbs-Duhem derived
            ll += -0.5 * ((aGa - E.A_GA[i, j]) / E.SIG_A) ** 2
    for xs, h in zip(E.ZIV_XSN, E.ZIV_HMIX):
        ll += -0.5 * ((E.hmix(xs, (a0, a1, a2)) - h) / E.SIG_H) ** 2
    xs = GASN_EUT[0]
    Te = F.liquidus((1 - xs, 0.0, xs), theta, 283.0, 400.0)
    if not np.isfinite(Te):
        return -np.inf
    ll += -0.5 * ((Te - GASN_EUT[1]) / E.SIG_EUT) ** 2
    for (x, TK), sig in zip(TERN, SIG_TERN):
        Tl = F.liquidus(x, theta, max(276.0, TK - 14), TK + 14)
        if not np.isfinite(Tl):
            Tl = F.liquidus(x, theta, 276.0, TK + 60)
        if not np.isfinite(Tl):
            return -np.inf
        ll += -0.5 * ((Tl - TK) / sig) ** 2
    return ll


def log_prob(theta):
    lp = E.log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    ll = log_like(theta)
    return lp + ll if np.isfinite(ll) else -np.inf


if __name__ == "__main__":
    from multiprocessing import Pool
    nsteps = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
    ndim, nwalk = 7, 24
    rng = np.random.RandomState(11)
    start = np.median(np.load(os.path.join(HERE, "emcee_chain.npy")), axis=0)
    p0 = start + np.array([150, 0.15, 200, 0.25, 400, 0.4, 600.]) * rng.randn(nwalk, ndim)
    out = os.path.join(HERE, "emcee_full_chain_raw.npy")
    with Pool(6) as pool:
        sampler = emcee.EnsembleSampler(nwalk, ndim, log_prob, pool=pool)
        t0 = time.time()
        for i, _ in enumerate(sampler.sample(p0, iterations=nsteps)):
            if (i + 1) % 100 == 0 or i + 1 == nsteps:
                chain = sampler.chain[:, : i + 1, :]          # emcee 2.x: (walkers, steps, dim)
                np.save(out, chain)
                np.save(os.path.join(HERE, "emcee_full_lnp_raw.npy"), sampler.lnprobability[:, : i + 1])
                try:
                    tau = np.array([emcee.autocorr.integrated_time(chain[:, :, k].mean(axis=0))
                                    for k in range(ndim)]).ravel()
                except Exception:
                    tau = np.full(ndim, np.nan)
                print(f"step {i+1}  {time.time()-t0:7.0f}s  acc {np.mean(sampler.acceptance_fraction):.2f}  "
                      f"tau {np.round(tau, 0)}", flush=True)
    print("done", flush=True)
