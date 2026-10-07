"""Ga-Sn-Sb (and InSb-Sn) against measurements, no ternary terms fitted.

Data (data/gasnsb/): Ga chemical potentials in Ga-Sb-Sn liquids at 1073 K from EMF (Katayama et al.
1998, Table 1), and the GaSb-Sn and InSb-Sn pseudobinary liquidus of Gerdes and Predel (1981), both
their Table 1 and the symbols of their Figs. 1 and 2, which disagree for GaSb-Sn. Writes
results/gasnsb_validation.json.
Usage: python validate_gasnsb.py [db.tdb]
"""
import csv
import json
import os
import sys
import warnings

import numpy as np
from pycalphad import Database, equilibrium, variables as v

from fast_eq import FastSystem

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data", "gasnsb")
TDB = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "databases", "GaInSnSb.tdb")
OUT = os.path.join(HERE, "..", "results", "gasnsb_validation.json")
F = 96485.33
SOL_GA = ["ORTHORHOMBIC_GA", "BCT_A5", "RHOMBOHEDRAL_A7", "ZINCBLENDE_B3", "SBSN", "SB3SN4"]
SOL_IN = ["TETRAGONAL_A6", "BCT_A5", "RHOMBOHEDRAL_A7", "ZINCBLENDE_B3", "SBSN", "SB3SN4", "BETA", "GAMMA"]


def read(name):
    with open(os.path.join(DATA, name), encoding="utf-8") as f:
        return list(csv.DictReader(l for l in f if not l.startswith("#")))


def rms(a):
    a = np.asarray(a, float)
    return float(np.sqrt(np.nanmean(a ** 2)))


def katayama(db):
    T = 1073.0
    g_ga = equilibrium(db, ["GA", "VA"], ["LIQUID"], {v.T: T, v.P: 101325, v.N: 1}).MU.sel(component="GA").values.item()
    rows = []
    for r in read("katayama1998_table1_emf_TABULATED.csv"):
        x = (float(r["x_Ga"]), float(r["x_Sb"]), float(r["x_Sn"]))
        e = equilibrium(db, ["GA", "SB", "SN", "VA"], ["LIQUID"],
                        {v.T: T, v.P: 101325, v.N: 1, v.X("GA"): x[0], v.X("SB"): x[1]})
        calc = (e.MU.sel(component="GA").values.item() - g_ga) / 1000
        meas = -3 * F * float(r["E_1073_mV"]) * 1e-3 / 1000
        rows.append({"y": float(r["section_y"]), "x": x, "meas_kJ": meas, "calc_kJ": calc})
    d = [r["calc_kJ"] - r["meas_kJ"] for r in rows]
    return rows, {"n": len(d), "rms_kJ": rms(d), "mean_kJ": float(np.mean(d)), "max_abs_kJ": float(np.max(np.abs(d)))}


def liquidus_set(S, rows, a_el):
    out = []
    for r in rows:
        xs = float(r["x_Sn"])
        if xs >= 0.995:
            continue
        x = {"GA": float(r["x_Ga"]), "IN": float(r["x_In"]), "SB": float(r["x_Sb"]), "SN": xs}
        TK = float(r["T_K"])
        Tl = S.liquidus(tuple(x[e] for e in S.el), [], 450.0, 1100.0, xtol=0.05)
        out.append({"x_Sn": xs, "meas_K": TK, "calc_K": Tl})
    return out


def trial_ternary(txt, res):
    """Fit L(Ga,Sb,Sn) = a + bT (and a constant) to the Katayama potentials and the Gerdes figure
    liquidus together, to judge whether a ternary term is justified. Not written to the database."""
    from scipy.optimize import least_squares
    t = txt.replace("$ Sb binary liquids (substitutional RK)",
                    "".join(f"PARAMETER G(LIQUID,GA,SB,SN;{o}) 298.15 VV0200+VV0201*T; 6000 N !\n" for o in range(3))
                    + "$ Sb binary liquids (substitutional RK)", 1)
    S = FastSystem(t, ["GA", "SB", "SN"], SOL_GA, ["VV0200", "VV0201"])
    kat = res["katayama1998_muGa_1073K"]["rows"]; fig = res["gerdes1981_GaSb-Sn_figure"]["points"]
    R, T = 8.31451, 1073.0
    eps = 1e-9

    def muga(x, th):
        mu = S.mu(np.array(x), T, th)
        mu0 = S.mu(np.array([1 - 2 * eps, eps, eps]), T, th)
        return (mu[0] - (mu0[0] - R * T * np.log(1 - 2 * eps))) / 1000

    def parts(th):
        a = np.array([muga(r["x"], th) - r["meas_kJ"] for r in kat])
        b = []
        for p in fig:
            xs = p["x_Sn"]
            Tl = S.liquidus((0.5 * (1 - xs), 0.5 * (1 - xs), xs), th, 450.0, 1100.0)
            b.append(Tl - p["meas_K"] if np.isfinite(Tl) else 100.0)
        return a, np.array(b)

    def resid(th):
        a, b = parts(th)
        return np.concatenate([a / 0.3, b / 5.0])

    out = {}
    c = least_squares(lambda t_: resid([t_[0], 0.0]), [-2000.0], x_scale=[3000.0])
    a, b = parts([c.x[0], 0.0])
    out["const"] = {"a_J": float(c.x[0]), "kat_rms_kJ": rms(a), "fig_rms_K": rms(b)}
    f = least_squares(resid, [-2000.0, 0.0], x_scale=[3000.0, 3.0])
    a, b = parts(f.x)
    cov = np.linalg.inv(f.jac.T @ f.jac) * (2 * f.cost / (len(f.fun) - 2))
    out["abT"] = {"a_J": float(f.x[0]), "b_J_K": float(f.x[1]), "a_sd": float(np.sqrt(cov[0, 0])),
                  "b_sd": float(np.sqrt(cov[1, 1])), "kat_rms_kJ": rms(a), "fig_rms_K": rms(b)}
    return out


def main():
    db = Database(TDB)
    txt = open(TDB).read()
    res = {}
    rows, stats = katayama(db)
    res["katayama1998_muGa_1073K"] = {**stats, "rows": rows}
    tab = read("gerdes1981_table1_liquidus_TABULATED.csv")
    fig = [r for r in read("gerdes1981_fig1_fig2_symbols_DIGITIZED.csv") if r["arrest"] == "liquidus"]
    for sysname, els, sol in [("GaSb-Sn", ["GA", "SB", "SN"], SOL_GA), ("InSb-Sn", ["IN", "SB", "SN"], SOL_IN)]:
        S = FastSystem(txt, els, sol, [])
        for tag, rr in [("table", [r for r in tab if r["system"] == sysname]),
                        ("figure", [r for r in fig if r["system"] == sysname])]:
            pts = liquidus_set(S, rr, els[0])
            d = [p["calc_K"] - p["meas_K"] for p in pts]
            res[f"gerdes1981_{sysname}_{tag}"] = {"n": len(d), "rms_K": rms(d), "mean_K": float(np.nanmean(d)), "points": pts}
        # degenerate eutectic on the Sn side
        xs_scan = np.linspace(0.90, 0.999, 60)
        Ts = [S.liquidus((0.5 * (1 - xs), 0.5 * (1 - xs), xs), [], 450.0, 700.0, xtol=0.02) for xs in xs_scan]
        res[f"model_{sysname}_liquidus_min_K"] = float(np.nanmin(Ts))
    res["trial_ternary_term"] = trial_ternary(txt, res)
    json.dump(res, open(OUT, "w"), indent=1)
    print(json.dumps({k: ({kk: vv for kk, vv in val.items() if kk not in ("rows", "points")} if isinstance(val, dict) else val)
                      for k, val in res.items()}, indent=1))


if __name__ == "__main__":
    main()
