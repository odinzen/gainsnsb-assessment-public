"""Ga-In-Sb against the Liao et al. (1982) data, in Liao's own metrics, from the written database.

Fitted: pseudobinary liquidus, ternary liquidus, enthalpies of mixing. Not fitted (predictions):
pseudobinary solidus, and the solid compositions on the 400 C (Miki 1975) and 500 C (Antypas 1972)
tie-lines. Writes results/gainsb_validation.json.
Usage: python validate_gainsb_liao.py [db.tdb]
"""
import json
import os
import sys
import warnings

import numpy as np
from pycalphad import Database, equilibrium, variables as v
from scipy.optimize import brentq

import fit_gainsb_liao as G
from fast_eq import FastSystem

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
TDB = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "databases", "GaInSnSb.tdb")
OUT = os.path.join(HERE, "..", "results", "gainsb_validation.json")
K = 273.15


_LGRID = None


def _liquid_grid(n=200):
    pts = [(a, b, 1 - a - b) for a in np.linspace(0, 1, n + 1) for b in np.linspace(0, 1 - a, max(int(round((1 - a) * n)) + 1, 1))]
    return np.clip(np.array(pts), 1e-9, 1)


def solidus(S, x):
    """Pseudobinary solidus from the solid's own tangent planes. (Ga,In)Sb at composition x fixes
    mu_Ga - mu_In but leaves mu_Sb free; the solid is liquid-free while some mu_Sb keeps every liquid
    composition above that plane. Solidus = the T where min over mu_Sb of the largest liquid driving
    force reaches zero."""
    from scipy.optimize import minimize_scalar
    global _LGRID
    if _LGRID is None:
        _LGRID = _liquid_grid()
    zb = next(f for ph, f, Y, X in S.sol if ph == "ZINCBLENDE_B3")
    yga = x[0] / (x[0] + x[1])
    order = [S.lc.index(e) for e in S.el]
    Ly = _LGRID[:, [S.el.index(c) for c in S.lc]]

    def gzb(y, T):   # J per formula unit (Ga,In)1Sb1
        return 2 * float(np.asarray(zb(S._cols(T, np.array([[y, 1 - y, 1.0]]), []))).ravel()[0])

    def worst(T):
        h = 1e-6
        d = (gzb(yga + h, T) - gzb(yga - h, T)) / (2 * h)          # mu_Ga - mu_In
        g = gzb(yga, T)
        GL = np.asarray(S._lg(S._cols(T, Ly, []))).ravel()

        def dfmax(musb):
            mu_in = g - musb - yga * d
            mu = np.array([mu_in + d, mu_in, musb])                 # Ga, In, Sb
            return float(np.max(_LGRID @ mu - GL))
        r = minimize_scalar(dfmax, bounds=(g - 2e5, g + 2e5), method="bounded", options=dict(xatol=0.5))
        return r.fun

    return brentq(worst, 780.0, 1000.0, xtol=0.05)


def tieline_solid(S, T, liq_coord, coord):
    """Liquid on the T isotherm with the given coordinate; returns the InSb fraction of the
    (Ga,In)Sb solid that crystallizes from it. coord 'xGa' (Fig. 14) or 'ga_ratio' (Fig. 15)."""
    def liquid(xsb):
        if coord == "xGa":
            xga = liq_coord
        else:
            xga = liq_coord * (1 - xsb)
        return np.array([xga, 1 - xga - xsb, xsb])

    f = lambda xsb: S.driving_force(liquid(xsb), T, [])[0]
    grid = np.linspace(1e-4, 0.45, 400)
    vals = [f(s) for s in grid]
    k = next((i for i in range(1, len(grid)) if vals[i - 1] < 0 <= vals[i]), None)
    if k is None:
        return np.nan, np.nan
    xsb = brentq(f, grid[k - 1], grid[k], xtol=1e-7)
    ph, xs = S.driving_force(liquid(xsb), T, [])[1]
    if ph != "ZINCBLENDE_B3":
        return np.nan, xsb
    return xs[1] / (xs[0] + xs[1]), xsb


def main():
    S = FastSystem(open(TDB).read(), ["GA", "IN", "SB"], G.SOLIDS, [])
    db = Database(TDB)
    phases = [p for p in db.phases if p != "GAS"]
    res = {}
    r = {"pb_liq": [], "tern_liq": [], "dH": []}
    for x, TK in G.PB_LIQ:
        r["pb_liq"].append(S.liquidus(x, [], TK - 80, TK + 80) - TK)
    for x, TK, _ in G.TL:
        r["tern_liq"].append(S.liquidus(x, [], TK - 120, TK + 120) - TK)
    for x, h in G.DH:
        r["dH"].append(S.hmix(x, 995.0, []) - h)
    r = {k: np.array(val) for k, val in r.items()}
    steep = np.array([abs(x[2] - 0.03) < 0.005 and abs(TK - 573.15) < 1 for x, TK, _ in G.TL])
    res["pseudobinary_liquidus_rms_C"] = [round(float(np.sqrt(np.mean(r["pb_liq"] ** 2))), 1), len(r["pb_liq"])]
    res["ternary_liquidus_rms_C"] = [round(float(np.sqrt(np.mean(r["tern_liq"] ** 2))), 1), len(r["tern_liq"])]
    res["ternary_liquidus_rms_C_without_dilute_300C_point"] = [round(float(np.sqrt(np.mean(r["tern_liq"][~steep] ** 2))), 1), int((~steep).sum())]
    res["dHmix_rms_J"] = [round(float(np.sqrt(np.mean(r["dH"] ** 2))), 0), len(r["dH"])]
    dh = np.array([h for _, h in G.DH]); big = np.abs(dh) > 100
    res["dHmix_fractional_rms_pct_above_100J"] = [round(float(100 * np.sqrt(np.mean((r["dH"][big] / dh[big]) ** 2))), 0), int(big.sum())]
    sol = [solidus(S, x) - TK for x, TK in G.PB_SOL]
    res["pseudobinary_solidus_rms_C_predicted"] = [round(float(np.sqrt(np.nanmean(np.square(sol)))), 1), len(sol)]
    tl = []
    for name, coord in [("gainsb_tieline_400C_fig14_DIGITIZED.csv", "xGa"), ("gainsb_tieline_500C_fig15_DIGITIZED.csv", "ga_ratio")]:
        for row in G.read(name):
            y_calc, _ = tieline_solid(S, float(row["T_C"]) + K, float(row["x_liq_axis"]), coord)
            tl.append((float(row["T_C"]), float(row["y_InSb_solid"]), y_calc))
    d = np.array([c - m for _, m, c in tl])
    res["tieline_solid_InSb_fraction_rms_predicted"] = [round(float(np.sqrt(np.nanmean(d ** 2))), 3), int(np.isfinite(d).sum())]
    res["tieline_points"] = tl
    json.dump(res, open(OUT, "w"), indent=1)
    print(json.dumps({k: val for k, val in res.items() if k != "tieline_points"}, indent=1))


if __name__ == "__main__":
    main()
