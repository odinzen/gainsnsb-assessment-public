"""Fit the Ga-In-Sb ternary liquid term to the measurements Liao et al. (1982) assessed.

Data (data/gainsb_liao/): GaSb-InSb pseudobinary liquidus (Woolley & Lees 1959, Blom & Plaskett
1971), ternary liquidus isotherm points (Blom & Plaskett 1971, Antypas 1972), and the ternary
liquid enthalpies of mixing at 995 K (Ansara, Gambino & Bros 1976). The (Ga,In)Sb solid carries
Liao's interaction and is not adjusted. The term is L = a + b*T, composition independent; the
enthalpies fix a, the liquidus fixes a + b*T near 600 C.
Usage: python fit_gainsb_liao.py [db.tdb]   (prints the fit; does not write the database)
"""
import csv
import os
import re
import sys

import numpy as np
from scipy.optimize import least_squares

from fast_eq import FastSystem

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data", "gainsb_liao")
TDB = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "databases", "GaInSnSb.tdb")
K = 273.15
SIG = {"pb_liq": 5.0, "tern_liq": 10.0, "dH": 150.0}   # degC, degC, J/mol
SOLIDS = ["ORTHORHOMBIC_GA", "TETRAGONAL_A6", "RHOMBOHEDRAL_A7", "ZINCBLENDE_B3", "BCT_A5"]


def read(name):
    with open(os.path.join(DATA, name), encoding="utf-8") as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


PB = [r for r in read("gainsb_pseudobinary_liquidus_solidus_fig12_DIGITIZED.csv")]
PB_LIQ = [((float(r["x_Ga"]), float(r["x_In"]), float(r["x_Sb"])), float(r["T_C"]) + K) for r in PB if r["boundary"] == "liquidus"]
PB_SOL = [((float(r["x_Ga"]), float(r["x_In"]), float(r["x_Sb"])), float(r["T_C"]) + K) for r in PB if r["boundary"] == "solidus"]
TL = [((float(r["x_Ga"]), float(r["x_In"]), float(r["x_Sb"])), float(r["T_C"]) + K, r["confidence"])
      for r in read("gainsb_ternary_liquidus_fig13_DIGITIZED.csv")]
DH = [((float(r["x_Ga"]), float(r["x_In"]), float(r["x_Sb"])), float(r["dH_mix_exp_J"]))
      for r in read("gainsb_ternary_liquid_dHmix_995K_ansara1976_TABULATED.csv")]


def symbolic_text(path):
    txt = open(path).read()
    for o in (0, 1, 2):
        txt, n = re.subn(rf"(PARAMETER G\(LIQUID,GA,IN,SB;{o}\) 298\.15 )[^;]*;", r"\g<1>VV0100+VV0101*T;", txt)
        assert n == 1
    return txt


def build(path=TDB):
    return FastSystem(symbolic_text(path), ["GA", "IN", "SB"], SOLIDS, ["VV0100", "VV0101"])


def residuals(S, th, parts=False):
    r = {"pb_liq": [], "tern_liq": [], "dH": []}
    for x, TK in PB_LIQ:
        Tl = S.liquidus(x, th, TK - 80, TK + 80)
        r["pb_liq"].append((Tl - TK) if np.isfinite(Tl) else 80.0)
    for x, TK, _ in TL:
        Tl = S.liquidus(x, th, TK - 120, TK + 120)
        r["tern_liq"].append((Tl - TK) if np.isfinite(Tl) else 120.0)
    for x, h in DH:
        r["dH"].append(S.hmix(x, 995.0, th) - h)
    if parts:
        return {k: np.array(v) for k, v in r.items()}
    return np.concatenate([np.array(r[k]) / SIG[k] for k in ("pb_liq", "tern_liq", "dH")])


def report(S, th, label):
    r = residuals(S, th, parts=True)
    dh = np.array([h for _, h in DH])
    frac = np.sqrt(np.mean((r["dH"] / dh) ** 2))
    print(f"{label}: a={th[0]:.0f} b={th[1]:.3f} | pseudobinary liquidus RMS {np.sqrt(np.mean(r['pb_liq']**2)):.1f} C "
          f"({len(PB_LIQ)} pts) | ternary liquidus RMS {np.sqrt(np.mean(r['tern_liq']**2)):.1f} C ({len(TL)} pts) | "
          f"dH RMS {np.sqrt(np.mean(r['dH']**2)):.0f} J/mol, fractional {100*frac:.0f}% ({len(DH)} pts)")
    return r


if __name__ == "__main__":
    S = build()
    report(S, [0.0, 0.0], "no ternary term")
    report(S, [-24000.0, 0.0], "previous -24000")
    fit = least_squares(lambda th: residuals(S, th), np.array([-10000.0, 0.0]),
                        x_scale=np.array([5000.0, 5.0]), diff_step=np.array([0.02, 0.2]))
    J = fit.jac
    cov = np.linalg.inv(J.T @ J) * (2 * fit.cost / max(len(fit.fun) - 2, 1))
    print("fit:", fit.x, "+/-", np.sqrt(np.diag(cov)))
    report(S, fit.x, "fitted")
    np.save(os.path.join(HERE, "..", "results", "gainsb_liao_fit.npy"), np.concatenate([fit.x, np.sqrt(np.diag(cov))]))
