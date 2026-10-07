"""Sensitivity of the Ga-In-Sb ternary term to the fit weights and to the one outlying liquidus point
(300 C, 3 at.% Sb). Writes results/gainsb_robustness.json; the database keeps the full-data fit."""
import json
import os

import numpy as np
from scipy.optimize import least_squares

import fit_gainsb_liao as G

HERE = os.path.dirname(os.path.abspath(__file__))
S = G.build()
BASE = dict(G.SIG)
OUTLIER = [i for i, (x, TK, _) in enumerate(G.TL) if abs(x[2] - 0.03) < 0.005 and abs(TK - 573.15) < 1]


def run(sig, drop):
    G.SIG.update(sig)
    keep = list(G.TL)
    if drop:
        G.TL[:] = [t for i, t in enumerate(keep) if i not in OUTLIER]
    f = least_squares(lambda th: G.residuals(S, th), np.array([3000.0, -24.0]),
                      x_scale=np.array([5000.0, 5.0]), diff_step=np.array([0.02, 0.2]))
    G.TL[:] = keep
    G.SIG.clear(); G.SIG.update(BASE)
    return {"a": float(f.x[0]), "b": float(f.x[1])}


cases = {"all_data": run({}, False), "no_outlier": run({}, True), "no_outlier_dH_x2": run({"dH": 75.0}, True),
         "no_outlier_dH_half": run({"dH": 300.0}, True), "no_outlier_tern_half": run({"tern_liq": 20.0}, True)}
a = [c["a"] for c in cases.values()]; b = [c["b"] for c in cases.values()]
out = {"cases": cases, "a_range": [min(a), max(a)], "b_range": [min(b), max(b)]}
json.dump(out, open(os.path.join(HERE, "..", "results", "gainsb_robustness.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
