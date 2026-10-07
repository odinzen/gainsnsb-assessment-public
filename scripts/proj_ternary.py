# -*- coding: utf-8 -*-
"""Liquidus projection for one bounding ternary of the Ga-In-Sn-Sb quaternary: liquidus T +
primary phase over the full triangle, cached to npz. Usage: proj_ternary.py <gasnsb|insnsb>.
Reads databases/GaInSnSb.tdb (with the fitted Ga-In-Sb term; In-Sn-Sb is Muggianu per Vassiliev,
Ga-Sn-Sb is Muggianu placeholder). calphad env."""
import os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np
from pycalphad import Database, equilibrium, variables as v

HERE = os.path.dirname(os.path.abspath(__file__))
db = Database(os.path.join(HERE, "..", "databases", "GaInSnSb.tdb"))
phases = [p for p in sorted(db.phases.keys()) if p != "GAS"]
K = 273.15

# system -> (three elements A,B,C with A=bottom-left, B=top, C=bottom-right; the two independent
# axes are X(B), X(C)); Tmax(K) covers the highest liquidus (compound/Sb rich corners).
SYS = {
    "gasnsb": dict(elems=("GA", "SN", "SB"), Tmax=1010.0, tag="GaSnSb"),
    "insnsb": dict(elems=("IN", "SN", "SB"), Tmax=1010.0, tag="InSnSb"),
}
key = sys.argv[1] if len(sys.argv) > 1 else "gasnsb"
S = SYS[key]
A, B, C = S["elems"]
comps = [A, B, C, "VA"]
STEP = 0.04
Ts = np.arange(S["Tmax"], 279.0, -3.0)

rows = []
for xB in np.round(np.arange(0.0, 1.0 + 1e-9, STEP), 4):
    for xC in np.round(np.arange(0.0, 1.0 + 1e-9, STEP), 4):
        if xB + xC > 1.0 + 1e-9:
            continue
        xA = max(1 - xB - xC, 0.0)
        xBc, xCc = xB, xC
        if xA < 2e-3:  # keep a sliver of A on the A=0 edge so pycalphad accepts it
            s = (1 - 2e-3) / max(xB + xC, 1e-9); xBc, xCc = xB * s, xC * s
        cond = {v.X(B): min(max(round(xBc, 4), 1e-4), 1 - 3e-4),
                v.X(C): min(max(round(xCc, 4), 1e-4), 1 - 3e-4),
                v.T: Ts, v.P: 101325, v.N: 1}
        e = equilibrium(db, comps, phases, cond)
        ph = e.Phase.values[0, 0]
        idx = [i for i in range(ph.shape[0]) if set(p for p in ph[i].ravel() if p) == {"LIQUID"}]
        if not idx:
            rows.append((xA, xB, xC, np.nan, "")); continue
        iL = max(idx); Tl = Ts[iL] - K
        prim = ""
        if iL + 1 < ph.shape[0]:
            s = set(p for p in ph[iL + 1].ravel() if p and p != "LIQUID")
            prim = sorted(s)[0] if s else ""
        rows.append((xA, xB, xC, Tl, prim))

a = np.array(rows, dtype=object)
out = os.path.join(HERE, "proj_%s.npz" % key)
np.savez(out, xA=a[:, 0].astype(float), xB=a[:, 1].astype(float), xC=a[:, 2].astype(float),
         T=a[:, 3].astype(float), prim=a[:, 4].astype(str),
         elems=np.array(S["elems"]), tag=np.array(S["tag"]))
u = {}
for r in rows:
    u[r[4]] = u.get(r[4], 0) + 1
print("saved proj_%s.npz  points:%d  fields:%s" % (key, len(rows), u))
