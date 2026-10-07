# -*- coding: utf-8 -*-
"""Quaternary Ga-In-Sn-Sb liquidus, shown as constant-Sb sections over the Ga-In-Sn triangle.
For each fixed x_Sb, the Ga:In:Sn ratios span the triangle (renormalized to 1-x_Sb); we record
the liquidus T and the primary solidifying phase. This is the central quaternary result: the
low-melting Galinstan region is progressively consumed by the (Ga,In)Sb zincblende field as Sb
rises. Caches one npz per level. calphad env."""
import os, warnings
warnings.filterwarnings("ignore")
import numpy as np
from pycalphad import Database, equilibrium, variables as v

HERE = os.path.dirname(os.path.abspath(__file__))
db = Database(os.path.join(HERE, "..", "databases", "GaInSnSb.tdb"))
phases = [p for p in sorted(db.phases.keys()) if p != "GAS"]
K = 273.15
STEP = 0.05
LEVELS = [0.0, 0.02, 0.05, 0.10]
# grid reaching all three edges/corners (values are floored off exact 0/1 so pycalphad is happy)
GRID = np.round(np.arange(0.0, 1.0 + 1e-9, STEP), 4)

for xSb in LEVELS:
    if xSb == 0.0:
        comps = ["GA", "IN", "SN", "VA"]; Tmax = 540.0  # base has no zincblende; max ~232 C (Sn)
    else:
        comps = ["GA", "IN", "SN", "SB", "VA"]; Tmax = 1010.0  # GaSb-rich corner runs hot
    Ts = np.arange(Tmax, 279.0, -3.0)
    rows = []
    for gIn in GRID:
        for gSn in GRID:
            if gIn + gSn > 1.0 + 1e-9:
                continue
            gGa = max(1 - gIn - gSn, 0.0)  # nominal barycentric position for the plot
            xIn = gIn * (1 - xSb); xSn = gSn * (1 - xSb)
            # keep a sliver of Ga on the Ga=0 edge so pycalphad has a valid nonzero fraction
            xGa = 1 - xIn - xSn - xSb
            if xGa < 2e-3:
                scale = (1 - xSb - 2e-3) / max(xIn + xSn, 1e-9)
                xIn *= scale; xSn *= scale
            cond = {v.X("IN"): max(round(xIn, 5), 1e-4), v.X("SN"): max(round(xSn, 5), 1e-4),
                    v.T: Ts, v.P: 101325, v.N: 1}
            if xSb > 0:
                cond[v.X("SB")] = round(xSb, 6)
            e = equilibrium(db, comps, phases, cond)
            ph = e.Phase.values[0, 0]
            idx = [i for i in range(ph.shape[0]) if set(p for p in ph[i].ravel() if p) == {"LIQUID"}]
            if not idx:
                rows.append((gGa, gIn, gSn, np.nan, "")); continue
            iL = max(idx); Tl = Ts[iL] - K
            prim = ""
            if iL + 1 < ph.shape[0]:
                s = set(p for p in ph[iL + 1].ravel() if p and p != "LIQUID")
                prim = sorted(s)[0] if s else ""
            rows.append((gGa, gIn, gSn, Tl, prim))
    a = np.array(rows, dtype=object)
    tag = "%02d" % round(xSb * 100)
    np.savez(os.path.join(HERE, "qsect_sb%s.npz" % tag),
             gGa=a[:, 0].astype(float), gIn=a[:, 1].astype(float), gSn=a[:, 2].astype(float),
             T=a[:, 3].astype(float), prim=a[:, 4].astype(str), xSb=np.array(xSb))
    u = {}
    for r in rows:
        u[r[4]] = u.get(r[4], 0) + 1
    print("xSb=%.2f -> qsect_sb%s.npz  points:%d  fields:%s" % (xSb, tag, len(rows), u))
