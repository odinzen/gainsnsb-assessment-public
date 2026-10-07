# -*- coding: utf-8 -*-
"""Ga-In-Sb ternary liquidus projection: liquidus T + primary phase over the full triangle.
Shows the (Ga,In)Sb zincblende compound field. Muggianu (pre-ternary-refinement). Caches to npz."""
import os, numpy as np, warnings
warnings.filterwarnings("ignore")
from pycalphad import Database, equilibrium, variables as v

HERE = os.path.dirname(os.path.abspath(__file__))
db = Database(os.path.join(HERE, "..", "databases", "GaInSnSb.tdb"))
comps = ["GA", "IN", "SB", "VA"]; phases = [p for p in sorted(db.phases.keys()) if p != "GAS"]; K = 273.15
STEP = 0.035
rows = []
for xIn in np.round(np.arange(0.0, 1.0 + 1e-9, STEP), 4):
    for xSb in np.round(np.arange(0.0, 1.0 + 1e-9, STEP), 4):
        if xIn + xSb > 1.0 + 1e-9:
            continue
        xGa = max(1 - xIn - xSb, 0.0)
        xInc, xSbc = xIn, xSb
        if xGa < 2e-3:  # keep a sliver of Ga on the Ga=0 edge so pycalphad accepts it
            s = (1 - 2e-3) / max(xIn + xSb, 1e-9); xInc, xSbc = xIn * s, xSb * s
        Ts = np.arange(1010.0, 279.0, -3.0)  # K; up to 737 C covers the GaSb corner (712 C)
        e = equilibrium(db, comps, phases, {v.X('IN'): min(max(round(xInc, 4), 1e-4), 1 - 3e-4),
                                            v.X('SB'): min(max(round(xSbc, 4), 1e-4), 1 - 3e-4), v.T: Ts, v.P: 101325, v.N: 1})
        ph = e.Phase.values[0, 0]
        idx = [i for i in range(ph.shape[0]) if set(p for p in ph[i].ravel() if p) == {'LIQUID'}]
        if not idx:
            rows.append((xGa, xIn, xSb, np.nan, "")); continue
        iL = max(idx); Tl = Ts[iL] - K
        prim = ""
        if iL + 1 < ph.shape[0]:
            s = set(p for p in ph[iL + 1].ravel() if p and p != 'LIQUID')
            prim = sorted(s)[0] if s else ""
        rows.append((xGa, xIn, xSb, Tl, prim))
a = np.array(rows, dtype=object)
np.savez(os.path.join(HERE, "gainsb_projection.npz"),
         xGa=a[:, 0].astype(float), xIn=a[:, 1].astype(float), xSb=a[:, 2].astype(float),
         T=a[:, 3].astype(float), prim=a[:, 4].astype(str))
u = {}
for r in rows:
    u[r[4]] = u.get(r[4], 0) + 1
print("saved gainsb_projection.npz  points:", len(rows), " fields:", u)
