# -*- coding: utf-8 -*-
"""Full-triangle Ga-In-Sn liquidus projection: at every composition, the liquidus temperature AND
the primary (first) solidifying phase. This gives the overall topology and the primary-phase fields
(Reviewer 1 point 5), not just the Ga-rich corner. Slow; caches to full_projection.npz."""
import os, numpy as np, warnings
warnings.filterwarnings("ignore")
from pycalphad import Database, equilibrium, variables as v

HERE = os.path.dirname(os.path.abspath(__file__))
db = Database(os.path.join(HERE, "..", "databases", "GaInSnSb.tdb"))
comps = ["GA", "IN", "SN", "VA"]
phases = [p for p in sorted(db.phases.keys()) if p != "GAS"]
K = 273.15
STEP = 0.03

rows = []
grid = [(round(xIn, 4), round(xSn, 4)) for xIn in np.arange(0.0, 1.0 + 1e-9, STEP)
        for xSn in np.arange(0.0, 1.0 + 1e-9, STEP) if xIn + xSn <= 1.0 + 1e-9]
for xIn, xSn in grid:
    xGa = max(1 - xIn - xSn, 0.0)
    xInc, xSnc = xIn, xSn
    if xGa < 2e-3:  # keep a sliver of Ga on the Ga=0 edge so pycalphad accepts it
        s = (1 - 2e-3) / max(xIn + xSn, 1e-9); xInc, xSnc = xIn * s, xSn * s
    # liquidus scan: melting points span ~10 C (Ga corner) to ~232 C (Sn). Scan wide.
    Ts = np.arange(510.0, 279.0, -1.0)
    e = equilibrium(db, comps, phases, {v.X('IN'): min(max(xInc, 1e-4), 1 - 3e-4),
                                        v.X('SN'): min(max(xSnc, 1e-4), 1 - 3e-4),
                                        v.T: Ts, v.P: 101325, v.N: 1})
    ph = e.Phase.values[0, 0]
    liq_only = [i for i in range(ph.shape[0]) if set(p for p in ph[i].ravel() if p) == {'LIQUID'}]
    if not liq_only:
        rows.append((xGa, xIn, xSn, np.nan, ""))
        continue
    iL = max(liq_only)
    Tliq = Ts[iL]
    # primary phase = the non-liquid phase present at the first temperature below the liquidus
    primary = ""
    if iL + 1 < ph.shape[0]:
        sset = set(p for p in ph[iL + 1].ravel() if p and p != 'LIQUID')
        primary = sorted(sset)[0] if sset else ""
    rows.append((xGa, xIn, xSn, Tliq - K, primary))

xGa = np.array([r[0] for r in rows]); xIn = np.array([r[1] for r in rows]); xSn = np.array([r[2] for r in rows])
T = np.array([r[3] for r in rows]); prim = np.array([r[4] for r in rows], dtype=object)
np.savez(os.path.join(HERE, "full_projection.npz"), xGa=xGa, xIn=xIn, xSn=xSn, T=T, prim=prim.astype(str))
uniq = {}
for p in prim:
    uniq[p] = uniq.get(p, 0) + 1
print("saved full_projection.npz  points:", len(rows))
print("primary-phase field sizes:", uniq)
