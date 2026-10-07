# -*- coding: utf-8 -*-
"""Quaternary result: liquidus isopleth from the Ga-In-Sn eutectic (Ga:In:Sn = 76:15:9) as Sb is
added. Shows the liquidus jump into the (Ga,In)Sb zincblende field. Muggianu (no ternary terms yet).
Greyscale, no baked caption. calphad env."""
import os, warnings, numpy as np
warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pycalphad import Database, equilibrium, variables as v

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "figures")
K = 273.15
plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif", "Times New Roman"],
                     "font.size": 13, "axes.linewidth": 1.1, "savefig.facecolor": "white",
                     "figure.facecolor": "white"})
db = Database(os.path.join(HERE, "..", "databases", "GaInSnSb.tdb"))
comps = ['GA', 'IN', 'SN', 'SB', 'VA']; phases = [p for p in sorted(db.phases.keys()) if p != 'GAS']

xsb = np.concatenate([[0.0, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 3e-3], np.arange(0.005, 0.1001, 0.0025)])
Tliq, prim = [], []
for xSb in xsb:
    xIn = 0.15 * (1 - xSb); xSn = 0.09 * (1 - xSb)
    Ts = np.arange(820.0, 279.0, -0.5)
    if xSb == 0:    # antimony-free ternary; any trace of Sb already puts GaSb on the liquidus
        e = equilibrium(db, comps[:3] + ['VA'], phases, {v.X('IN'): 0.15, v.X('SN'): 0.09, v.T: Ts, v.P: 101325, v.N: 1})
    else:
        e = equilibrium(db, comps, phases, {v.X('IN'): round(xIn, 6), v.X('SN'): round(xSn, 6),
                                            v.X('SB'): float(xSb), v.T: Ts, v.P: 101325, v.N: 1})
    ph = e.Phase.values[0, 0]
    idx = [i for i in range(ph.shape[0]) if set(p for p in ph[i].ravel() if p) == {'LIQUID'}]
    if idx:
        Tliq.append(Ts[max(idx)] - K)
        pr = [p for p in ph[max(idx) + 1].ravel() if p and p != 'LIQUID'] if max(idx) + 1 < ph.shape[0] else []
        prim.append(pr[0] if pr else '')
    else:
        Tliq.append(np.nan); prim.append('')
Tliq = np.array(Tliq); prim = np.array(prim, dtype=object)

fig, ax = plt.subplots(figsize=(7.2, 5.4))
fig.subplots_adjust(left=0.13, right=0.96, top=0.96, bottom=0.13)
m = np.isfinite(Tliq)
ax.plot(xsb[m] * 100, Tliq[m], color="0.0", lw=2.0, zorder=4)
import json
T0 = json.load(open(os.path.join(HERE, "..", "results", "headline_numbers.json")))["isopleth_76_15_9"]["0.000"]["T_C"]
ax.text(6.5, 575, "L", fontsize=14, fontweight="bold", ha="center")
ax.text(5.0, 250, "L + (Ga,In)Sb", fontsize=12, ha="center")
ax.set_xlabel("Antimony added, $x_{\\mathrm{Sb}}$ (at.%)  [Ga:In:Sn held at 76:15:9]")
ax.set_ylabel("Liquidus temperature ($^{\\circ}$C)")
ax.set_xlim(0, 10); ax.set_ylim(0, 620)
fig.savefig(os.path.join(OUT, "Fig_quaternary_sb_isopleth.png"), dpi=600)
plt.close(fig)
print("wrote Fig_quaternary_sb_isopleth.png")
print("liquidus at 0, 0.1, 1, 2, 5, 10 at% Sb (C):",
      [round(float(Tliq[np.argmin(abs(xsb-x))]), 1) for x in [0, 0.001, 0.01, 0.02, 0.05, 0.10]])
