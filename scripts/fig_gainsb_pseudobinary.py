# -*- coding: utf-8 -*-
"""Calculated InSb-GaSb pseudobinary liquidus and solidus against the measurements Liao et al. (1982)
assessed (Woolley & Lees 1959; Blom & Plaskett 1971), digitized from their Fig. 12
(data/gainsb_liao). The liquidus was fitted, the solidus is predicted. Greyscale, no caption on the image."""
import os, warnings
warnings.filterwarnings("ignore")
import numpy as np
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
comps = ["GA", "IN", "SB", "VA"]
phases = [p for p in sorted(db.phases.keys()) if p != "GAS"]  # condensed section; no vapor artifacts

# calculated liquidus + solidus across the join (x_Sb = 0.5), plotted vs x_GaSb = 1 - y_InSb.
# Liquidus from the exact driving force, solidus by bisection on the liquid fraction, both to 0.1 K.
from fast_eq import FastSystem
from validate_gainsb_liao import solidus
import fit_gainsb_liao as G
S = FastSystem(open(os.path.join(HERE, "..", "databases", "GaInSnSb.tdb")).read(), ["GA", "IN", "SB"], G.SOLIDS, [])
xgasb = np.round(np.arange(0.02, 0.981, 0.02), 3)
Tliq = np.array([S.liquidus((0.5 * x, 0.5 * (1 - x), 0.5), [], 790.0, 1000.0, xtol=0.1) - K for x in xgasb])
Tsol = np.array([solidus(S, (0.5 * x, 0.5 * (1 - x), 0.5)) - K for x in xgasb])
import csv
with open(os.path.join(HERE, "..", "data", "gainsb_liao", "gainsb_pseudobinary_liquidus_solidus_fig12_DIGITIZED.csv"), encoding="utf-8") as f:
    rows = list(csv.DictReader(l for l in f if not l.startswith("#")))
liao_liq = np.array([(float(r["x_GaSb"]), float(r["T_C"])) for r in rows if r["boundary"] == "liquidus"])
liao_sol = np.array([(float(r["x_GaSb"]), float(r["T_C"])) for r in rows if r["boundary"] == "solidus"])
import json
hn = json.load(open(os.path.join(HERE, "..", "results", "headline_numbers.json")))
T_gasb, T_insb = hn["GaSb_melting"]["T_C"], hn["InSb_melting"]["T_C"]

fig, ax = plt.subplots(figsize=(7.2, 5.6))
fig.subplots_adjust(left=0.12, right=0.96, top=0.96, bottom=0.12)
ax.plot(xgasb, Tliq, color="0.0", lw=2.0, zorder=3, label="Liquidus (this work)")
m = np.isfinite(Tsol)
ax.plot(xgasb[m], Tsol[m], color="0.0", lw=2.0, ls=(0, (5, 3)), zorder=3, label="Solidus (this work)")
ax.plot(liao_liq[:, 0], liao_liq[:, 1], "o", mfc="white", mec="0.0", mew=1.2, ms=6.5, ls="none",
        zorder=4, label="Liquidus, measured")
ax.plot(liao_sol[:, 0], liao_sol[:, 1], "^", mfc="0.35", mec="0.0", mew=0.9, ms=6.5, ls="none",
        zorder=4, label="Solidus, measured")
ax.set_xlabel("Mole fraction GaSb,  $x_{\\mathrm{GaSb}}$  (InSb $\\rightarrow$ GaSb)")
ax.set_ylabel("Temperature ($^{\\circ}$C)")
ax.set_xlim(0, 1); ax.set_ylim(500, 730)
# phase-region labels
ax.text(0.60, 722, "L", fontsize=13, ha="center", va="center")
ax.text(0.70, 636, "L + (Ga,In)Sb", fontsize=9, ha="center", va="center", rotation=18)
ax.text(0.40, 540, "(Ga,In)Sb", fontsize=12, ha="center", va="center")
ax.legend(loc="upper left", bbox_to_anchor=(0.02, 0.98), fontsize=10, frameon=False)
fig.savefig(os.path.join(OUT, "Fig_gainsb_pseudobinary.png"), dpi=600)
plt.close(fig)
o = np.argsort(liao_liq[:, 0])
Tmod = np.interp(liao_liq[o, 0], xgasb, Tliq)
print("wrote Fig_gainsb_pseudobinary.png ; liquidus RMS vs measured: %.1f C (n=%d)"
      % (float(np.sqrt(np.mean((Tmod - liao_liq[o, 1]) ** 2))), len(o)))
