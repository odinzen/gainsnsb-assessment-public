"""Ga-Sn-Sb and InSb-Sn without ternary terms against measurements. (a) Ga chemical potential in
Ga-Sb-Sn liquids at 1073 K, calculated against the EMF data of Katayama et al. (1998). (b) GaSb-Sn and
InSb-Sn pseudobinary liquidus against Gerdes and Predel (1981), figure symbols filled, their Table 1
open. Reads results/gasnsb_validation.json. Greyscale, no caption on the image."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fast_eq import FastSystem
import validate_gasnsb as V

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "figures", "Fig_gasnsb_validation.png")
plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif", "Times New Roman"],
                     "font.size": 12, "axes.linewidth": 1.1, "savefig.facecolor": "white",
                     "figure.facecolor": "white"})
K = 273.15
d = json.load(open(os.path.join(HERE, "..", "results", "gasnsb_validation.json")))
txt = open(V.TDB).read()

fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.0, 4.9))
rows = d["katayama1998_muGa_1073K"]["rows"]
for y, mk, fc in [(0.75, "o", "white"), (0.5, "s", "0.6"), (0.25, "^", "0.15")]:
    r = [q for q in rows if abs(q["y"] - y) < 1e-6]
    a1.plot([q["meas_kJ"] for q in r], [q["calc_kJ"] for q in r], mk, mfc=fc, mec="black", ms=7, ls="none",
            label=f"Sb:(Sb+Sn) = {y:.2f}")
lim = (-26, 0)
a1.plot(lim, lim, color="0.4", lw=1.0)
a1.set_xlim(lim); a1.set_ylim(lim); a1.set_aspect("equal")
a1.set_xlabel(r"Measured $\mu_{\mathrm{Ga}} - {}^{\circ}G_{\mathrm{Ga}}^{\mathrm{liq}}$ (kJ/mol)")
a1.set_ylabel(r"Calculated $\mu_{\mathrm{Ga}} - {}^{\circ}G_{\mathrm{Ga}}^{\mathrm{liq}}$ (kJ/mol)")
a1.legend(frameon=False, loc="upper left", fontsize=10)
a1.text(0.97, 0.03, "(a)", transform=a1.transAxes, ha="right", va="bottom", fontsize=13, fontweight="bold")

xs = np.linspace(0.0, 0.97, 60)
for sysname, els, sol, ls in [("GaSb-Sn", ["GA", "SB", "SN"], V.SOL_GA, "-"), ("InSb-Sn", ["IN", "SB", "SN"], V.SOL_IN, (0, (5, 3)))]:
    S = FastSystem(txt, els, sol, [])
    T = [S.liquidus((0.5 * (1 - x), 0.5 * (1 - x), x), [], 450.0, 1100.0, xtol=0.1) - K for x in xs]
    a2.plot(xs, T, color="black", lw=1.8, ls=ls, label=f"{sysname}, calculated")
    pf = d[f"gerdes1981_{sysname}_figure"]["points"]; pt = d[f"gerdes1981_{sysname}_table"]["points"]
    mk = "o" if sysname == "GaSb-Sn" else "s"
    a2.plot([p["x_Sn"] for p in pf], [p["meas_K"] - K for p in pf], mk, mfc="0.35", mec="black", ms=6, ls="none",
            label=f"{sysname}, Gerdes and Predel figure")
    a2.plot([p["x_Sn"] for p in pt], [p["meas_K"] - K for p in pt], mk, mfc="white", mec="0.35", ms=5, ls="none",
            label=f"{sysname}, Gerdes and Predel table")
a2.set_xlim(0, 1); a2.set_ylim(200, 750)
a2.set_xlabel(r"Mole fraction Sn, $x_{\mathrm{Sn}}$")
a2.set_ylabel(r"Liquidus temperature ($^{\circ}$C)")
a2.legend(frameon=False, loc="lower left", fontsize=8.5)
a2.text(0.97, 0.97, "(b)", transform=a2.transAxes, ha="right", va="top", fontsize=13, fontweight="bold")
fig.tight_layout()
fig.savefig(OUT, dpi=600)
print("saved", OUT)
