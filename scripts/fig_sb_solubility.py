"""Equilibrium Sb content of a Ga-In-Sn liquid at the Galinstan ratio (76:15:9) saturated with
(Ga,In)Sb, against temperature. Writes the curve to results/sb_solubility_curve.json and the figure.
Greyscale, no caption on the image."""
import json
import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pycalphad import Database, equilibrium, variables as v

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "figures", "Fig_sb_solubility.png")
plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif", "Times New Roman"],
                     "font.size": 13, "axes.linewidth": 1.1, "savefig.facecolor": "white",
                     "figure.facecolor": "white"})
db = Database(os.path.join(HERE, "..", "databases", "GaInSnSb.tdb"))
phases = [p for p in db.phases if p != "GAS"]
comps = ["GA", "IN", "SN", "SB", "VA"]
TC = np.arange(15.0, 401.0, 5.0)
xsb_tot = 0.03
r = equilibrium(db, comps, phases, {v.T: TC + 273.15, v.P: 101325, v.N: 1,
                                    v.X("GA"): 0.76 * (1 - xsb_tot), v.X("IN"): 0.15 * (1 - xsb_tot), v.X("SB"): xsb_tot})
n = len(TC)
PH = r.Phase.values.reshape(n, -1)
XX = r.X.values.reshape(n, PH.shape[1], -1)
isb = list(r.component.values).index("SB")
sol, ok = [], []
for i in range(n):
    ph = list(PH[i])
    good = "LIQUID" in ph and "ZINCBLENDE_B3" in ph
    sol.append(float(XX[i, ph.index("LIQUID"), isb]) if good else np.nan)
sol = np.array(sol)
json.dump({"T_C": TC.tolist(), "x_Sb_liquid": [None if not np.isfinite(s) else s for s in sol]},
          open(os.path.join(HERE, "..", "results", "sb_solubility_curve.json"), "w"), indent=1)
fig, ax = plt.subplots(figsize=(6.4, 5.0))
m = np.isfinite(sol)
ax.semilogy(TC[m], sol[m] * 1e6, color="black", lw=2.0)
ax.set_xlabel("Temperature ($^{\\circ}$C)")
ax.set_ylabel("Sb solubility in the liquid (ppm, atomic)")
ax.set_xlim(TC[m].min(), TC[m].max())
ax.grid(True, which="major", color="0.85", lw=0.6)
fig.tight_layout()
fig.savefig(OUT, dpi=600)
print("saved", OUT, "points", int(m.sum()))
