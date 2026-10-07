# -*- coding: utf-8 -*-
"""Plot the full-triangle Ga-In-Sn liquidus projection from full_projection.npz:
primary-solidification fields (greyscale) + liquidus isotherms + experimental points + eutectic.
Answers 'prove the overall liquidus topology'. House style, no baked caption. calphad env."""
import os, numpy as np, warnings
warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import ListedColormap
from proj_common import tri_xy, augment_edges, draw_field_labels

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "figures")
plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif", "Times New Roman"],
                     "font.size": 13, "savefig.facecolor": "white", "figure.facecolor": "white"})

d = np.load(os.path.join(HERE, "full_projection.npz"), allow_pickle=True)
xGa, xIn, xSn, T, prim = d["xGa"], d["xIn"], d["xSn"], d["T"], d["prim"].astype(str)
m = np.isfinite(T)
xGa, xIn, xSn, T, prim = xGa[m], xIn[m], xSn[m], T[m], prim[m]

X, Y = tri_xy(xGa, xIn, xSn)  # a=Ga bottom-left, b=In top, c=Sn bottom-right
tri = mtri.Triangulation(X, Y)  # isotherms + points on the raw data

# primary-phase fields -> integer codes, greyscale fills with readable labels
label = {"ORTHORHOMBIC_GA": "(Ga)", "TETRAGONAL_A6": "(In)", "BCT_A5": r"($\beta$-Sn)",
         "BETA": r"$\beta$", "GAMMA": r"$\gamma$"}
order = ["ORTHORHOMBIC_GA", "TETRAGONAL_A6", "BCT_A5", "BETA", "GAMMA"]
present = [p for p in order if p in set(prim)]
code = {p: i for i, p in enumerate(present)}
greys = ["0.97", "0.88", "0.80", "0.70", "0.62"][:len(present)]

# fill on an edge-augmented triangulation so the fields reach the frame
Xf, Yf, primf = augment_edges(xGa, xIn, xSn, prim)
Zf = np.array([code[p] for p in primf], dtype=float)
trif = mtri.Triangulation(Xf, Yf)

fig, ax = plt.subplots(figsize=(7.2, 6.4))
fig.subplots_adjust(left=0.03, right=0.97, top=0.93, bottom=0.05)
ax.tricontourf(trif, Zf, levels=np.arange(-0.5, len(present) + 0.5, 1), colors=greys)
cs = ax.tricontour(tri, T, levels=[15, 25, 50, 75, 100, 125, 150, 175, 200], colors="0.15", linewidths=0.8)
ax.clabel(cs, fmt="%d", fontsize=8, inline=True)

# experimental Evans-Prince liquidus points
ep = np.load(os.path.join(HERE, "..", "bayesian", "fig3_points.npy"))
ex, ey = tri_xy(ep[:, 0], ep[:, 1], ep[:, 2])
ax.plot(ex, ey, "o", mfc="white", mec="0.0", mew=1.1, ms=4.5, ls="none")
# calculated ternary eutectic, from results/headline_numbers.json
import json
eu = json.load(open(os.path.join(HERE, "..", "results", "headline_numbers.json")))["GaInSn_eutectic"]["liquid"]
eux, euy = tri_xy(eu["GA"], eu["IN"], eu["SN"])
ax.plot(eux, euy, marker="o", ms=9, mfc="0.0", mec="0.0")

# primary-field labels: interior for large fields, leader lines for small ones
draw_field_labels(ax, X, Y, prim, present, label,
                  nudge={"TETRAGONAL_A6": (0.05, 0.02), "BCT_A5": (-0.03, 0.08)})
# triangle frame + corners
c = np.array([tri_xy(1, 0, 0), tri_xy(0, 0, 1), tri_xy(0, 1, 0), tri_xy(1, 0, 0)])
ax.plot(c[:, 0], c[:, 1], color="0.0", lw=1.2)
ax.text(-0.03, -0.03, "Ga", ha="right", va="top", fontsize=14)
ax.text(1.03, -0.03, "Sn", ha="left", va="top", fontsize=14)
ax.text(0.5, np.sqrt(3) / 2 + 0.03, "In", ha="center", va="bottom", fontsize=14)
ax.set_aspect("equal"); ax.axis("off")
ax.set_xlim(-0.16, 1.16); ax.set_ylim(-0.12, 1.0)
fig.savefig(os.path.join(OUT, "Fig6_projection_full.png"), dpi=600)
plt.close(fig)
print("wrote Fig6_projection_full.png ; fields:", present)
