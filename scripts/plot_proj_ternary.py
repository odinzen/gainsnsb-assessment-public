# -*- coding: utf-8 -*-
"""Plot a bounding-ternary liquidus projection from proj_<key>.npz: greyscale primary-phase fields
+ liquidus isotherms. Usage: plot_proj_ternary.py <gasnsb|insnsb|gainsb>. House style, no baked
caption. For gainsb it reads the refreshed scripts/gainsb_projection.npz (with the fitted term).
calphad env."""
import os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from proj_common import tri_xy, augment_edges, draw_field_labels

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "figures")
plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif", "Times New Roman"],
                     "font.size": 13, "savefig.facecolor": "white", "figure.facecolor": "white"})

LABEL = {"ORTHORHOMBIC_GA": "(Ga)", "TETRAGONAL_A6": "(In)", "BCT_A5": r"($\beta$-Sn)",
         "RHOMBOHEDRAL_A7": "(Sb)", "ZINCBLENDE_B3": "(Ga,In)Sb", "SBSN": "SbSn",
         "SB3SN4": r"Sb$_3$Sn$_4$", "BETA": r"$\beta$", "GAMMA": r"$\gamma$"}
# per-system corner elements (A bottom-left, B top, C bottom-right) and the zincblende label
CFG = {
    "gainsb": dict(npz="gainsb_projection.npz", keys=("xGa", "xIn", "xSb"),
                   corners=("Ga", "In", "Sb"), out="Fig_gainsb_projection.png",
                   iso=[400, 450, 500, 550, 600, 650, 675], zb="(Ga,In)Sb"),
    "gasnsb": dict(npz="proj_gasnsb.npz", keys=("xA", "xB", "xC"),
                   corners=("Ga", "Sn", "Sb"), out="Fig_gasnsb_projection.png",
                   iso=[250, 350, 450, 550, 600], zb="GaSb",
                   force_leader=["RHOMBOHEDRAL_A7"]),  # (Sb) led out, matching panel (c)
    "insnsb": dict(npz="proj_insnsb.npz", keys=("xA", "xB", "xC"),
                   corners=("In", "Sn", "Sb"), out="Fig_insnsb_projection.png",
                   iso=[250, 350, 450, 500, 550, 600], zb="InSb",
                   force_leader=["RHOMBOHEDRAL_A7"]),  # (Sb) sits in a dark band -> lead it out
}
key = sys.argv[1] if len(sys.argv) > 1 else "gainsb"
c = CFG[key]
d = np.load(os.path.join(HERE, c["npz"]), allow_pickle=True)
A, Bv, C = (d[k].astype(float) for k in c["keys"])
T = d["T"].astype(float); prim = d["prim"].astype(str)
m = np.isfinite(T)
A, Bv, C, T, prim = A[m], Bv[m], C[m], T[m], prim[m]
lab = dict(LABEL); lab["ZINCBLENDE_B3"] = c["zb"]

X, Y = tri_xy(A, Bv, C)
tri = mtri.Triangulation(X, Y)  # isotherms on the raw data (edge positions untouched)

order = ["ORTHORHOMBIC_GA", "TETRAGONAL_A6", "BCT_A5", "BETA", "GAMMA",
         "ZINCBLENDE_B3", "SBSN", "SB3SN4", "RHOMBOHEDRAL_A7"]
present = [p for p in order if p in set(prim)]
code = {p: i for i, p in enumerate(present)}
greys = ["0.97", "0.90", "0.83", "0.76", "0.69", "0.62", "0.55", "0.48", "0.42"][:len(present)]

# fill on an edge-augmented triangulation so the fields reach the frame
Xf, Yf, primf = augment_edges(A, Bv, C, prim)
Zf = np.array([code[p] for p in primf], dtype=float)
trif = mtri.Triangulation(Xf, Yf)

fig, ax = plt.subplots(figsize=(7.2, 6.4))
fig.subplots_adjust(left=0.03, right=0.97, top=0.93, bottom=0.05)
ax.tricontourf(trif, Zf, levels=np.arange(-0.5, len(present) + 0.5, 1), colors=greys)
cs = ax.tricontour(tri, T, levels=c["iso"], colors="0.12", linewidths=0.8)
ax.clabel(cs, fmt="%d", fontsize=8, inline=True)

draw_field_labels(ax, X, Y, prim, present, lab, force_leader=c.get("force_leader", ()))
frame = np.array([tri_xy(1, 0, 0), tri_xy(0, 0, 1), tri_xy(0, 1, 0), tri_xy(1, 0, 0)])
ax.plot(frame[:, 0], frame[:, 1], color="0.0", lw=1.2)
ax.text(-0.03, -0.03, c["corners"][0], ha="right", va="top", fontsize=14)
ax.text(1.03, -0.03, c["corners"][2], ha="left", va="top", fontsize=14)
ax.text(0.5, np.sqrt(3) / 2 + 0.03, c["corners"][1], ha="center", va="bottom", fontsize=14)
ax.set_aspect("equal"); ax.axis("off")
ax.set_xlim(-0.16, 1.16); ax.set_ylim(-0.12, 1.0)
fig.savefig(os.path.join(OUT, c["out"]), dpi=600)
plt.close(fig)
print("wrote %s ; fields: %s" % (c["out"], present))
