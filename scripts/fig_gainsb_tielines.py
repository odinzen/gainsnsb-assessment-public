"""(Ga,In)Sb crystallized from In-rich Ga-In-Sb solutions: calculated InSb fraction of the solid against
the measured tie-lines (Miki et al. 1975 at 400 C, Antypas 1972 at 500 C, digitized in Liao et al.
1982). A prediction, not fitted. Reads results/gainsb_validation.json. Greyscale, no caption."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "figures", "Fig_gainsb_tielines.png")
plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif", "Times New Roman"],
                     "font.size": 13, "axes.linewidth": 1.1, "savefig.facecolor": "white",
                     "figure.facecolor": "white"})
pts = np.array(json.load(open(os.path.join(HERE, "..", "results", "gainsb_validation.json")))["tieline_points"], float)
fig, ax = plt.subplots(figsize=(5.2, 5.0))
ax.plot([0, 1], [0, 1], color="0.4", lw=1.0)
for T, mk, fc in [(400.0, "o", "white"), (500.0, "s", "0.35")]:
    m = pts[:, 0] == T
    ax.plot(pts[m, 1], pts[m, 2], mk, mfc=fc, mec="black", ms=6.5, mew=1.1, ls="none", label=f"{T:.0f} $^{{\\circ}}$C")
ax.set_xlim(0, 0.85); ax.set_ylim(0, 0.85); ax.set_aspect("equal")
ax.set_xlabel("Measured InSb fraction in (Ga,In)Sb")
ax.set_ylabel("Calculated InSb fraction in (Ga,In)Sb")
ax.legend(frameon=False, loc="upper left")
fig.tight_layout()
fig.savefig(OUT, dpi=600)
print("saved", OUT)
