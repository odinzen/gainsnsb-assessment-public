"""In-Sn-Sb liquid: calculated against measured chemical potential of In at 750 K (Vassiliev et al.
2001, Table 9). Binary extrapolation, no ternary term. Greyscale, no caption on the image."""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import validate_insnsb_vassiliev as V

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "figures", "Fig_insnsb_vassiliev.png")
plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif", "Times New Roman"],
                     "font.size": 13, "axes.linewidth": 1.1, "savefig.facecolor": "white",
                     "figure.facecolor": "white"})

rows = V.main()
meas = np.array([r[4] for r in rows])
calc = np.array([r[5] for r in rows])
fig, ax = plt.subplots(figsize=(5.2, 5.0))
lim = (-21.5, -5.0)
ax.plot(lim, lim, color="0.4", lw=1.0)
ax.plot(meas, calc, "o", mfc="white", mec="black", ms=7, mew=1.2)
ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel(r"Measured $\mu_{\mathrm{In}} - {}^{\circ}G_{\mathrm{In}}^{\mathrm{liq}}$ (kJ/mol)")
ax.set_ylabel(r"Calculated $\mu_{\mathrm{In}} - {}^{\circ}G_{\mathrm{In}}^{\mathrm{liq}}$ (kJ/mol)")
ax.set_aspect("equal")
fig.tight_layout()
fig.savefig(OUT, dpi=600)
print("saved", OUT)
