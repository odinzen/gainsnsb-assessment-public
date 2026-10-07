# -*- coding: utf-8 -*-
"""Plot the quaternary Ga-In-Sn-Sb result as four constant-Sb sections over the Ga-In-Sn triangle
(x_Sb = 0, 2, 5, 10 at%), each FILLED BY LIQUIDUS TEMPERATURE (shared greyscale) with isotherm
lines. Coloring by temperature (not primary phase) is the point: past ~2 at% Sb the (Ga,In)Sb
zincblende is primary across the whole triangle, so the story is the liquidus surface lifting and
the low-melting coolant pocket near the Galinstan corner vanishing. Reads qsect_sb*.npz. House
style, no baked caption. calphad env."""
import os, warnings
warnings.filterwarnings("ignore")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "figures")
plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif", "Times New Roman"],
                     "font.size": 12, "savefig.facecolor": "white", "figure.facecolor": "white"})

LEVELS = ["00", "02", "05", "10"]
# shared temperature scale across panels so the surface lifting is visually honest
TMIN, TMAX = 0.0, 600.0
FILL = np.arange(0, 601, 25)
ISO = [50, 100, 200, 300, 400, 500]


def tri_xy(a, b, cc):  # Ga bottom-left, In top, Sn bottom-right
    s = a + b + cc
    return 0.5 * (2 * cc + b) / s, (np.sqrt(3) / 2) * b / s


fig, axes = plt.subplots(2, 2, figsize=(10.6, 9.8))
fig.subplots_adjust(left=0.03, right=0.88, top=0.97, bottom=0.05, wspace=0.28, hspace=0.20)
cf = None
for ax, tag in zip(axes.ravel(), LEVELS):
    d = np.load(os.path.join(HERE, "qsect_sb%s.npz" % tag), allow_pickle=True)
    gGa, gIn, gSn = d["gGa"].astype(float), d["gIn"].astype(float), d["gSn"].astype(float)
    T = d["T"].astype(float); xSb = float(d["xSb"])
    m = np.isfinite(T)
    gGa, gIn, gSn, T = gGa[m], gIn[m], gSn[m], T[m]
    X, Y = tri_xy(gGa, gIn, gSn)
    tri = mtri.Triangulation(X, Y)
    Tc = np.clip(T, TMIN, TMAX)
    cf = ax.tricontourf(tri, Tc, levels=FILL, cmap="Greys", vmin=TMIN, vmax=TMAX, extend="max")
    cs = ax.tricontour(tri, T, levels=ISO, colors="0.05", linewidths=0.7)
    ax.clabel(cs, fmt="%d", fontsize=7, inline=True)
    frame = np.array([tri_xy(1, 0, 0), tri_xy(0, 0, 1), tri_xy(0, 1, 0), tri_xy(1, 0, 0)])
    ax.plot(frame[:, 0], frame[:, 1], color="0.0", lw=1.1)
    ax.text(-0.03, -0.05, "Ga", ha="right", va="top", fontsize=11)
    ax.text(1.03, -0.05, "Sn", ha="left", va="top", fontsize=11)
    ax.text(0.5, np.sqrt(3) / 2 + 0.015, "In", ha="center", va="bottom", fontsize=11)
    # mark the Galinstan eutectic composition (Ga:In:Sn = 76:15:9) on every panel
    ex, ey = tri_xy(0.76, 0.15, 0.09)
    ax.plot(ex, ey, marker="o", ms=7, mfc="none", mec="0.0", mew=1.4)
    # constant-Sb label in the upper-left, clear of the In apex
    ax.text(0.0, 1.02, r"$x_{\mathrm{Sb}}$ = %g at%%" % round(xSb * 100),
            transform=ax.transAxes, ha="left", va="bottom", fontsize=13, fontweight="bold")
    ax.set_xlim(-0.12, 1.12); ax.set_ylim(-0.14, np.sqrt(3) / 2 + 0.12)
    ax.set_aspect("equal"); ax.axis("off")

cax = fig.add_axes([0.92, 0.12, 0.02, 0.76])
cb = fig.colorbar(cf, cax=cax)
cb.set_label("Liquidus temperature ($^{\\circ}$C)", fontsize=12)
fig.savefig(os.path.join(OUT, "Fig_quaternary_sections.png"), dpi=600)
plt.close(fig)
print("wrote Fig_quaternary_sections.png")
