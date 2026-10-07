# -*- coding: utf-8 -*-
"""Compose two multi-panel manuscript figures from existing single PNGs:
Fig_montage_gasn (activity + enthalpy, 1x2) and Fig_montage_projections (the four bounding-ternary
liquidus projections, 2x2), each with (a)-(d) panel labels. No recomputation; just layout. calphad env."""
import os, warnings
warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "..", "figures")
plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif", "Times New Roman"],
                     "savefig.facecolor": "white", "figure.facecolor": "white"})


def panel(ax, fname, label):
    ax.imshow(mpimg.imread(os.path.join(FIG, fname)))
    ax.axis("off")
    ax.text(0.01, 0.99, label, transform=ax.transAxes, fontsize=16, fontweight="bold",
            va="top", ha="left")


# Fig: Ga-Sn reassessment, 1x2 (activity, enthalpy)
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, wspace=0.03)
panel(ax[0], "Fig4_activity_gasn.png", "(a)")
panel(ax[1], "Fig3_hmix_gasn.png", "(b)")
fig.savefig(os.path.join(FIG, "Fig_montage_gasn.png"), dpi=600)
plt.close(fig)
print("wrote Fig_montage_gasn.png")

# Fig: four bounding-ternary projections, 2x2
fig, ax = plt.subplots(2, 2, figsize=(11, 10))
fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, wspace=0.02, hspace=0.02)
panel(ax[0, 0], "Fig6_projection_full.png", "(a) Ga-In-Sn")
panel(ax[0, 1], "Fig_gainsb_projection.png", "(b) Ga-In-Sb")
panel(ax[1, 0], "Fig_insnsb_projection.png", "(c) In-Sn-Sb")
panel(ax[1, 1], "Fig_gasnsb_projection.png", "(d) Ga-Sn-Sb")
fig.savefig(os.path.join(FIG, "Fig_montage_projections.png"), dpi=600)
plt.close(fig)
print("wrote Fig_montage_projections.png")
