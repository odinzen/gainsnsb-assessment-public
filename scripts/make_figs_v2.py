# -*- coding: utf-8 -*-
"""Manuscript v2 figures, corrected-solid model. Stage 1: posterior corner plot and the
model-vs-Katayama activity figure (both fast, no pycalphad). House style: greyscale, serif,
no caption/description text baked on the image (only axis labels and data-feature labels).
Run from scripts/ with the calphad env."""
import sys, os
import numpy as np, warnings
warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "figures")
CHAIN = os.path.join(HERE, "..", "bayesian", "emcee_full_chain.npy")
MED = os.path.join(HERE, "..", "bayesian", "posterior_median.npy")
R = 8.3144626
plt.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif", "Times New Roman"],
    "font.size": 13, "axes.linewidth": 1.1, "axes.edgecolor": "black",
    "axes.facecolor": "white", "figure.facecolor": "white", "savefig.facecolor": "white",
})

# --- Katayama 1996 activities (the data fit, never previously plotted) ---
KAT_XGA = np.array([0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90])
KAT_T = np.array([1000.0, 1073.0, 1200.0])
A_GA = np.array([[0.140, 0.132, 0.121], [0.275, 0.258, 0.236], [0.400, 0.374, 0.340],
                 [0.514, 0.480, 0.434], [0.634, 0.585, 0.522], [0.701, 0.662, 0.611],
                 [0.765, 0.740, 0.706], [0.843, 0.824, 0.799], [0.908, 0.898, 0.884]])
A_SN = np.array([[0.901, 0.900, 0.900], [0.804, 0.804, 0.804], [0.711, 0.711, 0.713],
                 [0.622, 0.623, 0.625], [0.558, 0.530, 0.538], [0.464, 0.455, 0.443],
                 [0.394, 0.370, 0.338], [0.294, 0.267, 0.232], [0.187, 0.157, 0.124]])
med = np.load(MED)  # a0,b0,a1,b1,a2,b2,Lt


def activity(xGa, T, a):
    xSn = 1 - xGa; d = xGa - xSn
    L0, L1, L2 = a
    Gx = xGa * xSn * (L0 + d * L1 + d * d * L2)
    dG = (1 - 2 * xGa) * (L0 + d * L1 + d * d * L2) + xGa * xSn * (2 * L1 + 4 * d * L2)
    return xGa * np.exp((Gx + xSn * dG) / (R * T)), xSn * np.exp((Gx - xGa * dG) / (R * T))


def fig_activity():
    a0, b0, a1, b1, a2, b2 = med[:6]
    xg = np.linspace(0.02, 0.98, 200)
    shades = ["0.0", "0.40", "0.65"]     # 1000, 1073, 1200 K in greyscale
    markers = ["o", "s", "^"]
    fig, ax = plt.subplots(figsize=(7.0, 5.4))
    fig.subplots_adjust(left=0.12, right=0.96, top=0.96, bottom=0.12)
    for k, T in enumerate(KAT_T):
        aT = (a0 + b0 * T, a1 + b1 * T, a2 + b2 * T)
        aGa = np.array([activity(x, T, aT)[0] for x in xg])
        aSn = np.array([activity(x, T, aT)[1] for x in xg])
        ax.plot(xg, aGa, color=shades[k], lw=1.6)
        ax.plot(xg, aSn, color=shades[k], lw=1.6, ls=(0, (5, 2)))
        ax.plot(KAT_XGA, A_GA[:, k], markers[k], mfc=shades[k], mec=shades[k], ms=6, ls="none",
                label=f"{int(T)} K")
        ax.plot(KAT_XGA, A_SN[:, k], markers[k], mfc="white", mec=shades[k], ms=6, ls="none")
    ax.plot([0, 1], [0, 1], color="0.8", lw=0.8, zorder=0)
    ax.text(0.82, 0.78, "$a_{\\mathrm{Ga}}$", fontsize=14)
    ax.text(0.12, 0.80, "$a_{\\mathrm{Sn}}$", fontsize=14)
    ax.set_xlabel("Mole fraction Ga, $x_{\\mathrm{Ga}}$")
    ax.set_ylabel("Activity")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.legend(title="filled $a_{\\mathrm{Ga}}$, open $a_{\\mathrm{Sn}}$; lines this work",
              loc="lower center", fontsize=10, title_fontsize=9, frameon=False)
    for s in ax.spines.values():
        s.set_linewidth(1.1)
    fig.savefig(os.path.join(OUT, "Fig4_activity_gasn.png"), dpi=600)
    plt.close(fig)
    # report fit quality
    res = []
    for k, T in enumerate(KAT_T):
        aT = (a0 + b0 * T, a1 + b1 * T, a2 + b2 * T)
        for i, x in enumerate(KAT_XGA):
            g, s = activity(x, T, aT)
            res += [g - A_GA[i, k], s - A_SN[i, k]]
    res = np.array(res)
    print("Fig4 activity: RMS a_Ga (fitted, n=27) =", round(float(np.sqrt(np.mean(res[0::2] ** 2))), 4),
          "| a_Sn (Gibbs-Duhem derived, n=27) =", round(float(np.sqrt(np.mean(res[1::2] ** 2))), 4))


def fig_corner():
    import corner
    ch = np.load(CHAIN)
    labels = ["$a_0$", "$b_0$", "$a_1$", "$b_1$", "$a_2$", "$b_2$", "$L_t$"]
    fig = corner.corner(ch, labels=labels, quantiles=[0.16, 0.5, 0.84], show_titles=False,
                        label_kwargs={"fontsize": 12}, color="0.25", hist_kwargs={"color": "0.25"})
    fig.savefig(os.path.join(OUT, "Fig10_corner.png"), dpi=600)
    plt.close(fig)
    print("Fig10 corner: chain", ch.shape)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "activity"):
        fig_activity()
    if which in ("all", "corner"):
        fig_corner()
    print("stage-1 figures written to", os.path.normpath(OUT))
