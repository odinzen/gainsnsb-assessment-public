# -*- coding: utf-8 -*-
"""The three antimony binary phase diagrams as SEPARATE full-size figures for the Supplementary
Material (S1 Ga-Sb, S2 In-Sb, S3 Sn-Sb). Same mapping and styling as the 3-panel figure, one file
each so the detail is legible. calphad env."""
import os, warnings
warnings.filterwarnings("ignore")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pycalphad import Database, equilibrium, variables as v
from pycalphad.mapping import BinaryStrategy, plot_binary

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "figures")
K = 273.15
plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif", "Times New Roman"],
                     "font.size": 13, "axes.linewidth": 1.0, "savefig.facecolor": "white",
                     "figure.facecolor": "white"})
db = Database(os.path.join(HERE, "..", "databases", "GaInSnSb.tdb"))

PANELS = [
    dict(A="GA", B="SB", Tmax=1010.0, out="Fig_SbBin_GaSb.png",
         phases=["LIQUID", "ORTHORHOMBIC_GA", "ZINCBLENDE_B3", "RHOMBOHEDRAL_A7"], title="Ga-Sb",
         labels=[(0.13, 680, "L", None, None), (0.27, 250, "GaSb", 0.5, 250),
                 (0.82, 320, "(Sb)", None, None), (0.16, 130, "(Ga)", 0.01, 22)]),
    dict(A="IN", B="SB", Tmax=950.0, out="Fig_SbBin_InSb.png",
         phases=["LIQUID", "TETRAGONAL_A6", "ZINCBLENDE_B3", "RHOMBOHEDRAL_A7"], title="In-Sb",
         labels=[(0.20, 580, "L", None, None), (0.27, 250, "InSb", 0.5, 250),
                 (0.82, 340, "(Sb)", None, None), (0.17, 250, "(In)", 0.01, 95)]),
    dict(A="SN", B="SB", Tmax=950.0, out="Fig_SbBin_SnSb.png",
         phases=["LIQUID", "BCT_A5", "SBSN", "SB3SN4", "RHOMBOHEDRAL_A7"], title="Sn-Sb",
         labels=[(0.20, 560, "L", None, None), (0.28, 110, r"Sb$_3$Sn$_4$", 0.4286, 165),
                 (0.68, 110, "SbSn", 0.56, 165), (0.955, 330, "(Sb)", None, None),
                 (0.17, 130, r"($\beta$-Sn)", 0.012, 175)]),
]
yc = list(range(0, 701, 100))
for P in PANELS:
    strat = BinaryStrategy(db, [P["A"], P["B"], "VA"], P["phases"],
                           {v.X(P["B"]): (0.002, 0.998, 0.004), v.T: (250, P["Tmax"], 4),
                            v.P: 101325, v.N: 1})
    strat.do_map()
    fig, ax = plt.subplots(figsize=(6.4, 5.4))
    fig.subplots_adjust(left=0.12, right=0.97, top=0.95, bottom=0.12)
    plot_binary(strat, ax=ax, tielines=0)
    ax.set_title("")                      # plot_binary adds the system name as a title
    for ln in ax.lines:
        ln.set_color("0.0"); ln.set_linewidth(1.2)
    leg = ax.get_legend()
    if leg is not None:
        leg.remove()
    # the mapper misses invariants that sit within a step of a pure-element melting point (the Ga-side
    # eutectic of Ga-Sb); find each one directly and draw it, with the liquidus down to that point
    if P["A"] == "GA":
        Ts = np.arange(300.0, 306.0, 0.01)
        e = equilibrium(db, ["GA", "SB", "VA"], P["phases"], {v.X("SB"): 0.25, v.T: Ts, v.P: 101325, v.N: 1})
        ph = e.Phase.values.reshape(len(Ts), -1)
        Te = Ts[next(i for i in range(len(Ts)) if "LIQUID" in ph[i])]
        ax.plot([0, 0.5], [Te, Te], color="0.0", lw=1.2)
        xs = np.concatenate([[1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3], [0.002]])
        Tl = []
        for x in xs:
            Tg = np.arange(Te, 700.0, 0.25)
            r = equilibrium(db, ["GA", "SB", "VA"], P["phases"], {v.X("SB"): x, v.T: Tg, v.P: 101325, v.N: 1})
            p = r.Phase.values.reshape(len(Tg), -1)
            Tl.append(Tg[next(i for i in range(len(Tg)) if set(q for q in p[i] if q) == {"LIQUID"})])
        ax.plot(np.concatenate([[0], xs]), np.concatenate([[Te], Tl]), color="0.0", lw=1.2)
    for (lx, lT, lab, tx, tT) in P["labels"]:
        if tx is None:
            ax.text(lx, lT + K, lab, fontsize=12, ha="center", va="center")
        else:
            ax.annotate(lab, xy=(tx, tT + K), xytext=(lx, lT + K), fontsize=12, ha="center",
                        va="center", arrowprops=dict(arrowstyle="-", lw=0.8, color="0.0"))
    ax.set_yticks([c + K for c in yc]); ax.set_yticklabels([str(c) for c in yc])
    ax.set_xlabel(f"Mole fraction {P['B'].title()},  $x_{{\\mathrm{{{P['B'].title()}}}}}$")
    ax.set_ylabel("Temperature ($^{\\circ}$C)")
    ax.set_xlim(0, 1); ax.set_ylim(K, P["Tmax"])
    fig.savefig(os.path.join(OUT, P["out"]), dpi=600)
    plt.close(fig)
    print("wrote", P["out"])
