# -*- coding: utf-8 -*-
"""Manuscript v2, Stage 2: binary figures against the corrected-solid model.
Fig 3  Ga-Sn excess enthalpy of mixing: this work vs Anderson-Ansara 1992 vs Zivkovic data
Fig 5  Ga-Sn phase diagram: this work (corrected beta-Sn) vs Anderson-Ansara 1992 liquidus, exp eutectic
Fig 1  Ga-In phase diagram: this work (corrected In solid solubility)
House style: greyscale serif, no baked caption text. Run from scripts/ with the calphad env."""
import sys, os
import numpy as np, warnings
warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pycalphad import Database, equilibrium, variables as v

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "figures")
DB = os.path.join(HERE, "..", "databases")
SGTE = os.path.join(HERE, "..", "_work", "sgte", "mmc3", "datasets", "datasets")
K = 273.15
plt.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif", "Times New Roman"],
    "font.size": 13, "axes.linewidth": 1.1, "axes.edgecolor": "black",
    "axes.facecolor": "white", "figure.facecolor": "white", "savefig.facecolor": "white",
})
full = Database(os.path.join(DB, "GaInSnSb.tdb"))


def liquidus_curve(db, other, phases, xs, Thi, Tlo, step=0.5):
    comps = ["GA", other, "VA"]
    out = []
    for x in xs:
        Ts = np.arange(Thi, Tlo, -step)
        e = equilibrium(db, comps, phases, {v.X(other): round(float(x), 4), v.T: Ts, v.P: 101325, v.N: 1})
        ph = e.Phase.values[0, 0]
        idx = [i for i in range(ph.shape[0]) if set(p for p in ph[i].ravel() if p) == {'LIQUID'}]
        out.append(Ts[max(idx)] - K if idx else np.nan)
    return np.array(out)


def hmix_curve(xSn, a):
    xGa = 1 - xSn; d = xGa - xSn
    return xGa * xSn * (a[0] + d * a[1] + d * d * a[2])


def fig3_hmix():
    med = np.load(os.path.join(HERE, "..", "bayesian", "posterior_median.npy"))
    a_this = (med[0], med[2], med[4])          # temperature-independent (enthalpy) parts
    a_and = (3369.7, 528.9, 0.0)               # Anderson-Ansara 1992 GaSn liquid a-parts
    ZIV_XSN = np.array([0.050, 0.0742, 0.200, 0.400, 0.700])
    ZIV_H = np.array([95.0, 306.0, 673.0, 907.0, 723.0])
    x = np.linspace(0.001, 0.999, 300)
    fig, ax = plt.subplots(figsize=(7.0, 5.2))
    fig.subplots_adjust(left=0.14, right=0.96, top=0.96, bottom=0.13)
    ax.plot(x, hmix_curve(x, a_this), color="0.0", lw=1.8, label="this work")
    ax.plot(x, hmix_curve(x, a_and), color="0.55", lw=1.8, ls=(0, (5, 2)), label="Anderson and Ansara (1992)")
    ax.plot(ZIV_XSN, ZIV_H, "o", mfc="white", mec="0.0", mew=1.4, ms=7, ls="none",
            label="Zivkovic et al. (2003)")
    ax.axhline(0, color="0.8", lw=0.8, zorder=0)
    ax.set_xlabel("Mole fraction Sn, $x_{\\mathrm{Sn}}$")
    ax.set_ylabel("Excess enthalpy of mixing (J mol$^{-1}$)")
    ax.set_xlim(0, 1)
    ax.legend(frameon=False, fontsize=11, loc="lower center", bbox_to_anchor=(0.5, 0.03))
    fig.savefig(os.path.join(OUT, "Fig3_hmix_gasn.png"), dpi=600)
    plt.close(fig)
    r = hmix_curve(ZIV_XSN, a_this) - ZIV_H
    print("Fig3 hmix: this-work RMS vs Zivkovic =", round(float(np.sqrt(np.mean(r**2))), 0), "J/mol")


def fig5_gasn_pd():
    and92 = Database(os.path.join(SGTE, "GaSn-92And-LB.tdb"))
    xs = np.round(np.arange(0.01, 1.0, 0.02), 4)
    L_this = liquidus_curve(full, "SN", [p for p in sorted(full.phases.keys()) if p != "GAS"], xs, 505.0, 285.0)
    L_and = liquidus_curve(and92, "SN", sorted(and92.phases.keys()), xs, 505.0, 285.0)
    fig, ax = plt.subplots(figsize=(7.0, 5.2))
    fig.subplots_adjust(left=0.13, right=0.96, top=0.96, bottom=0.13)
    ax.plot(xs, L_this, color="0.0", lw=1.8, label="this work")
    ax.plot(xs, L_and, color="0.55", lw=1.8, ls=(0, (5, 2)), label="Anderson and Ansara (1992)")
    ax.plot([0.085], [20.5], "s", mfc="white", mec="0.0", mew=1.4, ms=8, ls="none",
            label="experiment (eutectic)")
    ax.plot([0], [29.8], marker="_", ms=12, color="0.0")
    ax.plot([1], [231.9], marker="_", ms=12, color="0.0")
    ax.text(0.5, 300, "L", fontsize=14, fontweight="bold", ha="center")
    ax.set_xlabel("Mole fraction Sn, $x_{\\mathrm{Sn}}$")
    ax.set_ylabel("Temperature ($^{\\circ}$C)")
    ax.set_xlim(0, 1); ax.set_ylim(0, 260)
    ax.legend(frameon=False, fontsize=11, loc="upper center")
    fig.savefig(os.path.join(OUT, "Fig5_gasn_pd.png"), dpi=600)
    plt.close(fig)
    print("Fig5 Ga-Sn PD: this-work liquidus min", round(float(np.nanmin(L_this)), 1), "C")


def fig1_gain_pd():
    xs = np.round(np.arange(0.01, 1.0, 0.02), 4)
    L_this = liquidus_curve(full, "IN", [p for p in sorted(full.phases.keys()) if p != "GAS"], xs, 445.0, 285.0)
    fig, ax = plt.subplots(figsize=(7.0, 5.2))
    fig.subplots_adjust(left=0.13, right=0.96, top=0.96, bottom=0.13)
    ax.plot(xs, L_this, color="0.0", lw=1.8, label="this work")
    ax.plot([0], [29.8], marker="_", ms=12, color="0.0")
    ax.plot([1], [156.6], marker="_", ms=12, color="0.0")
    ax.text(0.5, 180, "L", fontsize=14, fontweight="bold", ha="center")
    ax.set_xlabel("Mole fraction In, $x_{\\mathrm{In}}$")
    ax.set_ylabel("Temperature ($^{\\circ}$C)")
    ax.set_xlim(0, 1); ax.set_ylim(0, 180)
    ax.legend(frameon=False, fontsize=11, loc="upper center")
    fig.savefig(os.path.join(OUT, "Fig1_gain_pd.png"), dpi=600)
    plt.close(fig)
    print("Fig1 Ga-In PD: this-work liquidus min", round(float(np.nanmin(L_this)), 1), "C")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "3"):
        fig3_hmix()
    if which in ("all", "5"):
        fig5_gasn_pd()
    if which in ("all", "1"):
        fig1_gain_pd()
    print("stage-2 done")
