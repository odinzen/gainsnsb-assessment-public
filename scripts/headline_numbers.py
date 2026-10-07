"""Every number the manuscript quotes, recalculated from the shipped database.

Writes results/headline_numbers.json. Nothing in the text should be quoted from anywhere else.
Usage: python headline_numbers.py [path/to/db.tdb]
"""
import json
import os
import re
import sys
import warnings

import numpy as np
from pycalphad import Database, equilibrium, variables as v

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
TDB = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "databases", "GaInSnSb.tdb")
OUT = os.path.join(HERE, "..", "results", "headline_numbers.json")
K = 273.15
TXT = open(TDB).read()
DB = Database(TXT)
CONDENSED = [p for p in DB.phases if p != "GAS"]


def variant(**edits):
    """Database text with regex edits applied, for the counterfactual cases."""
    t = TXT
    for pat, rep in edits.values():
        t, n = re.subn(pat, rep, t)
        assert n >= 1, pat
    return Database(t)


def phases_and_liquid(db, comps, X, Ts, phases):
    cond = {v.T: Ts, v.P: 101325, v.N: 1, **{v.X(k): x for k, x in X.items()}}
    r = equilibrium(db, comps + ["VA"], phases, cond)
    n = len(Ts)
    PH = r.Phase.values.reshape(n, -1)
    NP = r.NP.values.reshape(n, -1)
    XX = r.X.values.reshape(n, PH.shape[1], -1)
    order = list(r.component.values)
    out = []
    for i in range(n):
        ph = list(PH[i])
        present = sorted({p for p, a in zip(ph, NP[i]) if p and a > 1e-9})
        liq = None
        if "LIQUID" in ph:
            m = ph.index("LIQUID")
            liq = {c: float(XX[i, m, order.index(c)]) for c in comps}
        out.append((float(Ts[i]), present, liq, (PH[i], XX[i], r.Y.values.reshape(n, PH.shape[1], -1)[i], order), i))
    return out


def first_liquid(db, comps, X, T0, T1, step=0.02, phases=None):
    """Lowest T at which liquid appears on heating: the invariant for an alloy inside its 3-phase
    (binary) or 4-phase (ternary) region. Returns T (C), assemblage below, liquid composition."""
    phases = phases or CONDENSED
    Ts = np.arange(T0, T1, 0.5)
    coarse = phases_and_liquid(db, comps, X, Ts, phases)
    for k, (T, pres, liq, *_) in enumerate(coarse):
        if "LIQUID" in pres:
            lo = coarse[k - 1][0]
            break
    else:
        return None
    fine = phases_and_liquid(db, comps, X, np.arange(lo, T + step, step), phases)
    below = fine[0][1]
    for T, pres, liq, *_ in fine:
        if "LIQUID" in pres:
            return {"T_C": round(T - K, 2), "below": below, "liquid": {c: round(x, 4) for c, x in liq.items()}}
        below = pres


def liquidus(db, comps, X, lo, hi, tol=0.05, phases=None):
    phases = phases or CONDENSED

    def state(T):
        return phases_and_liquid(db, comps, X, np.array([T]), phases)[0][1]

    if state(hi) != ["LIQUID"] or state(lo) == ["LIQUID"]:
        return None
    while hi - lo > tol:
        m = 0.5 * (lo + hi)
        if state(m) == ["LIQUID"]:
            hi = m
        else:
            lo = m
    return {"T_C": round(0.5 * (lo + hi) - K, 2), "primary": [p for p in state(lo) if p != "LIQUID"]}


def max_solubility(db, comps, solute, phase, X, Ts):
    best = 0.0
    for T, pres, liq, (ph, xx, yy, order), i in phases_and_liquid(db, comps, X, Ts, CONDENSED):
        ph = list(ph)
        if phase in ph:
            m = ph.index(phase)
            best = max(best, float(xx[m, order.index(solute)]))
    return round(best, 4)


def quat(ratio, xsb):
    s = 1 - xsb
    return {"GA": ratio[0] * s, "IN": ratio[1] * s, "SB": xsb}


def main():
    res = {"database": os.path.basename(TDB)}
    gis = ["GA", "IN", "SN"]
    res["GaInSn_eutectic"] = first_liquid(DB, gis, {"GA": 0.60, "IN": 0.25}, 278.0, 290.0)
    no_tern = variant(t=(r"(PARAMETER G\(LIQUID,GA,IN,SN;0\) 298\.15 )[^;]*;", r"\g<1>0;"))
    res["GaInSn_eutectic_no_ternary_term"] = first_liquid(no_tern, gis, {"GA": 0.60, "IN": 0.25}, 278.0, 292.0)
    pure = variant(a=(r"CONSTITUENT TETRAGONAL_A6 : GA,IN,SN : !", "CONSTITUENT TETRAGONAL_A6 : IN,SN : !"),
                   b=(r"CONSTITUENT BCT_A5 : GA,IN,SN,SB : !", "CONSTITUENT BCT_A5 : IN,SN,SB : !"))
    res["GaInSn_eutectic_pure_Ga_solids"] = first_liquid(pure, gis, {"GA": 0.60, "IN": 0.25}, 278.0, 290.0)
    res["GaSn_eutectic"] = first_liquid(DB, ["GA", "SN"], {"SN": 0.5}, 285.0, 300.0)
    res["GaSn_eutectic_pure_Ga_solids"] = first_liquid(pure, ["GA", "SN"], {"SN": 0.5}, 285.0, 300.0)
    res["GaIn_eutectic"] = first_liquid(DB, ["GA", "IN"], {"IN": 0.5}, 283.0, 295.0)
    res["InSn_eutectic"] = first_liquid(DB, ["IN", "SN"], {"SN": 0.5}, 385.0, 400.0)
    res["Ga_in_bSn_max"] = max_solubility(DB, ["GA", "SN"], "GA", "BCT_A5", {"SN": 0.5}, np.arange(285.0, 293.0, 0.25))
    res["Ga_in_In_max"] = max_solubility(DB, ["GA", "IN"], "GA", "TETRAGONAL_A6", {"IN": 0.5}, np.arange(280.0, 289.0, 0.25))
    # pure In: the In-Sn gamma phase of 04Dav becomes stable as pure In above about 844 K (an artifact
    # of the published assessment, noted by Hallstedt), so bracket below it
    res["melting"] = {el: liquidus(DB, [el], {}, lo, hi, tol=0.02)
                      for el, lo, hi in [("GA", 290.0, 1000.0), ("IN", 300.0, 800.0), ("SN", 300.0, 1000.0), ("SB", 300.0, 1000.0)]}
    gam = {}
    for T in [800.0, 850.0, 900.0]:
        gam[str(T)] = phases_and_liquid(DB, ["IN"], {}, np.array([T]), CONDENSED)[0][1]
    res["pure_In_stable_phase_vs_T_K"] = gam
    res["GaSb_melting"] = liquidus(DB, ["GA", "SB"], {"SB": 0.5}, 900.0, 1100.0)
    res["InSb_melting"] = liquidus(DB, ["IN", "SB"], {"SB": 0.5}, 700.0, 900.0)
    res["GaSb_eutectic_Sb_side"] = first_liquid(DB, ["GA", "SB"], {"SB": 0.95}, 840.0, 880.0)
    res["InSb_eutectic_Sb_side"] = first_liquid(DB, ["IN", "SB"], {"SB": 0.80}, 740.0, 790.0)
    res["InSb_eutectic_In_side"] = first_liquid(DB, ["IN", "SB"], {"SB": 0.10}, 415.0, 435.0)
    res["SnSb_invariants"] = {
        "L+Sb3Sn4->bSn": first_liquid(DB, ["SB", "SN"], {"SB": 0.15}, 505.0, 530.0),
        "L+SbSn->Sb3Sn4": first_liquid(DB, ["SB", "SN"], {"SB": 0.45}, 580.0, 610.0),
        "L+Sb->SbSn": first_liquid(DB, ["SB", "SN"], {"SB": 0.70}, 680.0, 710.0),
    }
    q = ["GA", "IN", "SN", "SB"]
    gal = (0.76, 0.15, 0.09)
    res["isopleth_76_15_9"] = {f"{x:.3f}": liquidus(DB, q, quat(gal, x), 270.0, 1050.0, tol=0.1)
                               for x in [0.0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.10]}
    sol = {}
    for Tc in [25, 50, 100, 200]:
        pts = phases_and_liquid(DB, q, quat(gal, 0.001), np.array([Tc + K]), CONDENSED)[0]
        _, pres, liq, (ph, xx, yy, order), i = pts
        ph = list(ph)
        zb = ph.index("ZINCBLENDE_B3") if "ZINCBLENDE_B3" in ph else None
        y_in = None
        if zb is not None:
            y = yy[zb]
            y_in = float(y[1] / (y[0] + y[1]))
        sol[str(Tc)] = {"phases": pres, "x_Sb_liquid": liq["SB"] if liq else None, "zincblende_In_fraction": y_in}
    res["Sb_solubility_galinstan"] = sol
    corner = {}
    for name, ratio in [("Sn", (0.0, 0.0, 1.0)), ("Ga2In2Sn96", (0.02, 0.02, 0.96)), ("Ga5In5Sn90", (0.05, 0.05, 0.90))]:
        corner[name] = {f"{x:.2f}": liquidus(DB, q, quat(ratio, x), 450.0, 1050.0, tol=0.1) for x in [0.02, 0.05, 0.10]}
    res["Sn_corner_liquidus"] = corner
    boil = None
    for TK in np.arange(1700.0, 2000.0, 0.5):
        r = equilibrium(DB, ["SB", "VA"], list(DB.phases), {v.T: TK, v.P: 101325, v.N: 1})
        if "GAS" in set(r.Phase.values.ravel()):
            boil = round(TK - K, 1)
            break
    res["Sb_boiling_1atm_C"] = boil
    r = equilibrium(DB, q + ["VA"], list(DB.phases), {v.T: 973.15, v.P: 101325, v.N: 1, **{v.X(k): x for k, x in quat(gal, 0.02).items()}})
    res["gas_at_700C_galinstan_2Sb"] = "GAS" in set(r.Phase.values.ravel())
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(res, open(OUT, "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
