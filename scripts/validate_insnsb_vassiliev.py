"""In-Sn-Sb liquid against the EMF chemical potentials of Vassiliev et al. (2001), Table 9, 750 K.

No ternary term is fitted; this tests the binary extrapolation. Alloys 14 and 15 are printed in
parentheses by Vassiliev (side reactions) and are excluded. For alloy 19 the Table 6 value
(-16.62 kJ/mol) is used because the printed Table 9 entry repeats alloy 12's.
Usage: python validate_insnsb_vassiliev.py [path/to/db.tdb]
"""
import csv
import os
import sys
import warnings

import numpy as np
from pycalphad import Database, equilibrium, variables as v

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
TDB = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "databases", "GaInSnSb.tdb")
DATA = os.path.join(HERE, "..", "data", "insnsb_vassiliev", "table9_mu_In_750K.csv")
T = 750.0
CORRECTED = {19: -16.62}


def rows():
    with open(DATA) as f:
        for r in csv.DictReader(line for line in f if not line.startswith("#")):
            yield r


def main():
    db = Database(TDB)
    g_in = equilibrium(db, ["IN", "VA"], ["LIQUID"], {v.T: T, v.P: 101325, v.N: 1}).MU.sel(component="IN").values.item()
    res, out = [], []
    for r in rows():
        n = int(r["alloy_no"])
        if r["measured_in_parentheses"] == "yes":
            continue
        meas = CORRECTED.get(n, float(r["mu_In_meas_kJ_mol"]))
        e = equilibrium(db, ["IN", "SB", "SN", "VA"], ["LIQUID"],
                        {v.T: T, v.P: 101325, v.N: 1, v.X("IN"): float(r["x_In"]), v.X("SB"): float(r["x_Sb"])})
        calc = (e.MU.sel(component="IN").values.item() - g_in) / 1000
        res.append(calc - meas)
        out.append((n, float(r["x_In"]), float(r["x_Sb"]), float(r["x_Sn"]), meas, calc))
    for n, xi, xb, xs, m, c in out:
        print(f"alloy {n:2d}  In{xi:.3f} Sb{xb:.3f} Sn{xs:.3f}  meas {m:7.2f}  calc {c:7.2f}  d {c - m:+.2f} kJ/mol")
    res = np.array(res)
    print(f"{len(res)} alloys: RMS {np.sqrt(np.mean(res ** 2)):.2f} kJ/mol, mean {res.mean():+.2f}, max |d| {np.abs(res).max():.2f}")
    import json
    calc_v = [float(r["mu_In_calc_kJ_mol"]) - CORRECTED.get(int(r["alloy_no"]), float(r["mu_In_meas_kJ_mol"]))
              for r in rows() if r["measured_in_parentheses"] != "yes"]
    json.dump({"n": len(res), "rms_kJ": float(np.sqrt(np.mean(res ** 2))), "mean_kJ": float(res.mean()),
               "max_abs_kJ": float(np.abs(res).max()),
               "vassiliev_own_extrapolation_rms_kJ": float(np.sqrt(np.mean(np.square(calc_v))))},
              open(os.path.join(HERE, "..", "results", "insnsb_validation.json"), "w"), indent=1)
    return out


if __name__ == "__main__":
    main()
