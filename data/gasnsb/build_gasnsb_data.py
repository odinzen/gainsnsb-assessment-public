"""Rebuild the Ga-Sb-Sn and In-Sb-Sn validation CSVs.

Sources
  Gerdes & Predel, J. Less-Common Met. 79 (1981) 281-288, doi:10.1016/0022-5088(81)90077-1
    Table 1 (p. 282) liquidus temperatures, Figs. 1 and 2 (p. 283) DTA symbols.
  Katayama, Fukuda, Hattori & Maruyama, Thermochim. Acta 314 (1998) 175-181
    Table 1 (p. 178) emf vs T of the zirconia cell, Table 2 (p. 179) binary Redlich-Kister fits.

Figure symbols were located in the 300 dpi page scan of Gerdes p. 283 (rendered with PyMuPDF at
300 dpi, which is the native CCITT scan resolution) by an open-circle detector: a pixel is a centre
candidate when dark pixels are found at radius 3.5-7 px in at least 92 % of 48 directions and the
5x5 box around it is at most 45 % dark; candidates are grouped and averaged. Every accepted centre
below was checked by eye on 3-8x zooms; detector hits on dashed-curve fragments were rejected.
Two symbols clipped by the left axis were placed by hand (noted). Axis calibration is a straight-line
least-squares fit to the printed gridlines (pixel rows/columns of the drawn lines).

Run with the calphad conda env:  python build_gasnsb_data.py
"""
import csv
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

GP = ("Gerdes & Predel 1981, J. Less-Common Met. 79, 281-288, doi:10.1016/0022-5088(81)90077-1")
KAT = ("Katayama, Fukuda, Hattori & Maruyama 1998, Thermochim. Acta 314, 175-181, "
       "PII S0040-6031(97)00464-4")

F = 96485.33212   # C/mol
R = 8.314462618   # J/(mol K)
K0 = 273.15

# ---------------------------------------------------------------- Gerdes & Predel, Table 1
# x_Sn and T(K) exactly as printed (decimal commas converted). Read from the page image; the
# PDF text layer has OCR errors here (InSb-Sn 797 -> '191', 777 -> '177', 547 -> '541').
TABLE1 = {
    "GaSb-Sn": [(0, 986), (0.04, 973), (0.05, 971), (0.06, 956), (0.10, 941), (0.15, 930),
                (0.20, 914), (0.25, 893), (0.30, 862), (0.40, 827), (0.50, 795), (0.60, 752),
                (0.70, 689), (0.80, 652), (0.85, 621), (0.90, 607)],
    "InSb-Sn": [(0, 797), (0.02, 790), (0.04, 783), (0.06, 779), (0.10, 777), (0.15, 769),
                (0.20, 753), (0.25, 742), (0.30, 728), (0.35, 716), (0.40, 692), (0.45, 675),
                (0.50, 660), (0.55, 648), (0.60, 628), (0.65, 611), (0.70, 606), (0.75, 547),
                (0.80, 546), (0.82, 531), (0.85, 529), (0.86, 521)],
}

# ---------------------------------------------------------------- Gerdes & Predel, Figs. 1, 2
# Gridline pixel positions in the 300 dpi scan of p. 283 (x = column, y = row).
CAL = {
    "Fig. 1": dict(xpx=[303.5, 460.5, 616.0, 773.0, 930.5, 1086.0], xv=[0, .2, .4, .6, .8, 1],
                   ypx=[692.5, 848.3, 1005.3, 1162.6, 1318.6, 1473.5],
                   yv=[1000, 900, 800, 700, 600, 500]),
    # The bottom border of Fig. 2 (400 K, row 2417) sits 5 px off the even spacing of the
    # 500-800 K gridlines; it is outside the data range and left out of the fit.
    "Fig. 2": dict(xpx=[304.5, 464.5, 623.5, 780.3, 939.5, 1095.5], xv=[0, .2, .4, .6, .8, 1],
                   ypx=[1796.5, 1953.5, 2109.5, 2266.0], yv=[800, 700, 600, 500]),
}

# (px, py, arrest, confidence, note). Arrest assignment: symbols on or near the horizontal line
# at the foot of the diagram are the eutectic halts; all others are liquidus points.
SYMBOLS = {
    ("GaSb-Sn", "Fig. 1"): [
        (302.5, 721.0, "liquidus", "medium", "clipped by the left axis (pure GaSb); centre by hand"),
        (337.6, 742.2, "liquidus", "high", ""),
        (377.7, 765.4, "liquidus", "high", ""),
        (417.6, 784.7, "liquidus", "high", ""),
        (456.2, 805.3, "liquidus", "high", ""),
        (495.7, 825.8, "liquidus", "high", ""),
        (536.9, 859.3, "liquidus", "high", ""),
        (617.2, 908.7, "liquidus", "high", ""),
        (693.7, 958.4, "liquidus", "high", ""),
        (773.6, 1018.6, "liquidus", "high", ""),
        (852.0, 1087.8, "liquidus", "high", ""),
        (931.0, 1169.8, "liquidus", "high", "on the calculated curves"),
        (969.6, 1234.9, "liquidus", "high", ""),
        (1008.5, 1283.5, "liquidus", "high", ""),
        (319.9, 1458.3, "eutectic", "high", ""),
        (351.1, 1461.5, "eutectic", "high", ""),
        (383.5, 1461.0, "eutectic", "high", ""),
        (462.5, 1462.6, "eutectic", "high", "on the x=0.2 gridline"),
        (617.5, 1463.5, "eutectic", "high", "on the x=0.4 gridline"),
        (695.5, 1464.5, "eutectic", "high", ""),
        (773.8, 1463.9, "eutectic", "high", "on the x=0.6 gridline"),
        (849.2, 1461.4, "eutectic", "high", ""),
        (930.9, 1467.7, "eutectic", "high", "on the x=0.8 gridline"),
        (1009.6, 1464.5, "eutectic", "high", ""),
        (1047.2, 1460.3, "eutectic", "high", ""),
        (1083.0, 1467.0, "eutectic", "medium", "clipped by the right axis (pure Sn); halt of Sn itself"),
    ],
    ("InSb-Sn", "Fig. 2"): [
        (308.5, 1802.5, "liquidus", "medium", "clipped by the left axis and top border (pure InSb); centre by hand"),
        (344.6, 1805.7, "liquidus", "high", ""),
        (382.9, 1824.7, "liquidus", "high", ""),
        (424.4, 1844.6, "liquidus", "high", ""),
        (462.3, 1863.3, "liquidus", "high", ""),
        (502.3, 1884.7, "liquidus", "high", ""),
        (541.7, 1908.4, "liquidus", "high", ""),
        (582.7, 1927.3, "liquidus", "high", ""),
        (622.3, 1967.3, "liquidus", "high", ""),
        (663.2, 1990.0, "liquidus", "high", ""),
        (702.9, 2012.7, "liquidus", "medium", "partly under the calculated curve"),
        (739.5, 2035.5, "liquidus", "high", ""),
        (778.5, 2059.0, "liquidus", "high", ""),
        (817.7, 2091.3, "liquidus", "high", ""),
        (856.3, 2104.6, "liquidus", "high", "on the 600 K gridline"),
        (936.7, 2172.3, "liquidus", "high", "on the x=0.8 gridline"),
        (976.5, 2219.5, "liquidus", "high", ""),
        (423.5, 2248.5, "eutectic", "high", ""),
        (464.6, 2247.2, "eutectic", "high", "on the x=0.2 gridline"),
        (502.5, 2248.5, "eutectic", "high", ""),
        (541.0, 2248.4, "eutectic", "high", ""),
        (581.5, 2247.5, "eutectic", "high", ""),
        (621.6, 2259.4, "eutectic", "high", "below the drawn eutectic line"),
        (662.7, 2258.3, "eutectic", "high", "below the drawn eutectic line"),
        (699.7, 2259.6, "eutectic", "high", "below the drawn eutectic line"),
        (779.3, 2249.7, "eutectic", "high", "on the x=0.6 gridline"),
        (858.1, 2257.0, "eutectic", "high", "below the drawn eutectic line"),
        (897.8, 2249.0, "eutectic", "high", ""),
        (937.5, 2255.5, "eutectic", "high", "below the drawn eutectic line"),
        (977.9, 2253.3, "eutectic", "high", "below the drawn eutectic line"),
        (1091.8, 2256.4, "eutectic", "medium", "clipped by the right axis (pure Sn); halt of Sn itself"),
    ],
}

# Drawn lines (not measurements): rows of the horizontal eutectic line and the composition where
# the drawn liquidus meets it, read on 4x zooms.
DRAWN = {
    "GaSb-Sn": dict(fig="Fig. 1", eut_row=1466.5, liq_end_px=1084.0),
    "InSb-Sn": dict(fig="Fig. 2", eut_row=2250.5, liq_end_px=1005.0),
}

# ---------------------------------------------------------------- Katayama et al., Table 1
# section y = x_Sb/(x_Sb+x_Sn); rows: x_Ga, a, b, c (E/mV = a + b*T +/- c), E(1073 K)/mV, +/-,
# a_Ga, +/-, gamma_Ga, alpha_Ga, all as printed.
KAT_T1 = {
    0.75: [(0.1, -4.54, 0.0819, 0.12, 83.34, 0.12, 0.0669, 0.0003, 0.669, -0.496),
           (0.3, 9.71, 0.0321, 0.20, 44.15, 0.20, 0.239, 0.002, 0.797, -0.463),
           (0.5, -29.57, 0.0487, 0.26, 22.69, 0.26, 0.479, 0.005, 0.958, -0.172),
           (0.75, -21.63, 0.0278, 0.20, 8.20, 0.20, 0.766, 0.006, 1.021, 0.333),
           (0.9, -17.84, 0.0194, 0.09, 2.96, 0.09, 0.908, 0.003, 1.009, 0.896)],
    0.50: [(0.1, -50.73, 0.1140, 0.35, 71.59, 0.35, 0.0979, 0.0012, 0.979, -0.026),
           (0.3, -21.12, 0.0564, 0.16, 39.40, 0.16, 0.278, 0.002, 0.927, -0.155),
           (0.5, -11.84, 0.0313, 0.19, 21.74, 0.19, 0.494, 0.004, 0.988, -0.048),
           (0.75, -36.10, 0.0401, 0.26, 6.93, 0.26, 0.799, 0.007, 1.065, 1.008),
           (0.9, -16.21, 0.0178, 0.11, 2.89, 0.11, 0.910, 0.004, 1.011, 1.094)],
    0.25: [(0.1, -56.03, 0.1136, 0.24, 65.86, 0.24, 0.118, 0.001, 1.180, 0.204),
           (0.3, -25.38, 0.0548, 0.18, 33.21, 0.18, 0.340, 0.003, 1.133, 0.255),
           (0.5, -5.75, 0.0239, 0.09, 19.89, 0.09, 0.524, 0.002, 1.048, 0.188),
           (0.75, -23.15, 0.0289, 0.19, 7.86, 0.19, 0.775, 0.005, 1.033, 0.519),
           (0.9, -18.09, 0.0190, 0.13, 2.30, 0.13, 0.928, 0.004, 1.031, 3.053)],
}
# Approximate measured temperature span per alloy, read from the symbols of Fig. 2 (p. 178).
# T = 1000 + (px - 356)/4.655 at 400 dpi; extreme-symbol centring +/-3 K.
KAT_TSPAN = {
    (0.75, 0.1): (1050, 1150), (0.50, 0.1): (1071, 1158), (0.25, 0.1): (1034, 1159),
    (0.75, 0.3): (1037, 1163), (0.50, 0.3): (1039, 1150), (0.25, 0.3): (1027, 1150),
    (0.75, 0.5): (1032, 1156), (0.50, 0.5): (1039, 1158), (0.25, 0.5): (1050, 1165),
    (0.75, 0.75): (1037, 1153), (0.50, 0.75): (1040, 1156), (0.25, 0.75): (1033, 1155),
    (0.75, 0.9): (1052, 1143), (0.50, 0.9): (1029, 1126), (0.25, 0.9): (1055, 1145),
}

# Table 2 (p. 179): G^xs = X_i X_j sum_v (X_i - X_j)^v L^(v), L^(v) = a + b T, 1000-1200 K.
# Units are not printed; J/mol is implied by the magnitudes.
KAT_T2 = [("Ga", "Sb", 0, -832.80, -6.6885), ("Ga", "Sb", 1, 5598.30, -6.6301),
          ("Ga", "Sb", 2, -120.56, -1.1081), ("Sn", "Ga", 0, 16792.00, -12.1380),
          ("Sn", "Ga", 1, -15241.00, 13.0650), ("Sn", "Ga", 2, 5651.50, -5.2132),
          ("Sb", "Sn", 0, -4855.50, -1.6688), ("Sb", "Sn", 1, 1092.70, -1.1787),
          ("Sb", "Sn", 2, 1752.80, -0.7129)]


def pseudo_comp(system, x_sn):
    """Atom fractions on AB-Sn: x_A = x_B = (1 - x_Sn)/2 (Gerdes Eqs. 12-14)."""
    half = (1 - x_sn) / 2
    xga, xin = (half, 0.0) if system == "GaSb-Sn" else (0.0, half)
    x_compound = (1 - x_sn) / (1 + x_sn)   # mole fraction of AB formula units
    return xga, xin, half, x_sn, x_compound


def calib(fig):
    c = CAL[fig]
    px = np.polyfit(c["xpx"], c["xv"], 1)
    py = np.polyfit(c["ypx"], c["yv"], 1)
    rx = np.polyval(px, c["xpx"]) - c["xv"]
    ry = np.polyval(py, c["ypx"]) - c["yv"]
    return px, py, abs(rx).max(), abs(ry).max()


def write(name, header, cols, rows):
    with open(HERE / name, "w", newline="", encoding="utf-8") as f:
        for line in header:
            f.write("# " + line + "\n")
        w = csv.writer(f)
        w.writerow(cols)
        w.writerows(rows)
    print(f"{name}: {len(rows)} rows")


def gerdes():
    fits = {fig: calib(fig) for fig in CAL}
    digit = []
    for (system, fig), syms in SYMBOLS.items():
        px, py, _, _ = fits[fig]
        for x_pix, y_pix, arrest, conf, note in syms:
            x = float(np.clip(np.polyval(px, x_pix), 0, 1))
            T = float(np.polyval(py, y_pix))
            digit.append((system, fig, arrest, x, T, x_pix, y_pix, conf, note))

    # Table 1 rows, each matched to a liquidus symbol at the same x (|dx| <= 0.006) if one exists.
    rows = []
    for system, data in TABLE1.items():
        liq = [d for d in digit if d[0] == system and d[2] == "liquidus"]
        for x_sn, T in data:
            xga, xin, xsb, xsn, xc = pseudo_comp(system, x_sn)
            near = [d for d in liq if abs(d[3] - x_sn) <= 0.006]
            if near:
                d = min(near, key=lambda d: abs(d[3] - x_sn))
                tf, dt = f"{d[4]:.1f}", f"{d[4] - T:+.1f}"
            else:
                tf, dt = "", ""
            if not near:
                flag = f"no symbol at this x in {'Fig. 1' if system == 'GaSb-Sn' else 'Fig. 2'}"
            elif system == "GaSb-Sn" and x_sn >= 0.10:
                flag = "suspect: Fig. 1 pairs this x with the T of the next row down (see notes.md)"
            elif system == "InSb-Sn" and x_sn in (0.10, 0.80):
                flag = "Fig. 2 symbol differs by more than 5 K (see notes.md)"
            else:
                flag = ""
            rows.append((system, f"{x_sn:.2f}", T, f"{T - K0:.2f}", f"{xga:.4f}", f"{xin:.4f}",
                         f"{xsb:.4f}", f"{xsn:.4f}", f"{xc:.4f}", "liquidus", tf, dt, flag))
    write("gerdes1981_table1_liquidus_TABULATED.csv", [
        "TABULATED. " + GP,
        "Table 1 (Tabelle 1, p. 282): liquidus temperatures by DTA on the pseudobinary sections GaSb-Sn and InSb-Sn.",
        "  Pt/Pt-Rh thermocouple calibrated on In, Sn, Pb, Bi, Sb, Al, Ag, Au melting points; corundum crucibles; Ar.",
        "  The table lists LIQUIDUS only. The eutectic halts are in the figures only (see the DIGITIZED file).",
        "  AlSb-Sn (third column of Table 1) is not transcribed here.",
        "Units: x_Sn_printed as printed (decimal commas -> points); T_K as printed (integer K); T_C = T_K - 273.15.",
        "Composition basis: x_Sn is the ATOM fraction of Sn. Gerdes Eqs. (12)-(14) set x_A = x_B = (1-x)/2, x_C = x",
        "  for section AB-C. So x_Ga (or x_In) = x_Sb = (1 - x_Sn)/2. X_AB_formula = (1-x_Sn)/(1+x_Sn) is the mole",
        "  fraction of AB formula units, given only for comparison with papers that use mol% GaSb or InSb.",
        "Read from the page image; the PDF text layer has OCR errors in this table (797->191, 777->177, 547->541).",
        "fig_symbol_T_K: the digitized Fig. 1/2 liquidus symbol at the same x (|dx|<=0.006), dT = figure - table.",
        "WARNING: for GaSb-Sn, x >= 0.10, Table 1 and Fig. 1 disagree by 11-58 K at the same x; see notes.md.",
    ], ["system", "x_Sn_printed", "T_K", "T_C", "x_Ga", "x_In", "x_Sb", "x_Sn", "X_AB_formula",
        "arrest", "fig_symbol_T_K", "dT_fig_minus_table_K", "flag"], rows)

    cal_lines = []
    for fig, (px, py, rx, ry) in fits.items():
        cal_lines.append(f"  {fig}: {0.1 / px[0]:.1f} px per 0.1 x_Sn, {-1 / py[0]:.3f} px per K; "
                         f"gridline residual max {rx:.4f} in x, {ry:.2f} K.")
    rows = []
    for system, fig, arrest, x, T, x_pix, y_pix, conf, note in digit:
        xga, xin, xsb, xsn, xc = pseudo_comp(system, x)
        rows.append((system, fig, arrest, f"{x:.4f}", f"{T:.1f}", f"{T - K0:.1f}", f"{xga:.4f}",
                     f"{xin:.4f}", f"{xsb:.4f}", f"{xsn:.4f}", f"{xc:.4f}", f"{x_pix:.1f}",
                     f"{y_pix:.1f}", conf, note))
    write("gerdes1981_fig1_fig2_symbols_DIGITIZED.csv", [
        "DIGITIZED. " + GP,
        "Fig. 1 (GaSb-Sn) and Fig. 2 (InSb-Sn), p. 283: experimental DTA symbols (open circles) only.",
        "  The solid (smoothed experimental), dashed (regular solution) and dash-dot (association model) curves were",
        "  NOT digitized as data. The authors state the eutectic is the only three-phase reaction in each section.",
        "arrest: 'liquidus' or 'eutectic' (halts on or just below the horizontal line at the foot of the diagram).",
        "Units: x_Sn atom fraction (same axis as Table 1); T_K read in kelvin (axis 'Temperatur in K'); T_C = T_K - 273.15.",
        "  x_Ga (or x_In) = x_Sb = (1 - x_Sn)/2; X_AB_formula = (1-x_Sn)/(1+x_Sn).",
        "Calibration: straight-line fit to the drawn gridlines in the 300 dpi scan (px, py = symbol centre).",
        *cal_lines,
        "Uncertainty: 'high' about +/-1.5 K and +/-0.003 in x_Sn; 'medium' (symbol clipped or under a curve) +/-3 K, +/-0.006.",
        "  This is reading error only; the paper gives no experimental uncertainty and notes scatter from undercooling.",
        "Symbols at x_Sn = 1 are the pure-Sn halt (Sn m.p. 505.08 K) and check the temperature scale.",
    ], ["system", "figure", "arrest", "x_Sn", "T_K", "T_C", "x_Ga", "x_In", "x_Sb", "x_Sn_atom",
        "X_AB_formula", "px", "py", "confidence", "note"], rows)

    rows = []
    for system, dr in DRAWN.items():
        px, py, _, _ = fits[dr["fig"]]
        eut = [d for d in digit if d[0] == system and d[2] == "eutectic" and d[3] < 0.99]
        Ts = np.array([d[4] for d in eut])
        t_line = float(np.polyval(py, dr["eut_row"]))
        x_end = float(np.polyval(px, dr["liq_end_px"]))
        sn = [d for d in digit if d[0] == system and d[2] == "eutectic" and d[3] >= 0.99][0]
        rows += [
            (system, dr["fig"], "eutectic T, drawn horizontal line", f"{t_line:.1f}", f"{t_line - K0:.1f}",
             "", "drawn line, authors' assignment; +/-1 K reading"),
            (system, dr["fig"], "eutectic T, mean of halt symbols x_Sn<1", f"{Ts.mean():.1f}",
             f"{Ts.mean() - K0:.1f}", "",
             f"n={len(Ts)}, sd={Ts.std(ddof=1):.1f} K, range {Ts.min():.1f}-{Ts.max():.1f} K"),
            (system, dr["fig"], "pure Sn halt symbol", f"{sn[4]:.1f}", f"{sn[4] - K0:.1f}", "1.0000",
             "clipped symbol at the Sn axis; Sn m.p. is 505.08 K"),
            (system, dr["fig"], "x_Sn where drawn liquidus meets eutectic line", "", "", f"{x_end:.3f}",
             "drawn curve, not a measurement; +/-0.01"),
        ]
    write("gerdes1981_invariants_DIGITIZED.csv", [
        "DIGITIZED. " + GP,
        "Figs. 1 and 2, p. 283. The paper prints no invariant temperature or eutectic composition; the text says",
        "  only that one eutectic is the sole three-phase reaction in each section. Values below are read from",
        "  the figures. Rows marked 'drawn' come from the authors' drawn lines, not from data symbols.",
        "Units: T in K and degC (T_C = T_K - 273.15); x_Sn atom fraction.",
    ], ["system", "figure", "quantity", "T_K", "T_C", "x_Sn", "basis"], rows)


def katayama():
    rows = []
    worst = dict(E=0.0, a=0.0, g=0.0, al=0.0)
    for y, data in KAT_T1.items():
        for xga, a, b, c, E, dE, aga, da, g, al in data:
            xsb, xsn = (1 - xga) * y, (1 - xga) * (1 - y)
            E_re = a + b * 1073
            a_re = math.exp(-3 * F * E * 1e-3 / (R * 1073))
            g_re = aga / xga
            al_re = math.log(g) / (1 - xga) ** 2
            worst["E"] = max(worst["E"], abs(E_re - E))
            worst["a"] = max(worst["a"], abs(a_re - aga) / aga)
            worst["g"] = max(worst["g"], abs(g_re - g))
            worst["al"] = max(worst["al"], abs(al_re - al))
            dG = -3 * F * E * 1e-3
            H = -3 * F * a * 1e-3
            S = 3 * F * b * 1e-3
            t0, t1 = KAT_TSPAN[(y, xga)]
            rows.append((f"{y:.2f}", xga, f"{xsb:.4f}", f"{xsn:.4f}", a, b, c, 1073, E, dE, aga, da,
                         g, al, f"{E_re:.2f}", f"{a_re:.4f}", f"{dG:.0f}", f"{dG - R * 1073 * math.log(xga):.0f}",
                         f"{H:.0f}", f"{S:.3f}", t0, t1))
    write("katayama1998_table1_emf_TABULATED.csv", [
        "TABULATED. " + KAT,
        "Table 1 (p. 178): emf of cell (-) W, Ga(l), Ga2O3 | ZrO2(+Y2O3) | Ga-Sb-Sn(l), Ga2O3 (+) W, and a_Ga at 1073 K.",
        "  Electrolyte 0.92 ZrO2 - 0.08 Y2O3; Ar atmosphere; stated range 1050-1150 K.",
        "Sections: Ga_x (Sb_y Sn_1-y)_1-x with y = 0.75, 0.50, 0.25 (section_y = x_Sb/(x_Sb+x_Sn)) and",
        "  x_Ga = 0.1, 0.3, 0.5, 0.75, 0.9. x_Sb = (1-x_Ga)*y, x_Sn = (1-x_Ga)*(1-y).",
        "Printed columns: emf_a_mV, emf_b_mV_per_K, emf_c_mV in E/mV = a + b*T +/- c (T in K, least squares);",
        "  E_1073_mV +/-; a_Ga +/-; gamma_Ga = a_Ga/x_Ga; alpha_Ga = ln(gamma_Ga)/(1-x_Ga)^2, all at 1073 K.",
        "Reference state: pure LIQUID Ga at the same T (reference electrode Ga(l) + Ga2O3). n = 3 electrons per Ga",
        "  (Ga2O3 + 6e- = 2 Ga + 3 O2-). Paper Eq. (1): -3EF = RT ln a_Ga = mu_Ga - G_Ga(liq).",
        "Recomputed (this file, not printed): E_1073_recalc = a + 1073 b; a_Ga_recalc = exp(-3 F E/(R T)) from printed E;",
        "  dmu_Ga_J = -3 F E (J/mol, = mu_Ga - G_Ga(liq) at 1073 K); dmu_Ga_xs_J = dmu_Ga - R T ln x_Ga;",
        "  dH_Ga_J = -3 F a and dS_Ga_J_per_K = +3 F b: the partial enthalpy and entropy of mixing of Ga implied by",
        "  the linear E(T), with a in V and b in V/K. F = 96485.33 C/mol, R = 8.314463 J/(mol K).",
        f"Transcription check vs printed values: max |E_recalc - E| = {worst['E']:.3f} mV; max rel. a_Ga error "
        f"{100 * worst['a']:.2f} %; max |a_Ga/x - gamma| = {worst['g']:.4f}; max |alpha recalc - alpha| = {worst['al']:.3f}.",
        "T_min_K, T_max_K: approximate span of the Fig. 2 data symbols for that alloy (DIGITIZED, +/-3 K).",
        "  Use E(T) only inside that span. Individual E,T points are plotted in Fig. 2 but not tabulated; not digitized.",
    ], ["section_y", "x_Ga", "x_Sb", "x_Sn", "emf_a_mV", "emf_b_mV_per_K", "emf_c_mV", "T_K", "E_1073_mV",
        "E_1073_err_mV", "a_Ga", "a_Ga_err", "gamma_Ga", "alpha_Ga", "E_1073_recalc_mV", "a_Ga_recalc",
        "dmu_Ga_J", "dmu_Ga_xs_J", "dH_Ga_J", "dS_Ga_J_per_K", "T_min_K", "T_max_K"], rows)

    write("katayama1998_table2_binary_RK_TABULATED.csv", [
        "TABULATED. " + KAT,
        "Table 2 (p. 179): Redlich-Kister fits of the authors' binary liquid G^xs, used for their Chou-model estimate.",
        "  G^xs = X_i X_j sum_v (X_i - X_j)^v L^(v)(T), L^(v) = a + b T. Fitted to G^xs computed at 50 K steps, 1000-1200 K.",
        "Units are not printed; J/mol of atoms is implied (b in J/(mol K)). Reference: pure liquids.",
        "Note the i-j order: Ga-Sb, Sn-Ga (not Ga-Sn), Sb-Sn. A TDB with Ga before Sn needs odd-order terms sign-flipped.",
        "Sources per the paper: Ga-Sb and Ga-Sn from Katayama's earlier emf work (ref. [1]); Sb-Sn from ref. [11]",
        "  (Tanaka et al., then 'to be published'). The text cites both [1] and [11] for Sb-Sn.",
    ], ["i", "j", "v", "a_J", "b_J_per_K", "T_range_K"],
        [(i, j, v, a, b, "1000-1200") for i, j, v, a, b in KAT_T2])


if __name__ == "__main__":
    gerdes()
    katayama()
