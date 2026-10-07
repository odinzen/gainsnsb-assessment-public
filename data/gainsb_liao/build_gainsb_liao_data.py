"""Rebuild the Ga-In-Sb data CSVs extracted from Liao et al. (1982) and Ansara et al. (1976).

Figure data are digitized from the 300 dpi scans of Liao et al., Calphad 6 (1982) 141.
Symbol centres were located by morphological opening of the binarized page (which removes
the thin calculated curves and keeps the filled symbols) and, where a symbol touches or is
partly hidden by a curve, by hand from 4x zooms with a 10 px grid. Those hand positions are
recorded below so the CSVs can be regenerated and audited. Axis calibration is a
least-squares affine fit to the tick marks, because the scans carry a small shear/rotation
(about 10 px over the plot height).

Run with the calphad conda env:  python build_gainsb_liao_data.py
"""
import csv
import os
from pathlib import Path

import numpy as np
import pymupdf

HERE = Path(__file__).resolve().parent
# The source PDF is not distributed; point GAINSNSB_PDF_DIR at a folder holding it.
LIAO_PDF = Path(os.environ.get("GAINSNSB_PDF_DIR", ".")) / "Liao_1982_Calphad_GaInSb_associated_solution.pdf"

LIAO = "Liao, Su, Tung & Brebrick, Calphad 6 (1982) 141-169, doi:10.1016/0364-5916(82)90009-8"
REF29 = "Wooley & Lees 1959, J. Less-Common Metals 1, 192 (Liao ref. 29)"
REF30 = "Blom & Plaskett 1971, J. Electrochem. Soc. 118, 1831 (Liao ref. 30)"
REF31 = "Antypas 1972, J. Crystal Growth 16, 181 (Liao ref. 31)"
REF32 = ("Miki, Segawa, Otsubo, Shirahata & Kiyibayashi [sic] 1975, Proc. 5th Int. Symp. GaAs "
         "and Related Compounds, Inst. Phys. Conf. Ser. 24 (Liao ref. 32)")


def affine(points):
    """Least-squares val = c0 + c1*px + c2*py from (px, py, val) control points."""
    a = np.array([[1.0, p[0], p[1]] for p in points])
    v = np.array([p[2] for p in points])
    c, *_ = np.linalg.lstsq(a, v, rcond=None)
    return c, v - a @ c


def apply(c, px, py):
    return c[0] + c[1] * px + c[2] * py


def write_csv(name, header_lines, columns, rows):
    path = HERE / name
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        for line in header_lines:
            fh.write(f"# {line}\n" if line else "#\n")
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(columns)
        w.writerows(["" if v is None else v for v in r] for r in rows)
    print(f"wrote {path.name}: {len(rows)} rows")


# ---------------------------------------------------------------- Fig. 12, p. 163
def fig12():
    top = [(800.0, .1), (1067.5, .3), (1201.5, .4), (1336.0, .5), (1471.0, .6),
           (1607.0, .7), (1742.0, .8), (1875.0, .9)]
    bot = [(790.5, .1), (1057.5, .3), (1192.0, .4), (1326.5, .5), (1462.0, .6),
           (1597.0, .7), (1732.0, .8), (1865.0, .9)]
    cx, rx = affine([(p, 420, v) for p, v in top] + [(p, 1760, v) for p, v in bot])
    left = [(402.5, 750), (676.5, 700), (949.5, 650), (1228.5, 600), (1772.5, 500)]
    right = [(403.5, 750), (678, 700), (951, 650), (1229.5, 600), (1503.5, 550), (1776, 500)]
    cy, ry = affine([(679, p, v) for p, v in left] + [(1985, p, v) for p, v in right])

    method = {"square": ("DTA", REF29), "circle": ("thermal analysis", REF29),
              "triangle": ("X-ray diffraction", REF29), "diamond": ("(not stated by Liao)", REF30)}
    # symbol, px, py, boundary, confidence, note
    pts = [
        ("square", 791.5, 1524.0, "liquidus", "high", ""),
        ("square", 792.0, 1601.0, "solidus", "high", ""),
        ("circle", 925.5, 1377.0, "liquidus", "high", ""),
        ("circle", 927.3, 1582.0, "solidus", "medium", "overlaps a triangle; centre placed by hand"),
        ("triangle", 952.0, 1565.0, "solidus", "medium", "overlaps a circle; centre placed by hand"),
        ("diamond", 981.5, 1227.0, "liquidus", "medium", "partly under calculated liquidus; by hand"),
        ("square", 1061.5, 1222.5, "liquidus", "high", ""),
        ("square", 1059.0, 1547.0, "solidus", "high", ""),
        ("triangle", 1160.5, 1521.0, "solidus", "high", "centre placed by hand"),
        ("circle", 1197.0, 1108.5, "liquidus", "high", ""),
        ("circle", 1193.5, 1505.5, "solidus", "high", "touches curve; centre placed by hand"),
        ("square", 1329.0, 983.0, "liquidus", "high", ""),
        ("square", 1328.5, 1438.0, "solidus", "high", ""),
        ("circle", 1465.0, 901.5, "liquidus", "high", ""),
        ("circle", 1463.5, 1383.0, "solidus", "high", ""),
        ("triangle", 1518.75, 1378.0, "solidus", "medium", "lower of two overlapping triangles"),
        ("triangle", 1531.25, 1355.25, "solidus", "medium", "upper of two overlapping triangles"),
        ("diamond", 1560.7, 825.0, "liquidus", "medium", "abuts a square; centre placed by hand"),
        ("square", 1601.0, 819.0, "liquidus", "high", ""),
        ("square", 1602.0, 1296.5, "solidus", "high", ""),
        ("square", 1735.5, 753.5, "liquidus", "high", ""),
        ("square", 1733.0, 1191.0, "solidus", "high", ""),
        ("square", 1872.5, 679.0, "liquidus", "high", ""),
        ("triangle", 1867.7, 1009.3, "solidus", "high", "centre placed by hand"),
        ("triangle", 1942.0, 847.5, "solidus", "high", "centre placed by hand"),
        ("square", 1951.5, 792.0, "solidus", "medium",
         "no paired liquidus point at this x; assigned to solidus because it lies on the solidus branch"),
    ]
    rows = []
    for sym, px, py, bnd, conf, note in pts:
        x = apply(cx, px, py)
        t = apply(cy, px, py)
        m, src = method[sym]
        rows.append((src, m, sym, bnd, f"{t:.1f}", f"{x:.4f}", f"{x / 2:.4f}", f"{(1 - x) / 2:.4f}",
                     "0.5000", px, py, conf, note))
    hdr = [
        "DIGITIZED from Liao et al. (1982) Fig. 12 (p. 163), GaSb-InSb pseudobinary liquidus and solidus.",
        f"Figure source: {LIAO}",
        "Caption: squares DTA, circles thermal analysis, triangles X-ray diffraction, all ref. 29; diamonds ref. 30.",
        "Only plotted SYMBOLS were digitized. The calculated liquidus/solidus curves were not.",
        "Columns: source = original data source as cited by Liao; method = per Liao caption; symbol;",
        "  boundary = liquidus or solidus (the higher of a same-method pair at one x is the liquidus; triangles are",
        "  XRD solid compositions); T_C in degC; x_GaSb = mole fraction GaSb on the pseudobinary;",
        "  x_Ga, x_In, x_Sb = atom fractions of that pseudobinary composition (liquid for liquidus rows, solid for",
        "  solidus rows); px, py = symbol centre in the 300 dpi page scan; confidence; note.",
        "Calibration: affine least squares on 16 x-ticks (top and bottom edges) and 11 T-ticks (both sides).",
        f"  x residual rms {np.sqrt((rx ** 2).mean()):.4f}, max {abs(rx).max():.4f};"
        f" T residual rms {np.sqrt((ry ** 2).mean()):.2f} degC, max {abs(ry).max():.2f} degC.",
        "  Scale: 134.5 px per 0.1 in x, 5.49 px per degC.",
        "Digitizing uncertainty (symbol centring + calibration): about +/-1 degC and +/-0.004 in x_GaSb for",
        "  'high' rows; +/-2 degC and +/-0.008 for 'medium' rows (overlapping or curve-obscured symbols).",
        "Liao does not state how many of these points entered sigma_L and sigma_S (Table V).",
    ]
    cols = ["source", "method", "symbol", "boundary", "T_C", "x_GaSb", "x_Ga", "x_In", "x_Sb",
            "px", "py", "confidence", "note"]
    write_csv("gainsb_pseudobinary_liquidus_solidus_fig12_DIGITIZED.csv", hdr, cols, rows)


# ---------------------------------------------------------------- Fig. 13, p. 164
def triangle_vertices():
    """Gibbs-triangle vertices from line fits of the three drawn edges (see module docstring)."""
    # Edge centre-line fits on the 300 dpi scan of p. 164 (robust linear fits, 194-216 samples each):
    #   bottom: row = 2025.64 + 0.0012993*col ; left: col = 1540.03 - 0.598739*row ;
    #   right: col = 896.861 + 0.587777*row
    b0, b1 = 2025.64065, 1.29925038e-3
    l0, l1 = 1540.03017, -0.598738893
    r0, r1 = 896.861063, 0.587776925

    def bottom_and(c0, c1):
        row = (b1 * c0 + b0) / (1 - b1 * c1)
        return np.array([c1 * row + c0, row])

    in_v = bottom_and(l0, l1)
    ga_v = bottom_and(r0, r1)
    row = (r0 - l0) / (l1 - r1)
    sb_v = np.array([l1 * row + l0, row])
    return in_v, ga_v, sb_v


def fig13():
    in_v, ga_v, sb_v = triangle_vertices()
    m = np.column_stack([ga_v - in_v, sb_v - in_v])
    minv = np.linalg.inv(m)
    # symbol, px, py, isotherm degC, confidence, note
    pts = [
        ("triangle", 1038.5, 1156.5, 600, "medium", "on the 600 degC curve, upper branch"),
        ("triangle", 1341.5, 1253.0, 675, "high", ""),
        ("triangle", 986.5, 1290.0, 600, "low",
         "at the crossing of the 600 degC curves with the x_Sb=0.5 line; shape mostly hidden"),
        ("triangle", 1009.0, 1431.0, 600, "high", ""),
        ("triangle", 1582.0, 1475.0, 675, "low", "base visible below the 675 degC curve, apex hidden"),
        ("triangle", 792.0, 1635.0, 500, "medium", "adjacent to a circle"),
        ("triangle", 1575.5, 1717.5, 600, "medium", "apex hidden in the curve"),
        ("triangle", 1066.5, 1762.0, 500, "low", "apex visible, base merged with curves"),
        ("triangle", 1327.0, 1843.0, 500, "low", "dark patch on the 500 degC curves; shape not resolvable"),
        ("triangle", 670.5, 1884.5, 400, "medium", "apex hidden in the curve"),
        ("triangle", 1241.6, 1982.25, 300, "high", "apex hidden in the curve, base clear"),
        ("circle", 677.5, 1484.0, 500, "high", ""),
        ("circle", 757.5, 1588.0, 500, "high", ""),
        ("circle", 814.5, 1646.5, 500, "medium", "adjacent to a triangle"),
        ("circle", 947.5, 1727.5, 500, "medium", "sits on the curve"),
        ("circle", 1205.6, 1802.75, 500, "low", "merged with curves; round top visible"),
        ("circle", 1266.0, 1819.5, 500, "medium", "sits on the curve"),
        ("circle", 1403.75, 1870.6, 500, "high", ""),
    ]
    rows = []
    for sym, px, py, t, conf, note in pts:
        x_ga, x_sb = minv @ (np.array([px, py]) - in_v)
        src = REF30 if sym == "triangle" else REF31
        rows.append((src, sym, t, f"{x_ga:.4f}", f"{1 - x_ga - x_sb:.4f}", f"{x_sb:.4f}", px, py, conf, note))
    hdr = [
        "DIGITIZED from Liao et al. (1982) Fig. 13 (p. 164), ternary liquidus isotherms.",
        f"Figure source: {LIAO}",
        "Caption: triangles ref. 30, circles ref. 31. Only plotted SYMBOLS were digitized; the solid (Liao), dotted",
        "  (ref. 9) and dashed (ref. 30) curves were not.",
        "T_C is NOT read from an axis: Liao gives no temperature per point. It is the labelled isotherm",
        "  (300, 400, 500, 600 or 675 degC) on which the symbol is drawn. Every symbol lies on one of these curves.",
        "Columns: source as cited by Liao; symbol; T_C (isotherm); x_Ga, x_In, x_Sb = liquid atom fractions;",
        "  px, py = symbol centre in the 300 dpi scan; confidence; note.",
        "Calibration: affine map from the three triangle vertices, each the intersection of straight-line fits to",
        f"  the drawn edges. In=({in_v[0]:.1f},{in_v[1]:.1f}) Ga=({ga_v[0]:.1f},{ga_v[1]:.1f})"
        f" Sb=({sb_v[0]:.1f},{sb_v[1]:.1f}) px. Edge-fit scatter 0.6-1.2 px.",
        "  Cross-check: the 0.05 tick marks on the In-Ga edge sit 3-6 px (0.002-0.003) Ga-ward of the",
        "  vertex-based positions, an unresolved drawing offset carried into the uncertainty below.",
        "Uncertainty: about +/-0.004 in each atom fraction for 'high', +/-0.008 for 'medium', +/-0.015 for 'low'.",
        "  'low' rows are blobs on the calculated curves whose shape cannot be fully resolved; treat with care.",
        "Not included: two dark blobs on the 300 degC curves at about (462,1975) px (x_Ga 0.06, x_In 0.91,",
        "  x_Sb 0.035) and (690,2008) px (x_Ga 0.20, x_In 0.79, x_Sb 0.012). Both widen downward like a partly",
        "  hidden triangle and may be two more ref. 30 points at 300 degC, but they cannot be separated from",
        "  overlapping dashes, so they are reported here only. Check against Blom & Plaskett (1971) directly.",
        "Liao fitted these as the ternary liquidus (sigma_L(TER) = 11.6 degC, Table V); number of points not stated.",
    ]
    cols = ["source", "symbol", "T_C", "x_Ga", "x_In", "x_Sb", "px", "py", "confidence", "note"]
    write_csv("gainsb_ternary_liquidus_fig13_DIGITIZED.csv", hdr, cols, rows)


# ---------------------------------------------------------------- Figs. 14, 15
def tieline(name, fig, page, temp, src, xlabel, xt, yt, pts, extra):
    cx, rx = affine(xt)
    cy, ry = affine(yt)
    rows = []
    for px, py, note in pts:
        xl = apply(cx, px, py)
        ys = apply(cy, px, py)
        rows.append((src, temp, f"{xl:.4f}", None, None, f"{ys:.4f}", f"{1 - ys:.4f}", px, py, note))
    hdr = [
        f"DIGITIZED from Liao et al. (1982) Fig. {fig} (p. {page}), liquid-solid tie-lines at {temp} degC.",
        f"Figure source: {LIAO}",
        f"Experimental points (filled triangles) from {src}. Only symbols were digitized,",
        "  not Liao's solid curve or the dotted/dashed ref. 8 curve.",
        f"x axis of the figure: {xlabel} in the LIQUID. y axis: mole fraction GaSb in the (Ga,In)Sb solid.",
        "The figure does NOT give the liquid Sb content, so x_In_liq and x_Sb_liq are left blank. They must be",
        "  taken from the original paper; they are not derivable from this figure without a liquidus model.",
        *extra,
        "Columns: source; T_C; x_liq_axis = the figure's x-axis quantity; x_In_liq, x_Sb_liq (blank);",
        "  y_GaSb_solid; y_InSb_solid = 1 - y_GaSb_solid; px, py = symbol centre (300 dpi scan); note.",
        "Symbol centre = bounding-box centre of the filled triangle (height about 36 px), placed from the apex or",
        "  base edge, whichever is not hidden by the calculated curve.",
        "Calibration: affine least squares on axis ticks (both opposite edges for each axis).",
        f"  x residual rms {np.sqrt((rx ** 2).mean()):.4f}, max {abs(rx).max():.4f};"
        f" y residual rms {np.sqrt((ry ** 2).mean()):.4f}, max {abs(ry).max():.4f}.",
        "  Scale about 134 px per 0.1 (x) and 140 px per 0.1 (y).",
        "Digitizing uncertainty: about +/-0.004 in both coordinates (+/-0.006 where two triangles overlap).",
    ]
    cols = ["source", "T_C", "x_liq_axis", "x_In_liq", "x_Sb_liq", "y_GaSb_solid", "y_InSb_solid",
            "px", "py", "note"]
    write_csv(name, hdr, cols, rows)


def fig14():
    xt = ([(p, 495, v) for p, v in [(637.5, 0), (772.5, .1), (906, .2), (1174, .4), (1308, .5), (1442, .6),
                                    (1575, .7), (1708.5, .8), (1842.5, .9), (1977.5, 1.0)]]
          + [(p, 1855, v) for p, v in [(649, 0), (1584.5, .7), (1718, .8), (1987, 1.0)]])
    left = [(482, 1.0), (619, .9), (761, .8), (902.5, .7), (1038.5, .6), (1177.5, .5), (1318, .4),
            (1456.5, .3), (1596.5, .2), (1737, .1), (1878, 0)]
    right = [(472, 1.0), (608.5, .9), (750.5, .8), (890.5, .7), (1027.5, .6), (1166.5, .5), (1307, .4),
             (1447, .3), (1587, .2), (1728, .1), (1868.5, 0)]
    yt = [(658, p, v) for p, v in left] + [(1965, p, v) for p, v in right]
    pts = [
        (691.75, 1540.5, ""), (699.0, 1498.0, ""), (718.25, 1450.0, ""),
        (736.0, 1123.75, "lower of two overlapping triangles"),
        (749.75, 1094.75, "upper of two overlapping triangles"),
        (762.25, 969.75, ""), (774.25, 873.5, ""), (841.5, 739.0, ""),
        (894.1, 689.0, "lower of two overlapping triangles"),
        (906.5, 648.75, "upper of two overlapping triangles"),
        (954.0, 631.25, ""), (1040.75, 666.25, "off the calculated curve"),
        (1167.5, 576.0, ""), (1202.4, 603.0, ""), (1416.0, 547.25, ""),
        (1441.0, 532.25, "upper of two triangles at the same x"),
        (1441.6, 563.75, "lower of two triangles at the same x"),
        (1574.0, 530.75, ""), (1708.75, 500.0, ""),
    ]
    tieline("gainsb_tieline_400C_fig14_DIGITIZED.csv", 14, 165, 400, REF32,
            "x_Ga = atom fraction Ga", xt, yt, pts,
            ["Liao gives the axis as 'atom fraction Ga in the liquid' (overall atom fraction, not Ga/(Ga+In))."])


def fig15():
    xt = ([(p, 1785, v) for p, v in [(628.5, 0), (764, .1), (1036, .3), (1306, .5), (1441.5, .6), (1579.5, .7),
                                     (1715, .8), (1850, .9), (1985, 1.0)]]
          + [(p, 438, v) for p, v in [(626, 0), (763.5, .1), (900, .2), (1036, .3), (1170, .4), (1305, .5),
                                      (1440, .6), (1577, .7), (1845.5, .9), (1980.5, 1.0)]])
    left = [(421, 1), (558.5, .9), (699, .8), (838, .7), (979.5, .6), (1118, .5), (1256.5, .4), (1393.5, .3),
            (1530.5, .2), (1670, .1), (1810, 0)]
    right = [(411, 1), (546.5, .9), (686, .8), (826.5, .7), (966, .6), (1106, .5), (1245, .4), (1381, .3),
             (1517.5, .2), (1655, .1), (1796, 0)]
    yt = [(642, p, v) for p, v in left] + [(1967, p, v) for p, v in right]
    pts = [(1061.5, 596.5, ""), (1210.5, 555.75, ""), (1278.0, 527.5, ""), (1372.0, 511.75, ""),
           (1399.0, 499.25, ""), (1468.0, 477.0, ""), (1521.0, 470.75, "")]
    tieline("gainsb_tieline_500C_fig15_DIGITIZED.csv", 15, 166, 500, REF31,
            "x_Ga/(x_Ga + x_In)", xt, yt, pts, [])


# ---------------------------------------------------------------- Ansara Table 4
def ansara():
    import importlib.util
    spec = importlib.util.spec_from_file_location("t4", HERE / "_ansara_table4.py")
    t4 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(t4)
    rows = []
    for key, data in t4.T4.items():
        r = t4.RATIO[key]
        for x_in, neg_exp, neg_eq3 in data:
            x_ga, _, x_sb = t4.comp(r, x_in)
            h_exp = -neg_exp
            mine = t4.eq3(x_ga, x_in, x_sb)
            rows.append(("Ansara, Gambino & Bros 1976 (Liao ref. 26)", 995, 721.85, key.replace("/", ":"),
                         f"{x_ga:.4f}", f"{x_in:.3f}", f"{x_sb:.4f}", h_exp, f"{h_exp * 4.184:.0f}",
                         -neg_eq3, f"{mine:.1f}"))
    hdr = [
        "TABULATED. Ansara, Gambino & Bros, J. Crystal Growth 32 (1976) 101-110, Tableau 4 (p. 106).",
        "Integral enthalpy of mixing of ternary Ga-In-Sb liquids at 995 K, relative to the pure LIQUID elements,",
        "  per mole of atoms (cal/mol as printed; the table prints -dH_M, signs flipped here so dH < 0 is exothermic).",
        "Method: Tian-Calvet high-temperature microcalorimeter, indirect drop; stated precision 6 %.",
        "Compositions: the paper gives x_In along four sections of fixed x_Ga/x_Sb = 3/1, 1/1, 1/2, 1/3.",
        "  x_Ga = (1-x_In)*r/(1+r), x_Sb = (1-x_In)/(1+r). x_In printed to 3 decimals (last digit subscripted).",
        "Transcription check: dH_eq3_recomputed is Ansara's Eq. (3) (Toop with In asymmetric, binary polynomials",
        "  and alpha=5552.64, beta=-5752.12, gamma=138.69 cal/mol from the paper). It reproduces the printed",
        "  'Eq.(3)' column within 1.5 cal/mol for 90 of 92 rows; the two exceptions (1:1 x_In=0.346, printed 302 vs",
        "  299.3; 1:3 x_In=0.301, printed 444 vs 440.6) were re-read from the scan and are as printed.",
        "The printed table has 92 rows; Liao (p. 162) says 'ninety' values were fitted. The 1:2 section ends with",
        "  two identical rows (x_In=0.613, 384, 416), probably a duplicated line; both are kept here, flag it.",
        "Liao plots these as Fig. 11 'at 722 degC' (995 K = 721.85 degC).",
        "Columns: source; T_K; T_C; ratio_Ga_Sb; x_Ga; x_In; x_Sb; dH_mix_exp_cal (cal/mol atoms);",
        "  dH_mix_exp_J (x4.184, J/mol atoms); dH_eq3_printed_cal; dH_eq3_recomputed_cal.",
    ]
    cols = ["source", "T_K", "T_C", "ratio_Ga_Sb", "x_Ga", "x_In", "x_Sb", "dH_mix_exp_cal", "dH_mix_exp_J",
            "dH_eq3_printed_cal", "dH_eq3_recomputed_cal"]
    write_csv("gainsb_ternary_liquid_dHmix_995K_ansara1976_TABULATED.csv", hdr, cols, rows)


if __name__ == "__main__":
    # Page scans are embedded at 300 dpi, so rendering at 300 dpi reproduces the scan pixel grid
    # the hand-placed coordinates above refer to.
    doc = pymupdf.open(LIAO_PDF)
    assert doc[22].rect.width == 540 and doc.page_count == 29
    fig12()
    fig13()
    fig14()
    fig15()
    ansara()
