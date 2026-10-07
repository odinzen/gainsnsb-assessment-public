# Vassiliev et al. 2001, In-Sb, Sb-Sn and In-Sb-Sn liquid EMF data

Source: V. Vassiliev, Y. Feutelais, M. Sghaier, B. Legendre, J. Alloys Compd. 314 (2001) 198-205,
doi:10.1016/S0925-8388(00)01243-3. The source PDF is not distributed.
Every value was read from the text layer and checked against a 220 dpi render of the page.

## Files

| File | Paper table | Page | Rows |
|---|---|---|---|
| table1_alloy_compositions.csv | 1 | 199 | 24 (2 In-Sb, 5 Sb-Sn, 17 ternary) |
| table3_InSb_emf.csv | 3 | 199 | 2 |
| table4_SbSn_emf.csv | 4 | 200 | 5 |
| table5_SbSn_dsc.csv | 5 | 201 | 18 |
| table6_ternary_emf.csv | 6 | 202 | 17 |
| table8_partial_H_S_In.csv | 8 | 203 | 16 printed rows (15 alloys plus a Hultgren end point) |
| table9_mu_In_750K.csv | 9 | 205 | 17 |

The paper prepared 24 alloys, but only 17 are ternary (nos. 8-24). Nos. 1-2 are In-Sb and 3-7 are Sb-Sn.
Table 6 prints only x_In, so x_Sb and x_Sn in the ternary CSVs come from Table 1.

## Cells and the E to mu conversion

Electrolyte is the LiCl-KCl eutectic (46.0 wt.% LiCl, T_eut = 625 K), dried and HCl-treated. Cells:

- (I) (-) In | In+ in electrolyte | In_x Sb_(1-x) (+), Table 3
- (II) (-) Sn | Sn2+ in electrolyte | Sn_x Sb_(1-x) (+), Table 4
- (III) (-) In | In+ in electrolyte | In_x Sb_y Sn_z (+), Table 6

The paper states mu_In = -F E = RT ln a_In and mu_Sn = -2 F E = RT ln a_Sn. So n = 1 for In (In+) and
n = 2 for Sn (Sn2+). The reference electrode is pure In (cells I, III) or pure Sn (cell II). The paper
says "pure indium" without the word liquid, but every cell temperature (624-835 K) is far above the
melting points of In (429.75 K) and Sn (505.08 K), so the reference state is pure liquid In or Sn.

For a ternary alloy at temperature T (K), with a in mV and b in mV/K:

    mu_In - G_In(liq) = -F (a + b T) * 1e-3      [J/mol],  F = 96485 C/mol
    partial H_In = -F a * 1e-3                   [J/mol]
    partial S_In = +F b * 1e-3                   [J/(mol K)]

Check, alloy 8 at 750 K. E = 28.39 + 0.10171*750 = 104.67 mV gives -10.10 kJ/mol, which is the printed
Table 9 measured value. Table 8 reproduces -F a and F b in the same way. The conversion was verified this
way for every row (see "Internal consistency" below).

Watch the b scaling. Tables 3 and 6 print b x 10^2, Table 4 prints b x 10^3. The CSVs carry both the
printed value and b in mV/K.

The a, b fits are only valid over each alloy's own [T_min, T_max] range, given per row.

## Binary liquid parameters used for the Muggianu calculation (Table 9 "calc")

Redlich-Kister, L_ij = sum_v vL_ij (x_i - x_j)^v, J/mol, element order as written.

- In-Sb, from Ansara et al., Calphad 18 (1994) 177 (ref [10])
  - 0L = -25631.2 + 102.9324 T - 13.45816 T lnT
  - 1L = -2115.4 - 1.31907 T
  - 2L = +2908.9
- In-Sn, the authors' own re-optimization from EMF, Vassiliev et al., Thermochim. Acta 315 (1998) 129
  (ref [3])
  - 0L = -783.19 - 0.59353 T
  - 1L = -149.39
- Sb-Sn, from Jonsson and Agren, Mater. Sci. Technol. 2 (1986) 913 (ref [13])
  - 0L = -5695.14 - 1.7090 T
  - 1L = +782.595
  - 2L = +1840.91

No ternary parameter was used. The authors conclude it is not needed (Sn behaves as a neutral solvent).
Unary data are SGTE (Dinsdale 1991). The calculations were done in Thermo-Calc.

## Ambiguities and flagged values (not corrected, kept as printed)

1. **Table 9 "calc" is not reproduced exactly by the printed parameters.** I recomputed mu_In at 750 K
   from the parameters above (Muggianu, SGTE ideal term, compositions from Table 1). It comes out
   0.05-0.30 kJ/mol less negative than the printed "calc" column for every alloy (mean 0.16 kJ/mol).
   The gap grows roughly with x_Sb*x_Sn, which points at the Sb-Sn description. Flipping the sign
   convention of the odd L terms does not fix it, and 770 K overshoots the other way. The cause is
   unknown. The parameters may be rounded in print, or the calculation may have used a slightly
   different description. Compare against "meas", not "calc", when validating a new model.
2. **Table 9, alloy 19, measured value.** Printed -16.77 kJ/mol. The Table 6 coefficients give -16.62.
   -16.77 is also the printed value for alloy 12 in the row above, so this may be a copy error in the
   table. Every other non-bracketed alloy agrees with -F E(750) to within 0.015 kJ/mol.
3. **Table 8, alloy 18, partial entropy.** Printed dS_meas = 4.43 J/(mol K), but F b from Table 6
   (b = 5.622e-2 mV/K) gives 5.42. Table 9 (meas -6.52) is consistent with b = 5.622, not with 4.43.
   Fig. 9 plots the point near 4.4, so the error, if it is one, is in Table 8 and Fig. 9 together, or
   in Table 6 b together with Table 9. Not resolvable from the paper.
4. **Table 6, alloy 18, s0^2 x 10^2.** Printed "37" with no decimal, unlike every other entry. Kept as
   printed. Alloy 18 also has T_min = 624 K, 1 K below the stated electrolyte eutectic (625 K).
5. **Table 6, Ebar for alloys 14 and 20.** For a least-squares line the mean EMF equals a + b*Tbar
   exactly. Both rows miss by about 0.45 mV (alloy 14: 228.04 vs 228.5; alloy 20: 198.35 vs 197.90).
   Small, but a or b may be off in the last digit.
6. **Table 8, alloy 20.** dH_meas printed -5685 but -F a gives -5697. dS_meas printed 17.50 but F b gives
   17.54. Small.
7. **Table 8 layout.** The table is two independent blocks printed side by side (left: alloy, x_Sb,
   dH sorted by x_Sb; right: alloy, x_In, dS sorted by x_In). Within a printed row the two alloy numbers
   differ. Join on alloy number, not on row. The last row is "[15]" (Hultgren 1973), x_Sb = 1,
   dH = -9194 J/mol (In infinitely dilute in Sb). Its right half prints x_In = 0, dS = 0, 0, which is
   meaningless because dS_In diverges as x_In -> 0. Do not use it.
8. **Table 8 "fit" column.** The fit values were evaluated at nominal compositions. Alloys 10, 11 and 12
   all show -5035 (x_Sb = 0.5) and 18 and 21 both show -1907 (x_Sb = 0.25). Eq. (7) at the actual x_Sb
   differs by up to 12 J/mol.
9. **Table 1, alloy 11.** Printed composition 0.3332 / 0.4991 / 0.1669 sums to 0.9992. Every other
   alloy sums to 1 within 2e-4. One of the three is probably mistyped. Kept as printed.
10. **Alloys 13-16 (x_In 0.05-0.11)** suffer side reactions with the pure-In and In-rich (no. 18)
    electrodes in the same cell (Section 3.1, Figs. 6-7). The authors excluded them from Eqs. (7)-(8)
    and bracket the measured mu_In of 14 and 15 in Table 9. They say only the first points of each E(T)
    run were used for these alloys. Treat 13-16 as low weight and 14-15 as unreliable. Alloy 18's EMF
    drifts upward by the same mechanism (Fig. 7, about 4 mV over 80 days at 754 K).
11. **Table 5 reaction assignment.** The paper does not assign a reaction to any of T1-T4, neither in
    the table nor in the text. Fig. 3 shows them only as "invariant" points on the Jonsson-Agren
    calculated diagram, at about 242-246 degC, about 323-326 degC and about 424-426 degC. Whether the
    T1 (about 242 degC) and T2 (about 245 degC) arrests are two separate reactions is not stated.
    The `reaction_assigned` column is "not stated". Take assignments from the Sb-Sn solid-state paper
    (ref [1], Vassiliev, Lelaurain, Hertz, J. Alloys Compd. 247 (1997) 223) rather than inferring them
    here. T_liq values are printed as integers (degC).
12. **Text typo.** Section 2.2.4 says "our measured mu_Sb data" for Sb-Sn. Cell II measures Sn, and the
    Fig. 2 caption says mu_Sn. Table 4 is mu_Sn data.
13. **Table 7** ("minimum of the Gibbs function of mixing at 900 K": Sn-Sb -2340, In-Sn -4310, In-Sb
    -9070 J/mol) is not consistent with the parameters above. The ideal term alone at x = 0.5 is
    -5187 J/mol, so the Sn-Sb and In-Sn entries cannot be total Delta G_mix. In-Sb (-9070) does match the
    total. Not extracted, mentioned only as a caution.
14. The +/- on a and b in Table 6 are not defined as 1 sigma or otherwise. Tables 3 and 4 print no
    uncertainty on a or b.

## Katayama et al. 1996

`Katayama_1996_MaterTransJIM_GaSn_EMF.pdf` (Mater. Trans. JIM 37 (1996) 988-990) is binary Ga-Sn only.
It uses a zirconia cell at 1000-1200 K. There are no Ga-Sn-Sb or other Sb data. Its reference list cites
Ga-Sb, Ga-Sb-Te, Ga-Sb-In and Ga-Sb-Bi work by the same group, but not Ga-Sb-Sn.
