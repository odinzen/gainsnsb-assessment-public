# Liao et al. (1982) model for the (Ga,In)Sb zincblende solid solution

Source: P.-K. Liao, C.-H. Su, T. Tung, R. F. Brebrick, "Quantitative simultaneous fit to the liquidus
surface and thermodynamic data for the Ga-In-Sb system using an associated solution model for the
liquid", Calphad 6 (1982) 141-169, doi:10.1016/0364-5916(82)90009-8.
Read from the page scans (pp. 148-149, 162, 164), not from OCR text.

## Model as written

The solid is written (Ga_u In_1-u)_2-y Sb_y with y close to 1, and treated as a pseudobinary solution of the
components GaSb and InSb (p. 148, Eq. 56). The liquidus equations are (p. 149)

    mu1 + mu3 = RT ln(u Gamma13) + dG_f(GaSb(s))          (Eq. 64)
    mu2 + mu3 = RT ln((1-u) Gamma23) + dG_f(InSb(s))      (Eq. 65)

where u is the mole fraction of GaSb in the solid and dG_f is the standard Gibbs energy of formation of the
pure compound from its pure liquid elements (Eqs. 62-63). The solid is "quasiregular" (p. 149):

    RT ln Gamma13 = (W_S - V_S T)(1 - u)^2                 (Eq. 66)
    RT ln Gamma23 = (W_S - V_S T) u^2                      (Eq. 67)

"where W_S and V_S are adjustable constants independent of T and u." The equivalent integral excess Gibbs
energy (Gibbs-Duhem integration of Eqs. 66-67) is

    G^ex = u (1 - u) (W_S - V_S T)

Fitted values, Table V (p. 164), "W_S in cal, V_S in cal/K":

| W_S | V_S |
|---|---|
| 871.7 cal | -1 cal/K |

V_S is printed as the integer -1. W_S and V_S were fitted to the GaSb-InSb pseudobinary liquidus and solidus
and the ternary liquidus points together (p. 162). The paper reports that the best-fit pairs lie along a line,

    W_S = 1750.91 + 879.18 V_S       (Eq. 79, p. 162)

and compares it with an earlier pseudobinary-only quasiregular fit (Brebrick & Panlener 1974, Liao ref. 33),

    W_S = 1338 + 863 V_S             (Eq. 80)

Consistency check: Eq. 79 at V_S = -1 gives W_S = 1750.91 - 879.18 = 871.73 cal, which matches Table V.

## Basis (per mole of what)

The basis is **one mole of formula unit (Ga,In)Sb**, i.e. one mole of GaSb + InSb (2 g-atoms). The paper does
not state this in words. The evidence:

1. Eqs. 64-67 write the solid chemical potentials as mu_Ga + mu_Sb = mu_GaSb (Eq. 56), which is per mole of
   compound, and Gamma13, Gamma23 are activity coefficients of the components GaSb and InSb.
2. Table I (p. 151, "kcal/mol") is per mole of compound. For GaSb, dH_f(T_m.p.) = -16.30 kcal/mol, which
   equals -H^m + 2 x dH_mix(x=1/2) = -15.80 + 2 x (-0.25) kcal, using the liquid enthalpy of mixing per g-atom
   of -250 cal/g-atom quoted on p. 159. The factor 2 only works on a per-formula-unit basis.

By contrast, the liquid excess quantities are per g-atom (Eq. 13 converts from per mole of species), and the
liquid interaction coefficients (alpha, omega) are per mole of species.

## Conversion to J per mole of formula unit (Ga,In)1Sb1

Using the thermochemical calorie, 1 cal = 4.184 J (the IT calorie, 4.1868 J, would change the values by 0.07 %):

    W_S = 871.7 cal x 4.184 J/cal = 3647.2 J/mol
    V_S = -1 cal/K x 4.184 J/cal = -4.184 J/(mol K)

    W_S - V_S T = 3647.2 - (-4.184) T = 3647.2 + 4.184 T   J/mol

so for a two-sublattice description (Ga,In)1(Sb)1 with G^ex = y_Ga y_In L0,

    L0(Ga,In:Sb) = +3647.2 + 4.184*T   J/mol of formula unit

Examples: 6994 J/mol at 800 K, 7413 J/mol at 900 K. The positive T coefficient means a negative excess entropy,
S^ex = -4.184 u(1-u) J/(mol K). If the solid is instead modelled per mole of atoms, (Ga,In)0.5(Sb)0.5, halve both
numbers: L0 = 1823.6 + 2.092*T J/mol of atoms.

## What the data actually pin down

Eq. 79 says the fits cannot separate W_S from V_S. Writing W_S - V_S T = 1750.91 + V_S (879.18 - T) shows that
every pair on the line gives the same interaction energy at T = 879.18 K (606.0 degC):

    W_S - V_S T at 879.18 K = 1750.91 cal = 7325.8 J/mol of formula unit, independent of V_S

That temperature sits inside the pseudobinary melting range (525-709 degC), as expected. The choice V_S = -1
fixes the temperature slope; the data only weakly constrain it. The older pseudobinary-only fit (Eq. 80) gives
1338 cal = 5598 J/mol at its own pivot of 863 K, so the two analyses differ by about 1.7 kJ/mol in the pinned
value. This matters when choosing priors for the reassessment.

## Liao's liquid model, for reference

Associated solution with species Ga(1), In(2), Sb(3), GaSb(4), InSb(5). Excess Gibbs energy per mole of species
(Eq. 19, p. 144) with quadratic coefficients alpha_ij and only two cubic terms kept (beta_13, beta_23, p. 145).
Binary quadratic coefficients are re-parameterized (Eqs. 28-30, p. 146; In-Sb analogue implied):

    alpha13 = 2 Omega_S,  alpha14 = (Omega_S - Omega_R + Omega_A)/2,  alpha34 = (Omega_S - Omega_R - Omega_A)/2
    Omega_i = omega_i - nu_i T  (i = S, R; Eq. 48),  Omega_A constant (p. 151)
    beta13 constant (p. 151)

Associate dissociation: dG_i = dH_i - T dS_i (Eq. 9), dissociation of GaSb (i=4) and InSb (i=5) to pure liquids.

In-Sb (Table II, p. 152): z* = 0.6092; omega_R = 1902.99 cal; nu_R = 1.5858 cal/K; omega_S = 1720.62 cal;
nu_S = 1.1627 cal/K; Omega_A = -460.68 cal; beta23 = 1321.78 cal; dH5 = 2817 cal; dS5 = -0.45723 cal/K.

Ga-Sb (Table IV, p. 158, entry 3, "the set chosen for further calculations"): z* = 0.2233; omega_R = -200.7;
nu_R = 4.185; Omega_A = -302.1 (units not printed in Table IV; by analogy with Table II, cal and cal/K).
dH4 = 1,668 cal/mol and dS4 = -1.464 cal/(K mol) (p. 158). omega_S, nu_S and beta13 for Ga-Sb are fixed by the
auxiliary constraints (Eqs. 68, 70, 71, 74) and are not printed.

Ga-In (pp. 160-161): dH_M = 2 omega12 x_In x_Ga with omega12 = 530 cal/g-atom; dS_M^xs = 2 nu12 x_In x_Ga with
nu12 = -0.5746 cal/(g-atom K).

Ternary (Table V, p. 164, cal): alpha15 (Ga-InSb) = 525; alpha24 (In-GaSb) = -93.9; alpha45 (GaSb-InSb) = 1176.8;
all temperature independent (p. 151).

Compound data (Table I, p. 151; kcal/mol, cal/(K mol), per mole of compound):

| | T_m.p. (degC) | H^m | dH_f(T_m.p.) | dS_f(T_m.p.) | dH_f (eff.) | dS_f (eff.) | dC_p(1/2, T_m.p.) |
|---|---|---|---|---|---|---|---|
| InSb | 525 | 11.41 | -13.27 | -11.53 | -13.02 | -11.22 | 1.00 |
| GaSb | 709.2 | 15.80 | -16.30 | -13.10 | -16.23 | -13.03 | 0.575 |

Formation quantities are relative to the pure liquid elements. Sb terminal solid for the eutectic valley used
T_m = 631 degC and H^m = 4760 cal/g-atom (p. 163).
