# Liao et al. (1982): reported fit quality

Source: Liao, Su, Tung & Brebrick, Calphad 6 (1982) 141-169. Values read from the page scans.

## Metric definitions

- Phase-boundary temperatures (Eq. 76, p. 151), RMS deviation in degC:
  sigma_T^2 = sum_j (T_j,cal - T_j,obs)^2 / M
- Thermochemical quantities Z (activity, enthalpy of mixing, partial molar enthalpy) (Eq. 77, p. 153), RMS of the
  fractional deviation: sigma(Z)^2 = sum_j [(Z_j,obs - Z_j,cal)/Z_j,obs]^2 / M, reported in %.
  Exception: Table IV gives sigma(dH_M) for Ga-Sb as an absolute value in cal/g-atom.

## Ternary and pseudobinary (Table V, p. 164; text p. 162)

| Dataset | Metric | Value | Points | Fitted to |
|---|---|---|---|---|
| GaSb-InSb pseudobinary liquidus | sigma_L (RMS, degC) | 5.88 | not stated (11 symbols in Fig. 12) | W_S, V_S |
| GaSb-InSb pseudobinary solidus | sigma_S (RMS, degC) | 2.7 | not stated (15 symbols in Fig. 12) | W_S, V_S |
| Ternary liquidus | sigma_L(TER) (RMS, degC) | 11.6 | not stated (18 clear symbols in Fig. 13, 2 more unresolved) | W_S, V_S |
| Ternary liquid enthalpy of mixing, 995 K | sigma(dH_M), fractional RMS | 15 % | 90 (Ansara et al. 1976 Table 4 prints 92 rows) | alpha15, alpha24, alpha45 |
| Tie-lines, 400 degC (Miki 1975) and 500 degC (Antypas 1972) | none | "excellent" agreement (p. 163) | 19 and 7 symbols (Figs. 14, 15) | not fitted, compared only |

The three ternary alpha parameters were fitted to the enthalpies alone; W_S and V_S were then fitted to the
pseudobinary liquidus, solidus and ternary liquidus together (p. 162). The point counts in parentheses are what
I digitized from the figures, not numbers stated by Liao.

## In-Sb binary (Tables II, III; pp. 151-156)

| Dataset | Metric | Value | Points (source) |
|---|---|---|---|
| InSb liquidus | sigma_T (degC) | 8.43 | 16 (13 thermal analysis, Liu & Peretti 1952; 3 dissolution, Shunk 1969) |
| Sb liquidus | sigma_V (degC) | 1.99 | 12 (Liu & Peretti 1952) |
| Calculated InSb-Sb eutectic | T, x_Sb | 495.2 degC, 0.691 | vs. 494 +/- 0.5 degC, 0.69 (exp.) |
| a_In at 627, 700, 800 degC | sigma(a_In), % | 6.8, 10.6, 10.4 | Hoshino et al. 1965; Chatterji & Smith 1973 |
| dH_M at 627, 680, 684, 713, 778, 860, 911 degC | sigma(dH_M), % | 5.0, 4.2, 3.7, 3.2, 4.0, 4.1, 7.2 | Predel & Oehme 1976 (680 degC, tabulated); Predel (others, scaled from a graph) |
| Partial enthalpies at 627 degC | sigma(dH_In), sigma(dH_Sb), % | 15.2, 15.3 | Hultgren et al. 1973 |

Only the liquidus was minimized (sigma_T); the thermochemical sigmas are a posteriori comparisons.

## Ga-Sb binary (Table IV, p. 158)

Liquidus and enthalpy of mixing fitted simultaneously. Data: 6 liquidus points scaled from a graph of
dissolution experiments (Hall 1963), 22 tabulated DTA points (Maglione & Potier 1968), 28 enthalpies of mixing
(Ansara et al. 1976; Gambino & Bros 1975).

| Entry | sigma_T (degC), GaSb liquidus | sigma_V (degC), Sb liquidus | sigma(dH_M) cal/g-atom, 727 / 750 degC | Eutectic degC, x |
|---|---|---|---|---|
| 1 (best liquidus) | 7.6 (1.6) | 0.91 | 65.6 / 103.5 | 588.6, 0.883 |
| 2 (best enthalpy) | 16.1 (7.7) | 2.4 | 9.1 / 26.9 | 578.7, 0.8742 |
| **3 (adopted)** | **9.5 (7.4)** | **0.92** | **18.8 / 45.0** | **589.2, 0.8881** |

The bracketed sigma_T excludes the 6 points scaled from the dissolution experiments (p. 157). The two
sigma(dH_M) columns are headed "727 degC" and "750" as printed, although the Ga-Sb enthalpy data are at 995 K
(721.85 degC); the headers are reported as printed.

## Ga-In binary (pp. 160-161)

Not fitted to phase-diagram data. omega12 = 530 cal/g-atom reproduces the selected enthalpy of mixing
(maximum 265 cal/g-atom at x = 0.5) "over the whole composition range"; nu12 from the entropy of mixing. No
numerical metric is given.
