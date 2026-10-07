"""Data, priors and closed-form excess properties shared by the Ga-In-Sn Bayesian runs."""
import numpy as np

R = 8.3144626

# ---- data ----
KAT_XGA = np.array([0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90])
KAT_T = np.array([1000.0, 1073.0, 1200.0])
A_GA = np.array([[0.140, 0.132, 0.121], [0.275, 0.258, 0.236], [0.400, 0.374, 0.340],
                 [0.514, 0.480, 0.434], [0.634, 0.585, 0.522], [0.701, 0.662, 0.611],
                 [0.765, 0.740, 0.706], [0.843, 0.824, 0.799], [0.908, 0.898, 0.884]])
A_SN = np.array([[0.901, 0.900, 0.900], [0.804, 0.804, 0.804], [0.711, 0.711, 0.713],
                 [0.622, 0.623, 0.625], [0.558, 0.530, 0.538], [0.464, 0.455, 0.443],
                 [0.394, 0.370, 0.338], [0.294, 0.267, 0.232], [0.187, 0.157, 0.124]])
ZIV_XSN = np.array([0.050, 0.0742, 0.200, 0.400, 0.700])
ZIV_HMIX = np.array([95.0, 306.0, 673.0, 907.0, 723.0])

SIG_A, SIG_H, SIG_TL, SIG_EUT = 0.02, 60.0, 2.5, 2.0


def gasn_activity(xGa, T, a):
    """a_Ga, a_Sn in binary liquid Ga-Sn from excess-only closed form. a=(L0,L1,L2)."""
    xSn = 1 - xGa
    d = xGa - xSn
    L0, L1, L2 = a
    Gx = xGa * xSn * (L0 + d * L1 + d * d * L2)
    dG = (1 - 2 * xGa) * (L0 + d * L1 + d * d * L2) + xGa * xSn * (2 * L1 + 4 * d * L2)
    muGa_x = Gx + xSn * dG
    muSn_x = Gx - xGa * dG
    return xGa * np.exp(muGa_x / (R * T)), xSn * np.exp(muSn_x / (R * T))


def hmix(xSn, a):
    xGa = 1 - xSn; d = xGa - xSn
    return xGa * xSn * (a[0] + d * a[1] + d * d * a[2])


def log_prior(theta):
    a0, b0, a1, b1, a2, b2, Lt = theta
    if not (-5000 < a0 < 15000 and -20000 < a1 < 20000 and -20000 < a2 < 20000):
        return -np.inf
    if not (abs(b0) < 15 and abs(b1) < 15 and abs(b2) < 15 and -25000 < Lt < 15000):
        return -np.inf
    # weakly-informative: enthalpy near Zivkovic scale, entropy small
    lp = -0.5 * ((b0 / 4) ** 2 + (b1 / 4) ** 2 + (b2 / 4) ** 2)
    return lp


# ---- Ga-In-Sn liquidus (Evans and Prince 1978) ----
# Fig. 2 thermal arrests on the Ga-17 wt% Sn to In section, digitized: (In wt%, T degC); at In = w wt%,
# Ga = 0.83(100-w) and Sn = 0.17(100-w) wt%. Reading uncertainty about +/-2.5 C, +/-1.5 wt% In.
_M = {"GA": 69.723, "IN": 114.818, "SN": 118.710}
FIG2 = [
    (0.5, 34.5), (5.0, 27.5), (9.0, 24.0), (13.0, 19.5), (16.0, 15.5),   # liq+(Sn) branch
    (30.0, 26.0), (40.0, 38.0), (49.0, 51.0), (60.0, 64.0), (68.0, 76.0),  # liq+(In) branch
    (77.0, 84.0), (85.0, 92.0), (90.0, 103.0), (93.0, 114.0),
]
EUT_WT, EUT_T = (66.0, 20.5, 13.5), 10.7     # chemically analysed ternary eutectic, Ga/In/Sn wt%
# liquidus points of Yatsenko et al. and of van Ingen, as quoted by Evans and Prince (at.%)
EXTRA = [((0.731, 0.179, 0.090), 16.6), ((0.711, 0.188, 0.101), 14.6)]


def wt_to_at(ga, ing, sn):
    n = np.array([ga / _M["GA"], ing / _M["IN"], sn / _M["SN"]])
    n = n / n.sum()
    return (n[0], n[1], n[2])


TERN_LIQUIDUS = ([(wt_to_at(0.83 * (100 - w), w, 0.17 * (100 - w)), Tc + 273.15) for w, Tc in FIG2]
                 + [(wt_to_at(*EUT_WT), EUT_T + 273.15)]
                 + [(x, Tc + 273.15) for x, Tc in EXTRA])
