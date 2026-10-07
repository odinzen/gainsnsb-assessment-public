"""Forward model for the Ga-In-Sn Bayesian assessment, built from the shipped database.

Every phase energy comes from pycalphad's Model of databases/GaInSnSb.tdb itself, with the five
sampled liquid parameters (Ga-Sn a+bT per order, Ga-In-Sn ternary ;0) swapped for free symbols.
Nothing is re-typed, so the fit and the shipped file cannot drift apart; tests/test_fit_matches_tdb.py
checks that anyway. The solids carry their full Ga, In and Sn solubility.
"""
import os
import re
import warnings

import numpy as np
import symengine as se
from pycalphad import Database, Model, variables as v

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
TDB = os.path.join(HERE, "..", "databases", "GaInSnSb.tdb")
COMPS = ["GA", "IN", "SN", "VA"]
SOLIDS = ["ORTHORHOMBIC_GA", "TETRAGONAL_A6", "BCT_A5", "BETA", "GAMMA"]
PNAMES = ["VV0000", "VV0001", "VV0002", "VV0003", "VV0004", "VV0005", "VV0006"]  # a0 b0 a1 b1 a2 b2 Lt
R = 8.31451


def _symbolic_db():
    txt = open(TDB).read()
    subs = {
        r"G\(LIQUID,GA,SN;0\)": "VV0000+VV0001*T",
        r"G\(LIQUID,GA,SN;1\)": "VV0002+VV0003*T",
        r"G\(LIQUID,GA,SN;2\)": "VV0004+VV0005*T",
        r"G\(LIQUID,GA,IN,SN;0\)": "VV0006",
    }
    for key, expr in subs.items():
        pat = rf"(PARAMETER {key} 298\.15 )[^;]*;"
        txt, n = re.subn(pat, rf"\g<1>{expr};", txt)
        assert n == 1, key
    funcs = "".join(f"FUNCTION {p} 298.15 0; 6000 N !\n" for p in PNAMES)
    txt = txt.replace("TYPE_DEFINITION", funcs + "TYPE_DEFINITION", 1)
    return Database(txt)


DBS = _symbolic_db()
_T = v.T


class _Phase:
    """GM(T, y, theta) and, for the liquid, dGM/dy, as vectorized callables."""

    def __init__(self, name):
        m = Model(DBS, COMPS, name, parameters=PNAMES)
        self.name = name
        self.y = [s for s in m.site_fractions]
        self.consts = [str(s.species.name) for s in self.y]
        gm = m.GM.xreplace({v.P: 101325})
        args = [_T] + self.y + [se.Symbol(p) for p in PNAMES]
        self._g = se.Lambdify(args, [gm], backend="llvm")
        self._dg = se.Lambdify(args, [gm.diff(s) for s in self.y], backend="llvm") if name == "LIQUID" else None

    def _args(self, T, Y, theta):
        Y = np.atleast_2d(Y)
        n = Y.shape[0]
        cols = [np.full(n, float(T))] + [Y[:, i] for i in range(Y.shape[1])] + [np.full(n, t) for t in theta]
        return np.column_stack(cols)

    def G(self, T, Y, theta):
        return np.asarray(self._g(self._args(T, Y, theta))).reshape(-1)

    def dG(self, T, Y, theta):
        return np.asarray(self._dg(self._args(T, Y, theta))).reshape(np.atleast_2d(Y).shape[0], -1)


LIQ = _Phase("LIQUID")
SOL = {p: _Phase(p) for p in SOLIDS}


def _grid(consts):
    """Sampling grid in each solid's own constituent order. Fine near the In and Sn corners
    (where the Ga solubility lives), coarse elsewhere."""
    k = len(consts)
    if k == 1:
        return np.ones((1, 1))
    if k == 2:
        return np.column_stack([np.linspace(1e-6, 1 - 1e-6, 801), 1 - np.linspace(1e-6, 1 - 1e-6, 801)])
    pts = []
    for a in np.linspace(0, 1, 51):
        for b in np.linspace(0, 1 - a, max(int(round((1 - a) * 50)) + 1, 1)):
            pts.append((a, b, 1 - a - b))
    iga = consts.index("GA")
    others = [i for i in range(3) if i != iga]
    for xg in np.linspace(0, 0.12, 25):
        for f in np.linspace(0, 1, 201):
            p = [0.0, 0.0, 0.0]
            p[iga] = xg
            p[others[0]] = (1 - xg) * f
            p[others[1]] = (1 - xg) * (1 - f)
            pts.append(tuple(p))
    g = np.clip(np.array(pts), 1e-9, 1)
    return g / g.sum(axis=1, keepdims=True)


GRIDS = {p: _grid(SOL[p].consts) for p in SOLIDS}
_ORDER = ["GA", "IN", "SN"]


def mu_liquid(x, T, theta):
    """Liquid chemical potentials (Ga, In, Sn) at overall x = (xGa, xIn, xSn)."""
    y = np.array([[x[_ORDER.index(c)] for c in LIQ.consts]])
    g = LIQ.G(T, y, theta)[0]
    d = LIQ.dG(T, y, theta)[0]
    mu = g + d - np.dot(y[0], d)
    return {c: mu[i] for i, c in enumerate(LIQ.consts)}


def max_driving_force(x, T, theta):
    mu = mu_liquid(x, T, theta)
    best = -np.inf
    for p in SOLIDS:
        ph = SOL[p]
        Y = GRIDS[p]
        mus = np.array([mu[c] for c in ph.consts])
        df = Y @ mus - ph.G(T, Y, theta)
        best = max(best, float(df.max()))
    return best


def liquidus(x, theta, lo, hi):
    from scipy.optimize import brentq
    x = np.clip(np.asarray(x, float), 1e-9, 1)
    f = lambda T: max_driving_force(x, T, theta)
    flo, fhi = f(lo), f(hi)
    if flo < 0 or fhi > 0:
        return np.nan
    return brentq(f, lo, hi, xtol=0.02)


def gm(phase, T, x_by_const, theta):
    """Molar G of a phase at one composition, for the consistency test."""
    ph = LIQ if phase == "LIQUID" else SOL[phase]
    y = np.array([[x_by_const[c] for c in ph.consts]])
    return float(ph.G(T, y, theta)[0])
