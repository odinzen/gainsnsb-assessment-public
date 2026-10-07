"""Fast liquidus and mixing-enthalpy evaluation straight from a TDB, for parameter fitting.

Phase energies are pycalphad Models of the given database text, compiled once, with chosen
parameters left as free symbols (FUNCTION VVxxxx placeholders). Solids are sampled on pycalphad's
own point grid. Used by fit_gainsb_liao.py; the result is always re-checked with pycalphad's
equilibrium on the written database.
"""
import warnings

import numpy as np
import symengine as se
from pycalphad import Database, Model, calculate, variables as v

warnings.filterwarnings("ignore")


class FastSystem:
    def __init__(self, tdb_text, elements, solids, pnames, pdens=3000):
        funcs = "".join(f"FUNCTION {p} 298.15 0; 6000 N !\n" for p in pnames)
        txt = tdb_text.replace("TYPE_DEFINITION", funcs + "TYPE_DEFINITION", 1)
        self.db = Database(txt)
        self.el = list(elements)
        comps = self.el + ["VA"]
        self.pn = list(pnames)
        psym = [se.Symbol(p) for p in pnames]
        m = Model(self.db, comps, "LIQUID", parameters=self.pn)
        self.ly = list(m.site_fractions)
        self.lc = [str(s.species.name) for s in self.ly]
        g = m.GM.xreplace({v.P: 101325})
        args = [v.T] + self.ly + psym
        self._lg = se.Lambdify(args, [g], backend="llvm")
        self._ldg = se.Lambdify(args, [g.diff(s) for s in self.ly], backend="llvm")
        self._lh = se.Lambdify(args, [g - v.T * g.diff(v.T)], backend="llvm")
        self.sol = []
        for ph in solids:
            ms = Model(self.db, comps, ph, parameters=self.pn)
            ys = list(ms.site_fractions)
            gs = ms.GM.xreplace({v.P: 101325})
            f = se.Lambdify([v.T] + ys + psym, [gs], backend="llvm")
            pts = calculate(self.db, comps, ph, T=500, P=101325, pdens=pdens)
            Y = pts.Y.values.reshape(-1, pts.Y.shape[-1])[:, : len(ys)]
            X = pts.X.values.reshape(-1, pts.X.shape[-1])
            order = list(pts.component.values)
            ok = ~np.isnan(Y).any(axis=1)
            X = np.column_stack([X[ok, order.index(e)] for e in self.el])
            self.sol.append((ph, f, Y[ok], X))

    @staticmethod
    def _cols(T, Y, theta):
        n = Y.shape[0]
        return np.column_stack([np.full(n, float(T))] + [Y[:, i] for i in range(Y.shape[1])]
                               + [np.full(n, float(t)) for t in theta])

    def _ly(self, x):
        return np.array([[max(x[self.el.index(c)], 1e-12) for c in self.lc]])

    def mu(self, x, T, theta):
        y = self._ly(x)
        a = self._cols(T, y, theta)
        g = np.asarray(self._lg(a)).ravel()[0]
        d = np.asarray(self._ldg(a)).ravel()
        mu = g + d - np.dot(y[0], d)
        return np.array([mu[self.lc.index(e)] for e in self.el])

    def driving_force(self, x, T, theta):
        mu = self.mu(x, T, theta)
        best, who = -np.inf, None
        for ph, f, Y, X in self.sol:
            df = X @ mu - np.asarray(f(self._cols(T, Y, theta))).ravel()
            k = int(np.argmax(df))
            if df[k] > best:
                best, who = float(df[k]), (ph, X[k])
        return best, who

    def liquidus(self, x, theta, lo, hi, xtol=0.05):
        from scipy.optimize import brentq
        f = lambda T: self.driving_force(x, T, theta)[0]
        if f(lo) < 0 or f(hi) > 0:
            return np.nan
        return brentq(f, lo, hi, xtol=xtol)

    def hmix(self, x, T, theta):
        """Integral enthalpy of mixing of the liquid relative to the pure liquid elements."""
        y = self._ly(x)
        h = np.asarray(self._lh(self._cols(T, y, theta))).ravel()[0]
        pure = 0.0
        for i, e in enumerate(self.el):
            if x[i] > 0:
                xp = np.zeros(len(self.el)); xp[i] = 1.0
                pure += x[i] * np.asarray(self._lh(self._cols(T, self._ly(xp), theta))).ravel()[0]
        return h - pure
