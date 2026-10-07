"""The Bayesian forward model and the shipped database must describe the same phases.

Rule: a fitting tool never carries its own copy of a phase energy. forward_full builds every phase
from databases/GaInSnSb.tdb, so these tests guard against that link being broken (for example by a
future hand-written shortcut) and check that the shipped parameters are the posterior medians.
"""
import os
import re
import sys

import numpy as np
import pytest
from pycalphad import Database, calculate

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "bayesian"))
import forward_full as F  # noqa: E402

TDB = os.path.join(ROOT, "databases", "GaInSnSb.tdb")
CHAIN = os.path.join(ROOT, "bayesian", "emcee_full_chain.npy")
NAMES = [("GA,SN;0", 0, 1), ("GA,SN;1", 2, 3), ("GA,SN;2", 4, 5)]


def shipped_theta():
    txt = open(TDB).read()
    th = []
    for key, _, _ in NAMES:
        m = re.search(rf"PARAMETER G\(LIQUID,{key}\) 298\.15 ([+-]?[\d.]+)([+-][\d.]+)\*T;", txt)
        th += [float(m.group(1)), float(m.group(2))]
    m = re.search(r"PARAMETER G\(LIQUID,GA,IN,SN;0\) 298\.15 ([+-]?[\d.]+);", txt)
    return th + [float(m.group(1))]


@pytest.mark.parametrize("phase", ["LIQUID"] + F.SOLIDS)
@pytest.mark.parametrize("T", [283.15, 350.0, 500.0, 900.0])
def test_phase_energy_matches_file(phase, T):
    db = Database(TDB)
    th = shipped_theta()
    cons = (F.LIQ if phase == "LIQUID" else F.SOL[phase]).consts
    rng = np.random.default_rng(0)
    for _ in range(6):
        y = rng.dirichlet(np.ones(len(cons))) if len(cons) > 1 else np.ones(1)
        fit = F.gm(phase, T, dict(zip(cons, y)), th)
        ship = float(calculate(db, F.COMPS, phase, T=T, P=101325, points=np.array([y])).GM.values.squeeze())
        assert abs(fit - ship) < 1.0, (phase, T, y, fit, ship)


def test_shipped_parameters_are_posterior_medians():
    if not os.path.exists(CHAIN):
        pytest.skip("posterior chain not present")
    med = np.median(np.load(CHAIN), axis=0)
    th = np.array(shipped_theta())
    # shipped values are rounded to 0.1 J/mol and 1e-4 J/mol/K
    tol = np.array([0.06, 6e-5, 0.06, 6e-5, 0.06, 6e-5, 0.06])
    assert np.all(np.abs(med - th) <= tol), (med, th)
