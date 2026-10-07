# GaInSnSb-CALPHAD

Reproducibility package for:

**Primary crystallization of III-V (Ga,In)Sb from low-melting Ga-In-Sn melts: a thermodynamic assessment of the Ga-In-Sn-Sb system**
Michael E. Bustamante, Kristina Lilova (submitted to the *Journal of Alloys and Compounds*)

Everything quoted in the paper is calculated from `databases/GaInSnSb.tdb` by the scripts here; the
numbers land in `results/*.json` and the figures in `figures/`.

## Contents

| Path | Description |
|------|-------------|
| `databases/GaInSnSb.tdb` | The Ga-In-Sn-Sb database: liquid, terminal solid solutions with their solubilities, (Ga,In)Sb zincblende, SbSn and Sb3Sn4, the In-Sn beta and gamma phases, and the gas. Sources are listed in the file header. |
| `databases/LICENSE.md` | The database license (CC BY-NC 4.0) |
| `bayesian/emcee_full.py` | Bayesian joint optimization of the Ga-Sn liquid and the Ga-In-Sn ternary term (emcee) |
| `bayesian/forward_full.py` | Forward model; every phase is built from the database itself |
| `bayesian/assess_data.py` | The data, uncertainties and priors of the fit |
| `bayesian/finalize_full.py` | Convergence diagnostics, posterior medians written into the database, posterior-predictive eutectics, the three-way transferability test (Table 2) and the eutectic diagnostic |
| `bayesian/emcee_full_chain.npy` | Posterior samples after burn-in (36,000 x 7: a0, b0, a1, b1, a2, b2, Lt) |
| `bayesian/emcee_full_chain_raw.npy`, `emcee_full_lnp_raw.npy` | Full walker chains and log-probabilities (24 walkers x 3000 steps) |
| `bayesian/posterior_median.npy` | Posterior medians, the values in the database |
| `bayesian/emcee_chain.npy` | Chain of the earlier pure-solid run, used only as the sampler's starting point |
| `bayesian/fig3_points.npy` | Evans and Prince liquidus points plotted in Fig. 7a |
| `scripts/headline_numbers.py` | Invariants, melting points, isopleth, antimony solubility, Sn-corner liquidus and vapor checks |
| `scripts/fit_gainsb_liao.py`, `gainsb_robustness.py`, `validate_gainsb_liao.py` | Ga-In-Sb ternary term: fit, sensitivity, and comparison with the data in Liao's metrics (Table 3) |
| `scripts/validate_insnsb_vassiliev.py` | In-Sn-Sb against Vassiliev et al. (Fig. 5) |
| `scripts/validate_gasnsb.py` | Ga-Sn-Sb against Katayama et al. (1998) and Gerdes and Predel (1981), with the trial ternary term (Fig. 6) |
| `scripts/fast_eq.py` | Fast liquidus and enthalpy evaluation from the database, used by the fits |
| `scripts/fig_*.py`, `make_figs_v2*.py`, `*projection*.py`, `*sections*.py`, `make_montages.py` | Figures |
| `data/` | Digitized and transcribed datasets with their provenance (Liao et al. 1982, Vassiliev et al. 2001, Katayama et al. 1998, Gerdes and Predel 1981). The `build_*.py` files record the digitizing; the source PDFs are not distributed. |
| `results/` | Every calculated number quoted in the paper |
| `tests/test_fit_matches_tdb.py` | Checks that the fitting code and the database give the same Gibbs energy for every phase and that the database carries the posterior medians |
| `figures/Figure1.png` ... `Figure11.png`, `FigureS1.png` ... `FigureS3.png` | Figures as submitted |
| `references.json` | The bibliography (CSL-JSON from Crossref) |

## Requirements

Python 3.10+, pycalphad 0.11, numpy, scipy, matplotlib, symengine, emcee 2.2, corner, pytest
(conda-forge).

```
conda create -n gainsnsb -c conda-forge python=3.12 pycalphad=0.11 numpy scipy matplotlib symengine emcee=2.2 corner pytest
conda activate gainsnsb
bash reproduce.sh
```

`reproduce.sh` recalculates all results and figures from the database and runs the tests. The
projections and sections compute equilibrium grids and take about an hour; the scripts write
working file names in `figures/`, of which `FigureN.png` are the submitted copies. The Markov chain
itself is rerun with `python bayesian/emcee_full.py 3000` (several hours on six cores), and
`bash reproduce.sh --posterior` repeats the post-processing that writes the medians into the
database.

## Notes

- Calculations use pycalphad. Liquidus temperatures are located by bisection to 0.05 C, invariants
  by stepping in 0.02 C.
- The In-Sn gamma endmember of David et al. is modified above 600 K only, so that pure indium does
  not become solid again above about 844 K; below 600 K it is as published.

## License

- **Database** (`databases/GaInSnSb.tdb`): CC BY-NC 4.0 (`databases/LICENSE.md`). Free to use, share and adapt
  for non-commercial purposes with credit to the paper and authors. Any commercial use, including use in paid
  products or services, requires a license from Odinzen LLC.
- **Code** (scripts, the fitting and analysis code, tests): MIT (`LICENSE`).
- **Figures and data** (`figures/`, `data/`, `results/`, the Markov chain): CC BY 4.0.
