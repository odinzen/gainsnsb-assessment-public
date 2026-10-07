#!/usr/bin/env bash
# Recalculate every result and figure from databases/GaInSnSb.tdb. Run from the repository root.
# The Markov chain is not resampled here; bayesian/emcee_full.py reruns it (hours), and
# `--posterior` redoes the post-processing that writes the posterior medians into the database.
set -euo pipefail
cd "$(dirname "$0")"
if [[ "${1:-}" == "--posterior" ]]; then
  (cd bayesian && python finalize_full.py diag tdb ppc three eutectic)
fi
cd scripts
python headline_numbers.py > /dev/null
python validate_gainsb_liao.py
python gainsb_robustness.py > /dev/null
python validate_insnsb_vassiliev.py | tail -1
python validate_gasnsb.py > /dev/null
python make_figs_v2.py
python make_figs_v2_binaries.py 3
python fig_gainsb_pseudobinary.py
python fig_gainsb_tielines.py
python fig_insnsb_vassiliev.py | tail -1
python fig_gasnsb_validation.py
python fig_sb_separate.py
python fig_sb_isopleth.py
python fig_sb_solubility.py
python quaternary_sections.py
python plot_quaternary_sections.py
python full_projection_compute.py
python full_projection_plot.py
python gainsb_projection_compute.py
python proj_ternary.py gasnsb
python proj_ternary.py insnsb
for k in gainsb gasnsb insnsb; do python plot_proj_ternary.py "$k"; done
python make_montages.py
cd ..
python -m pytest tests -q
