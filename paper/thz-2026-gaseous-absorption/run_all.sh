#!/bin/sh
# Reproduce every number, table and figure of the T-TST paper (paper 4).
# Needs: ns-3-dev built with contrib/thz-ntn, the Python in .venv (numpy,
# scipy, matplotlib, openpyxl, itur==0.4.0, pyrtlib==1.2.0; see
# requirements.txt) and am 14.0 (set AM_BIN).
# SKIP_CPP=1 reuses out/cpp_*.csv. SKIP_AM=1 reuses the am 14.0 outputs held in
# am_check/am_compare.json, am_slant.json, meas_check.json and
# slant_cfg/meas2_am.json, and recomputes everything else in those stages.
# Every macro file is deleted first and rewritten from scratch by this run, so
# no value from an earlier run can survive; check_macros.py then fails if the
# manuscript (main_v6.tex, sections_v6/) uses a macro that this run did not
# produce or if a macro file predates the run. PAR caps the parallel stages.
set -e
cd "$(dirname "$0")"
PY=${PY:-.venv/bin/python}
AM_BIN=${AM_BIN:-am}
PAR=${PAR:-2}
OMP_NUM_THREADS=1; OPENBLAS_NUM_THREADS=1; MKL_NUM_THREADS=1; AM_WORKERS=$PAR
export AM_BIN OMP_NUM_THREADS OPENBLAS_NUM_THREADS MKL_NUM_THREADS AM_WORKERS
RUN_START=$(date +%s); export RUN_START
[ "${SKIP_CPP:-0}" = 1 ] || ./cpp_sweeps.sh
$PY test_annex2_rev13.py
rm -f ../macros.tex ../am_macros.tex ../am_slant_macros.tex ../meas_macros.tex ../meas2_macros.tex \
      ../heff_macros.tex ../oosci_macros.tex ../budget_macros.tex ../sonde_rate_macros.tex ../p835v7_macros.tex ../arctic_macros.tex \
      out/fig_extra.json out/heff_model.json out/sonde.json out/sonde_rate.npz
$PY horiz.py
# caches first (atm defines the realizable atmospheres the later stages use),
# then the model grids, then one invocation of every summary, which writes
# macros.tex from scratch
$PY analysis.py --no-macros ref verify short atm lines prof avail sonde
$PY run_models2.py out/models.npz
$PY analysis.py verify sshort satm slines sprof models savail ssonde
$PY heff_model.py
$PY oos_block_bootstrap.py
$PY arctic.py
$PY p835v7_check.py
$PY make_figs.py
$PY make_fig1.py
$PY am_check/am_compare.py
$PY am_check/am_slant.py
$PY am_check/meas_check.py
MODELS="P676 R98 R17 R22 R24 MWL24 P676np P676loc"
[ "${SKIP_AM:-0}" = 1 ] || MODELS="$MODELS am"
for m in $MODELS; do rm -f am_check/slant_cfg/meas2_$m.json; done
echo $MODELS | tr ' ' '\n' | xargs -P "$PAR" -I{} $PY am_check/meas_check2.py --model {}
$PY am_check/meas_check2.py
$PY budget.py
$PY sonde_rate.py
$PY check_macros.py
echo "RUN_ALL_DONE $(date '+%F %T')"
