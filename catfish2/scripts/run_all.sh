#!/usr/bin/env bash
# Run every simulation, then the figures. Pass "quick" for a fast smoke test.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p ../logs

if [ "${1:-}" = "quick" ]; then
  A1="50 200 50"; A2="80 80"; A3="80 80"; A4="30 200"; A8="80 80"; A9="80 80"
  echo "QUICK MODE -- results are not publication grade"
else
  # 02 and 03 now take B_CAL B_EVAL (calibration draws, held-out draws) over a
  # shared overlapping-pathway universe, not NREP B.
  A1="500 1000 200"; A2="500 500"; A3="500 500"; A4="200 1000"
  A8="400 400"; A9="500 500"
fi

run () { echo "=== $1 ==="; Rscript "$1" $2 2>&1 | tee "../logs/${1%.R}.log"; }
run 01_null_calibration.R "$A1"
run 02_breakage_stress.R  "$A2"
run 03_adaptive_loo.R     "$A3"
run 04_wy_vs_mvn.R        "$A4"
run 08_missing_correlation.R "$A8"
run 09_component_calibration.R "$A9"
run 05_figures.R          ""
run 07_tables_and_numbers.R ""
echo
echo "done. results/, figures/ and tex/auto/ populated."
echo "Now rebuild the manuscript:"
echo "  cd ../tex && pdflatex catfish_wy && bibtex catfish_wy && pdflatex catfish_wy && pdflatex catfish_wy"
echo
echo "06_maize_candidate_genes.R is NOT run here: it needs the MAGMA outputs"
echo "for the maize grain-nitrogen GWAS, which live outside this repo."
echo "Run it separately once those paths are available."
