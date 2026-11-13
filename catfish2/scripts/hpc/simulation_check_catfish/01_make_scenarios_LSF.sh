#!/bin/bash
# =============================================================================
# Step 01: design, causal SNPs and per-person vectors u (01_make_scenarios.R).
#   bsub < 01_make_scenarios_LSF.sh                  (pilot)
#   MODE=full EFFECT_SCALE=<from pilot> bsub < 01_make_scenarios_LSF.sh
# =============================================================================
#BSUB -J sim_design
##BSUB -q sara
#BSUB -n 1
#BSUB -W 240
#BSUB -o logs/design.%J.out
#BSUB -e logs/design.%J.err
set -euo pipefail
source /share/maize/ntanduk/CATFISH/simulation_check_catfish/config.sh
for f in "${SIM_PANEL}.bed" "$SIM_ANNOT" "$SIM_SETS" "${IN}/genes.tsv"; do
  [ -s "$f" ] || { echo "missing $f: run step 00" >&2; exit 1; }
done
echo "[design] MODE=$MODE EFFECT_SCALE=$EFFECT_SCALE -> $RUN"
Rscript "${SIM_ROOT}/01_make_scenarios.R"
