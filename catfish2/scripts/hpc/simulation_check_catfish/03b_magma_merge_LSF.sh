#!/bin/bash
# =============================================================================
# Step 03b: merge each replicate's 10 chromosome batches (as 15_magma_merge),
# giving RUN/magma/rep_<r>/rep_<r>.genes.out and .genes.raw. One element per
# replicate. Submit when step 03 has finished.
#   pilot: bsub -J "sim_merge[1-20]" < 03b_magma_merge_LSF.sh
#   full : MODE=full bsub -J "sim_merge[1-850]" < 03b_magma_merge_LSF.sh
# =============================================================================
#BSUB -J "sim_merge[1-20]"
##BSUB -q sara
#BSUB -n 1
#BSUB -W 60
#BSUB -o logs/merge.%J.%I.out
#BSUB -e logs/merge.%J.%I.err
set -euo pipefail
source /share/maize/ntanduk/CATFISH/simulation_check_catfish/config.sh
RP=$(printf "rep_%04d" "$LSB_JOBINDEX"); PREFIX="${RUN}/magma/${RP}/${RP}"
missing=0
for c in $CHROMS; do
  [ -s "${PREFIX}.batch${c}_chr.genes.raw" ] || { echo "[merge] $RP MISSING chr $c"; missing=$((missing + 1)); }
done
[ "$missing" -eq 0 ] || { echo "[merge] $RP: $missing chromosome(s) missing; rerun step 03 for them" >&2; exit 1; }
"$MAGMA" --merge "$PREFIX" --out "$PREFIX"
echo "[merge] $RP: $(( $(wc -l < "${PREFIX}.genes.out") - 1 )) genes"
