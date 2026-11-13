#!/bin/bash
# =============================================================================
# Step 06: MAGMA gene-set analysis on every replicate, competitive and
# --model self-contained (as compare/11_magma_gsa_real_LSF.sh), on SIM_SETS.
# After step 03b.
#   pilot: bsub -J "sim_gsa[1-20]" < 06_magma_gsa_LSF.sh
#   full : MODE=full bsub -J "sim_gsa[1-850]" < 06_magma_gsa_LSF.sh
# =============================================================================
#BSUB -J "sim_gsa[1-20]"
##BSUB -q sara
#BSUB -n 1
#BSUB -W 60
#BSUB -o logs/gsa.%J.%I.out
#BSUB -e logs/gsa.%J.%I.err
set -euo pipefail
source /share/maize/ntanduk/CATFISH/simulation_check_catfish/config.sh
RP=$(printf "rep_%04d" "$LSB_JOBINDEX"); RAW="${RUN}/magma/${RP}/${RP}.genes.raw"
[ -s "$RAW" ] || { echo "missing $RAW: run step 03b" >&2; exit 1; }
mkdir -p "${RUN}/magma_gsa"
"$MAGMA" --gene-results "$RAW" --set-annot "$SIM_SETS" --model self-contained \
         --out "${RUN}/magma_gsa/${RP}"
[ -s "${RUN}/magma_gsa/${RP}.gsa.out" ] || { echo "no gsa.out for $RP" >&2; exit 1; }
