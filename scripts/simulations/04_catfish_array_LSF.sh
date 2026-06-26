#!/bin/bash
# =============================================================================
# 04_catfish_array_LSF.sh
# Single production CATFISH array entry point for phenotype permutations.
#
# Defaults match the final chosen setup:
#   - mvn_calibrate_components = FALSE  (direct omnibus MVN calibration)
#   - TFisher tau grid         = 0.01,0.05,0.5,1
#   - tail_mode                = hybrid_gpd
#   - B_perm                   = 10000
# =============================================================================
#BSUB -J "catfish[1-100]%100"
#BSUB -n 1
#BSUB -W 3000
##BSUB -q sara
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=4GB]"
#BSUB -o logs/catfish.%J.%I.out
#BSUB -e logs/catfish.%J.%I.err

set -euo pipefail
PIPE_DIR=/share/maize/ntanduk/CATFISH/simulation/permutation_pipeline
source "${PIPE_CONFIG:-$PIPE_DIR/config.sh}"
mkdir -p logs

export MVN_CALIBRATE_COMPONENTS="${MVN_CALIBRATE_COMPONENTS:-false}"
export PATCH_EMP_NULL_LOWER="${PATCH_EMP_NULL_LOWER:-false}"
export CATFISH_THREADS=1
export TAU_GRID="${TAU_GRID:-0.01,0.05,0.5,1}"
export CATFISH_B_PERM="${CATFISH_B_PERM:-10000}"
export TAIL_MODE="${TAIL_MODE:-hybrid_gpd}"
export TAIL_SWITCH_EXCEED="${TAIL_SWITCH_EXCEED:-10}"
export TAIL_GPD_K="${TAIL_GPD_K:-250}"
export TAIL_MIN_B="${TAIL_MIN_B:-10000}"
export TAIL_MIN_TAIL="${TAIL_MIN_TAIL:-50}"

B=${LSB_JOBINDEX}
Rscript "$PIPE_DIR/04_catfish_perm.R" "$B"
