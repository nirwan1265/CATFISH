#!/bin/bash
#BSUB -J "sap_catfish[1-100]%100"
#BSUB -n 1
#BSUB -W 3000
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=4GB]"
#BSUB -o logs/catfish.%J.%I.out
#BSUB -e logs/catfish.%J.%I.err

set -euo pipefail
WORKING_DIR=/share/maize/ntanduk/CATFISH/simulation/sap_stem_volume_permutation
source "$WORKING_DIR/config.sap_stem_volume.sh"
mkdir -p "$WORKING_DIR/logs"

export MVN_CALIBRATE_COMPONENTS="${MVN_CALIBRATE_COMPONENTS:-false}"
export PATCH_EMP_NULL_LOWER="${PATCH_EMP_NULL_LOWER:-false}"
export CATFISH_THREADS="${CATFISH_THREADS:-1}"
export TAU_GRID="${TAU_GRID:-0.1,0.05,0.01}"
export CATFISH_B_PERM="${CATFISH_B_PERM:-10000}"
export TAIL_MODE="${TAIL_MODE:-hybrid_gpd}"
export TAIL_SWITCH_EXCEED="${TAIL_SWITCH_EXCEED:-10}"
export TAIL_GPD_K="${TAIL_GPD_K:-250}"
export TAIL_MIN_B="${TAIL_MIN_B:-10000}"
export TAIL_MIN_TAIL="${TAIL_MIN_TAIL:-50}"

B=${LSB_JOBINDEX}
"${RSCRIPT:-Rscript}" "$WORKING_DIR/05_catfish_perm.R" "$B"
