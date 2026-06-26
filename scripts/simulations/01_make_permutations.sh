#!/bin/bash
#BSUB -n 1
#BSUB -W 14400
#BSUB -q sara
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=20GB]"
#BSUB -J make_permutations
#BSUB -o stdout.%J
#BSUB -e stderr.%J
# =============================================================================
# 01_make_permutations.sh  --  Submit wrapper for 01_make_permutations.R
# Submit:  bsub < 01_make_permutations.sh
# =============================================================================

#set -euo pipefail
#source "$(dirname "$0")/config.sh"

#Rscript "$(dirname "$0")/01_make_permutations.R"



set -euo pipefail
WORKING_DIR=/share/maize/ntanduk/CATFISH/simulation/permutation_pipeline
source "$WORKING_DIR/config.sh"

Rscript "$WORKING_DIR/01_make_permutations.R"
