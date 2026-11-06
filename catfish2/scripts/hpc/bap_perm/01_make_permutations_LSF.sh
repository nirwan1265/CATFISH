#!/bin/bash
# Step 1: N_PERM permuted phenotype files.   bsub < 01_make_permutations_LSF.sh
#BSUB -J bap_perm01
#BSUB -q sara
#BSUB -n 1
#BSUB -W 60
#BSUB -o logs/perm01.%J.out
#BSUB -e logs/perm01.%J.err
set -euo pipefail
source /share/maize/ntanduk/CATFISH/simulation/bap_perm_wy/config.sh
Rscript "$WORKING_DIR/01_make_permutations.R"
