#!/bin/bash
# Step 5: Westfall-Young per permutation, one array element per permutation.
#   bsub < 05_wy_array_LSF.sh
#BSUB -J "bap_wy[1-100]%25"
#BSUB -q sara
#BSUB -n 1
#BSUB -W 720
#BSUB -o logs/wy.%J.%I.out
#BSUB -e logs/wy.%J.%I.err
set -euo pipefail
source /share/maize/ntanduk/CATFISH/simulation/bap_perm_wy/config.sh
Rscript "$WORKING_DIR/05_wy_perm.R" "$LSB_JOBINDEX"
