#!/bin/bash
# Step 9: pool MAGMA's gene-set results over the permutations.   bsub < 09_aggregate_magma_gsa_LSF.sh
#BSUB -J bap_magma_aggregate
#BSUB -q sara
#BSUB -n 1
#BSUB -W 60
#BSUB -o logs/magma_aggregate.%J.out
#BSUB -e logs/magma_aggregate.%J.err
set -euo pipefail
source /share/maize/ntanduk/CATFISH/simulation/bap_perm_wy/config.sh
Rscript "$WORKING_DIR/09_aggregate_magma_gsa.R"
