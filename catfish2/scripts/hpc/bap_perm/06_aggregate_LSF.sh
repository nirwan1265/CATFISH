#!/bin/bash
# Step 6: pool the permutations.   bsub < 06_aggregate_LSF.sh
#BSUB -J bap_aggregate
#BSUB -q sara
#BSUB -n 1
#BSUB -W 60
#BSUB -o logs/aggregate.%J.out
#BSUB -e logs/aggregate.%J.err
set -euo pipefail
source /share/maize/ntanduk/CATFISH/simulation/bap_perm_wy/config.sh
Rscript "$WORKING_DIR/06_aggregate.R"
