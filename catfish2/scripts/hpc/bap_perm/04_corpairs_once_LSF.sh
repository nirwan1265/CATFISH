#!/bin/bash
# Step 4: gene-gene correlation pairs, once (genotype-only, identical across
# permutations).   bsub < 04_corpairs_once_LSF.sh
#BSUB -J bap_corpairs
#BSUB -q sara
#BSUB -n 1
#BSUB -W 240
#BSUB -o logs/corpairs.%J.out
#BSUB -e logs/corpairs.%J.err
set -euo pipefail
source /share/maize/ntanduk/CATFISH/simulation/bap_perm_wy/config.sh
Rscript "$WORKING_DIR/04_corpairs_once.R"
