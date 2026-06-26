#!/bin/bash
#BSUB -n 1
#BSUB -W 5000
#BSUB -q sara
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=5GB]"
#BSUB -J make_permutations
#BSUB -o stdout.%J
#BSUB -e stderr.%J

#conda activate /rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/env_RTIGER

set -euo pipefail
WORKING_DIR=/share/maize/ntanduk/CATFISH/simulation/permutation_pipeline
source "$WORKING_DIR/config.sh"

Rscript "$WORKING_DIR/00_corpairs_once.R"
