#!/bin/bash
#BSUB -J make_sap_perms
#BSUB -n 1
#BSUB -W 120
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=4GB]"
#BSUB -o logs/make_perms.%J.out
#BSUB -e logs/make_perms.%J.err

set -euo pipefail
WORKING_DIR=/share/maize/ntanduk/CATFISH/simulation/sap_stem_volume_permutation
source "$WORKING_DIR/config.sap_stem_volume.sh"
mkdir -p "$WORKING_DIR/logs"
"${RSCRIPT:-Rscript}" "$WORKING_DIR/01_make_permutations.R"
