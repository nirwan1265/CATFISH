#!/bin/bash
#BSUB -J sap_corpairs
#BSUB -n 1
#BSUB -W 180
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=4GB]"
#BSUB -o logs/corpairs.%J.out
#BSUB -e logs/corpairs.%J.err

set -euo pipefail
WORKING_DIR=/share/maize/ntanduk/CATFISH/simulation/sap_stem_volume_permutation
source "$WORKING_DIR/config.sap_stem_volume.sh"
mkdir -p "$WORKING_DIR/logs"
"${RSCRIPT:-Rscript}" "$WORKING_DIR/04_corpairs_once.R"
