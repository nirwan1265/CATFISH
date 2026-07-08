#!/bin/bash
#BSUB -J "sap_magma[1-1000]%100"
#BSUB -n 1
#BSUB -W 3000
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=6GB]"
#BSUB -o logs/magma.%J.%I.out
#BSUB -e logs/magma.%J.%I.err

set -euo pipefail
WORKING_DIR=/share/maize/ntanduk/CATFISH/simulation/sap_stem_volume_permutation
source "$WORKING_DIR/config.sap_stem_volume.sh"
mkdir -p "$WORKING_DIR/logs"

CHR_ARR=($CHROMS)
NCHR=${#CHR_ARR[@]}
TID=$((LSB_JOBINDEX - 1))
B=$(( TID / NCHR + 1 ))
CHR=${CHR_ARR[$(( TID % NCHR ))]}

"${RSCRIPT:-Rscript}" "$WORKING_DIR/03_magma_catfish.R" "$B" "$CHR"
