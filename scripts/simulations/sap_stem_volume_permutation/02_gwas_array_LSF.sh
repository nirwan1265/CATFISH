#!/bin/bash
#BSUB -J "sap_gwas[1-1000]%80"
#BSUB -n 1
#BSUB -W 5000
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=8GB]"
#BSUB -o logs/gwas.%J.%I.out
#BSUB -e logs/gwas.%J.%I.err

set -euo pipefail
WORKING_DIR=/share/maize/ntanduk/CATFISH/simulation/sap_stem_volume_permutation
source "$WORKING_DIR/config.sap_stem_volume.sh"
mkdir -p "$WORKING_DIR/logs"

CHR_ARR=($CHROMS)
NCHR=${#CHR_ARR[@]}
TID=$((LSB_JOBINDEX - 1))
B=$(( TID / NCHR + 1 ))
CHR=${CHR_ARR[$(( TID % NCHR ))]}

bash "$WORKING_DIR/02_run_gwas_single.sh" "$B" "$CHR"
