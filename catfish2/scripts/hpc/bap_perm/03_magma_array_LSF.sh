#!/bin/bash
# Step 3: MAGMA gene analysis, one array element per (permutation, chromosome),
# same index mapping as step 2. Plain MAGMA (03_magma_single.sh), no R, with the
# genotype-only annotation from step 2b.   bsub < 03_magma_array_LSF.sh
#BSUB -J "bap_magma[1-1000]%25"
#BSUB -q sara
#BSUB -n 1
#BSUB -W 720
#BSUB -o logs/magma.%J.%I.out
#BSUB -e logs/magma.%J.%I.err
set -euo pipefail
source /share/maize/ntanduk/CATFISH/simulation/bap_perm_wy/config.sh
CHR_ARR=($CHROMS); NCHR=${#CHR_ARR[@]}
TID=$((LSB_JOBINDEX - 1)); B=$(( TID / NCHR + 1 )); CHR=${CHR_ARR[$(( TID % NCHR ))]}
bash "$WORKING_DIR/03_magma_single.sh" "$B" "$CHR"
