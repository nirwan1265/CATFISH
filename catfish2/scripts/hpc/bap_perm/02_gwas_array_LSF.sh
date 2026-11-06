#!/bin/bash
# Step 2: GWAS, one array element per (permutation, chromosome).
# Index i -> perm = (i-1) %/% 10 + 1, chr = CHROMS[(i-1) %% 10]. 100 perms x 10 chr = 1000.
# Submit with the vcf2gwas conda environment active:  bsub < 02_gwas_array_LSF.sh
#BSUB -J "bap_gwas[1-1000]%25"
#BSUB -q sara
#BSUB -n 1
#BSUB -W 720
#BSUB -o logs/gwas.%J.%I.out
#BSUB -e logs/gwas.%J.%I.err
set -euo pipefail
source /share/maize/ntanduk/CATFISH/simulation/bap_perm_wy/config.sh
CHR_ARR=($CHROMS); NCHR=${#CHR_ARR[@]}
TID=$((LSB_JOBINDEX - 1)); B=$(( TID / NCHR + 1 )); CHR=${CHR_ARR[$(( TID % NCHR ))]}
bash "$WORKING_DIR/02_run_gwas_single.sh" "$B" "$CHR"
