#!/bin/bash
# =============================================================================
# 03_magma_array_LSF.sh  --  Permuted/synthetic MAGMA gene analysis (LSF array)
#                            via CATFISH::magma_gene (see 03_magma_catfish.R).
#                            One array element = one (phenotype index, chr).
#                            THE bottleneck step (~8h/chr like your real run).
# Submit AFTER step 02:  bsub < 03_magma_array_LSF.sh
# Array size = N_PHENO * 10.
# =============================================================================
#BSUB -J "magma[1-1000]%50"
#BSUB -n 1
#BSUB -W 5000
##BSUB -q sara
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=20GB]"
#BSUB -o logs/magma.%J.%I.out
#BSUB -e logs/magma.%J.%I.err

#set -euo pipefail
#source "${PIPE_CONFIG:-$(dirname "$0")/config.sh}"

set -euo pipefail
WORKING_DIR=/share/maize/ntanduk/CATFISH/simulation/permutation_pipeline
source "$WORKING_DIR/config.sh"

mkdir -p logs

CHR_ARR=($CHROMS); NCHR=${#CHR_ARR[@]}
TID=$((LSB_JOBINDEX - 1))
B=$(( TID / NCHR + 1 ))
CHR=${CHR_ARR[$(( TID % NCHR ))]}
echo "[magma] pheno=$B chr=$CHR (LSB_JOBINDEX=$LSB_JOBINDEX)"

Rscript "$WORKING_DIR/03_magma_catfish.R" "$B" "$CHR"
