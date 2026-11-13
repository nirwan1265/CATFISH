#!/bin/bash
# =============================================================================
# Step 07a: each replicate's SNP p-values for PascalX (SNP P, chr 1-10), from
# its simulated z (03_write_pval.R). Submit with the R env active. Step 07b
# (PascalX env) reads these files and removes them when it has scored them.
#   pilot: bsub -J "sim_pxpval[1-20]" < 07a_pascalx_pval_LSF.sh
#   full : MODE=full bsub -J "sim_pxpval[1-850]" < 07a_pascalx_pval_LSF.sh
# =============================================================================
#BSUB -J "sim_pxpval[1-20]"
##BSUB -q sara
#BSUB -n 1
#BSUB -W 60
#BSUB -o logs/pxpval.%J.%I.out
#BSUB -e logs/pxpval.%J.%I.err
set -euo pipefail
source /share/maize/ntanduk/CATFISH/simulation_check_catfish/config.sh
R="$LSB_JOBINDEX"; RP=$(printf "rep_%04d" "$R")
mkdir -p "${RUN}/pascalx"
Rscript "${SIM_ROOT}/03_write_pval.R" "$R" all pascalx "${RUN}/pascalx/${RP}.pval.tmp"
echo "[pxpval] $RP written"
