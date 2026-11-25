#!/bin/bash
# =============================================================================
# Step 07d: preranked GSEA (07d_gsea.R) on the MAGMA gene z-scores, one array
# element = one design cell (all its replicates). Needs step 03b (MAGMA genes.out).
# Submit with the R env active (this script activates nothing); fgsea must be
# installed once in FGSEA_RLIB (README, "GSEA").
#   full:  MODE=full bsub -J "sim_gsea[1-17]" < 07d_gsea_LSF.sh
# =============================================================================
#BSUB -J "sim_gsea[1-17]"
##BSUB -q sara
#BSUB -n 1
#BSUB -W 480
#BSUB -o logs/gsea.%J.%I.out
#BSUB -e logs/gsea.%J.%I.err
set -euo pipefail
source /share/maize/ntanduk/CATFISH/simulation_check_catfish/config.sh
Rscript "${SIM_ROOT}/07d_gsea.R" "$LSB_JOBINDEX"
