#!/bin/bash
# =============================================================================
# Step 22: preranked GSEA (22_gsea_lipids.R) on the ten lipid sets. Needs the
# merged .genes.raw files in inputs/lipids/ (as step 11). Submit with the R env
# active (this script activates nothing); fgsea installed once in FGSEA_RLIB.
#   bsub < 22_gsea_lipids_LSF.sh
# =============================================================================
#BSUB -J "cmp_gsea_lipids"
##BSUB -q sara
#BSUB -n 1
#BSUB -W 480
#BSUB -o logs/gsea_lipids.%J.out
#BSUB -e logs/gsea_lipids.%J.err
set -euo pipefail
source /share/maize/ntanduk/CATFISH/compare/config.sh
Rscript "${ROOT}/22_gsea_lipids.R"
