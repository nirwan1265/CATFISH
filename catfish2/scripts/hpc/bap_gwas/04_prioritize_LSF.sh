#!/bin/bash
#BSUB -n 1
#BSUB -W 720
#BSUB -q sara
#BSUB -J bap_prioritize
#BSUB -o /share/maize/ntanduk/BAP_gwas/logs/prioritize.%J.out
#BSUB -e /share/maize/ntanduk/BAP_gwas/logs/prioritize.%J.err
############################################################################
## Step 04: CATFISH gene prioritisation for ADF, NDF, NFC, Brix (milk,
## maturity) and the genes they share (04_prioritize_bap.R). One job, traits in
## turn. Submit with the R env that has CATFISH active (R_plot).
##   bsub < 04_prioritize_LSF.sh
## Other traits: TRAITS=ADF_pct_DM,Lignin_pct_DM bsub < 04_prioritize_LSF.sh
############################################################################
set -eu
Rscript /share/maize/ntanduk/BAP_gwas/scripts/04_prioritize_bap.R
