#!/bin/bash
#BSUB -n 1
#BSUB -W 720
#BSUB -q sara
#BSUB -J "bap_logo[4,10,11]"
#BSUB -o /share/maize/ntanduk/BAP_gwas/logs/logo.%J.%I.out
#BSUB -e /share/maize/ntanduk/BAP_gwas/logs/logo.%J.%I.err
############################################################################
## Step 06: leave-one-gene-out pathway evidence for the stem-composition gene
## score (06_logo_bap.R), one array element per trait, same trait numbering as
## 01-03 (4 = Brix_maturity, 10 = NDF_pct_DM, 11 = NFC_pct_DM). Needs steps
## 03 and 04 done for the trait. Submit with the R env that has CATFISH (R_plot).
##   bsub < 06_logo_LSF.sh
############################################################################
set -eu
TRAITS=(Anthesis_date_days Harvest_date_days Total_fresh_weight_kg Brix_maturity Brix_milk
        Dry_weight_kg Stalk_height_cm Dry_tons_per_acre ADF_pct_DM NDF_pct_DM NFC_pct_DM Lignin_pct_DM)
TRAIT=${TRAITS[$((LSB_JOBINDEX - 1))]}
Rscript /share/maize/ntanduk/BAP_gwas/scripts/06_logo_bap.R "$TRAIT"
