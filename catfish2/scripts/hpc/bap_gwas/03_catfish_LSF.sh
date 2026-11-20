#!/bin/bash
#BSUB -n 1
#BSUB -W 720
#BSUB -q sara
#BSUB -J "bap_catfish[1-12]"
#BSUB -o /share/maize/ntanduk/BAP_gwas/logs/catfish.%J.%I.out
#BSUB -e /share/maize/ntanduk/BAP_gwas/logs/catfish.%J.%I.err
############################################################################
## Step 03: CATFISH Westfall-Young (self-contained + competitive) per BAP
## trait, one array element per trait (03_catfish_bap.R). Submit with the R env
## that has CATFISH active (R_plot).
##   all:   bsub < 03_catfish_LSF.sh
##   one:   bsub -J "bap_catfish[7]" < 03_catfish_LSF.sh     (7 = Stalk_height_cm)
############################################################################
set -eu
TRAITS=(Anthesis_date_days Harvest_date_days Total_fresh_weight_kg Brix_maturity Brix_milk
        Dry_weight_kg Stalk_height_cm Dry_tons_per_acre ADF_pct_DM NDF_pct_DM NFC_pct_DM Lignin_pct_DM)
TRAIT=${TRAITS[$((LSB_JOBINDEX - 1))]}
Rscript /share/maize/ntanduk/BAP_gwas/scripts/03_catfish_bap.R "$TRAIT"
