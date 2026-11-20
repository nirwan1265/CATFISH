#!/bin/bash
# Step 32: sARTP inputs for the ten lipid sets (32_artp_prep.R), one job each.
# Needs ARTP_LAMBDA and ARTP_N set in config.sh.
#   bsub < 32_artp_prep_LSF.sh
#BSUB -J "cmp_artp_prep[1-10]"
#BSUB -q sara
#BSUB -n 1
#BSUB -W 600
#BSUB -o logs/artp_prep.%J.%I.out
#BSUB -e logs/artp_prep.%J.%I.err
set -eo pipefail
source /share/maize/ntanduk/CATFISH/compare/config.sh
source activate "$ARTP_CONDA" 2>/dev/null || conda activate "$ARTP_CONDA"
export R_LIBS_USER="$ARTP_RLIB"
A=($LIPIDS); TA=${A[$((LSB_JOBINDEX - 1))]}
Rscript "${ROOT}/32_artp_prep.R" "$TA"
