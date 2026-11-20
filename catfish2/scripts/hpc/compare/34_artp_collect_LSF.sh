#!/bin/bash
# Step 34: gather the sARTP chunks into out/artp/<ta>.csv.
#   bsub < 34_artp_collect_LSF.sh
#BSUB -J cmp_artp_collect
#BSUB -q sara
#BSUB -n 1
#BSUB -W 60
#BSUB -o logs/artp_collect.%J.out
#BSUB -e logs/artp_collect.%J.err
set -eo pipefail
source /share/maize/ntanduk/CATFISH/compare/config.sh
source activate "$ARTP_CONDA" 2>/dev/null || conda activate "$ARTP_CONDA"
export R_LIBS_USER="$ARTP_RLIB"
Rscript "${ROOT}/34_artp_collect.R"
