#!/bin/bash
# =============================================================================
# Step 33: sARTP for the ten lipid sets, each split into ARTP_NCHUNK chunks of
# pathways: job i runs set ceil(i / ARTP_NCHUNK), chunk ((i - 1) % ARTP_NCHUNK) + 1.
# With ARTP_NCHUNK = 30 the array is 1-300. Test one chunk first:
#   bsub -J "cmp_artp[1]" < 33_artp_run_LSF.sh
#   bsub < 33_artp_run_LSF.sh
# A rerun skips pathways already written.
# =============================================================================
#BSUB -J "cmp_artp[1-300]"
#BSUB -q sara
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -W 2880
#BSUB -o logs/artp.%J.%I.out
#BSUB -e logs/artp.%J.%I.err
set -eo pipefail
source /share/maize/ntanduk/CATFISH/compare/config.sh
source activate "$ARTP_CONDA" 2>/dev/null || conda activate "$ARTP_CONDA"
export R_LIBS_USER="$ARTP_RLIB"
A=($LIPIDS); i=$((LSB_JOBINDEX - 1))
TA=${A[$((i / ARTP_NCHUNK))]}; K=$(( (i % ARTP_NCHUNK) + 1 ))
[ -s "${ARTP_IN}/${TA}.summary.rds" ] || { echo "no inputs for $TA: run step 32" >&2; exit 1; }
Rscript "${ROOT}/33_artp_run.R" "$TA" "$K"
