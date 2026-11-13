#!/bin/bash
# =============================================================================
# Step 07b: PascalX gene and pathway scoring on every replicate, with the
# settings of compare/21_pascalx_lipids_LSF.sh (+/- 10 kb, MAF 0.01, NCBI37.3
# genes with Entrez IDs) on the simulation panel (step 00b) and SIM_SETS.
# Reads the p-value file from step 07a and removes it once scored.
# Submit with the PascalX env active.
#   pilot: bsub -J "sim_pascalx[1-20]" < 07b_pascalx_LSF.sh
#   full : MODE=full bsub -J "sim_pascalx[1-850]" < 07b_pascalx_LSF.sh
# =============================================================================
#BSUB -J "sim_pascalx[1-20]"
##BSUB -q sara
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -W 480
#BSUB -o logs/pascalx.%J.%I.out
#BSUB -e logs/pascalx.%J.%I.err
set -eo pipefail
source /share/maize/ntanduk/CATFISH/simulation_check_catfish/config.sh
R="$LSB_JOBINDEX"; RP=$(printf "rep_%04d" "$R")
TMP="${RUN}/pascalx/${RP}.pval.tmp"
[ -s "${PASCALX_REF}/sim_eur.chr1.db" ] || { echo "PascalX panel not imported: run step 00b" >&2; exit 1; }
[ -s "$PASCALX_GENOME" ] || { echo "missing $PASCALX_GENOME (written by compare/21)" >&2; exit 1; }
[ -s "$TMP" ] || { echo "missing $TMP: run step 07a" >&2; exit 1; }
export LD_LIBRARY_PATH="${PASCALX_DIR}/build/lib:${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
## -c: the panel has chromosomes 1-10 only (PascalX reads a Python list, e.g. [1,2,3])
## -p: cores, the -n requested above
PX_CHR="[${CHROMS// /,}]"
"${PASCALX_DIR}/pascalx" -p "${LSB_DJOB_NUMPROC:-1}" -c "$PX_CHR" -n True -m "$PASCALX_MAF" -w "$PASCALX_WINDOW" \
  -pw "$SIM_SETS" -cn 0 -cs 1 \
  "$PASCALX_GENOME" "${PASCALX_REF}/sim_eur" "${RUN}/pascalx/${RP}" \
  genescoring -sh True -cr 0 -cp 1 "$TMP"
[ -s "${RUN}/pascalx/${RP}.tsv" ] || { echo "no pathway output for $RP (p-value file kept: $TMP)" >&2; exit 1; }
rm -f "$TMP"                                   # only once PascalX has written its output
