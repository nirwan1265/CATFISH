#!/bin/bash
# =============================================================================
# Step 03r: the same MAGMA gene analysis as step 03 (MAGMA --batch <chr> chr,
# default gene model, the replicate's p-values written per chromosome to a
# temporary file), but one array element = one replicate, doing its 10
# chromosomes one after the other. Chromosomes already done are skipped.
# Ten times fewer job starts than step 03, for when many short jobs fail at
# start-up ("~/.lsbatch/...shell: No such file or directory").
# Element i = replicate i. Submit with the R env active, e.g. only the
# replicates still missing a chromosome:
#   MODE=full bsub -q sara -J "sim_magma_rep[10,67,88]%20" < 03r_magma_by_replicate_LSF.sh
# =============================================================================
#BSUB -J "sim_magma_rep[1-20]"
##BSUB -q sara
#BSUB -n 1
#BSUB -W 240
#BSUB -o logs/magma_rep.%J.%I.out
#BSUB -e logs/magma_rep.%J.%I.err
set -euo pipefail
source /share/maize/ntanduk/CATFISH/simulation_check_catfish/config.sh
R="$LSB_JOBINDEX"; RP=$(printf "rep_%04d" "$R")
OUTD="${RUN}/magma/${RP}"; mkdir -p "$OUTD"
for CHR in $CHROMS; do
  if [ -s "${OUTD}/${RP}.batch${CHR}_chr.genes.raw" ]; then echo "$RP chr $CHR already done"; continue; fi
  TMP="${OUTD}/${RP}.chr${CHR}.pval.tmp"
  Rscript "${SIM_ROOT}/03_write_pval.R" "$R" "$CHR" magma "$TMP"
  "$MAGMA" --bfile "$SIM_PANEL" --pval "$TMP" ncol=N --gene-annot "$SIM_ANNOT" \
           --batch "$CHR" chr --out "${OUTD}/${RP}"
  rm -f "$TMP"
  echo "[magma] $RP chr $CHR done"
done
n=$(ls "${OUTD}"/*batch*_chr.genes.raw 2>/dev/null | wc -l)
echo "[magma] $RP: $n of 10 chromosomes"
