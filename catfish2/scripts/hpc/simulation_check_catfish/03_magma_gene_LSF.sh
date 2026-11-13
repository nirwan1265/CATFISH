#!/bin/bash
# =============================================================================
# Step 03: MAGMA gene analysis, one array element = one replicate x one
# chromosome (MAGMA --batch <chr> chr, default gene model, as the lipid runs in
# 14_magma_trait_LSF.sh). Element i is replicate REP_FROM + (i - 1) / 10,
# chromosome (i - 1) % 10 + 1. The replicate's SNP p-values for that chromosome
# are written to a temporary file and removed when MAGMA has finished.
#   pilot (20 replicates):  bsub -J "sim_magma[1-200]" < 03_magma_gene_LSF.sh
#   full, in blocks of 100 replicates (1,000 elements each):
#     MODE=full REP_FROM=1   bsub -J "sim_magma[1-1000]" < 03_magma_gene_LSF.sh
#     MODE=full REP_FROM=101 bsub -J "sim_magma[1-1000]" < 03_magma_gene_LSF.sh  ...
# =============================================================================
#BSUB -J "sim_magma[1-200]"
##BSUB -q sara
#BSUB -n 1
#BSUB -W 240
#BSUB -o logs/magma.%J.%I.out
#BSUB -e logs/magma.%J.%I.err
set -euo pipefail
source /share/maize/ntanduk/CATFISH/simulation_check_catfish/config.sh
REP_FROM="${REP_FROM:-1}"
i=$((LSB_JOBINDEX - 1)); R=$((REP_FROM + i / 10)); CHR=$((i % 10 + 1))
NREP=$(( $(wc -l < "${RUN}/design.tsv") - 1 ))
[ "$R" -le "$NREP" ] || { echo "replicate $R beyond the design ($NREP); nothing to do"; exit 0; }
RP=$(printf "rep_%04d" "$R")
OUTD="${RUN}/magma/${RP}"; mkdir -p "$OUTD"
[ -s "${OUTD}/${RP}.batch${CHR}_chr.genes.raw" ] && { echo "$RP chr $CHR already done"; exit 0; }
TMP="${OUTD}/${RP}.chr${CHR}.pval.tmp"
Rscript "${SIM_ROOT}/03_write_pval.R" "$R" "$CHR" magma "$TMP"
"$MAGMA" --bfile "$SIM_PANEL" --pval "$TMP" ncol=N --gene-annot "$SIM_ANNOT" \
         --batch "$CHR" chr --out "${OUTD}/${RP}"
rm -f "$TMP"
echo "[magma] $RP chr $CHR done"
