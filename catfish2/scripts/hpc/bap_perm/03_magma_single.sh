#!/bin/bash
# =============================================================================
# 03_magma_single.sh -- plain MAGMA gene analysis for ONE (permutation, chr),
# no R: the same MAGMA call the real dry-tons run made through magma_gene()
#   SNP p-value file : SNP, P = p_wald, N = N_TOTAL - n_miss  (from GEMMA assoc)
#   magma --bfile <chr bfile> --gene-annot <chr annot> --pval <file> ncol=N use=SNP,P
#         --gene-model $GENE_MODEL --out perm_XXX_chrNN.<model>
#   bash 03_magma_single.sh <perm> <chr, e.g. 01>
# BATCH=1: the same call in MAGMA chromosome-batch mode (--batch <chr> chr),
# written to ${MAGMA_OUT}_batch/perm_XXX/perm_XXX.<model>.batch<k>_chr.*, so
# that magma --merge can combine the ten chromosomes into one .genes.raw for
# MAGMA's own gene-set analysis (07, 08). The default run is unchanged.
# =============================================================================
set -euo pipefail
source /share/maize/ntanduk/CATFISH/simulation/bap_perm_wy/config.sh
B="$1"; CHR="$2"; BPAD=$(printf "%03d" "$B")
ASSOC="${GWAS_DIR}/perm_${BPAD}/chr${CHR}/assoc.txt"
ANNOT="$(echo "$ANNOT_PATH_TEMPLATE" | sed "s/@CHR@/${CHR}/g")"
BFILE="${GENO_DIR}/$(echo "$GENO_PREFIX_PATTERN" | sed "s/@CHR@/${CHR}/g")"
N_TOTAL=$(head -1 "$NTOTAL_FILE")
TAG=$(echo "$GENE_MODEL" | sed 's/[^A-Za-z0-9]\{1,\}/_/g')
if [ "${BATCH:-0}" = "1" ]; then
  K=$((10#$CHR))
  OUTDIR="${MAGMA_OUT}_batch/perm_${BPAD}"
  PREFIX="${OUTDIR}/perm_${BPAD}.${TAG}"
  BATCH_ARGS=(--batch "$K" chr)
  OUT_STEM="${PREFIX}.batch${K}_chr"
else
  OUTDIR="${MAGMA_OUT}/perm_${BPAD}/chr${CHR}"
  PREFIX="${OUTDIR}/perm_${BPAD}_chr${CHR}.${TAG}"
  BATCH_ARGS=()
  OUT_STEM="$PREFIX"
fi
PVAL="${OUTDIR}/perm_${BPAD}_chr${CHR}.pval.txt"
for f in "$ASSOC" "$ANNOT" "${BFILE}.bed" "$MAGMA"; do [ -e "$f" ] || { echo "missing $f" >&2; exit 1; }; done
mkdir -p "$OUTDIR"
echo "[magma] perm=$B chr=$CHR N_TOTAL=$N_TOTAL model=$GENE_MODEL batch=${BATCH:-0}"

awk -v NT="$N_TOTAL" '
  NR == 1 { for (i = 1; i <= NF; i++) { if ($i == "rs") s = i; if ($i == "p_wald") p = i; if ($i == "n_miss") m = i }
            if (!s || !p || !m) { print "assoc.txt has no rs / p_wald / n_miss column" > "/dev/stderr"; exit 1 }
            print "SNP\tP\tN"; next }
  { print $s "\t" $p "\t" (NT - $m) }' "$ASSOC" > "$PVAL"

"$MAGMA" --bfile "$BFILE" \
         --gene-annot "$ANNOT" \
         --pval "$PVAL" ncol=N use=SNP,P \
         --gene-model "$GENE_MODEL" \
         ${BATCH_ARGS[@]+"${BATCH_ARGS[@]}"} \
         --out "$PREFIX"

[ -s "${OUT_STEM}.genes.out" ] && [ -s "${OUT_STEM}.genes.raw" ] || { echo "[magma] no genes.out/genes.raw for $OUT_STEM; see ${OUT_STEM}.log" >&2; exit 1; }
rm -f "$PVAL"
echo "[magma] -> ${OUT_STEM}.genes.out ($(($(grep -vc '^#' "${OUT_STEM}.genes.out") - 1)) genes)"
