#!/bin/bash
# =============================================================================
# Step 10: MAGMA gene analysis of the real dry-biomass GWAS with the SNP-wise
# mean model (the model CATFISH used; the original run was multi=snp-wise),
# one chromosome per array element, in chromosome-batch mode so 11 can merge
# them. Same SNP p-values, annotation and genotypes as the original run.
#   bsub < 10_magma_dry_mean_LSF.sh
# =============================================================================
#BSUB -J "cmp_magma_dry[1-10]"
#BSUB -q sara
#BSUB -n 1
#BSUB -W 720
#BSUB -o logs/magma_dry.%J.%I.out
#BSUB -e logs/magma_dry.%J.%I.err
set -euo pipefail
source /share/maize/ntanduk/CATFISH/compare/config.sh
CHR_ARR=($CHROMS); CHR=${CHR_ARR[$((LSB_JOBINDEX - 1))]}; K=$((10#$CHR))
PVAL=$(echo "$DRY_PVAL_TEMPLATE" | sed "s/@CHR@/${CHR}/g")
ANNOT=$(echo "$DRY_ANNOT_TEMPLATE" | sed "s/@CHR@/${CHR}/g")
BFILE=$(echo "$DRY_BFILE_TEMPLATE" | sed "s/@CHR@/${CHR}/g")
for f in "$PVAL" "$ANNOT" "${BFILE}.bed" "$MAGMA"; do [ -e "$f" ] || { echo "missing $f" >&2; exit 1; }; done
PREFIX="${OUT}/magma_dry_mean/Dry_tons_per_acre.snp_wise_mean"
"$MAGMA" --bfile "$BFILE" --gene-annot "$ANNOT" --pval "$PVAL" ncol=N use=SNP,P \
         --gene-model "$GENE_MODEL" --batch "$K" chr --out "$PREFIX"
[ -s "${PREFIX}.batch${K}_chr.genes.raw" ] || { echo "no genes.raw for chr $CHR" >&2; exit 1; }
