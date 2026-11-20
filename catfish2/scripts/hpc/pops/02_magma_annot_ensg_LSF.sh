#!/bin/bash
# =============================================================================
# Step 02: MAGMA SNP-to-gene annotation in Ensembl IDs, for PoPS. Same panel
# and window as the paper's annotation (g1000_<panel>.genes.annot, +/- 10 kb);
# only the gene set and IDs differ (PoPS's gene_annot_jun10.txt, 18,383 genes).
#   PANEL=eur bsub < 02_magma_annot_ensg_LSF.sh       (default eur)
# Writes ${ROOT}/g1000_<panel>_ensg.genes.annot and ${ROOT}/pops_ensg.genes.loc
# Then the gene analysis with the paper's scripts, from scripts/hpc:
#   ANNOT_OVERRIDE=${ROOT}/g1000_eur_ensg.genes.annot OUT_SUFFIX=_ensg TRAIT=ldl bsub < 14_magma_trait_LSF.sh
#   OUT_SUFFIX=_ensg TRAIT=ldl bsub < 15_magma_merge_trait_LSF.sh
# =============================================================================
#BSUB -J pops_annot
##BSUB -q sara
#BSUB -n 1
#BSUB -W 120
#BSUB -o logs/pops_annot.%J.out
#BSUB -e logs/pops_annot.%J.err
set -euo pipefail
source /share/maize/ntanduk/CATFISH/pops/scripts/config.sh
PANEL="${PANEL:-eur}"
LOC=${ROOT}/pops_ensg.genes.loc
# MAGMA gene location: ID CHR START STOP (from ENSGID NAME CHR START END TSS)
awk -F'\t' 'NR > 1 { print $1 "\t" $3 "\t" $4 "\t" $5 }' "$GENE_ANNOT" > "$LOC"
echo "[annot] $(wc -l < "$LOC") genes in $LOC"
"$MAGMA" --annotate window=${WINDOW} \
         --snp-loc "${ROOT}/g1000_${PANEL}.bim" \
         --gene-loc "$LOC" \
         --out "${ROOT}/g1000_${PANEL}_ensg"
echo "[annot] wrote ${ROOT}/g1000_${PANEL}_ensg.genes.annot"
