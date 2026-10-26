#!/usr/bin/env bash
# =============================================================================
# 02_annotate_local.sh -- SNP-to-gene annotation, run on the LAPTOP.
# MAGMA here is a macOS binary, so this cannot run inside the Linux sandbox.
#   bash scripts/human_ldl/02_annotate_local.sh
# Takes a few minutes. Produces ldl_eur.genes.annot, which you copy to the
# server along with everything else -- the cluster then only runs the gene
# analysis.
# =============================================================================
set -euo pipefail

MAGMA="${MAGMA:-$HOME/Documents/software/magma_v1.10_mac/magma}"
DATA="${DATA:-$HOME/Documents/Github/MAGCAT_data}"

BFILE="$DATA/g1000_eur"
GENELOC="$DATA/NCBI37.3.gene.loc"
SNPLOC="$DATA/g1000_eur.snploc"
PREFIX="$DATA/ldl_eur"
WINDOW="10,10"          # matches the sorghum analysis

[ -x "$MAGMA" ] || { echo "MAGMA not found at $MAGMA" >&2; exit 1; }

if [ ! -s "$SNPLOC" ]; then
  echo "[1/2] building SNP location file from ${BFILE}.bim ..."
  awk -F'\t' -v OFS='\t' '{print $2, $1, $4}' "${BFILE}.bim" > "$SNPLOC"
fi
echo "      snploc rows: $(wc -l < "$SNPLOC")"

echo "[2/2] annotating (window=$WINDOW) ..."
"$MAGMA" --annotate window="$WINDOW" \
         --snp-loc  "$SNPLOC" \
         --gene-loc "$GENELOC" \
         --out      "$PREFIX"

echo
echo "done -> ${PREFIX}.genes.annot"
echo "genes annotated: $(grep -vc '^#' "${PREFIX}.genes.annot")"
