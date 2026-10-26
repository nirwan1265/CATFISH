#!/usr/bin/env bash
# MAGMA gene analysis for GLGC LDL-C (EUR), build 37.
# Produces the gene-level p-values and the .genes.raw gene-gene correlation
# matrix that CATFISH consumes.
#
# Run 01_munge_glgc_ldl.sh first.
set -euo pipefail

MAGMA="${MAGMA:-$HOME/Documents/software/magma_v1.10_mac/magma}"
DATA="${DATA:-$HOME/Documents/Github/MAGCAT_data}"
OUTDIR="$DATA/ldl/magma"
PREFIX="$OUTDIR/ldl_eur"

GENELOC="$DATA/NCBI37.3.gene.loc"
BFILE="$DATA/g1000_eur"
PVAL="$DATA/ldl/ldl_eur.magma.pval"
SNPLOC="$DATA/ldl/g1000_eur.snploc"

mkdir -p "$OUTDIR"
[ -x "$MAGMA" ] || { echo "MAGMA not executable at $MAGMA" >&2; exit 1; }
[ -f "$PVAL" ]  || { echo "missing $PVAL -- run 01_munge_glgc_ldl.sh first" >&2; exit 1; }

# --- SNP location file: rsID, chromosome, bp (from the reference .bim) ------
if [ ! -s "$SNPLOC" ]; then
  echo "[1/3] building SNP location file ..."
  awk -F'\t' -v OFS='\t' '{print $2, $1, $4}' "$BFILE.bim" > "$SNPLOC"
fi
echo "      snploc rows: $(wc -l < "$SNPLOC")"

# --- Annotation: map SNPs to genes, 35 kb upstream / 10 kb downstream -------
# 35,10 is the conventional window in the MAGMA literature; it captures most
# promoter and proximal-regulatory variation without merging adjacent genes.
echo "[2/3] annotating SNPs to genes ..."
"$MAGMA" --annotate window=35,10 \
         --snp-loc "$SNPLOC" \
         --gene-loc "$GENELOC" \
         --out "$PREFIX"

# --- Gene analysis ----------------------------------------------------------
# ncol=N uses the per-variant sample size from the GLGC file.
# This is the slow step: expect roughly 1-3 hours.
echo "[3/3] gene analysis (slow -- 1-3 hours) ..."
"$MAGMA" --bfile "$BFILE" \
         --pval "$PVAL" ncol=N \
         --gene-annot "$PREFIX.genes.annot" \
         --out "$PREFIX"

echo
echo "done. key outputs:"
ls -lh "$PREFIX".genes.out "$PREFIX".genes.raw 2>/dev/null
echo
echo "  .genes.out -> gene-level p-values (CATFISH input)"
echo "  .genes.raw -> gene-gene correlations (CATFISH MVN calibration)"
