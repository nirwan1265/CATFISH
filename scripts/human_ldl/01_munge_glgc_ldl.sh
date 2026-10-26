#!/usr/bin/env bash
# Convert GLGC 2021 LDL-C (EUR) summary statistics into MAGMA --pval format.
#
# Two things this handles that otherwise break the analysis:
#
# 1. GLGC's `pvalue` column underflows to exactly 0 at the strongest loci
#    (473 such variants on chr1 alone, including rs12740374 at SORT1, where
#    -log10(p) = 1725). We rebuild p from `pvalue_neg_log10` and floor it at
#    1e-300, the smallest value a double reliably represents.
# 2. MAGMA can only use SNPs present in the 1000G EUR reference panel, so we
#    filter to those here: ~91M rows down to ~22M, saving about 2 GB.
#
# Input  cols: 1 rsID 2 CHROM 3 POS_b37 4 REF 5 ALT 6 N 7 N_studies
#              8 POOLED_ALT_AF 9 EFFECT_SIZE 10 SE 11 pvalue_neg_log10
#              12 pvalue 13 pvalue_neg_log10_GC 14 pvalue_GC
# Output cols: SNP P N          (overwrites OUT each run)
set -euo pipefail

DATA="${DATA:-$HOME/Documents/Github/MAGCAT_data}"
IN="$DATA/ldl/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"
BIM="$DATA/g1000_eur.bim"
OUT="$DATA/ldl/ldl_eur.magma.pval"

[ -f "$IN" ]  || { echo "missing: $IN"  >&2; exit 1; }
[ -f "$BIM" ] || { echo "missing: $BIM" >&2; exit 1; }

if command -v gzcat >/dev/null 2>&1; then ZCAT=gzcat; else ZCAT=zcat; fi
echo "streaming $IN ..."

"$ZCAT" "$IN" | awk -F'\t' -v OFS='\t' -v BIM="$BIM" '
  BEGIN {
    while ((getline line < BIM) > 0) { split(line, b, "\t"); ref[b[2]] = 1; nref++ }
    close(BIM)
    printf("reference SNPs loaded: %d\n", nref) > "/dev/stderr"
    print "SNP", "P", "N"
  }
  NR == 1 { next }
  $1 !~ /^rs/  { nonrs++;  next }
  !($1 in ref) { notref++; next }
  {
    nl = $11 + 0
    if (nl != nl || $11 == "") { bad++; next }
    if (nl > 300) { p = 1e-300; floored++ } else { p = 10 ^ (-nl) }
    if (p <= 0)   { p = 1e-300; floored++ }
    printf "%s\t%.6g\t%d\n", $1, p, $6
    kept++
  }
  END {
    printf("kept=%d  not_in_panel=%d  non_rs=%d  bad=%d  floored=%d\n",
           kept, notref, nonrs, bad, floored) > "/dev/stderr"
  }
' > "$OUT"

echo "wrote $OUT"
ls -lh "$OUT"
