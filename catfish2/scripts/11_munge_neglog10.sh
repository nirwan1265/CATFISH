#!/usr/bin/env bash
# =============================================================================
# 11_munge_neglog10.sh -- write SNP + -log10(p) for the GLGC LDL-C summary
# statistics, with no floor.
#
# WHY
#   01_munge_glgc_ldl.sh produces MAGMA's --pval format, which takes a
#   probability. GLGC's strongest variants have -log10(p) up to 1725, far below
#   the smallest representable double (~2.2e-308), so that script floors them at
#   1e-300. MAGMA needs that floor and is unaffected by it -- but the
#   candidate-gene GWAS layer in 06_candidate_genes.R ranks genes by their
#   minimum SNP p-value, and the floor collapses the 53 strongest genes
#   (LDLR, PCSK9, APOB, APOE, LPA, HMGCR among them) into a single tie.
#
#   -log10(p) has no such limit. This script carries it through unfloored so the
#   GWAS layer can order the top of the list.
#
# Output: SNP  NEGLOG10P   (tab separated, with header)
# Usage:  bash 11_munge_neglog10.sh
# =============================================================================
set -euo pipefail

DATA="${DATA:-$HOME/Documents/Github/MAGCAT_data}"
IN="${IN:-$DATA/ldl/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz}"
BIM="${BIM:-$DATA/g1000_eur.bim}"
OUT="${OUT:-$DATA/ldl/ldl_eur.neglog10.txt}"

[ -f "$IN" ]  || { echo "missing: $IN"  >&2; exit 1; }
[ -f "$BIM" ] || { echo "missing: $BIM" >&2; exit 1; }

if command -v gzcat >/dev/null 2>&1; then ZCAT=gzcat; else ZCAT=zcat; fi
echo "streaming $IN ..."

"$ZCAT" "$IN" | awk -F'\t' -v OFS='\t' -v BIM="$BIM" '
  BEGIN {
    while ((getline line < BIM) > 0) { split(line, b, "\t"); ref[b[2]] = 1; nref++ }
    close(BIM)
    printf("reference SNPs loaded: %d\n", nref) > "/dev/stderr"
    print "SNP", "NEGLOG10P"
  }
  NR == 1 { next }
  $1 !~ /^rs/  { nonrs++;  next }
  !($1 in ref) { notref++; next }
  {
    nl = $11 + 0
    if (nl != nl || $11 == "") { bad++; next }
    if (nl > maxnl) { maxnl = nl }
    if (nl > 300) { above300++ }
    printf "%s\t%.6f\n", $1, nl
    kept++
  }
  END {
    printf("kept=%d  not_in_panel=%d  non_rs=%d  bad=%d\n",
           kept, notref, nonrs, bad) > "/dev/stderr"
    printf("max -log10(p) = %.1f ; variants above 300 (floored by 01) = %d\n",
           maxnl, above300) > "/dev/stderr"
  }
' > "$OUT"

echo "wrote $OUT"
ls -lh "$OUT"
