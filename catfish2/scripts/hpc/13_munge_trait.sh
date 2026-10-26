#!/usr/bin/env bash
# =============================================================================
# 13_munge_trait.sh -- trait-parameterised version of 01_munge_glgc_ldl.sh.
#
# Writes BOTH files the pipeline needs, in one pass over the 2.3 GB input:
#   <trait>_eur.magma.pval    SNP P N        for MAGMA (p floored at 1e-300)
#   <trait>_eur.neglog10.txt  SNP NEGLOG10P  for the candidate-gene GWAS layer
#                                            (unfloored -- see 11_munge_neglog10.sh)
#
#   TRAIT=HDL bash 13_munge_trait.sh
#   for T in HDL logTG TC; do TRAIT=$T bash 13_munge_trait.sh; done
# =============================================================================
set -euo pipefail

TRAIT="${TRAIT:?set TRAIT, e.g. TRAIT=HDL}"
. "$(dirname "$0")/env.sh"

t="$(echo "$TRAIT" | tr '[:upper:]' '[:lower:]')"
TDIR="$DATA/$t"                      # mirrors the existing $DATA/ldl/
IN="${IN:-$TDIR/${TRAIT}_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz}"
PVAL="$TDIR/${t}_eur.magma.pval"
NLOG="$TDIR/${t}_eur.neglog10.txt"
mkdir -p "$TDIR"

[ -f "$IN" ]  || { echo "missing: $IN  (run 12_fetch_glgc.sh first)" >&2; exit 1; }
[ -f "$BIM" ] || { echo "missing: $BIM" >&2; exit 1; }

if command -v gzcat >/dev/null 2>&1; then ZCAT=gzcat; else ZCAT=zcat; fi
echo "streaming $IN ..."

"$ZCAT" "$IN" | awk -F'\t' -v OFS='\t' -v BIM="$BIM" -v PVAL="$PVAL" -v NLOG="$NLOG" '
  BEGIN {
    while ((getline line < BIM) > 0) { split(line, b, "\t"); ref[b[2]] = 1; nref++ }
    close(BIM)
    printf("reference SNPs loaded: %d\n", nref) > "/dev/stderr"
    print "SNP", "P", "N"        > PVAL
    print "SNP", "NEGLOG10P"     > NLOG
  }
  NR == 1 { next }
  $1 !~ /^rs/  { nonrs++;  next }
  !($1 in ref) { notref++; next }
  {
    nl = $11 + 0
    if (nl != nl || $11 == "") { bad++; next }
    if (nl > maxnl) maxnl = nl
    if (nl > 300) { p = 1e-300; floored++ } else { p = 10 ^ (-nl) }
    if (p <= 0)   { p = 1e-300; floored++ }
    printf "%s\t%.6g\t%d\n", $1, p, $6 > PVAL
    printf "%s\t%.6f\n",     $1, nl     > NLOG
    kept++
  }
  END {
    printf("kept=%d not_in_panel=%d non_rs=%d bad=%d floored=%d max_neglog10=%.1f\n",
           kept, notref, nonrs, bad, floored, maxnl) > "/dev/stderr"
  }
'
echo "wrote $PVAL"
echo "wrote $NLOG"
ls -lh "$PVAL" "$NLOG"
