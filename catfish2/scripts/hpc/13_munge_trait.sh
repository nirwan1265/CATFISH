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

#   ANC=AFR TRAIT=LDL bash 13_munge_trait.sh   # other ancestries (default EUR);
#   SNPs are filtered to that ancestry's 1000G panel (HIS uses the AMR panel)
TRAIT="${TRAIT:?set TRAIT, e.g. TRAIT=HDL}"
USER_BIM="${BIM:-}"
. "$(dirname "$0")/env.sh"
# ancestry -> GLGC file suffix and MAGMA 1000G panel (default EUR = unchanged)
ANC="$(echo "${ANC:-EUR}" | tr '[:lower:]' '[:upper:]')"
a="$(echo "$ANC" | tr '[:upper:]' '[:lower:]')"
case "$ANC" in
  EUR) SUFFIX="HRC_1KGP3_others_ALL"; PANEL=eur ;;
  AFR) SUFFIX="HRC_1KGP3_others_ALL"; PANEL=afr ;;
  SAS) SUFFIX="HRC_1KGP3_others_ALL"; PANEL=sas ;;
  EAS) SUFFIX="1KGP3_ALL";            PANEL=eas ;;
  HIS) SUFFIX="1KGP3_ALL";            PANEL=amr ;;   # Hispanic -> 1000G AMR panel
  *) echo "ANC must be EUR, AFR, EAS, SAS or HIS (got $ANC)" >&2; exit 1 ;;
esac
[ -n "$USER_BIM" ] || BIM="$DATA/g1000_${PANEL}.bim"

t="$(echo "$TRAIT" | tr '[:upper:]' '[:lower:]')"
TDIR="$DATA/$t"                      # mirrors the existing $DATA/ldl/
IN="${IN:-$TDIR/${TRAIT}_INV_${ANC}_${SUFFIX}.meta.singlevar.results.gz}"
PVAL="$TDIR/${t}_${a}.magma.pval"
NLOG="$TDIR/${t}_${a}.neglog10.txt"
echo "trait=$TRAIT ancestry=$ANC panel=g1000_$PANEL"
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
