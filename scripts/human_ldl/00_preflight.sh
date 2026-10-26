#!/bin/bash
# Check everything is in place before submitting 22 jobs.  bash 00_preflight.sh
ROOT=/share/maize/ntanduk/CATFISH
MAGMA=/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/magma_v1.10/magma
ok=0

chk () { if [ -s "$2" ]; then printf "  OK    %-22s %s\n" "$1" "$(du -h "$2"|cut -f1)"; else printf "  FAIL  %-22s MISSING: %s\n" "$1" "$2"; ok=1; fi; }

echo "MAGMA binary:"
if [ -x "$MAGMA" ]; then echo "  OK    $MAGMA"; else echo "  FAIL  not executable: $MAGMA"; ok=1; fi

echo "inputs:"
chk "reference .bed"  "$ROOT/g1000_eur.bed"
chk "reference .bim"  "$ROOT/g1000_eur.bim"
chk "reference .fam"  "$ROOT/g1000_eur.fam"
chk "gene annotation" "$ROOT/ldl_eur.genes.annot"

PVAL=$ROOT/ldl/ldl_eur.magma.pval
[ -f "$PVAL" ] || PVAL=$ROOT/ldl_eur.magma.pval
chk "p-value file"    "$PVAL"

if [ -s "$PVAL" ]; then
  rows=$(wc -l < "$PVAL")
  echo "  pval rows: $rows  (expect 44,664,015)"
  [ "$rows" -lt 40000000 ] && { echo "  FAIL  pval file looks truncated"; ok=1; }
  [ "$(tail -1 "$PVAL" | awk '{print NF}')" -ne 3 ] && { echo "  FAIL  last pval line truncated"; ok=1; }
  head -1 "$PVAL" | grep -q '^SNP' || { echo "  FAIL  pval header is not 'SNP P N'"; ok=1; }
fi

if [ -s "$ROOT/ldl_eur.genes.annot" ]; then
  echo "  annotated genes: $(grep -vc '^#' "$ROOT/ldl_eur.genes.annot")"
fi

echo
[ "$ok" -eq 0 ] && echo "ALL CHECKS PASSED -- bsub < 04_magma_array_LSF.sh" \
               || echo "FIX THE ABOVE BEFORE SUBMITTING"
exit $ok
