#!/bin/bash
# Check everything is in place before submitting 22 jobs.  bash 00_preflight.sh
ROOT=/share/maize/ntanduk/CATFISH
MAGMA=/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/magma_v1.10/magma
ok=0

chk () { if [ -s "$2" ]; then printf "  OK    %-22s %s\n" "$1" "$(du -h "$2"|cut -f1)"; else printf "  FAIL  %-22s MISSING: %s\n" "$1" "$2"; ok=1; fi; }

echo "MAGMA binary:"
if [ -x "$MAGMA" ]; then echo "  OK    $MAGMA"; else echo "  FAIL  not executable: $MAGMA"; ok=1; fi

# Same TRAIT / ANC / REF / ANNOT / PVAL as the 04 run (defaults: LDL-C EUR), e.g.
#   TRAIT=hdl ANC=afr REF=g1000_afr ANNOT=$ROOT/g1000_afr.genes.annot bash 00_preflight.sh
TRAIT=${TRAIT:-ldl}; ANC=${ANC:-eur}; REF=${REF:-g1000_eur}
ANNOT=${ANNOT:-$ROOT/ldl_eur.genes.annot}
echo "inputs (trait=$TRAIT anc=$ANC ref=$REF):"
chk "reference .bed"  "$ROOT/$REF.bed"
chk "reference .bim"  "$ROOT/$REF.bim"
chk "reference .fam"  "$ROOT/$REF.fam"
chk "gene annotation" "$ANNOT"

if [ -z "${PVAL:-}" ]; then
  if [ "${TRAIT}_${ANC}" = "ldl_eur" ]; then
    PVAL=$ROOT/ldl/ldl_eur.magma.pval
    [ -f "$PVAL" ] || PVAL=$ROOT/ldl_eur.magma.pval
  else
    PVAL=$ROOT/lipids/${TRAIT}_${ANC}.magma.pval
  fi
fi
chk "p-value file"    "$PVAL"

if [ -s "$PVAL" ]; then
  rows=$(wc -l < "$PVAL")
  if [ "${TRAIT}_${ANC}" = "ldl_eur" ]; then
    echo "  pval rows: $rows  (expect 44,664,015)"
    [ "$rows" -lt 40000000 ] && { echo "  FAIL  pval file looks truncated"; ok=1; }
  else
    echo "  pval rows: $rows  (no fixed expectation for this trait/ancestry; compare with the munge log's kept= count)"
  fi
  [ "$(tail -1 "$PVAL" | awk '{print NF}')" -ne 3 ] && { echo "  FAIL  last pval line truncated"; ok=1; }
  head -1 "$PVAL" | grep -q '^SNP' || { echo "  FAIL  pval header is not 'SNP P N'"; ok=1; }
fi

if [ -s "$ANNOT" ]; then
  echo "  annotated genes: $(grep -vc '^#' "$ANNOT")"
fi

echo
[ "$ok" -eq 0 ] && echo "ALL CHECKS PASSED -- bsub < 04_magma_array_LSF.sh" \
               || echo "FIX THE ABOVE BEFORE SUBMITTING"
exit $ok
