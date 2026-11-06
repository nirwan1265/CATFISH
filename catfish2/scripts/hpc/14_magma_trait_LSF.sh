#!/bin/bash
# =============================================================================
# 14_magma_trait_LSF.sh -- MAGMA gene analysis for one GLGC lipid trait.
#   Identical to 04_magma_array_LSF.sh except that the trait is a variable.
#   Every MAGMA flag and every LSF resource request is unchanged, so the other
#   lipid traits are analysed exactly as LDL-C was and the results are
#   comparable. Do not add gene settings here that the LDL run did not use.
#
#   One array element = one autosome, via MAGMA's native --batch N chr.
#
#   TRAIT=hdl bsub < 14_magma_trait_LSF.sh
#   When all 22 finish:  TRAIT=hdl bsub < 15_magma_merge_trait_LSF.sh
#
#   Other ancestries (default EUR, unchanged): same flags, that ancestry's
#   1000G panel and its own annotation g1000_<panel>.genes.annot
#   (window 10,10 as for EUR); output in magma_<trait>_<anc>/.
#     ANC=AFR TRAIT=ldl bsub -J "magma_ldl_afr[1-22]%22" < 14_magma_trait_LSF.sh
# =============================================================================
#BSUB -J "lipid_magma[1-22]%22"
#BSUB -n 1
#BSUB -W 5000
##BSUB -q sara
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=20GB]"
#BSUB -o logs/magma.%J.%I.out
#BSUB -e logs/magma.%J.%I.err

set -euo pipefail

TRAIT="${TRAIT:?set TRAIT, e.g. TRAIT=hdl}"

ROOT=/share/maize/ntanduk/CATFISH
MAGMA=/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/magma_v1.10/magma
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

BFILE=$ROOT/g1000_$PANEL
if [ "$ANC" = "EUR" ]; then
  # the annotation is trait-independent: it maps genes to reference-panel
  # variants and never sees the GWAS, so the LDL one is reused as is
  ANNOT=$ROOT/g1000_eur.genes.annot
  [ -f "$ANNOT" ] || ANNOT=$ROOT/ldl_eur.genes.annot
else
  ANNOT=$ROOT/g1000_$PANEL.genes.annot
  [ -f "$ANNOT" ] || { echo "[magma] missing $ANNOT -- annotate the $PANEL panel first" >&2; exit 1; }
fi
# tolerate any of the locations the munge step has used
for cand in "$ROOT/lipids/${TRAIT}_${a}.magma.pval" \
            "$ROOT/$TRAIT/${TRAIT}_${a}.magma.pval" \
            "$ROOT/${TRAIT}_${a}.magma.pval"; do
  [ -f "$cand" ] && PVAL="$cand" && break
done
if [ -z "${PVAL:-}" ]; then
  echo "[magma] no .magma.pval found for trait=$TRAIT. Looked in:" >&2
  echo "  $ROOT/lipids/ , $ROOT/$TRAIT/ , $ROOT/" >&2
  exit 1
fi

if [ "$ANC" = "EUR" ]; then OUTDIR=$ROOT/magma_$TRAIT; else OUTDIR=$ROOT/magma_${TRAIT}_$a; fi
PREFIX=$OUTDIR/${TRAIT}_${a}

mkdir -p "$OUTDIR" logs

CHR=$LSB_JOBINDEX
echo "[magma] trait=$TRAIT anc=$ANC chr=$CHR  panel=$BFILE  annot=$ANNOT  pval=$PVAL"

# ncol=N uses GLGC's per-variant sample size column.
"$MAGMA" --bfile      "$BFILE" \
         --pval       "$PVAL" ncol=N \
         --gene-annot "$ANNOT" \
         --batch      "$CHR" chr \
         --out        "$PREFIX"

echo "[magma] trait=$TRAIT chr $CHR done"
