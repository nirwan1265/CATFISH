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

BFILE=$ROOT/g1000_eur
# the annotation is trait-independent: it maps genes to reference-panel
# variants and never sees the GWAS, so the LDL one is reused as is
ANNOT=$ROOT/g1000_eur.genes.annot
[ -f "$ANNOT" ] || ANNOT=$ROOT/ldl_eur.genes.annot
# tolerate any of the locations the munge step has used
for cand in "$ROOT/lipids/${TRAIT}_eur.magma.pval" \
            "$ROOT/$TRAIT/${TRAIT}_eur.magma.pval" \
            "$ROOT/${TRAIT}_eur.magma.pval"; do
  [ -f "$cand" ] && PVAL="$cand" && break
done
if [ -z "${PVAL:-}" ]; then
  echo "[magma] no .magma.pval found for trait=$TRAIT. Looked in:" >&2
  echo "  $ROOT/lipids/ , $ROOT/$TRAIT/ , $ROOT/" >&2
  exit 1
fi

OUTDIR=$ROOT/magma_$TRAIT
PREFIX=$OUTDIR/${TRAIT}_eur

mkdir -p "$OUTDIR" logs

CHR=$LSB_JOBINDEX
echo "[magma] trait=$TRAIT chr=$CHR  pval=$PVAL"

# ncol=N uses GLGC's per-variant sample size column.
"$MAGMA" --bfile      "$BFILE" \
         --pval       "$PVAL" ncol=N \
         --gene-annot "$ANNOT" \
         --batch      "$CHR" chr \
         --out        "$PREFIX"

echo "[magma] trait=$TRAIT chr $CHR done"
