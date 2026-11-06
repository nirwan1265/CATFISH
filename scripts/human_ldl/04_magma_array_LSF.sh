#!/bin/bash
# =============================================================================
# 04_magma_array_LSF.sh -- MAGMA gene analysis for human LDL-C (GLGC 2021 EUR).
#   One array element = one autosome, via MAGMA's native --batch N chr.
#   Self-contained: no config file to source.
#
#   Run 00_preflight.sh first, then:  bsub < 04_magma_array_LSF.sh
#   When all 22 finish:               bsub < 05_magma_merge_LSF.sh
# =============================================================================
#BSUB -J "ldl_magma[1-22]%11"
#BSUB -n 1
#BSUB -W 5000
##BSUB -q sara
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=20GB]"
#BSUB -o logs/magma.%J.%I.out
#BSUB -e logs/magma.%J.%I.err

set -euo pipefail

ROOT=/share/maize/ntanduk/CATFISH
MAGMA=/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/magma_v1.10/magma

# Parameters (env, passed through bsub; defaults reproduce the LDL-C EUR run):
#   TRAIT ldl|hdl  ANC eur|afr|eas|sas|his  REF g1000_<panel>
#   ANNOT, PVAL, OUTDIR  explicit paths
# e.g.  TRAIT=hdl ANC=afr REF=g1000_afr ANNOT=$ROOT/g1000_afr.genes.annot \
#         OUTDIR=$ROOT/magma_hdl_afr bsub -J "magma_hdl_afr[1-22]%11" < 04_magma_array_LSF.sh
TRAIT=${TRAIT:-ldl}
ANC=${ANC:-eur}
REF=${REF:-g1000_eur}
BFILE=$ROOT/$REF
ANNOT=${ANNOT:-$ROOT/ldl_eur.genes.annot}
if [ -z "${PVAL:-}" ]; then
  if [ "${TRAIT}_${ANC}" = "ldl_eur" ]; then
    PVAL=$ROOT/ldl/ldl_eur.magma.pval
    [ -f "$PVAL" ] || PVAL=$ROOT/ldl_eur.magma.pval      # tolerate either location
  else
    PVAL=$ROOT/lipids/${TRAIT}_${ANC}.magma.pval
  fi
fi

OUTDIR=${OUTDIR:-$ROOT/magma}
PREFIX=$OUTDIR/${TRAIT}_${ANC}

mkdir -p "$OUTDIR" logs

CHR=$LSB_JOBINDEX
echo "[magma] chr=$CHR  trait=$TRAIT  anc=$ANC  ref=$BFILE  annot=$ANNOT  pval=$PVAL"

# ncol=N uses GLGC's per-variant sample size column.
"$MAGMA" --bfile      "$BFILE" \
         --pval       "$PVAL" ncol=N \
         --gene-annot "$ANNOT" \
         --batch      "$CHR" chr \
         --out        "$PREFIX"

echo "[magma] chr $CHR done"
