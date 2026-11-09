#!/bin/bash
# =============================================================================
# Step 11: MAGMA's own gene-set analysis (competitive, plus --model
# self-contained) on the real analyses: elements 1-10 the lipid sets (Reactome),
# element 11 sorghum dry biomass (SorghumBicolorCyc; merges step 10 first).
#   bsub < 11_magma_gsa_real_LSF.sh
# =============================================================================
#BSUB -J "cmp_magma_gsa[1-11]"
#BSUB -q sara
#BSUB -n 1
#BSUB -W 120
#BSUB -o logs/magma_gsa.%J.%I.out
#BSUB -e logs/magma_gsa.%J.%I.err
set -euo pipefail
source /share/maize/ntanduk/CATFISH/compare/config.sh
A=($LIPIDS Dry_tons_per_acre); TA=${A[$((LSB_JOBINDEX - 1))]}
if [ "$TA" = "Dry_tons_per_acre" ]; then
  PREFIX="${OUT}/magma_dry_mean/Dry_tons_per_acre.snp_wise_mean"
  for c in $CHROMS; do k=$((10#$c)); [ -s "${PREFIX}.batch${k}_chr.genes.raw" ] || { echo "step 10 missing chr $c" >&2; exit 1; }; done
  "$MAGMA" --merge "$PREFIX" --out "$PREFIX"
  RAW="${PREFIX}.genes.raw"; SETS="$SORGHUM_SETS"
else
  RAW="${ROOT}/inputs/lipids/${TA}.genes.raw"; SETS="$REACTOME_SETS"
fi
for f in "$RAW" "$SETS"; do [ -s "$f" ] || { echo "missing $f" >&2; exit 1; }; done
"$MAGMA" --gene-results "$RAW" --set-annot "$SETS" --model self-contained --out "${OUT}/magma_gsa/${TA}"
[ -s "${OUT}/magma_gsa/${TA}.gsa.out" ] || { echo "no gsa.out for $TA" >&2; exit 1; }
echo "[magma gsa] $TA done"
