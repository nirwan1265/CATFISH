#!/bin/bash
# =============================================================================
# Step 03: PoPS for one trait, default (recommended) settings: ridge, MAGMA
# covariates and error covariance, feature selection p < 0.05, HLA removed,
# control features always included, seed 42.
# MAGMA input: the Ensembl run of ../14 and ../15 (magma_<trait>[_<anc>]_ensg/).
#   TRAIT=ldl bsub < 03_pops_LSF.sh                  (EUR)
#   TRAIT=ldl ANC=afr bsub -J pops_ldl_afr < 03_pops_LSF.sh
# PoPS loads the whole munged feature matrix into memory; no memory request is
# set here.
# Writes ${POPS_OUT}/<trait>_<anc>.preds (PoPS_Score per ENSGID), .coefs, .marginals
# =============================================================================
#BSUB -J pops
##BSUB -q sara
#BSUB -n 1
#BSUB -W 720
#BSUB -o logs/pops.%J.out
#BSUB -e logs/pops.%J.err
set -euo pipefail
source /share/maize/ntanduk/CATFISH/pops/scripts/config.sh
TRAIT="${TRAIT:?set TRAIT, e.g. TRAIT=ldl}"
a="$(echo "${ANC:-eur}" | tr '[:upper:]' '[:lower:]')"
if [ "$a" = "eur" ]; then MP=${ROOT}/magma_${TRAIT}_ensg/${TRAIT}_eur
else MP=${ROOT}/magma_${TRAIT}_${a}_ensg/${TRAIT}_${a}; fi
for f in "${MP}.genes.out" "${MP}.genes.raw"; do
  [ -s "$f" ] || { echo "[pops] missing $f (run ../14 and ../15 with OUT_SUFFIX=_ensg)" >&2; exit 1; }
done
NCH=$(ls ${FEATURE_PREFIX}.mat.*.npy 2>/dev/null | wc -l)
[ "$NCH" -ge 1 ] || { echo "[pops] no munged features at ${FEATURE_PREFIX} (run step 01)" >&2; exit 1; }
echo "[pops] trait=$TRAIT anc=$a magma=$MP chunks=$NCH"
python "${POPS_DIR}/pops.py" \
  --gene_annot_path "$GENE_ANNOT" \
  --feature_mat_prefix "$FEATURE_PREFIX" \
  --num_feature_chunks "$NCH" \
  --magma_prefix "$MP" \
  --control_features_path "$CONTROL_FEATURES" \
  --out_prefix "${POPS_OUT}/${TRAIT}_${a}" \
  --verbose
echo "[pops] wrote ${POPS_OUT}/${TRAIT}_${a}.preds"
