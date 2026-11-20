#!/bin/bash
# =============================================================================
# Step 01: munge the PoPS feature matrix once (munge_feature_directory.py).
# Reads every file in FEATURE_RAW_DIR (only PoPS.features.txt.gz may be there;
# pandas reads .gz directly), aligns rows to GENE_ANNOT, and writes
#   ${FEATURE_PREFIX}.mat.<i>.npy, .cols.<i>.txt, .rows.txt
# The whole feature file is loaded at once as float64 (18,383 genes x ~57,000
# features is ~8.5 GB before parsing overhead), so this job needs a lot of
# memory; no memory request is set here.
#   bsub < 01_munge_features_LSF.sh
# =============================================================================
#BSUB -J pops_munge
##BSUB -q sara
#BSUB -n 1
#BSUB -W 720
#BSUB -o logs/pops_munge.%J.out
#BSUB -e logs/pops_munge.%J.err
set -euo pipefail
source /share/maize/ntanduk/CATFISH/pops/scripts/config.sh
n_files=$(find "$FEATURE_RAW_DIR" -maxdepth 1 -type f | wc -l)
[ "$n_files" -ge 1 ] || { echo "no feature files in $FEATURE_RAW_DIR" >&2; exit 1; }
echo "[munge] $n_files file(s) in $FEATURE_RAW_DIR -> ${FEATURE_PREFIX}"
python "${POPS_DIR}/munge_feature_directory.py" \
  --gene_annot_path "$GENE_ANNOT" \
  --feature_dir "$FEATURE_RAW_DIR" \
  --save_prefix "$FEATURE_PREFIX" \
  --max_cols "$MAX_COLS"
echo "[munge] chunks written: $(ls ${FEATURE_PREFIX}.mat.*.npy | wc -l)"
