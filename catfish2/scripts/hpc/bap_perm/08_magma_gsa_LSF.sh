#!/bin/bash
# =============================================================================
# Step 8: MAGMA's own gene-set analysis (competitive, plus --model
# self-contained) for each permutation, on the same SNP-wise mean gene
# analysis that CATFISH used.
#   MERGE_MODE=renamed (default)  merge the existing per-chromosome runs,
#                                 linked under batch names (only after 07 says
#                                 identical)
#   MERGE_MODE=batch              merge chromosome-batch reruns (BATCH=1 in 03)
# Needs INPUT_DIR/sorghumbicolorcyc.sets.annot (08a_write_set_annot.R).
#   bsub < 08_magma_gsa_LSF.sh
# =============================================================================
#BSUB -J "bap_magma_gsa[1-100]%25"
#BSUB -q sara
#BSUB -n 1
#BSUB -W 120
#BSUB -o logs/magma_gsa.%J.%I.out
#BSUB -e logs/magma_gsa.%J.%I.err
set -euo pipefail
source /share/maize/ntanduk/CATFISH/simulation/bap_perm_wy/config.sh
B=$LSB_JOBINDEX; P=$(printf "perm_%03d" "$B")
TAG=$(echo "$GENE_MODEL" | sed 's/[^A-Za-z0-9]\{1,\}/_/g')
SETS="${INPUT_DIR}/sorghumbicolorcyc.sets.annot"
GSA="${RUN_ROOT}/magma_gsa"; mkdir -p "$GSA"
[ -s "$SETS" ] || { echo "missing $SETS: run 08a_write_set_annot.R" >&2; exit 1; }
case "${MERGE_MODE:-renamed}" in
  renamed)
    MDIR="${RUN_ROOT}/magma_merged/${P}"; mkdir -p "$MDIR"
    for c in $CHROMS; do
      k=$((10#$c))
      for ext in genes.out genes.raw log log.suppl; do
        src="${MAGMA_OUT}/${P}/chr${c}/${P}_chr${c}.${TAG}.${ext}"
        [ -e "$src" ] && ln -sf "$src" "${MDIR}/${P}.${TAG}.batch${k}_chr.${ext}"
      done
      [ -e "${MDIR}/${P}.${TAG}.batch${k}_chr.genes.raw" ] || { echo "missing chr $c raw for $P" >&2; exit 1; }
    done ;;
  batch) MDIR="${MAGMA_OUT}_batch/${P}" ;;
  *) echo "MERGE_MODE must be renamed or batch" >&2; exit 1 ;;
esac
"$MAGMA" --merge "${MDIR}/${P}.${TAG}" --out "${MDIR}/${P}.${TAG}"
[ -s "${MDIR}/${P}.${TAG}.genes.raw" ] || { echo "merge produced no genes.raw for $P" >&2; exit 1; }
"$MAGMA" --gene-results "${MDIR}/${P}.${TAG}.genes.raw" \
         --set-annot "$SETS" \
         --model self-contained \
         --out "${GSA}/${P}"
[ -s "${GSA}/${P}.gsa.out" ] || { echo "no gsa.out for $P; see ${GSA}/${P}.log" >&2; exit 1; }
echo "[magma gsa] $P -> ${GSA}/${P}.gsa.out"
