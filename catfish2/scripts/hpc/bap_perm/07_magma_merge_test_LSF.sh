#!/bin/bash
# =============================================================================
# Step 7: can the existing per-chromosome MAGMA runs be merged for MAGMA's own
# gene-set analysis? MAGMA documents --merge only for chromosome-batch runs
# (--batch <k> chr). Permutation 1 is therefore run both ways:
#   batch     BATCH=1 bsub -J "bap_magma[1-10]" < 03_magma_array_LSF.sh   (run first)
#   renamed   the existing per-chromosome outputs, linked under batch names
# and the two merged .genes.raw / .genes.out are compared. If they are
# identical, 08 merges every permutation from the existing runs (MERGE_MODE=renamed);
# if not, every permutation has to be rerun with BATCH=1 (MERGE_MODE=batch).
#   bsub < 07_magma_merge_test_LSF.sh
# =============================================================================
#BSUB -J bap_magma_merge_test
#BSUB -q sara
#BSUB -n 1
#BSUB -W 60
#BSUB -o logs/magma_merge_test.%J.out
#BSUB -e logs/magma_merge_test.%J.err
set -uo pipefail
source /share/maize/ntanduk/CATFISH/simulation/bap_perm_wy/config.sh
P=perm_001
TAG=$(echo "$GENE_MODEL" | sed 's/[^A-Za-z0-9]\{1,\}/_/g')
BDIR="${MAGMA_OUT}_batch/${P}"; RDIR="${RUN_ROOT}/magma_merge_test/${P}_renamed"
mkdir -p "$RDIR"
for c in $CHROMS; do
  k=$((10#$c))
  [ -s "${BDIR}/${P}.${TAG}.batch${k}_chr.genes.raw" ] || { echo "batch run missing for chr $c: run BATCH=1 first" >&2; exit 1; }
  for ext in genes.out genes.raw log log.suppl; do
    src="${MAGMA_OUT}/${P}/chr${c}/${P}_chr${c}.${TAG}.${ext}"
    [ -e "$src" ] && ln -sf "$src" "${RDIR}/${P}.${TAG}.batch${k}_chr.${ext}"
  done
done
echo "== merge, batch runs";   "$MAGMA" --merge "${BDIR}/${P}.${TAG}" --out "${BDIR}/${P}.${TAG}"; echo "exit $?"
echo "== merge, renamed runs"; "$MAGMA" --merge "${RDIR}/${P}.${TAG}" --out "${RDIR}/${P}.${TAG}"; echo "exit $?"

same=1
echo "== per chromosome: batch run vs existing run (data lines only)"
for c in $CHROMS; do
  k=$((10#$c))
  for ext in genes.out genes.raw; do
    a="${BDIR}/${P}.${TAG}.batch${k}_chr.${ext}"; b="${MAGMA_OUT}/${P}/chr${c}/${P}_chr${c}.${TAG}.${ext}"
    if cmp -s <(grep -v '^#' "$a") <(grep -v '^#' "$b"); then echo "chr$c $ext identical"; else echo "chr$c $ext DIFFERENT"; same=0; fi
  done
done
echo "== merged: batch vs renamed"
for ext in genes.out genes.raw; do
  a="${BDIR}/${P}.${TAG}.${ext}"; b="${RDIR}/${P}.${TAG}.${ext}"
  if [ ! -s "$a" ] || [ ! -s "$b" ]; then echo "merged $ext MISSING (batch: $(ls -la "$a" 2>&1 | wc -l), renamed: $(ls -la "$b" 2>&1 | wc -l))"; same=0
  elif cmp -s <(grep -v '^#' "$a") <(grep -v '^#' "$b"); then echo "merged $ext identical ($(grep -vc '^#' "$a") lines)"
  else echo "merged $ext DIFFERENT"; same=0; fi
done
echo "== header lines of the merged .genes.raw (batch, then renamed)"
grep '^#' "${BDIR}/${P}.${TAG}.genes.raw" 2>/dev/null; echo "--"; grep '^#' "${RDIR}/${P}.${TAG}.genes.raw" 2>/dev/null
if [ "$same" = "1" ]; then echo "RESULT: identical -> use MERGE_MODE=renamed in 08"; else echo "RESULT: not identical -> do not merge the renamed runs; send this log"; fi
