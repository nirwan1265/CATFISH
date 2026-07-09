#!/bin/bash
set -euo pipefail
WORKING_DIR=/share/maize/ntanduk/CATFISH/simulation/sap_stem_volume_permutation
source "$WORKING_DIR/config.sap_stem_volume.sh"

MAXPERM=${1:-$N_PERM}
tag=$(echo "$GENE_MODEL" | sed 's/[^A-Za-z0-9]\+/_/g')
OUT="$WORKING_DIR/missing_magma_jobs.txt"
: > "$OUT"
complete=0

for b in $(seq 1 "$MAXPERM"); do
  bp=$(printf %03d "$b")
  miss=0
  for chr in $CHROMS; do
    f="$MAGMA_OUT/perm_${bp}/chr${chr}/perm_${bp}_chr${chr}.${tag}.genes.out"
    if [ ! -s "$f" ]; then
      echo "$bp $chr" >> "$OUT"
      miss=1
    fi
  done
  [ "$miss" -eq 0 ] && complete=$((complete+1))
done

K=$(wc -l < "$OUT")
echo "MAGMA-complete permutations: $complete / $MAXPERM"
echo "missing/incomplete MAGMA jobs: $K -> $OUT"
echo "First few missing jobs:"
head -n 20 "$OUT" || true
echo ">>> set the recover array header to [1-$K]"

