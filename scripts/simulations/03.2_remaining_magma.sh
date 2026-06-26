#!/bin/bash
set -euo pipefail
PIPE_DIR=/share/maize/ntanduk/CATFISH/simulation/permutation_pipeline
source "$PIPE_DIR/config.sh"
MAXPERM=${1:-84}                                  # we only did 1..84
tag=$(echo "$GENE_MODEL" | sed 's/[^A-Za-z0-9]\+/_/g')   # -> multi_snp_wise
OUT="$PIPE_DIR/missing_jobs.txt"; : > "$OUT"
complete=0
for b in $(seq 1 "$MAXPERM"); do
  bp=$(printf %03d "$b"); miss=0
  for chr in $CHROMS; do
    f="$MAGMA_OUT/perm_${bp}/chr${chr}/perm_${bp}_chr${chr}.${tag}.genes.out"
    if [ ! -s "$f" ]; then echo "$bp $chr" >> "$OUT"; miss=1; fi
  done
  [ "$miss" -eq 0 ] && complete=$((complete+1))
done
K=$(wc -l < "$OUT")
echo "complete permutations: $complete / $MAXPERM"
echo "missing (perm,chr) jobs: $K   -> $OUT"
echo ">>> set the recover array header to [1-$K]"
