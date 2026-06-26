#!/bin/bash
set -euo pipefail
PIPE_DIR=/share/maize/ntanduk/CATFISH/simulation/permutation_pipeline
source "$PIPE_DIR/config.sh"
MAXPERM=${1:-84}
OUT="$PIPE_DIR/missing_gwas.txt"; : > "$OUT"
complete=0
for b in $(seq 1 "$MAXPERM"); do
  bp=$(printf %03d "$b"); miss=0
  for chr in $CHROMS; do
    [ -s "$GWAS_DIR/perm_${bp}/chr${chr}/assoc.txt" ] || { echo "$bp $chr" >> "$OUT"; miss=1; }
  done
  [ "$miss" -eq 0 ] && complete=$((complete+1))
done
K=$(wc -l < "$OUT")
echo "GWAS-complete permutations: $complete / $MAXPERM"
echo "missing GWAS (perm,chr): $K   -> $OUT"
echo ">>> set the gwas recover array to [1-$K]"
