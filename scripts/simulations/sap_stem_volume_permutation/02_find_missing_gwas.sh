#!/bin/bash
set -euo pipefail
WORKING_DIR=/share/maize/ntanduk/CATFISH/simulation/sap_stem_volume_permutation
source "$WORKING_DIR/config.sap_stem_volume.sh"

MAXPERM=${1:-$N_PERM}
OUT="$WORKING_DIR/missing_gwas.txt"
: > "$OUT"
complete=0

for b in $(seq 1 "$MAXPERM"); do
  bp=$(printf %03d "$b")
  miss=0
  for chr in $CHROMS; do
    f="$GWAS_DIR/perm_${bp}/chr${chr}/assoc.txt"
    if [ ! -s "$f" ]; then
      echo "$bp $chr" >> "$OUT"
      miss=1
    fi
  done
  [ "$miss" -eq 0 ] && complete=$((complete+1))
done

K=$(wc -l < "$OUT")
echo "GWAS-complete permutations: $complete / $MAXPERM"
echo "missing GWAS (perm,chr): $K -> $OUT"
echo "First few missing:"
head -n 20 "$OUT" || true

