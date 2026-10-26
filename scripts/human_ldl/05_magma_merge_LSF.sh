#!/bin/bash
# =============================================================================
# 05_magma_merge_LSF.sh -- merge the 22 per-chromosome batches.
#   Produces ldl_eur.genes.out (gene p-values) and ldl_eur.genes.raw
#   (gene-gene correlations -- the input CATFISH's MVN calibration needs).
#   Submit only after all of 04 has finished:  bsub < 05_magma_merge_LSF.sh
# =============================================================================
#BSUB -J "ldl_merge"
#BSUB -n 1
#BSUB -W 240
##BSUB -q sara
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=16GB]"
#BSUB -o logs/merge.%J.out
#BSUB -e logs/merge.%J.err

set -euo pipefail

ROOT=/share/maize/ntanduk/CATFISH
MAGMA=/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/magma_v1.10/magma
PREFIX=$ROOT/magma/ldl_eur

mkdir -p logs

missing=0
for c in $(seq 1 22); do
  if ! ls "${PREFIX}.batch${c}_chr.genes.raw" >/dev/null 2>&1; then
    echo "[merge] MISSING batch $c"; missing=$((missing+1))
  fi
done
if [ "$missing" -gt 0 ]; then
  echo "[merge] $missing batch(es) missing -- rerun those array indices first." >&2
  exit 1
fi

"$MAGMA" --merge "$PREFIX" --out "$PREFIX"

echo "[merge] done"
ls -lh "${PREFIX}.genes.out" "${PREFIX}.genes.raw"
echo "[merge] genes analysed: $(( $(wc -l < "${PREFIX}.genes.out") - 1 ))"
