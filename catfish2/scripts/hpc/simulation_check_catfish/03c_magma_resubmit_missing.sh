#!/bin/bash
# =============================================================================
# Step 03c: resubmit every MAGMA chromosome job (step 03) whose .genes.raw is
# missing, for all replicates of the design. Run it on the login node (it only
# submits jobs), with the R env active:
#   MODE=full bash 03c_magma_resubmit_missing.sh
# Array element i = replicate (i - 1) / 10 + 1, chromosome (i - 1) % 10 + 1.
# Submits in groups (index lists under 200 characters), at most 200 running per group, and waits
# while more than 800 of your jobs are pending (the cluster's pending limit).
# Run it again until it reports 0 missing.
# =============================================================================
set -euo pipefail
source /share/maize/ntanduk/CATFISH/simulation_check_catfish/config.sh
cd "$SIM_ROOT"
NREP=$(( $(wc -l < "${RUN}/design.tsv") - 1 ))
TODO="${RUN}/magma_todo.txt"
: > "$TODO"
for r in $(seq 1 "$NREP"); do
  for c in $CHROMS; do
    f=$(printf "%s/magma/rep_%04d/rep_%04d.batch%d_chr.genes.raw" "$RUN" "$r" "$r" "$c")
    [ -s "$f" ] || echo $(( (r - 1) * 10 + c )) >> "$TODO"
  done
done
N=$(wc -l < "$TODO")
echo "[03c] MODE=$MODE: $NREP replicates, $N chromosome jobs missing"
[ "$N" -gt 0 ] || { echo "[03c] nothing to submit; run the merge (03b)"; exit 0; }
## LSF rejects long job names ("Bad job name"), so consecutive elements are
## written as ranges (2001-3000) and each submission's index list is kept under
## 200 characters
IDX_GROUPS=$(awk '
  function flush() { if (s == "") return; r = (s == e) ? s : s "-" e;
                     if (length(cur) + length(r) + 1 > 200) { print cur; cur = "" }
                     cur = (cur == "") ? r : cur "," r }
  NR == 1 { s = $1; e = $1; next }
  $1 == e + 1 { e = $1; next }
  { flush(); s = $1; e = $1 }
  END { flush(); if (cur != "") print cur }' "$TODO")
for IDX in $IDX_GROUPS; do
  while [ "$(bjobs 2>/dev/null | grep -c PEND || true)" -gt 800 ]; do sleep 60; done
  echo "[03c] submitting elements ${IDX:0:70}..."
  bsub -J "sim_magma[${IDX}]%200" < 03_magma_gene_LSF.sh
done
echo "[03c] all groups submitted"
