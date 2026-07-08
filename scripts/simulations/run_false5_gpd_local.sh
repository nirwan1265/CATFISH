#!/bin/zsh
set -euo pipefail
setopt no_bg_nice

REPO="/Users/nirwantandukar/Documents/Github/MAGCAT"
OUT_ROOT="$REPO/Final/final_main_fig_false5_gpd"
PER_DIR="$OUT_ROOT/per_perm"
LOG_DIR="$OUT_ROOT/logs"
STATUS_FILE="$OUT_ROOT/status.txt"

mkdir -p "$PER_DIR" "$LOG_DIR"
print -r -- "RUNNING" > "$STATUS_FILE"

export MAGMA_OUT="/Users/nirwantandukar/Documents/Research/data/BAP/simulations/permutation_runs/magma"

perms=(1 5 10 15 20)
pids=()

for b in "${perms[@]}"; do
  log_file="$LOG_DIR/perm_$(printf "%03d" "$b").log"
  Rscript "$REPO/scripts/simulations/run_false5_gpd_local_wrapper.R" "$b" > "$log_file" 2>&1 &
  pids+=($!)
done

print -r -- "${(j: :)pids}" > "$OUT_ROOT/pids.txt"

for pid in "${pids[@]}"; do
  wait "$pid"
done

Rscript "$REPO/scripts/simulations/05_aggregate_and_plot.R" "$PER_DIR" "$OUT_ROOT" > "$LOG_DIR/aggregate.log" 2>&1
print -r -- "DONE" > "$STATUS_FILE"
