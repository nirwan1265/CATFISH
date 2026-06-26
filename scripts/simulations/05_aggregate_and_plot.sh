#!/bin/bash
#BSUB -n 1
#BSUB -W 60
#BSUB -q sara
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=8GB]"
#BSUB -J catfish_aggregate
#BSUB -o stdout.%J
#BSUB -e stderr.%J


set -euo pipefail
PIPE_DIR=/share/maize/ntanduk/CATFISH/simulation/permutation_pipeline
source "${PIPE_CONFIG:-$PIPE_DIR/config.sh}"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

Rscript "$PIPE_DIR/05_aggregate_and_plot.R" "$CATFISH_OUT" "$FINAL_DIR"
