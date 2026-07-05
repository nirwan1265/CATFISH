#!/bin/bash
#BSUB -J spike_eval
#BSUB -n 1
#BSUB -W 120
##BSUB -q sara
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=8GB]"
#BSUB -o logs/spike_eval.%J.out
#BSUB -e logs/spike_eval.%J.err

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SUBMIT_DIR="${LS_SUBCWD:-${LSB_SUB_PWD:-$PWD}}"
PIPE_DIR="$SCRIPT_DIR"
if [[ ! -f "$PIPE_DIR/config.spikein.sh" && -f "$SUBMIT_DIR/config.spikein.sh" ]]; then
  PIPE_DIR="$SUBMIT_DIR"
fi
CONFIG_PATH="${PIPE_CONFIG:-$PIPE_DIR/config.spikein.sh}"
source "$CONFIG_PATH"
mkdir -p "$PIPE_DIR/logs"

"${RSCRIPT:-Rscript}" "$PIPE_DIR/09_spikein_eval.R" "$CATFISH_OUT" "$PHENO_DIR" "$SPIKE_SUMMARY_DIR"
