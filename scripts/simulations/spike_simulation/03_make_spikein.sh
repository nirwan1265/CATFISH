#!/bin/bash
#BSUB -J make_spikein_stage12
#BSUB -n 1
#BSUB -W 240
##BSUB -q sara
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=16GB]"
#BSUB -o logs/make_spikein_stage12.%J.out
#BSUB -e logs/make_spikein_stage12.%J.err

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

"${RSCRIPT:-Rscript}" "$PIPE_DIR/03_make_spikein.R"
