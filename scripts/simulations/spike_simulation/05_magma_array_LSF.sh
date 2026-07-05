#!/bin/bash
#BSUB -J "magma[1-1000]%50"
#BSUB -n 1
#BSUB -W 5000
##BSUB -q sara
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=20GB]"
#BSUB -o logs/magma.%J.%I.out
#BSUB -e logs/magma.%J.%I.err

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

CHR_ARR=($CHROMS)
NCHR=${#CHR_ARR[@]}
TID=$((LSB_JOBINDEX - 1))
B=$(( TID / NCHR + 1 ))
CHR=${CHR_ARR[$(( TID % NCHR ))]}

echo "[magma] spike=$B chr=$CHR (LSB_JOBINDEX=$LSB_JOBINDEX)"
Rscript "$PIPE_DIR/05_magma_catfish.R" "$B" "$CHR"
