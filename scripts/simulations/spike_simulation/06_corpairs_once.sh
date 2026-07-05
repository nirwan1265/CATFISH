#!/bin/bash

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

Rscript "$PIPE_DIR/06_corpairs_once.R"
