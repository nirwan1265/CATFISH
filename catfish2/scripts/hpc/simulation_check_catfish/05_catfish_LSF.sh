#!/bin/bash
# Step 05: CATFISH on every replicate (05_catfish.R). After steps 03b and 04.
#   pilot: bsub -J "sim_catfish[1-20]" < 05_catfish_LSF.sh
#   full : MODE=full bsub -J "sim_catfish[1-850]" < 05_catfish_LSF.sh
#BSUB -J "sim_catfish[1-20]"
##BSUB -q sara
#BSUB -n 1
#BSUB -W 480
#BSUB -o logs/catfish.%J.%I.out
#BSUB -e logs/catfish.%J.%I.err
set -euo pipefail
source /share/maize/ntanduk/CATFISH/simulation_check_catfish/config.sh
Rscript "${SIM_ROOT}/05_catfish.R" "$LSB_JOBINDEX"
