#!/bin/bash
# Step 04: gene-gene correlations, once (04_corpairs_once.R). After step 03b.
#   bsub < 04_corpairs_once_LSF.sh
#BSUB -J sim_corpairs
##BSUB -q sara
#BSUB -n 1
#BSUB -W 120
#BSUB -o logs/corpairs.%J.out
#BSUB -e logs/corpairs.%J.err
set -euo pipefail
source /share/maize/ntanduk/CATFISH/simulation_check_catfish/config.sh
Rscript "${SIM_ROOT}/04_corpairs_once.R"
