#!/bin/bash
# Step 08: power and off-target rates (08_aggregate.R). After steps 05-07.
#   bsub < 08_aggregate_LSF.sh            (MODE=full bsub ... for the full run)
#BSUB -J sim_aggregate
##BSUB -q sara
#BSUB -n 1
#BSUB -W 120
#BSUB -o logs/aggregate.%J.out
#BSUB -e logs/aggregate.%J.err
set -euo pipefail
source /share/maize/ntanduk/CATFISH/simulation_check_catfish/config.sh
Rscript "${SIM_ROOT}/08_aggregate.R"
