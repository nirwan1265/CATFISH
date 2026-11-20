#!/bin/bash
# =============================================================================
# Step 07c: two-stage ACAT gene-set test (07c_acat.R), one array element = one
# design cell (all its replicates). Needs only steps 01-02 (the simulated z).
#   full:  MODE=full bsub -J "sim_acat[1-17]" < 07c_acat_LSF.sh
# =============================================================================
#BSUB -J "sim_acat[1-17]"
##BSUB -q sara
#BSUB -n 1
#BSUB -W 240
#BSUB -o logs/acat.%J.%I.out
#BSUB -e logs/acat.%J.%I.err
set -euo pipefail
source /share/maize/ntanduk/CATFISH/simulation_check_catfish/config.sh
Rscript "${SIM_ROOT}/07c_acat.R" "$LSB_JOBINDEX"
