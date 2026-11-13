#!/bin/bash
# =============================================================================
# Step 02: GWAS statistics for every replicate of one design cell: plink2
# --variant-score gives each panel SNP's sum of genotype x u over the 503
# people for every replicate column of RUN/u/u_cell<k>.tsv, and
# 02_simulate_z.R standardises it to z (and checks it against the exact z of
# the causal SNPs from step 01). One array element per cell.
#   pilot: bsub -J "sim_z[1-4]" < 02_simulate_z_LSF.sh
#   full : MODE=full bsub -J "sim_z[1-17]" < 02_simulate_z_LSF.sh
# =============================================================================
#BSUB -J "sim_z[1-4]"
##BSUB -q sara
#BSUB -n 1
#BSUB -W 240
#BSUB -o logs/simz.%J.%I.out
#BSUB -e logs/simz.%J.%I.err
set -euo pipefail
source /share/maize/ntanduk/CATFISH/simulation_check_catfish/config.sh
K=$(printf "%02d" "$LSB_JOBINDEX")
U="${RUN}/u/u_cell${K}.tsv"
[ -s "$U" ] || { echo "missing $U: run step 01" >&2; exit 1; }
mkdir -p "${RUN}/z"
"$PLINK2" --bfile "$SIM_PANEL" --variant-score "$U" \
          --out "${RUN}/z/vscore_cell${K}"
Rscript "${SIM_ROOT}/02_simulate_z.R" "$K"
