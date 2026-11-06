#!/bin/bash
# Step 2b: MAGMA annotation, one array element per chromosome (built once;
# genotype-only, so every permutation reuses it).   bsub < 02b_annotate_LSF.sh
#BSUB -J "bap_annot[1-10]"
#BSUB -q sara
#BSUB -n 1
#BSUB -W 120
#BSUB -o logs/annot.%J.%I.out
#BSUB -e logs/annot.%J.%I.err
set -euo pipefail
source /share/maize/ntanduk/CATFISH/simulation/bap_perm_wy/config.sh
CHR_ARR=($CHROMS)
CHR=${CHR_ARR[$((LSB_JOBINDEX - 1))]}
Rscript "$WORKING_DIR/02b_annotate.R" "$CHR"
