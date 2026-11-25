#!/bin/bash
# Leave-one-gene-out pathway evidence for the positive controls, one job per
# analysis (same inputs and settings as wy_gc_real_LSF.sh). SETS (space-
# separated) picks the analyses; default the ten lipid sets. CONTROLS_FILE and
# LOGO_OUT pass through to wy_logo.R. Submit with the R env that has CATFISH.
#   bsub < wy_logo_LSF.sh
#   benchmark traits (gold-standard genes):
#   SETS="cad_eur brca_eur cd_eur uc_eur af_eur sle_eur asthma_eur" \
#     CONTROLS_FILE=/share/maize/ntanduk/CATFISH/gc_real/inputs/45_gold_genes.csv \
#     LOGO_OUT=/share/maize/ntanduk/CATFISH/gc_real/out_logo_gold \
#     bsub -J "wy_logo_bench[1-7]" < wy_logo_LSF.sh
#BSUB -J "wy_logo[1-10]"
##BSUB -q sara
#BSUB -n 1
#BSUB -W 720
#BSUB -o logs/wy_logo.%J.%I.out
#BSUB -e logs/wy_logo.%J.%I.err
set -euo pipefail
ROOT=/share/maize/ntanduk/CATFISH/gc_real
A=(${SETS:-ldl_eur hdl_eur ldl_afr hdl_afr ldl_eas hdl_eas ldl_sas hdl_sas ldl_his hdl_his})
cd "$ROOT"
GC_ROOT="$ROOT" Rscript "$ROOT/wy_logo.R" "${A[$((LSB_JOBINDEX - 1))]}"
