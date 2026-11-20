#!/bin/bash
# Leave-one-gene-out pathway evidence for the lipid positive controls, one job
# per lipid analysis (same inputs and settings as wy_gc_real_LSF.sh).
#   bsub < wy_logo_LSF.sh
#BSUB -J "wy_logo[1-10]"
##BSUB -q sara
#BSUB -n 1
#BSUB -W 720
#BSUB -o logs/wy_logo.%J.%I.out
#BSUB -e logs/wy_logo.%J.%I.err
set -euo pipefail
ROOT=/share/maize/ntanduk/CATFISH/gc_real
A=(ldl_eur hdl_eur ldl_afr hdl_afr ldl_eas hdl_eas ldl_sas hdl_sas ldl_his hdl_his)
cd "$ROOT"
GC_ROOT="$ROOT" Rscript "$ROOT/wy_logo.R" "${A[$((LSB_JOBINDEX - 1))]}"
