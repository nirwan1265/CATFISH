#!/bin/bash
# Westfall-Young with gene-level GC for the 12 real analyses, one job each.
#   bsub < wy_gc_real_LSF.sh            (GC=0 bsub ... reruns without GC)
#BSUB -J "wy_gc_real[1-12]"
#BSUB -q sara
#BSUB -n 1
#BSUB -W 720
#BSUB -o logs/wy_gc_real.%J.%I.out
#BSUB -e logs/wy_gc_real.%J.%I.err
set -euo pipefail
ROOT=/share/maize/ntanduk/CATFISH/gc_real
A=(ldl_eur hdl_eur ldl_afr hdl_afr ldl_eas hdl_eas ldl_sas hdl_sas ldl_his hdl_his Dry_tons_per_acre Lignin_DM_perc)
cd "$ROOT"
GC_ROOT="$ROOT" Rscript "$ROOT/wy_gc_real.R" "${A[$((LSB_JOBINDEX - 1))]}"
