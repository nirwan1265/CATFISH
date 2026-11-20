#!/bin/bash
# =============================================================================
# config.sh -- PoPS (Weeks et al., Polygenic Priority Score) for the gene-
# prioritisation benchmark. Sourced by every step in this folder.
# Steps:
#   01_munge_features_LSF.sh   feature matrix -> numpy chunks (once)
#   02_magma_annot_ensg_LSF.sh MAGMA annotation in Ensembl IDs (once per panel)
#   then MAGMA gene analysis with ../14 and ../15 (ANNOT_OVERRIDE, OUT_SUFFIX=_ensg)
#   03_pops_LSF.sh             PoPS, one job per trait
# The conda env (POPS_ENV) is activated by you before submitting.
# =============================================================================
export POPS_DIR=/usr/local/usrapps/maize/ntanduk/pops                    # git clone of FinucaneLab/pops
export POPS_ENV=/usr/local/usrapps/maize/ntanduk/pops_env                # python env for PoPS
export POPS_DATA=/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/pops_data
export GENE_ANNOT=${POPS_DATA}/gene_annot_jun10.txt                       # ENSGID NAME CHR START END TSS (18,383 genes)
export CONTROL_FEATURES=${POPS_DATA}/features_jul17_control.txt
export FEATURE_RAW_DIR=${POPS_DATA}/features_raw                          # only PoPS.features.txt.gz in here
export FEATURE_PREFIX=${POPS_DATA}/features_munged/pops_features
export MAX_COLS=5000                                                      # columns per munged chunk (script default)

export ROOT=/share/maize/ntanduk/CATFISH                                  # panel, MAGMA runs (as ../14, ../15)
export MAGMA=/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/magma_v1.10/magma
export WINDOW="10,10"                                                     # kb up/down, as the paper's MAGMA annotation
export POPS_OUT=${ROOT}/pops                                              # PoPS results
mkdir -p "$POPS_OUT" "${POPS_DATA}/features_munged" logs 2>/dev/null || true
