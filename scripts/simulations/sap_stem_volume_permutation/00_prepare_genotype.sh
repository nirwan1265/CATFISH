#!/bin/bash
#BSUB -J prep_sap_genotype
#BSUB -n 1
#BSUB -W 240
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=16GB]"
#BSUB -o logs/prep_genotype.%J.out
#BSUB -e logs/prep_genotype.%J.err

set -euo pipefail
WORKING_DIR=/share/maize/ntanduk/CATFISH/simulation/sap_stem_volume_permutation
source "$WORKING_DIR/config.sap_stem_volume.sh"
mkdir -p "$WORKING_DIR/logs"

mkdir -p \
  "$PIPE_ROOT" "$INPUT_DIR" "$PREP_ROOT" "$PHENO_DIR" "$GWAS_DIR" "$MAGMA_OUT" \
  "$CATFISH_OUT" "$FINAL_DIR" "$GENO_DIR" "$ANNOT_DIR"

echo "[prep] PIPE_ROOT=$PIPE_ROOT"
echo "[prep] input dir:   $INPUT_DIR"
echo "[prep] prep dir:    $PREP_ROOT"
echo "[prep] gwas dir:    $GWAS_DIR"
echo "[prep] magma dir:   $MAGMA_OUT"
echo "[prep] catfish dir: $CATFISH_OUT"
echo "[prep] final dir:   $FINAL_DIR"
echo "[prep] config:      $WORKING_DIR/config.sap_stem_volume.sh"

"${RSCRIPT:-Rscript}" "$WORKING_DIR/00_prepare_genotype.R"
