#!/bin/bash
#BSUB -J "redosapmagma[1-1000]%100"
#BSUB -n 1
#BSUB -W 3000
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=6GB]"
#BSUB -o logs/redo_magma.%J.%I.out
#BSUB -e logs/redo_magma.%J.%I.err

set -euo pipefail
WORKING_DIR=/share/maize/ntanduk/CATFISH/simulation/sap_stem_volume_permutation
source "$WORKING_DIR/config.sap_stem_volume.sh"

LINE=$(sed -n "${LSB_JOBINDEX}p" "$WORKING_DIR/missing_magma_jobs.txt")
B=$(echo "$LINE" | awk '{print $1}')
CHR=$(echo "$LINE" | awk '{print $2}')

if [[ -z "${B:-}" || -z "${CHR:-}" ]]; then
  echo "ERROR: no entry found for LSB_JOBINDEX=$LSB_JOBINDEX in $WORKING_DIR/missing_magma_jobs.txt" >&2
  exit 1
fi

echo "[redosapmagma] perm=$B chr=$CHR"

GWAS_FILE="${GWAS_DIR}/perm_${B}/chr${CHR}/assoc.txt"
ANNOT_FILE="${ANNOT_DIR}/chr${CHR}/stem_volume_chr${CHR}.genes.annot"
BFILE_PREFIX="${GENO_DIR}/sap_chr${CHR}"

[ -s "$GWAS_FILE" ] || { echo "ERROR: no assoc.txt for perm $B chr $CHR -- rerun GWAS first" >&2; exit 1; }
[ -s "$ANNOT_FILE" ] || { echo "ERROR: missing genes.annot for chr $CHR -- rerun 00_prepare_genotype.sh first" >&2; exit 1; }
[ -s "${BFILE_PREFIX}.bed" ] || { echo "ERROR: missing PLINK bed for chr $CHR -- rerun 00_prepare_genotype.sh first" >&2; exit 1; }

"${RSCRIPT:-Rscript}" "$WORKING_DIR/03_magma_catfish.R" "$((10#$B))" "$CHR"

