#!/bin/bash
#BSUB -J "redogwas[16-54]%50"
#BSUB -n 1
#BSUB -W 5000
##BSUB -q sara
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=4GB]"
#BSUB -o stdout.%J.%I
#BSUB -e stderr.%J.%I
set -euo pipefail
PIPE_DIR=/share/maize/ntanduk/CATFISH/simulation/permutation_pipeline
source "${PIPE_CONFIG:-$PIPE_DIR/config.sh}"

LINE=$(sed -n "${LSB_JOBINDEX}p" "$PIPE_DIR/missing_gwas.txt")
B=$(echo "$LINE" | awk '{print $1}'); CHR=$(echo "$LINE" | awk '{print $2}')
echo "[redogwas] perm=$B chr=$CHR"

PHENO="${PHENO_DIR}/pheno_perm_${B}.csv"
[ -s "$PHENO" ] || { echo "ERROR: no phenotype $PHENO" >&2; exit 1; }
VCF="${VCF_DIR}/$(echo "$VCF_PREFIX_PATTERN" | sed "s/@CHR@/${CHR}/g")"
[ -s "$VCF" ] || { echo "ERROR: no VCF $VCF" >&2; exit 1; }

# --- private per-job tmp: copy the VCF here so its .csi index is fully isolated.
#     This folder is the ONLY thing this script ever deletes (its own tmp).
TMP="${OUT_ROOT}/tmp_gwas/perm_${B}_chr_${CHR}_${LSB_JOBID}"
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT
cp "$VCF" "$TMP/g.vcf.gz"
cd "$TMP"

"$VCF2GWAS" -v "$TMP/g.vcf.gz" -pf "$PHENO" -ap -cf 'PCA' -c 3 -lmm 1 -np -nl -T 1

# --- keep ONLY the assoc result; never remove or overwrite existing output ---
WORKDIR="${GWAS_DIR}/perm_${B}/chr${CHR}"
mkdir -p "$WORKDIR"
mapfile -t ASSOC_FILES < <(find "$TMP/Output" -name '*.assoc.txt' -print | sort)
[ "${#ASSOC_FILES[@]}" -gt 0 ] || { echo "ERROR: no assoc produced" >&2; exit 1; }
[ "${#ASSOC_FILES[@]}" -eq 1 ] || {
  echo "ERROR: expected exactly 1 assoc file, found ${#ASSOC_FILES[@]}" >&2
  printf '  %s\n' "${ASSOC_FILES[@]}" >&2
  exit 1
}
cp "${ASSOC_FILES[0]}" "$WORKDIR/assoc.txt"
echo "[redogwas] assoc -> $WORKDIR/assoc.txt"
