#!/bin/bash
#BSUB -J "gwas[1-1000]%80"
#BSUB -n 1
#BSUB -W 5000
##BSUB -q sara
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=5GB]"
#BSUB -o logs/gwas.%J.%I.out
#BSUB -e logs/gwas.%J.%I.err

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SUBMIT_DIR="${LS_SUBCWD:-${LSB_SUB_PWD:-$PWD}}"
PIPE_DIR="$SCRIPT_DIR"
if [[ ! -f "$PIPE_DIR/config.spikein.sh" && -f "$SUBMIT_DIR/config.spikein.sh" ]]; then
  PIPE_DIR="$SUBMIT_DIR"
fi
CONFIG_PATH="${PIPE_CONFIG:-$PIPE_DIR/config.spikein.sh}"
source "$CONFIG_PATH"
mkdir -p "$PIPE_DIR/logs"

CHR_ARR=($CHROMS)
NCHR=${#CHR_ARR[@]}
TID=$((LSB_JOBINDEX - 1))
B=$(( TID / NCHR + 1 ))
CHR=${CHR_ARR[$(( TID % NCHR ))]}
BPAD=$(printf "%03d" "$B")

echo "[gwas] spike=$B chr=$CHR (LSB_JOBINDEX=$LSB_JOBINDEX)"

VCF_SRC="${VCF_DIR}/$(echo "$VCF_PREFIX_PATTERN" | sed "s/@CHR@/${CHR}/g")"
PHENO="${PHENO_DIR}/pheno_perm_${BPAD}.csv"
WORKDIR="${GWAS_DIR}/perm_${BPAD}/chr${CHR}"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

rm -rf Output

TMPROOT="${TMPDIR:-/tmp}"
JOB_TMP=$(mktemp -d "${TMPROOT}/spike_gwas_${BPAD}_chr${CHR}_XXXXXX")
cleanup() {
  rm -rf "$JOB_TMP"
}
trap cleanup EXIT

VCF_LOCAL="${JOB_TMP}/$(basename "$VCF_SRC")"
cp "$VCF_SRC" "$VCF_LOCAL"
if [[ -f "${VCF_SRC}.csi" ]]; then
  cp "${VCF_SRC}.csi" "${VCF_LOCAL}.csi"
fi
if [[ -f "${VCF_SRC}.tbi" ]]; then
  cp "${VCF_SRC}.tbi" "${VCF_LOCAL}.tbi"
fi

echo "[gwas] using local VCF copy: $VCF_LOCAL"
"$VCF2GWAS" -v "$VCF_LOCAL" -pf "$PHENO" -ap -cf 'PCA' -c 3 -lmm 1 -np -nl -T 1

mapfile -t ASSOC_FILES < <(find Output -name '*.assoc.txt' -print | sort)
if [[ "${#ASSOC_FILES[@]}" -eq 0 ]]; then
  echo "[gwas] ERROR: no .assoc.txt produced by vcf2gwas in $WORKDIR/Output" >&2
  exit 1
fi
if [[ "${#ASSOC_FILES[@]}" -ne 1 ]]; then
  echo "[gwas] ERROR: expected exactly 1 .assoc.txt, found ${#ASSOC_FILES[@]}" >&2
  printf '  %s\n' "${ASSOC_FILES[@]}" >&2
  exit 1
fi

cp "${ASSOC_FILES[0]}" "${WORKDIR}/assoc.txt"
echo "[gwas] assoc -> ${WORKDIR}/assoc.txt"
