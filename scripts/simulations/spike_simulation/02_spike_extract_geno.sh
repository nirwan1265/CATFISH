#!/bin/bash
#BSUB -J spike_extract
#BSUB -n 1
#BSUB -W 14400
##BSUB -q sara
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=20GB]"
#BSUB -o logs/spike_extract.%J.out
#BSUB -e logs/spike_extract.%J.err

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

SNPLIST="${SPIKE_WORK}/causal_snps.txt"
TMP="${SPIKE_WORK}/extract_tmp"
mkdir -p "$TMP"
: > "${TMP}/merge_list.txt"
FIRST=""

for CHR in $CHROMS; do
  GENO_PFX="${GENO_DIR}/$(echo "$GENO_PREFIX_PATTERN" | sed "s/@CHR@/${CHR}/g")"
  OUT="${TMP}/causal_chr${CHR}"
  if "$PLINK" --bfile "$GENO_PFX" --extract "$SNPLIST" --make-bed \
              --out "$OUT" --allow-no-sex --allow-extra-chr 2>/dev/null; then
    if [[ -f "${OUT}.bed" ]]; then
      if [[ -z "$FIRST" ]]; then
        FIRST="$OUT"
      else
        echo "$OUT" >> "${TMP}/merge_list.txt"
      fi
    fi
  fi
done

if [[ -z "$FIRST" ]]; then
  echo "[spike] ERROR: no chromosomes produced extracted causal SNPs" >&2
  exit 1
fi

MERGED="${SPIKE_WORK}/causal_geno"
if [[ -s "${TMP}/merge_list.txt" ]]; then
  "$PLINK" --bfile "$FIRST" --merge-list "${TMP}/merge_list.txt" \
           --make-bed --out "$MERGED" --allow-no-sex --allow-extra-chr
else
  cp "${FIRST}.bed" "${MERGED}.bed"
  cp "${FIRST}.bim" "${MERGED}.bim"
  cp "${FIRST}.fam" "${MERGED}.fam"
fi

"$PLINK" --bfile "$MERGED" --recode A --out "${SPIKE_WORK}/causal_geno" \
         --allow-no-sex --allow-extra-chr
echo "[spike] dosage matrix -> ${SPIKE_WORK}/causal_geno.raw"
