#!/bin/bash

set -euo pipefail
WORKING_DIR=/share/maize/ntanduk/CATFISH/simulation/sap_stem_volume_permutation
source "$WORKING_DIR/config.sap_stem_volume.sh"

B="$1"
CHR="$2"
BPAD=$(printf "%03d" "$B")
CHR_NUM="${CHR#0}"
PHENO="${PHENO_DIR}/pheno_perm_${BPAD}.csv"
OUTDIR="${GWAS_DIR}/perm_${BPAD}/chr${CHR}"
ASSOC="${OUTDIR}/assoc.txt"
VCF_SRC="${VCF_DIR}/$(echo "$VCF_PREFIX_PATTERN" | sed "s/@CHR@/${CHR}/g" | sed "s/@CHR_NUM@/${CHR_NUM}/g")"

mkdir -p "$OUTDIR"
echo "[gwas] perm=$B chr=$CHR"
echo "[gwas] config=$WORKING_DIR/config.sap_stem_volume.sh"
rm -rf "$OUTDIR/Output"
rm -f "$ASSOC"

TMPROOT="${TMPDIR:-/tmp}"
JOB_TMP=$(mktemp -d "${TMPROOT}/sap_gwas_${BPAD}_chr${CHR}_XXXXXX")
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

cd "$OUTDIR"
"$VCF2GWAS" -v "$VCF_LOCAL" -pf "$PHENO" -ap -cf 'PCA' -c 3 -lmm 1 -np -nl -T 1

mapfile -t ASSOC_FILES < <(find Output -name '*.assoc.txt' -print | sort)
if [[ "${#ASSOC_FILES[@]}" -eq 0 ]]; then
  echo "[gwas] ERROR: no .assoc.txt produced by vcf2gwas in $OUTDIR/Output" >&2
  exit 1
fi
if [[ "${#ASSOC_FILES[@]}" -ne 1 ]]; then
  echo "[gwas] ERROR: expected exactly 1 .assoc.txt, found ${#ASSOC_FILES[@]}" >&2
  printf '  %s\n' "${ASSOC_FILES[@]}" >&2
  exit 1
fi
cp "${ASSOC_FILES[0]}" "$ASSOC"
echo "[gwas] assoc -> $ASSOC"
