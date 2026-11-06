#!/bin/bash
# One (permutation, chromosome) GWAS with vcf2gwas, using the genome-wide
# kinship (-k) and PCs (-cf) from 00_prepare_genome_LSF.sh: nothing is
# recomputed per job, and the kinship does not come from the chromosome being
# tested. Runs on a private copy of the chromosome VCF: vcf2gwas writes and
# deletes a .csi next to the VCF, which breaks other jobs reading the same file.
#   bash 02_run_gwas_single.sh <perm> <chr, e.g. 01>
set -euo pipefail
source /share/maize/ntanduk/CATFISH/simulation/bap_perm_wy/config.sh
B="$1"; CHR="$2"; BPAD=$(printf "%03d" "$B")
PHENO="${PHENO_DIR}/pheno_perm_${BPAD}.csv"
VCF_SRC="$(echo "$CHR_VCF_PATTERN" | sed "s/@CHR@/${CHR}/g")"
OUTDIR="${GWAS_DIR}/perm_${BPAD}/chr${CHR}"; ASSOC="${OUTDIR}/assoc.txt"
for f in "$PHENO" "$VCF_SRC" "$KINSHIP" "$PCS_CSV"; do [ -f "$f" ] || { echo "missing $f" >&2; exit 1; }; done
mkdir -p "$OUTDIR"; rm -rf "$OUTDIR/Output"; rm -f "$ASSOC"
JOB_TMP=$(mktemp -d "${TMPDIR:-/tmp}/bap_gwas_${BPAD}_chr${CHR}_XXXXXX")
trap 'rm -rf "$JOB_TMP"' EXIT
VCF_LOCAL="${JOB_TMP}/$(basename "$VCF_SRC")"
cp "$VCF_SRC" "$VCF_LOCAL"
cd "$OUTDIR"
echo "[gwas] perm=$B chr=$CHR"
"$VCF2GWAS" -v "$VCF_LOCAL" -pf "$PHENO" -ap -k "$KINSHIP" -cf "$PCS_CSV" -ac -lmm 1 -np -nl -T 1
mapfile -t F < <(find Output -name '*.assoc.txt' -print | sort)
[ "${#F[@]}" -eq 1 ] || { echo "[gwas] expected 1 .assoc.txt, found ${#F[@]}" >&2; printf '  %s\n' "${F[@]}" >&2; exit 1; }
cp "${F[0]}" "$ASSOC"
rm -rf "$OUTDIR/Output"
echo "[gwas] -> $ASSOC ($(($(wc -l < "$ASSOC") - 1)) SNPs)"
