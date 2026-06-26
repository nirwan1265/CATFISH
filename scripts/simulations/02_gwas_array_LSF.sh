#!/bin/bash
# =============================================================================
# 02_gwas_array_LSF.sh  --  Permuted/synthetic GWAS via vcf2gwas (LSF array),
#                           matching the REAL command:
#     vcf2gwas -v <chr.vcf.gz> -pf <pheno.csv> -ap -cf 'PCA' -c 3 -lmm 1 -np -nl -T 3
#
# One array element = one (phenotype index, chromosome). vcf2gwas handles
# relatedness + PCA internally, so there is NO separate kinship step.
# Submit:  bsub < 02_gwas_array_LSF.sh   (export PIPE_CONFIG first for spike-in)
# Array size = N_PHENO * 10  (e.g. 100 -> [1-1000]).
# =============================================================================

#BSUB -J "gwas[1-1000]%80"
#BSUB -n 1
#BSUB -W 5000
##BSUB -q sara
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=5GB]"
#BSUB -o logs/gwas.%J.%I.out
#BSUB -e logs/gwas.%J.%I.err



set -euo pipefail
WORKING_DIR=/share/maize/ntanduk/CATFISH/simulation/permutation_pipeline
source "$WORKING_DIR/config.sh"


#set -euo pipefail
#source "${PIPE_CONFIG:-$(dirname "$0")/config.sh}"
mkdir -p logs

CHR_ARR=($CHROMS); NCHR=${#CHR_ARR[@]}
TID=$((LSB_JOBINDEX - 1))
B=$(( TID / NCHR + 1 ))
CHR=${CHR_ARR[$(( TID % NCHR ))]}
BPAD=$(printf "%03d" "$B")
echo "[gwas] pheno=$B chr=$CHR (LSB_JOBINDEX=$LSB_JOBINDEX)"

VCF="${VCF_DIR}/$(echo "$VCF_PREFIX_PATTERN" | sed "s/@CHR@/${CHR}/g")"
PHENO="${PHENO_DIR}/pheno_perm_${BPAD}.csv"
WORKDIR="${GWAS_DIR}/perm_${BPAD}/chr${CHR}"
mkdir -p "$WORKDIR"
cd "$WORKDIR"     # vcf2gwas writes results into ./Output here

# --- avoid reusing stale vcf2gwas outputs from prior attempts ----------------
rm -rf Output

# --- run vcf2gwas exactly as the real pipeline did ---------------------------
"$VCF2GWAS" -v "$VCF" -pf "$PHENO" -ap -cf 'PCA' -c 3 -lmm 1 -np -nl -T 1

# --- locate the GEMMA assoc output deterministically -------------------------
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
echo "[gwas] assoc -> ${WORKDIR}/assoc.txt  (from ${ASSOC_FILES[0]})"
