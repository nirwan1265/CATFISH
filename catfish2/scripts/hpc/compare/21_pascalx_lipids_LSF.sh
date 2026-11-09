#!/bin/bash
# =============================================================================
# Step 21: PascalX gene scoring and pathway scoring for the ten lipid sets,
# on the same SNP p-values MAGMA analysed, the matching 1000 Genomes panel,
# NCBI37.3 gene locations and the same Reactome sets (Entrez IDs). Writes
# out/pascalx/<ta>.tsv (pathway, p). Needs step 20 for that set's panel.
#   bsub < 21_pascalx_lipids_LSF.sh
# =============================================================================
#BSUB -J "cmp_pascalx_lipids[1-10]"
#BSUB -q sara
#BSUB -n 1
#BSUB -W 2880
#BSUB -o logs/pascalx_lipids.%J.%I.out
#BSUB -e logs/pascalx_lipids.%J.%I.err
set -eo pipefail
source /share/maize/ntanduk/CATFISH/compare/config.sh
source activate "$PASCALX_ENV" 2>/dev/null || conda activate "$PASCALX_ENV"
export LD_LIBRARY_PATH="${PASCALX_DIR}/build/lib:${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
A=($LIPIDS); TA=${A[$((LSB_JOBINDEX - 1))]}; t=${TA%_*}; a=${TA#*_}
case "$a" in eur) P=eur ;; afr) P=afr ;; eas) P=eas ;; sas) P=sas ;; his) P=amr ;; esac
## gene file: NCBI37.3 autosomal genes, with the Entrez ID also as the symbol so
## that the Entrez-ID Reactome sets match (written once, atomically)
if [ ! -s "$PASCALX_GENOME" ]; then
  [ -s "$GENE_LOC_37" ] || { echo "missing $GENE_LOC_37" >&2; exit 1; }
  tmp="${PASCALX_GENOME}.$$"
  awk -F'\t' 'BEGIN{OFS="\t"} $2 ~ /^([1-9]|1[0-9]|2[0-2])$/ {print $1,$2,$3,$4,$5,$1}' "$GENE_LOC_37" > "$tmp"
  mv -n "$tmp" "$PASCALX_GENOME"; rm -f "$tmp"
fi
PVAL=""
for cand in "${LIPID_ROOT}/lipids/${TA}.magma.pval" "${LIPID_ROOT}/${t}/${TA}.magma.pval" "${LIPID_ROOT}/${TA}.magma.pval"; do
  [ -s "$cand" ] && PVAL="$cand" && break
done
[ -n "$PVAL" ] || { echo "no ${TA}.magma.pval under ${LIPID_ROOT}" >&2; exit 1; }
[ -s "${PASCALX_REF}/g1000_${P}.chr1.db" ] || { echo "panel g1000_${P} not imported: run step 20" >&2; exit 1; }
echo "[pascalx] $TA  panel=$P  pval=$PVAL  window=$PASCALX_WINDOW  maf=$PASCALX_MAF"
"${PASCALX_DIR}/pascalx" -p 1 -n True -m "$PASCALX_MAF" -w "$PASCALX_WINDOW" \
  -pw "$REACTOME_SETS" -cn 0 -cs 1 \
  "$PASCALX_GENOME" "${PASCALX_REF}/g1000_${P}" "${OUT}/pascalx/${TA}" \
  genescoring -sh True -cr 0 -cp 1 "$PVAL"
[ -s "${OUT}/pascalx/${TA}.tsv" ] || { echo "no pathway output for $TA" >&2; exit 1; }
echo "[pascalx] $TA -> ${OUT}/pascalx/${TA}.tsv"
