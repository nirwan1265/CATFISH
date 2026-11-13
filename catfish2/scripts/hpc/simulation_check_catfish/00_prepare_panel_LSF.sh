#!/bin/bash
# =============================================================================
# Step 00: the simulation panel and its gene and pathway definitions, once.
#   - g1000_eur (the lipid panel), chromosomes 1-10, HapMap3 SNPs, MAF >= 0.01
#     (a SNP with no variation in the 503 people cannot be standardised)
#   - genotype counts of every panel SNP (to standardise genotypes in step 02)
#   - MAGMA's EUR annotation (+/- 10 kb) restricted to the panel's SNPs
#   - Reactome sets restricted to genes with panel SNPs, 10-200 genes
#   bsub < 00_prepare_panel_LSF.sh
# =============================================================================
#BSUB -J sim_prepare
##BSUB -q sara
#BSUB -n 1
#BSUB -W 240
#BSUB -o logs/prepare.%J.out
#BSUB -e logs/prepare.%J.err
set -euo pipefail
source /share/maize/ntanduk/CATFISH/simulation_check_catfish/config.sh
[ -s "$HM3" ] || { echo "missing $HM3 (rsync w_hm3.snplist from the Mac)" >&2; exit 1; }
[ -s "${PANEL_FULL}.bed" ] || { echo "missing ${PANEL_FULL}.bed" >&2; exit 1; }
[ -s "$ANNOT_FULL" ] || ANNOT_FULL="${LIPID_ROOT}/ldl_eur.genes.annot"
[ -s "$ANNOT_FULL" ] || { echo "no EUR gene annotation in ${LIPID_ROOT}" >&2; exit 1; }
mkdir -p "$(dirname "$SIM_PANEL")"

awk 'NR > 1 {print $1}' "$HM3" > "${IN}/hm3.ids"
"$PLINK2" --bfile "$PANEL_FULL" --chr 1-10 --extract "${IN}/hm3.ids" --maf 0.01 \
          --make-bed --out "$SIM_PANEL"
"$PLINK2" --bfile "$SIM_PANEL" --geno-counts --out "$SIM_PANEL"
echo "[prepare] panel: $(wc -l < "${SIM_PANEL}.bim") SNPs, $(wc -l < "${SIM_PANEL}.fam") people"

ANNOT_FULL="$ANNOT_FULL" Rscript "${SIM_ROOT}/00_prepare_sets.R"
