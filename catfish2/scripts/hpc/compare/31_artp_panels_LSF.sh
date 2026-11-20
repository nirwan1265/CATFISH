#!/bin/bash
# =============================================================================
# Step 31: split each 1000 Genomes panel the lipid MAGMA runs used
# (g1000_{eur,afr,eas,sas,amr}) into chromosomes, so that each sARTP call reads
# only the chromosomes its pathway touches. No SNP or sample is filtered.
#   bsub < 31_artp_panels_LSF.sh
# =============================================================================
#BSUB -J "cmp_artp_panels[1-5]"
#BSUB -q sara
#BSUB -n 1
#BSUB -W 600
#BSUB -o logs/artp_panels.%J.%I.out
#BSUB -e logs/artp_panels.%J.%I.err
set -eo pipefail
source /share/maize/ntanduk/CATFISH/compare/config.sh
P=(eur afr eas sas amr); P=${P[$((LSB_JOBINDEX - 1))]}
IN="${LIPID_ROOT}/g1000_${P}"
[ -s "${IN}.bed" ] || { echo "missing ${IN}.bed" >&2; exit 1; }
for c in $(seq 1 22); do
  out="${ARTP_REF}/g1000_${P}.chr${c}"
  [ -s "${out}.bed" ] && continue
  "$PLINK" --bfile "$IN" --chr "$c" --make-bed --out "$out" > /dev/null
done
ls "${ARTP_REF}"/g1000_${P}.chr*.bed | wc -l
