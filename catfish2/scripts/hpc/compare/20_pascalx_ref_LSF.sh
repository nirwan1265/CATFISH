#!/bin/bash
# =============================================================================
# Step 20: the 1000 Genomes reference panels MAGMA used for the lipids
# (g1000_<panel>), written per chromosome as .tped.gz and imported into
# PascalX's database once, so that the scoring jobs (21) only read them.
# PascalX imports PLINK data from .tped (its VCF import requires a QUAL value
# that PLINK-written VCFs do not have). One array element per panel.
#   bsub < 20_pascalx_ref_LSF.sh
# =============================================================================
#BSUB -J "cmp_pascalx_ref[1-5]"
#BSUB -q sara
#BSUB -n 1
#BSUB -W 720
#BSUB -o logs/pascalx_ref.%J.%I.out
#BSUB -e logs/pascalx_ref.%J.%I.err
set -eo pipefail
source /share/maize/ntanduk/CATFISH/compare/config.sh
source activate "$PASCALX_ENV" 2>/dev/null || conda activate "$PASCALX_ENV"
export LD_LIBRARY_PATH="${PASCALX_DIR}/build/lib:${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
PANELS=(eur afr eas sas amr); P=${PANELS[$((LSB_JOBINDEX - 1))]}
BFILE="${LIPID_ROOT}/g1000_${P}"
[ -s "${BFILE}.bed" ] || { echo "missing ${BFILE}.bed" >&2; exit 1; }
command -v "$PLINK" >/dev/null || { echo "plink not found (set PLINK=/path/to/plink)" >&2; exit 1; }
for c in $(seq 1 22); do
  out="${PASCALX_REF}/g1000_${P}.chr${c}"
  [ -s "${out}.tped.gz" ] && continue
  "$PLINK" --bfile "$BFILE" --chr "$c" --recode transpose 12 --out "$out"
  gzip -f "${out}.tped"
done
python - <<PY
from PascalX import refpanel
R = refpanel.refpanel()
R.set_refpanel(filename="${PASCALX_REF}/g1000_${P}", parallel=1, chrlist=list(range(1, 23)))
print("imported g1000_${P}")
PY
