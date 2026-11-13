#!/bin/bash
# =============================================================================
# Step 00b: the simulation panel written per chromosome as .tped.gz and imported
# into PascalX's database, as compare/20_pascalx_ref_LSF.sh did for the lipids.
# Run after step 00.
#   bsub < 00b_pascalx_ref_LSF.sh
# =============================================================================
#BSUB -J sim_pascalx_ref
##BSUB -q sara
#BSUB -n 1
#BSUB -W 720
#BSUB -o logs/pascalx_ref.%J.out
#BSUB -e logs/pascalx_ref.%J.err
set -eo pipefail
source /share/maize/ntanduk/CATFISH/simulation_check_catfish/config.sh
export LD_LIBRARY_PATH="${PASCALX_DIR}/build/lib:${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
[ -s "${SIM_PANEL}.bed" ] || { echo "missing ${SIM_PANEL}.bed: run step 00" >&2; exit 1; }
command -v "$PLINK" >/dev/null || { echo "plink not found (set PLINK=/path/to/plink)" >&2; exit 1; }
mkdir -p "$PASCALX_REF"
for c in $CHROMS; do
  out="${PASCALX_REF}/sim_eur.chr${c}"
  [ -s "${out}.tped.gz" ] && continue
  "$PLINK" --bfile "$SIM_PANEL" --chr "$c" --recode transpose 12 --out "$out"
  gzip -f "${out}.tped"
done
python - <<PY
from PascalX import refpanel
R = refpanel.refpanel()
R.set_refpanel(filename="${PASCALX_REF}/sim_eur", parallel=1, chrlist=[${CHROMS// /, }])
print("imported sim_eur, chromosomes ${CHROMS}")
PY
