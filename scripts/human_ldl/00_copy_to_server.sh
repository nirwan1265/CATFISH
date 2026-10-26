#!/usr/bin/env bash
# =============================================================================
# 00_copy_to_server.sh -- stage the human LDL MAGMA inputs onto the HPC.
# Run on the LAPTOP, after 01_munge_glgc_ldl.sh has finished.
#   bash scripts/human_ldl/00_copy_to_server.sh
# =============================================================================
set -euo pipefail

LOCAL_DATA="$HOME/Documents/Github/MAGCAT_data"
LOCAL_REPO="$HOME/Documents/Github/MAGCAT"

SERVER="${SERVER:-ntanduk@login.hpc.ncsu.edu}"          # <-- set your login
REMOTE="${REMOTE:-/share/maize/ntanduk/CATFISH}"

PVAL="$LOCAL_DATA/ldl/ldl_eur.magma.pval"

# --- refuse to ship a truncated pval file ------------------------------------
if [ ! -s "$PVAL" ]; then
  echo "ERROR: $PVAL missing. Run 01_munge_glgc_ldl.sh first." >&2; exit 1
fi
rows=$(wc -l < "$PVAL")
if [ "$rows" -lt 15000000 ]; then
  echo "ERROR: $PVAL has only $rows rows -- expected ~22M." >&2
  echo "       This looks like a partial run. Re-run 01_munge_glgc_ldl.sh." >&2
  exit 1
fi
if [ "$(tail -1 "$PVAL" | awk '{print NF}')" -ne 3 ]; then
  echo "ERROR: last line of $PVAL is truncated. Re-run 01_munge_glgc_ldl.sh." >&2
  exit 1
fi
echo "[check] pval OK: $rows rows"

ssh "$SERVER" "mkdir -p '$REMOTE/data' '$REMOTE/scripts' '$REMOTE/logs' '$REMOTE/magma'"

echo "[copy] reference panel + gene locations + pathway sets (~3.6 GB)"
rsync -avh --progress \
  "$LOCAL_DATA/g1000_eur.bed" \
  "$LOCAL_DATA/g1000_eur.bim" \
  "$LOCAL_DATA/g1000_eur.fam" \
  "$LOCAL_DATA/g1000_eur.synonyms" \
  "$LOCAL_DATA/NCBI37.3.gene.loc" \
  "$LOCAL_DATA/c2.cp.v2026.2.Hs.entrez.sets.annot" \
  "$LOCAL_DATA/humancyc.entrez.sets.annot" \
  "$PVAL" \
  "$SERVER:$REMOTE/data/"

echo "[copy] LSF scripts"
rsync -avh --progress \
  "$LOCAL_REPO/scripts/human_ldl/config_human_ldl.sh" \
  "$LOCAL_REPO/scripts/human_ldl/03_magma_annot_LSF.sh" \
  "$LOCAL_REPO/scripts/human_ldl/04_magma_array_LSF.sh" \
  "$LOCAL_REPO/scripts/human_ldl/05_magma_merge_LSF.sh" \
  "$SERVER:$REMOTE/scripts/"

echo
echo "done. On the server:"
echo "  cd $REMOTE"
echo "  bsub < scripts/03_magma_annot_LSF.sh"
