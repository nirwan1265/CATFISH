#!/usr/bin/env bash
# =============================================================================
# 12_fetch_glgc.sh -- download GLGC 2021 ancestry-specific (EUR) summary stats
# for the lipid traits used in the multi-trait ablation.
#
# Source: https://csg.sph.umich.edu/willer/public/glgc-lipids2021/
#         results/ancestry_specific/
#
# All five traits share one filename pattern, one column layout and one
# reference panel, so the entire downstream pipeline is unchanged per trait.
#
# SIZE WARNING: each file is roughly 2.3 GB compressed. Four traits is ~9 GB
# plus the LDL file you already have.
#
# RUN THIS WHERE THE PIPELINE RUNS -- normally the cluster, since the files are
# ~2.3 GB each and MAGMA reads them there. Paths come from env.sh.
#
#   bash 12_fetch_glgc.sh            # HDL logTG TC
#   TRAITS="HDL logTG TC nonHDL" bash 12_fetch_glgc.sh
# =============================================================================
set -euo pipefail

. "$(dirname "$0")/env.sh"

BASE="https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/ancestry_specific"
TRAITS="${TRAITS:-HDL logTG TC}"


# The cluster login node may not have outbound HTTPS. Check before downloading
# 7 GB and failing, and say what to do instead.
if ! curl -fsSI --max-time 20 "$BASE/" >/dev/null 2>&1; then
  echo "cannot reach $BASE from this host." >&2
  echo "  Either run this on a host with outbound HTTPS and copy the files to" >&2
  echo "  $OUT, or use a data-transfer node if the cluster has one." >&2
  exit 1
fi

for T in $TRAITS; do
  t="$(echo "$T" | tr '[:upper:]' '[:lower:]')"
  OUT="$DATA/$t"                    # mirrors the existing $DATA/ldl/
  mkdir -p "$OUT"
  F="${T}_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"
  if [ -s "$OUT/$F" ]; then
    echo "already have $F"; continue
  fi
  echo "downloading $F ..."
  curl -fL --retry 3 --retry-delay 5 -C - -o "$OUT/$F.part" "$BASE/$F"
  mv "$OUT/$F.part" "$OUT/$F"
  echo "  $(ls -lh "$OUT/$F" | awk '{print $5}')"
done

echo
echo "done."
for T in $TRAITS; do
  t="$(echo "$T" | tr '[:upper:]' '[:lower:]')"; ls -lh "$DATA/$t" 2>/dev/null
done
