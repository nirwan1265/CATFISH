#!/usr/bin/env bash
# 01_fetch_sumstats.sh -- download each trait's GWAS Catalog harmonised summary
# statistics (and its meta.yaml) listed in traits.tsv. LOGIN NODE ONLY (compute
# nodes have no internet). wget -c resumes a partial download if rerun.
#   nohup bash 01_fetch_sumstats.sh > fetch.log 2>&1 &
#   TRAITS="t2d cad" bash 01_fetch_sumstats.sh     # only some traits
set -u
DIR="$(cd "$(dirname "$0")" && pwd)"
OUT=/share/maize/ntanduk/CATFISH/benchmark/sumstats
FTP=https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics
mkdir -p "$OUT"
awk -F'\t' 'NR > 1 { print $1 "\t" $6 }' "$DIR/traits.tsv" | while IFS=$'\t' read -r t f; do
  if [ -n "${TRAITS:-}" ] && ! echo " $TRAITS " | grep -q " $t "; then continue; fi
  [ -n "$f" ] || { echo "[$t] no file in traits.tsv, skipped"; continue; }
  mkdir -p "$OUT/$t"
  echo "[$t] $(date '+%H:%M') downloading $(basename "$f")"
  wget -c -q -O "$OUT/$t/$(basename "$f")" "$FTP/$f" || echo "[$t] FAILED"
  wget -c -q -O "$OUT/$t/$(basename "$f")-meta.yaml" "$FTP/$f-meta.yaml" || echo "[$t] no meta.yaml"
  ls -lh "$OUT/$t"
done
echo "done $(date)"
