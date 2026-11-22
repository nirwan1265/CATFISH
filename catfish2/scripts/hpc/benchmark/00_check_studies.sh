#!/usr/bin/env bash
# =============================================================================
# 00_check_studies.sh -- before downloading anything, check each candidate GWAS
# Catalog study for the gene-prioritisation benchmark: author, year, trait,
# sample size, whether full summary statistics exist, and the file names in its
# summary-statistics folder (harmonised/ included). Light: a few small web
# requests per study, no downloads.
#   bash 00_check_studies.sh                 # the 10 studies in traits.tsv
#   bash 00_check_studies.sh GCST006862 ...  # extra candidate accessions
# =============================================================================
set -u
DIR="$(cd "$(dirname "$0")" && pwd)"
if [ "$#" -gt 0 ]; then IDS="$*"; else IDS=$(awk -F'\t' 'NR > 1 { print $2 }' "$DIR/traits.tsv"); fi
API=https://www.ebi.ac.uk/gwas/rest/api/studies
FTP=https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics

range_of() {   # GCST006414 -> GCST006001-GCST007000 ; GCST90128471 -> GCST90128001-GCST90129000
  local id=$1 num w lo hi
  num=${id#GCST}; w=${#num}; num=$((10#$num))
  lo=$(( (num - 1) / 1000 * 1000 + 1 )); hi=$(( lo + 999 ))
  printf "GCST%0${w}d-GCST%0${w}d" "$lo" "$hi"
}

for id in $IDS; do
  echo "================ $id ================"
  curl -s --max-time 60 -H "Accept: application/json" "$API/$id" | python3 -c '
import json, sys
try:
    d = json.load(sys.stdin)
except Exception as e:
    print("  API: no JSON (", e, ")"); sys.exit()
p = d.get("publicationInfo", {}) or {}
a = (p.get("author") or {}).get("fullname", "?")
print("  publication :", a, p.get("publicationDate", "?"), "|", (p.get("title") or "")[:90])
print("  trait       :", (d.get("diseaseTrait") or {}).get("trait", "?"))
print("  sample      :", d.get("initialSampleSize", "?"))
print("  full p set  :", d.get("fullPvalueSet", "?"))
'
  url="$FTP/$(range_of "$id")/$id"
  for sub in "" "harmonised/"; do
    echo "  files in $url/$sub"
    curl -s --max-time 60 "$url/$sub" | grep -o 'href="[^"?/][^"]*"' | sed 's/href="//; s/"$//' | sed 's/^/    /'
  done
  sleep 2
done
