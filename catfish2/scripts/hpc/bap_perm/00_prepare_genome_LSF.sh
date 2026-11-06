#!/bin/bash
# =============================================================================
# 00_prepare_genome_LSF.sh -- built ONCE; the genotypes never change between
# permutations, so neither do the kinship matrix and the PCs.
#   1. keep.txt  = lines with a TRAIT value, in .fam order (the order
#                  01_make_permutations.R writes the phenotype files in)
#   2. per-chromosome VCFs of exactly those lines, in that order, plus their
#      concatenation
#   3. kinship   = vcf2gwas -gk on the genome-wide VCF (GEMMA centred relatedness)
#   4. 3 PCs     = plink --pca on the genome-wide VCF -> pcs.csv
# Every per-chromosome GWAS job then uses the genome-wide kinship and PCs. The
# chromosome VCFs hold only phenotyped lines in the kinship's order, so the -k
# matrix lines up row for row with every job.
# Submit with the vcf2gwas conda environment active (bcftools, plink):
#   bsub < 00_prepare_genome_LSF.sh
# =============================================================================
#BSUB -J bap_prep_genome
#BSUB -q sara
#BSUB -n 1
#BSUB -W 1440
#BSUB -o logs/prep_genome.%J.out
#BSUB -e logs/prep_genome.%J.err
set -euo pipefail
source /share/maize/ntanduk/CATFISH/simulation/bap_perm_wy/config.sh
for b in bcftools plink vcf2gwas python3; do command -v $b >/dev/null || { echo "not on PATH: $b (activate the vcf2gwas env)" >&2; exit 1; }; done
cd "$REF_DIR"

## 1. phenotyped lines, in .fam order
FAM="${GENO_DIR}/$(echo "$GENO_PREFIX_PATTERN" | sed 's/@CHR@/01/g').fam"
python3 - "$FAM" "$PHENO_CSV" "$TRAIT" keep.txt <<'PY'
import sys, csv
fam, pheno, trait, out = sys.argv[1:]
order = [l.split()[1] for l in open(fam)]
with open(pheno) as fh:
    r = csv.reader(fh); h = next(r); j = h.index(trait)
    val = {row[0]: row[j] for row in r}
keep = [i for i in order if val.get(i, "NA") not in ("NA", "", "nan", "NaN")]
open(out, "w").write("\n".join(keep) + "\n")
print(f"[prep] {len(order)} genotyped lines, {len(keep)} with {trait}")
PY

## 2. per-chromosome VCFs of those lines, in that order (used by every GWAS
##    job), and their concatenation (used only for the kinship and PCs)
LIST=""
for c in $CHROMS; do
  SRC="${VCF_DIR}/$(echo "$VCF_PREFIX_PATTERN" | sed "s/@CHR@/${c}/g")"
  DST="$(echo "$CHR_VCF_PATTERN" | sed "s/@CHR@/${c}/g")"
  echo "##contig=<ID=${c}>" > contig_${c}.hdr      # the BAP VCF headers declare no contigs
  bcftools annotate -h contig_${c}.hdr "$SRC" -Ou | bcftools view -S keep.txt -Oz -o "$DST"
  cmp -s <(bcftools query -l "$DST") keep.txt || { echo "[prep] chr$c sample order differs from keep.txt" >&2; exit 1; }
  echo "[prep] chr$c: $DST"
  LIST="$LIST $DST"
done
bcftools concat $LIST -Oz -o "$GENOME_VCF"
bcftools index -f "$GENOME_VCF"
NS=$(bcftools query -l "$GENOME_VCF" | wc -l)
cmp -s <(bcftools query -l "$GENOME_VCF") keep.txt || { echo "[prep] VCF sample order differs from keep.txt" >&2; exit 1; }
echo "[prep] $GENOME_VCF: $NS lines, $(bcftools index -n "$GENOME_VCF") variants"

## 3. kinship, once
rm -rf kin && mkdir kin && cd kin
vcf2gwas -v "$GENOME_VCF" -gk -T 1
K=$(find . -name '*.cXX.txt' | head -1)
[ -n "$K" ] || { echo "[prep] vcf2gwas -gk wrote no .cXX.txt" >&2; exit 1; }
cp "$K" "$KINSHIP"; cd "$REF_DIR"
NK=$(wc -l < "$KINSHIP")
[ "$NK" -eq "$NS" ] || { echo "[prep] kinship has $NK rows for $NS lines" >&2; exit 1; }
echo "[prep] kinship: $KINSHIP ($NK x $NK)"

## 4. PCs, once
plink --vcf "$GENOME_VCF" --double-id --allow-extra-chr --pca "$N_PC" --out pca > pca.plink.log 2>&1
awk -v k="$N_PC" 'BEGIN{printf "ID"; for(i=1;i<=k;i++) printf ",PC%d", i; print ""}
     {printf "%s", $2; for(i=3;i<3+k;i++) printf ",%s", $i; print ""}' pca.eigenvec > "$PCS_CSV"
[ "$(($(wc -l < "$PCS_CSV") - 1))" -eq "$NS" ] || { echo "[prep] pcs.csv rows != lines" >&2; exit 1; }
echo "[prep] PCs: $PCS_CSV"; head -3 "$PCS_CSV"
echo "[prep] done"
