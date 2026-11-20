#!/bin/bash
#BSUB -n 4
#BSUB -W 1440
#BSUB -q sara
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=64GB]"
#BSUB -J bap_merge_vcf
#BSUB -o /share/maize/ntanduk/BAP_gwas/logs/merge_vcf.%J.out
#BSUB -e /share/maize/ntanduk/BAP_gwas/logs/merge_vcf.%J.err
############################################################################
## Step 00 (once): join the 10 per-chromosome BAP VCFs (MAF >= 0.01) into one
## genome-wide VCF, so that GEMMA's kinship in step 01 uses all chromosomes.
## Also reports how many phenotyped lines are in the VCF.
## Submit (from /share/maize/ntanduk/BAP_gwas/scripts):  bsub < 00_merge_vcf_LSF.sh
############################################################################
ENV=/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/vcf2gwas
export PATH="$ENV/bin:$PATH"

ROOT=/share/maize/ntanduk/BAP_gwas
IN=/rsstu/users/r/rrellan/DOE_CAREER/BAP/maf_01_perc
VCF=$ROOT/geno/BAP_allchr.maf01.vcf.gz
PHENO=$ROOT/pheno/BAP_pheno_vcf2gwas.csv

mkdir -p "$ROOT/geno" "$ROOT/logs"
echo "=== merge on $(hostname) at $(date) ==="
files=()
for c in 01 02 03 04 05 06 07 08 09 10; do
  f=$IN/NEW_TERRA_HA_fixed_Genotyped.recalibrated.filtered.snps_only.Chr${c}.maf01.vcf.gz
  [ -f "$f" ] || { echo "missing: $f"; exit 1; }
  files+=("$f")
done
bcftools concat --threads 4 -Oz -o "$VCF" "${files[@]}"
bcftools index -t "$VCF"
echo "variants: $(bcftools index -n "$VCF")"

## phenotyped lines present in the VCF
bcftools query -l "$VCF" | sort > "$ROOT/geno/vcf_samples.txt"
awk -F, 'NR > 1 { print $1 }' "$PHENO" | sort > "$ROOT/geno/pheno_samples.txt"
echo "VCF samples: $(wc -l < "$ROOT/geno/vcf_samples.txt")"
echo "phenotyped lines: $(wc -l < "$ROOT/geno/pheno_samples.txt")"
echo "phenotyped lines in the VCF: $(comm -12 "$ROOT/geno/vcf_samples.txt" "$ROOT/geno/pheno_samples.txt" | wc -l)"
echo "=== done $(date) ==="
