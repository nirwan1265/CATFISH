#!/bin/bash
#BSUB -n 4
#BSUB -W 14400
#BSUB -q sara
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=64GB]"
#BSUB -J "bap_gwas[1-12]"
#BSUB -o /share/maize/ntanduk/BAP_gwas/logs/gwas.%J.%I.out
#BSUB -e /share/maize/ntanduk/BAP_gwas/logs/gwas.%J.%I.err
############################################################################
## Step 01: vcf2gwas (GEMMA LMM, Wald test) on the BAP, one phenotype per
## array element, genome-wide (kinship from all chromosomes), MAF >= 0.05,
## 3 genotype PCs as covariates (-cf PCA -c 3, as the paper's dry-tons GWAS
## and run_sorghum.sh), no plots (-np -nq -nl).
## vcf2gwas names its working files after the phenotype file, the VCF and a
## start timestamp, so parallel runs collide unless these differ: each trait
## gets its own folder, its own copy of the phenotype file, its own timestamp,
## and the jobs start 20 s apart. vcf2gwas also writes working files next to
## the VCF, so each trait runs on its own copy of the VCF (runs/<trait>/vcf/),
## removed when the trait finishes without error.
## Trait 6 (Dry_weight_kg) is Dry_tons_per_acre / 32, the same GWAS: skip it.
##   bsub -J "bap_gwas[1-5,7-12]" < 01_gwas_trait_LSF.sh
## Activate the vcf2gwas env before submitting.
############################################################################
# ENV=/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/vcf2gwas
# export PATH="$ENV/bin:$PATH"

ROOT=/share/maize/ntanduk/BAP_gwas
VCF=$ROOT/geno/BAP_allchr.maf01.vcf.gz
PHENO=$ROOT/pheno/BAP_pheno_vcf2gwas.csv
TRAITS=(Anthesis_date_days Harvest_date_days Total_fresh_weight_kg Brix_maturity Brix_milk
        Dry_weight_kg Stalk_height_cm Dry_tons_per_acre ADF_pct_DM NDF_pct_DM NFC_pct_DM Lignin_pct_DM)
TRAIT=${TRAITS[$((LSB_JOBINDEX - 1))]}
OUT=$ROOT/runs/$TRAIT

for f in "$VCF" "$PHENO"; do [ -f "$f" ] || { echo "missing: $f"; exit 1; }; done
mkdir -p "$OUT"
cd "$OUT" || exit 1
sleep $(( (LSB_JOBINDEX - 1) * 20 ))
cp "$PHENO" "$OUT/pheno_${TRAIT}.csv"
mkdir -p "$OUT/vcf"
VCF_RUN=$OUT/vcf/$(basename "$VCF")
cp "$VCF" "$VCF_RUN"
[ -f "${VCF}.tbi" ] && cp "${VCF}.tbi" "${VCF_RUN}.tbi"
[ -f "${VCF}.csi" ] && cp "${VCF}.csi" "${VCF_RUN}.csi"
echo "=== $TRAIT on $(hostname) at $(date) ==="
echo "VCF:   $VCF_RUN (copy of $VCF)"
echo "PHENO: $OUT/pheno_${TRAIT}.csv ($(( $(wc -l < "$PHENO") - 1 )) lines)"

vcf2gwas -v "$VCF_RUN" -pf "$OUT/pheno_${TRAIT}.csv" -p "$TRAIT" --timestamp "$TRAIT" \
         -lmm 1 -cf PCA -c 3 -q 0.05 -np -nq -nl -T 4 -M 60000 -o "$OUT" \
  || { echo "vcf2gwas failed for $TRAIT (VCF copy kept in $OUT/vcf)"; exit 1; }
rm -rf "$OUT/vcf"

echo "=== done $TRAIT $(date) ==="
