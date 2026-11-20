#!/bin/bash
#BSUB -n 1
#BSUB -W 60
#BSUB -q sara
#BSUB -R "rusage[mem=8GB]"
#BSUB -J bap_chr4_locus
#BSUB -o /share/maize/ntanduk/BAP_gwas/logs/chr4_locus.%J.out
#BSUB -e /share/maize/ntanduk/BAP_gwas/logs/chr4_locus.%J.err
############################################################################
## Step 05: inputs for the chr4 sugar/fibre locus figure
## (catfish2/scripts/final/43_bap_chr4_locus.R). Window chr4 63.85-64.10 Mb.
##   region_<trait>.tsv  every SNP in the window: chr rs ps allele1 allele0 af beta se p_wald
##                       (Brix_maturity, Brix_milk, NDF, NFC, ADF)
##   leads.tsv           each of those traits' strongest SNP in the window
##   lead_effects.tsv    the lead SNPs' GEMMA effects on all 11 traits
##   leads_geno.raw      genotype dosage of the lead SNPs per line (plink --recode A),
##                       from the chr4 bfile MAGMA used (same SNP IDs as GEMMA)
##   leads_ld.ld/.snplist  r2 between the lead SNPs
## Output in /share/maize/ntanduk/BAP_gwas/chr4_locus
##   bsub < 05_chr4_locus_LSF.sh
############################################################################
set -eu
ROOT=/share/maize/ntanduk/BAP_gwas
PLINK=/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/plink
BFILE=/rsstu/users/r/rrellan/DOE_CAREER/BAP/maf_01_perc/plink_bfiles/NEW_TERRA_HA_fixed_Genotyped.recalibrated.filtered.snps_only.Chr04.maf01
CHR=4; FROM=63850000; TO=64100000
LEAD_TRAITS="Brix_maturity Brix_milk NDF_pct_DM NFC_pct_DM ADF_pct_DM"
ALL_TRAITS="Anthesis_date_days Harvest_date_days Total_fresh_weight_kg Brix_maturity Brix_milk Stalk_height_cm Dry_tons_per_acre ADF_pct_DM NDF_pct_DM NFC_pct_DM Lignin_pct_DM"
OUT=$ROOT/chr4_locus
mkdir -p "$OUT"; cd "$OUT"

assoc_of() { find "$ROOT/runs/$1" -name "*.assoc.txt" | head -1; }
COLS='chr rs ps allele1 allele0 af beta se p_wald'

## window SNPs and the lead SNP per trait (columns by header name)
: > leads.tsv
for t in $LEAD_TRAITS; do
  f=$(assoc_of "$t"); [ -n "$f" ] || { echo "no .assoc.txt for $t"; exit 1; }
  awk -v c=$CHR -v a=$FROM -v b=$TO -v cols="$COLS" 'BEGIN{OFS="\t"; n=split(cols,k," ")}
    NR==1{for(i=1;i<=NF;i++) h[$i]=i; for(j=1;j<=n;j++) printf "%s%s", k[j], (j<n?OFS:ORS); next}
    $h["chr"]+0==c && $h["ps"]>=a && $h["ps"]<=b {for(j=1;j<=n;j++) printf "%s%s", $h[k[j]], (j<n?OFS:ORS)}' \
    "$f" > "region_${t}.tsv"
  awk 'NR>1' "region_${t}.tsv" | sort -t$'\t' -g -k9,9 | head -1 | awk -v t="$t" 'BEGIN{OFS="\t"}{print t,$2,$3,$9}' >> leads.tsv
  echo "$t: $(( $(wc -l < "region_${t}.tsv") - 1 )) SNPs in window; lead $(tail -1 leads.tsv | cut -f2-4)"
done
cut -f2 leads.tsv | sort -u > leads.snps

## lead SNPs' effects on every trait
printf "trait\t%s\n" "$(echo $COLS | tr ' ' '\t')" > lead_effects.tsv
for t in $ALL_TRAITS; do
  f=$(assoc_of "$t"); [ -n "$f" ] || { echo "no .assoc.txt for $t (skipped)"; continue; }
  awk -v t="$t" -v cols="$COLS" 'BEGIN{OFS="\t"; n=split(cols,k," ")}
    NR==FNR{s[$1]=1; next}
    FNR==1{for(i=1;i<=NF;i++) h[$i]=i; next}
    ($h["rs"] in s){printf "%s", t; for(j=1;j<=n;j++) printf "%s%s", OFS, $h[k[j]]; printf "\n"}' \
    leads.snps "$f" >> lead_effects.tsv
done

## genotypes and LD of the lead SNPs
$PLINK --bfile "$BFILE" --allow-extra-chr --extract leads.snps --recode A --out leads_geno
$PLINK --bfile "$BFILE" --allow-extra-chr --extract leads.snps --r2 square --write-snplist --out leads_ld

echo "=== leads ==="; cat leads.tsv
echo "=== r2 ==="; cat leads_ld.snplist; cat leads_ld.ld
echo "done $(date)"
