#!/bin/bash
#BSUB -n 1
#BSUB -W 720
#BSUB -q sara
#BSUB -J "bap_magma[1-120]%40"
#BSUB -o /share/maize/ntanduk/BAP_gwas/logs/magma.%J.%I.out
#BSUB -e /share/maize/ntanduk/BAP_gwas/logs/magma.%J.%I.err
############################################################################
## Step 02: MAGMA gene analysis for the BAP traits, one array element per
## (trait, chromosome): element i -> trait (i-1)/10 + 1, chromosome (i-1)%10 + 1.
## Same MAGMA call as the paper's dry-tons run and the BAP permutation null
## (scripts/hpc/bap_perm/03_magma_single.sh):
##   SNP p-values  SNP = rs, P = p_wald, N = lines with the trait - n_miss
##                 (from the trait's GEMMA .assoc.txt, step 01)
##   annotation    the genotype-only per-chromosome annotation already built
##                 for dry tons (+/- 25 kb, sorghum.genes.loc)
##   model         --gene-model snp-wise=mean
## Writes magma/<trait>/<trait>_chrNN.snp_wise_mean.genes.out / .genes.raw
##   all:   bsub < 02_magma_LSF.sh
##   one trait (e.g. Stalk_height_cm = trait 7 -> elements 61-70):
##          bsub -J "bap_magma[61-70]" < 02_magma_LSF.sh
############################################################################
set -eu
ROOT=/share/maize/ntanduk/BAP_gwas
MAGMA=/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/magma_v1.10/magma
BFILE_DIR=/rsstu/users/r/rrellan/DOE_CAREER/BAP/maf_01_perc/plink_bfiles
ANNOT_DIR=/share/maize/ntanduk/CATFISH/simulation/bap_perm_wy/inputs/annot
PHENO=$ROOT/pheno/BAP_pheno_vcf2gwas.csv
TRAITS=(Anthesis_date_days Harvest_date_days Total_fresh_weight_kg Brix_maturity Brix_milk
        Dry_weight_kg Stalk_height_cm Dry_tons_per_acre ADF_pct_DM NDF_pct_DM NFC_pct_DM Lignin_pct_DM)

i=$((LSB_JOBINDEX - 1))
TRAIT=${TRAITS[$((i / 10))]}
K=$((i % 10 + 1)); CHR=$(printf "%02d" "$K")
BFILE=$BFILE_DIR/NEW_TERRA_HA_fixed_Genotyped.recalibrated.filtered.snps_only.Chr${CHR}.maf01
ANNOT=$ANNOT_DIR/Dry_tons_per_acre_chr${CHR}.genes.annot
OUTDIR=$ROOT/magma/$TRAIT
PREFIX=$OUTDIR/${TRAIT}_chr${CHR}.snp_wise_mean
mkdir -p "$OUTDIR"

[ -s "${PREFIX}.genes.out" ] && [ -s "${PREFIX}.genes.raw" ] && { echo "$TRAIT chr$CHR already done"; exit 0; }
for f in "${BFILE}.bed" "${BFILE}.bim" "$ANNOT" "$MAGMA" "$PHENO" "$ROOT/geno/vcf_samples.txt"; do
  [ -e "$f" ] || { echo "missing: $f" >&2; exit 1; }
done

## the trait's GEMMA output from step 01 (exactly one .assoc.txt)
mapfile -t A < <(find "$ROOT/runs/$TRAIT" -name "*.assoc.txt" -type f)
[ "${#A[@]}" -eq 1 ] || { echo "expected one .assoc.txt under runs/$TRAIT, found ${#A[@]}" >&2; exit 1; }
ASSOC=${A[0]}

## N_TOTAL: genotyped lines with a value for this trait
N_TOTAL=$(awk -F, -v t="$TRAIT" '
  FNR == NR { g[$1] = 1; next }
  FNR == 1  { for (j = 2; j <= NF; j++) if ($j == t) c = j; next }
  ($1 in g) && $c != "NA" && $c != "" { n++ }
  END { print n + 0 }' "$ROOT/geno/vcf_samples.txt" "$PHENO")
echo "[magma] $TRAIT chr$CHR  N_TOTAL=$N_TOTAL  assoc=$ASSOC"

## SNP P N for this chromosome (chr column compared as a number: 1 == 01)
PVAL=$OUTDIR/${TRAIT}_chr${CHR}.pval.txt
awk -v NT="$N_TOTAL" -v K="$K" '
  NR == 1 { for (j = 1; j <= NF; j++) { if ($j == "chr") c = j; if ($j == "rs") s = j; if ($j == "p_wald") p = j; if ($j == "n_miss") m = j }
            if (!c || !s || !p || !m) { print "assoc has no chr / rs / p_wald / n_miss column" > "/dev/stderr"; exit 1 }
            print "SNP\tP\tN"; next }
  ($c + 0) == K { print $s "\t" $p "\t" (NT - $m) }' "$ASSOC" > "$PVAL"

## the SNP IDs must be the .bim IDs, or MAGMA silently uses nothing
n_p=$(($(wc -l < "$PVAL") - 1))
n_in=$(awk 'FNR == NR { b[$2] = 1; next } FNR > 1 && ($1 in b) { n++ } END { print n + 0 }' "${BFILE}.bim" "$PVAL")
echo "[magma] chr$CHR: $n_p SNPs in the GWAS, $n_in of them in the .bim"
[ "$n_in" -gt 0 ] || { echo "no GWAS SNP ID matches the .bim IDs: check the rs column of $ASSOC" >&2; exit 1; }

"$MAGMA" --bfile "$BFILE" \
         --gene-annot "$ANNOT" \
         --pval "$PVAL" ncol=N use=SNP,P \
         --gene-model snp-wise=mean \
         --out "$PREFIX"

[ -s "${PREFIX}.genes.out" ] && [ -s "${PREFIX}.genes.raw" ] || { echo "no genes.out/genes.raw; see ${PREFIX}.log" >&2; exit 1; }
rm -f "$PVAL"
echo "[magma] -> ${PREFIX}.genes.out ($(($(grep -vc '^#' "${PREFIX}.genes.out") - 1)) genes)"
