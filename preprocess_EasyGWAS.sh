
# # 1) Clean CSV: drop comment lines, strip leading '>' on line 1
# awk '
#   BEGIN{FS=OFS=","}
#   /^[[:space:]]*#/ {next}                 # skip metadata/comment lines
#   !done {sub(/^[[:space:]]*>/, "", $0); done=1}   # first real line: remove leading >
#   {print}
# ' Bio6_coldest_temp.csv > Bio6_coldest_temp_clean.csv



# # ARABIDOPSIS SAMPLE: 1001
# # 93

# # https://easygwas.biochem.mpg.de/gwas/results/summary/b9e4f418-e5b1-4af1-af29-c427d8c292d3/



# # here is the genotype:
# # https://1001genomes.org/data/GMI-MPI/releases/v3.1/SNP_matrix_imputed_hdf5/

# # Take just the snps
# bcftools view -v snps \
#   -Oz \
#   -o 1001genomes_snps.vcf.gz \
#   1001genomes_snp-short-indel_only_ACGTN.vcf.gz

# bcftools index 1001genomes_snps.vcf.gz

# # Conver to ped map file for MAGMA
# plink \
#   --vcf 1001genomes_snps.vcf.gz \
#   --allow-extra-chr \
#   --recode \
#   --out 1001genomes