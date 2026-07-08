#!/bin/bash

export TRAIT="Stem_volume"
export CHROMS="01 02 03 04 05 06 07 08 09 10"
export N_PERM=100
export MASTER_SEED=20260706

# Core tools
export CATFISH_REPO="/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/CATFISH"
export RSCRIPT="Rscript"
export PLINK="/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/plink"
export MAGMA="/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/magma_v1.10/magma"
export VCF2GWAS="vcf2gwas"
export GENE_MODEL="multi=snp-wise"

# Pipeline output root
export PIPE_ROOT="/share/maize/ntanduk/CATFISH/simulation/SAP_results_simulation"
export INPUT_DIR="${PIPE_ROOT}/inputs"

# Input data
# Put the phenotype CSV directly in the script folder on HPC as:
#   /share/maize/ntanduk/CATFISH/simulation/sap_stem_volume_permutation/stem_volume.csv
# source on laptop:
#   /Users/nirwantandukar/Documents/Research/data/SAP/GWAS_phenotype/Stem_volume/stem_volume.csv
export PHENO_CSV="/share/maize/ntanduk/CATFISH/simulation/sap_stem_volume_permutation/stem_volume.csv"
# Put the real stem-volume GWAS assoc file directly in the script folder on HPC as:
#   /share/maize/ntanduk/CATFISH/simulation/sap_stem_volume_permutation/stem_volume_real_gwas.assoc.txt
# source on laptop:
#   /Users/nirwantandukar/Documents/Research/results/GWAS/SAP/Stem_diameter/Stem_volume_mod_sub_stem_volume_SAP_bialleles_MAF_0.05_11.assoc.txt
export REAL_GWAS_FILE="/share/maize/ntanduk/CATFISH/simulation/sap_stem_volume_permutation/stem_volume_real_gwas.assoc.txt"

# SAP genotype VCFs. Adjust this pattern if your HPC files use a different name.
export VCF_DIR="/rsstu/users/r/rrellan/DOE_CAREER/SorghumGEA/data/SAP/chromosome"
export VCF_PREFIX_PATTERN="SAP_only_samples_bialleles_MAF_0.05_chr@CHR_NUM@.vcf.gz"

# CATFISH reference files from the real HPC CATFISH checkout
export GENE_LOC_FILE="${CATFISH_REPO}/inst/extdata/sorghum.genes.loc"
export PATHWAY_FILE="${CATFISH_REPO}/inst/extdata/pathway/sorghumbicolorcyc_pathways.20230103.SORBI"

export PREP_ROOT="${PIPE_ROOT}/prep"
export PHENO_DIR="${PIPE_ROOT}/permuted_phenotypes"
export GWAS_DIR="${PIPE_ROOT}/gwas"
export MAGMA_OUT="${PIPE_ROOT}/magma"
export CATFISH_OUT="${PIPE_ROOT}/catfish"
export FINAL_DIR="${PIPE_ROOT}/final"
export NTOTAL_FILE="${PIPE_ROOT}/N_TOTAL.txt"
export CORPAIRS="${PIPE_ROOT}/stem_volume_gene_cor_pairs.txt"

# Prepared genotype / annotation locations
export GENO_DIR="${PREP_ROOT}/plink_bfiles"
export GENO_PREFIX_PATTERN="sap_chr@CHR@"
export MERGED_BFILE="${PREP_ROOT}/sap_merged"
export FILTERED_GWAS_FILE="${PREP_ROOT}/stem_volume_gwas_filtered.txt"
export ANNOT_DIR="${PREP_ROOT}/annot"
export ANNOT_PATH_TEMPLATE="${ANNOT_DIR}/chr@CHR@/stem_volume_chr@CHR@.genes.annot"

# CATFISH settings
export MVN_CALIBRATE_COMPONENTS="false"
export PATCH_EMP_NULL_LOWER="false"
export CATFISH_THREADS=1
export TAU_GRID="0.1,0.05,0.01"
export CATFISH_B_PERM=10000
export TAIL_MODE="hybrid_gpd"
export TAIL_SWITCH_EXCEED=10
export TAIL_GPD_K=250
export TAIL_MIN_B=10000
export TAIL_MIN_TAIL=50

mkdir -p \
  "$PIPE_ROOT" "$INPUT_DIR" "$PREP_ROOT" "$PHENO_DIR" "$GWAS_DIR" "$MAGMA_OUT" \
  "$CATFISH_OUT" "$FINAL_DIR" "$GENO_DIR" "$ANNOT_DIR"
