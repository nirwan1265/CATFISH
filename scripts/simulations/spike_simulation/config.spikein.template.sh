#!/bin/bash

# Prefilled from the existing HPC permutation_pipeline config.
# Copy to config.spikein.sh if you want a separate editable copy.

export TRAIT="Dry_tons_per_acre"
export PHENO_CSV="/share/maize/ntanduk/landadapt/convert/BAP_energy_traits_renamed_for_genotype.csv"

export CHROMS="01 02 03 04 05 06 07 08 09 10"
export MASTER_SEED=20260627

HPC_TOOLS=/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk
export CATFISH_REPO="${HPC_TOOLS}/CATFISH"
export PATHWAY_FILE="${CATFISH_REPO}/inst/extdata/pathway/sorghumbicolorcyc_pathways.20230103.SORBI"
export GENE_LOC_FILE="${CATFISH_REPO}/inst/extdata/sorghum.genes.loc"

export VCF_DIR="/rsstu/users/r/rrellan/DOE_CAREER/BAP/maf_01_perc"
export VCF_PREFIX_PATTERN="NEW_TERRA_HA_fixed_Genotyped.recalibrated.filtered.snps_only.Chr@CHR@.maf01.vcf.gz"

export GENO_DIR="/rsstu/users/r/rrellan/DOE_CAREER/BAP/maf_01_perc/plink_bfiles"
export GENO_PREFIX_PATTERN="NEW_TERRA_HA_fixed_Genotyped.recalibrated.filtered.snps_only.Chr@CHR@.maf01"

export PLINK_BIN="${HPC_TOOLS}/plink"
export PLINK="${PLINK_BIN}"
export VCF2GWAS="vcf2gwas"
export MAGMA="${HPC_TOOLS}/magma_v1.10/magma"
export RSCRIPT="Rscript"
export GENE_MODEL="multi=snp-wise"
export SPECIES="sorghum"

# Reuse the same genotype-only MAGMA annotation strategy as the old pipeline.
export ANNOT_PATH_TEMPLATE="/share/maize/ntanduk/landadapt/convert/Lignin_magma_single_chr/chr@CHR@/annot/Dry_tons_per_acre_chr@CHR@.genes.annot"

# New spike-in run root. This is separate from the phenotype-permutation run.
export SPIKE_ROOT="/share/maize/ntanduk/CATFISH/simulation/permutation_spike_run"
export PHENO_DIR="${SPIKE_ROOT}/permuted_phenotypes"
export GWAS_DIR="${SPIKE_ROOT}/gwas"
export MAGMA_OUT="${SPIKE_ROOT}/magma"
export CATFISH_OUT="${SPIKE_ROOT}/catfish"
export FINAL_DIR="${SPIKE_ROOT}/final"
export SPIKE_SUMMARY_DIR="${SPIKE_ROOT}/summary"
export SPIKE_WORK="${SPIKE_ROOT}/spike_build"
export CORPAIRS="${SPIKE_ROOT}/gene_cor_pairs_CONSTANT.txt"
export NTOTAL_FILE="${SPIKE_ROOT}/N_TOTAL.txt"

# Single-pathway staged spike controls (preferred stable route: 10 -> 11 -> 12).
# Set SPIKE_PATHWAY to a real sorghum pathway_id before running staged spike-ins.
export SPIKE_PATHWAY=""
export SPIKE_GENE_COL=""
export SPIKE_ARCH="MDS"
export SPIKE_H2=0.10
export SPIKE_N_STRONG=3
export SPIKE_MSA_FRAC=0.6
export SPIKE_N_SIM=100
export SPIKE_SEED=99001

# CATFISH settings
export MVN_CALIBRATE_COMPONENTS=false
export PATCH_EMP_NULL_LOWER=false
export CATFISH_THREADS=1
export TAU_GRID="0.1,0.05,0.01"
export CATFISH_B_PERM=10000
export TAIL_MODE="empirical"
export TAIL_SWITCH_EXCEED=10
export TAIL_GPD_K=250
export TAIL_MIN_B=10000
export TAIL_MIN_TAIL=50

mkdir -p "$PHENO_DIR" "$GWAS_DIR" "$MAGMA_OUT" "$CATFISH_OUT" "$FINAL_DIR" "$SPIKE_SUMMARY_DIR" "$SPIKE_WORK"
