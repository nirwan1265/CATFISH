#!/bin/bash
# =============================================================================
# config.sh -- BAP phenotype-permutation null for the Westfall-Young paper.
# Sourced by every step. Paths are the ones the earlier permutation_pipeline
# used (scripts/simulations/spike_simulation/config.spikein.template.sh).
# =============================================================================
export WORKING_DIR=/share/maize/ntanduk/CATFISH/simulation/bap_perm_wy   # scripts live here
export RUN_ROOT="${WORKING_DIR}/run"                                     # all outputs
export INPUT_DIR="${WORKING_DIR}/inputs"                                  # phenotype, annotations, pathways

## trait and permutations
export TRAIT="Dry_tons_per_acre"
export PHENO_CSV="${INPUT_DIR}/BAP_energy_traits_renamed_for_genotype.csv"   # first column = line IDs as in the VCF
export N_PERM=100
export MASTER_SEED=20260627        # same seed as the earlier run: permutation b uses seed MASTER_SEED + b
export CHROMS="01 02 03 04 05 06 07 08 09 10"

## genotypes (BAP, MAF >= 1%)
export VCF_DIR="/rsstu/users/r/rrellan/DOE_CAREER/BAP/maf_01_perc"
export VCF_PREFIX_PATTERN="NEW_TERRA_HA_fixed_Genotyped.recalibrated.filtered.snps_only.Chr@CHR@.maf01.vcf.gz"
export GENO_DIR="/rsstu/users/r/rrellan/DOE_CAREER/BAP/maf_01_perc/plink_bfiles"
export GENO_PREFIX_PATTERN="NEW_TERRA_HA_fixed_Genotyped.recalibrated.filtered.snps_only.Chr@CHR@.maf01"

## tools
HPC_TOOLS=/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk
export CATFISH_REPO="${HPC_TOOLS}/CATFISH"    # git pull first: needs catfish_wy_pathways()
export VCF2GWAS="vcf2gwas"
export MAGMA="${HPC_TOOLS}/magma_v1.10/magma"

## MAGMA: the genotype-only per-chromosome annotation the real analysis used
## (+/- 25 kb), and the gene model.
##   snp-wise=mean   the paper's gene p-value (P_SNPWISE_MEAN), fast
##   multi=snp-wise  also gives P_MULTI / P_SNPWISE_TOP1, but ~8 h per chromosome
export ANNOT_PATH_TEMPLATE="${INPUT_DIR}/annot/Dry_tons_per_acre_chr@CHR@.genes.annot"
## built by 02b_annotate_LSF.sh exactly as scripts/run_magma_single_chr_dry_tons_server.R
## did: SNP locations from the PLINK .bim (all genotyped lines), the package's
## sorghum gene locations, window 25,25 kb
export GENE_LOC_FILE="${CATFISH_REPO}/inst/extdata/sorghum.genes.loc"
export WINDOW_KB="25,25"
export GENE_MODEL="snp-wise=mean"

## pathways and Westfall-Young settings: identical to 23_westfall_young_sorghum.R
export PATHWAY_FILE="${INPUT_DIR}/sorghumbicolorcyc_pathways.20230103.SORBI"
export MIN_GENES=5
export MAX_GENES=200
export WY_B=10000

## built once by 00_prepare_genome_LSF.sh: per-chromosome VCFs of the
## phenotyped lines (in .fam order), and from all ten together the genome-wide
## kinship matrix and 3 PCs that every per-chromosome GWAS reuses
export REF_DIR="${RUN_ROOT}/ref"
export CHR_VCF_PATTERN="${REF_DIR}/bap_phenotyped_chr@CHR@.vcf.gz"
export GENOME_VCF="${REF_DIR}/bap_genome_phenotyped.vcf.gz"
export KINSHIP="${REF_DIR}/kinship.cXX.txt"
export PCS_CSV="${REF_DIR}/pcs.csv"
export N_PC=3

## derived
export PHENO_DIR="${RUN_ROOT}/permuted_phenotypes"
export GWAS_DIR="${RUN_ROOT}/gwas"
export MAGMA_OUT="${RUN_ROOT}/magma"
export WY_OUT="${RUN_ROOT}/wy"
export FINAL_DIR="${RUN_ROOT}/final"
export NTOTAL_FILE="${RUN_ROOT}/N_TOTAL.txt"
export CORPAIRS="${RUN_ROOT}/gene_cor_pairs_CONSTANT.txt"
mkdir -p "$INPUT_DIR" "$REF_DIR" "$PHENO_DIR" "$GWAS_DIR" "$MAGMA_OUT" "$WY_OUT" "$FINAL_DIR" "$WORKING_DIR/logs"
