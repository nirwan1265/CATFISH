#!/bin/bash
# =============================================================================
# config.sh -- simulation_check_catfish: CATFISH, MAGMA and PascalX on simulated
# GWAS statistics with known truth (real 1000 Genomes EUR LD, chromosomes 1-10,
# HapMap3 SNPs). Sourced by every step. See README.md.
#   MODE=pilot (default) or MODE=full selects the design and the output folder.
# =============================================================================
export SIM_ROOT=/share/maize/ntanduk/CATFISH/simulation_check_catfish   # scripts live here
export MODE="${MODE:-pilot}"
export IN="${SIM_ROOT}/inputs"                 # panel, annotation, pathways (shared by pilot and full)
export RUN="${SIM_ROOT}/run_${MODE}"           # design, simulated z, MAGMA, CATFISH, MAGMA gsa, PascalX
mkdir -p "$IN" "$RUN" "${SIM_ROOT}/logs" 2>/dev/null || true

## ---- tools ---------------------------------------------------------------------
export HPC_TOOLS=/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk
export MAGMA="${HPC_TOOLS}/magma_v1.10/magma"
export PLINK2="${HPC_TOOLS}/plink2"
export PLINK="${HPC_TOOLS}/plink"            # PLINK 1.9, only to write .tped for PascalX (as compare/20)
export PASCALX_ENV=/usr/local/usrapps/maize/ntanduk/pascalx
export PASCALX_DIR=/usr/local/usrapps/maize/ntanduk/PascalX
## DEV_WY=<path to westfall_young.R> uses that file's catfish_wy_pathways() on top
## of the installed CATFISH (as gc_real/wy_gc_real.R); leave empty once the
## installed package has the competitive test
export DEV_WY="${DEV_WY:-/share/maize/ntanduk/CATFISH/gc_real/dev/westfall_young.R}"

## ---- inputs from the lipid analyses (unchanged) ---------------------------------
export LIPID_ROOT=/share/maize/ntanduk/CATFISH
export PANEL_FULL="${LIPID_ROOT}/g1000_eur"                 # the panel MAGMA used for the lipids
export ANNOT_FULL="${LIPID_ROOT}/g1000_eur.genes.annot"     # +/- 10 kb; ldl_eur.genes.annot if absent
export REACTOME_SETS=/share/maize/ntanduk/CATFISH/gc_real/inputs/c2.cp.reactome.v2026.2.Hs.entrez.sets.annot
export PASCALX_GENOME=/share/maize/ntanduk/CATFISH/compare/inputs/pascalx_genome_entrez.tsv
export HM3="${IN}/w_hm3.snplist"                            # LDSC HapMap3 list (rsync from the Mac)

## ---- simulation panel and pathways (written by step 00) ----------------------
export CHROMS="1 2 3 4 5 6 7 8 9 10"
export SIM_PANEL="${IN}/panel/g1000_eur_hm3_chr1_10"        # PLINK files
export SIM_COUNTS="${IN}/panel/g1000_eur_hm3_chr1_10.gcount" # genotype counts (plink2 --geno-counts)
export SIM_ANNOT="${IN}/sim.genes.annot"                    # MAGMA annotation restricted to the panel
export SIM_SETS="${IN}/reactome_chr1_10.sets.annot"         # Reactome, genes on chr 1-10, 10-200 genes
export PASCALX_REF="${IN}/pascalx_ref"                      # PascalX import of the simulation panel
export MIN_GENES=10
export MAX_GENES=200

## ---- design --------------------------------------------------------------------
export MASTER_SEED=20261109                # replicate r uses MASTER_SEED + r
export N_REP_FULL=50                       # replicates per scenario x target size (full run)
export N_REP_PILOT=5
export SIM_N=10000                         # value of the N column in the SNP p-value files
## signal: one causal SNP per causal gene, expected z = lambda (random sign)
##   archetype layers (17_archetype_power_wy.R): lambda = EFFECT_SCALE x delta
export EFFECT_SCALE="${EFFECT_SCALE:-4}"   # full run: set from the pilot
export PILOT_SCALES="2 4 6"                # pilot: S3 at these scales
export BG_FRAC=0.02                        # background: share of genes with one causal SNP (S1-S6)
export BG_LAMBDA=3                         # background causal SNP expected z (random sign)
export LD_WINDOW_KB=1000                  # a causal SNP's signal reaches the SNPs within this, through their panel LD
export BUFFER_KB=1000                      # a pathway with no causal SNP within this of its genes is "clean"

## ---- methods -------------------------------------------------------------------
export WY_B=2000                           # CATFISH null replicates (self-contained and competitive)
export WY_SEED=20260706                    # as in the paper's real and permutation analyses
export PASCALX_WINDOW=10000                # +/- 10 kb, as MAGMA
export PASCALX_MAF=0.01
