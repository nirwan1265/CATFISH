#!/bin/bash
# =============================================================================
# config.sh -- comparison of CATFISH with MAGMA's own gene-set analysis (and,
# next, PascalX) on the real analyses. Sourced by every step in this folder.
# =============================================================================
export ROOT=/share/maize/ntanduk/CATFISH/compare          # scripts, inputs/, out/, logs/
export HPC_TOOLS=/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk
export MAGMA="${HPC_TOOLS}/magma_v1.10/magma"
export GC_ROOT=/share/maize/ntanduk/CATFISH/gc_real        # Reactome set file used by CATFISH
export BAP_ROOT=/share/maize/ntanduk/CATFISH/simulation/bap_perm_wy

## lipids: MAGMA's merged .genes.raw of each trait-ancestry set (rsync from the
## Mac into inputs/lipids/<ta>.genes.raw) and the Reactome sets CATFISH used
export LIPIDS="ldl_eur hdl_eur ldl_afr hdl_afr ldl_eas hdl_eas ldl_sas hdl_sas ldl_his hdl_his"
export REACTOME_SETS="${GC_ROOT}/inputs/c2.cp.reactome.v2026.2.Hs.entrez.sets.annot"

## sorghum dry biomass, rerun here with snp-wise=mean in chromosome-batch mode
## so that the ten chromosomes can be merged:
##   SNP p-values  the real run's own per-chromosome .pval.txt files (rsync from
##                 the Mac: Research/results/CATFISH/MAGMA/Dry_tons_per_acre)
##   annotation    the per-chromosome annotation built for the permutation null
##                 exactly as the real run built it (bap_perm/02b, +/- 25 kb)
export DRY_PVAL_TEMPLATE="${ROOT}/inputs/sorghum/Dry_tons_per_acre_chr@CHR@.pval.txt"
export DRY_ANNOT_TEMPLATE="${BAP_ROOT}/inputs/annot/Dry_tons_per_acre_chr@CHR@.genes.annot"
export DRY_BFILE_TEMPLATE="/rsstu/users/r/rrellan/DOE_CAREER/BAP/maf_01_perc/plink_bfiles/NEW_TERRA_HA_fixed_Genotyped.recalibrated.filtered.snps_only.Chr@CHR@.maf01"
export SORGHUM_SETS="${BAP_ROOT}/inputs/sorghumbicolorcyc.sets.annot"    # written by bap_perm/08a
export CHROMS="01 02 03 04 05 06 07 08 09 10"
export GENE_MODEL="snp-wise=mean"

export OUT="${ROOT}/out"
mkdir -p "${ROOT}/inputs/lipids" "${ROOT}/inputs/sorghum" "${OUT}/magma_dry_mean" "${OUT}/magma_gsa" "${ROOT}/logs"

## ---- PascalX (lipids) ---------------------------------------------------------
## installed under usrapps; every job activates this env and adds the C++ libs
export PASCALX_ENV=/usr/local/usrapps/maize/ntanduk/pascalx
export PASCALX_DIR=/usr/local/usrapps/maize/ntanduk/PascalX
export PLINK="${PLINK:-plink}"                     # PLINK 1.9, used once to write the reference panels as .tped
export LIPID_ROOT=/share/maize/ntanduk/CATFISH       # g1000_<panel> PLINK panels and lipids/<ta>.magma.pval, as in 14_magma_trait_LSF.sh
export GENE_LOC_37="${LIPID_ROOT}/NCBI37.3.gene.loc" # the gene locations MAGMA used for the lipids
export PASCALX_REF="${ROOT}/inputs/pascalx_ref"      # g1000_<panel>.chr<N>.tped.gz and PascalX's imported .db files
export PASCALX_GENOME="${ROOT}/inputs/pascalx_genome_entrez.tsv"   # NCBI37.3 genes, Entrez ID as gene and symbol
## settings to match the MAGMA lipid analysis: +/- 10 kb window, all SNPs the
## panel holds (PascalX's own defaults are 50 kb and MAF >= 0.05)
export PASCALX_WINDOW=10000
export PASCALX_MAF=0.01
mkdir -p "$PASCALX_REF" "${OUT}/pascalx"

## ---- sARTP (ARTP2 package; github.com/zhangh12/ARTP3 now serves ARTP2) -------
## the ten lipid sets, the same GLGC summary statistics MAGMA was given (raw p,
## not GLGC's GC p), the same 1000 Genomes panel, MAGMA's SNP-to-gene map
## (g1000_<panel>.genes.annot, +/- 10 kb, Entrez IDs) and the same Reactome
## pathways CATFISH tested for that set (10-200 of the set's genes)
export ARTP_RLIB=/usr/local/usrapps/maize/ntanduk/R_libs        # ARTP2 installed here, never in ~
export ARTP_CONDA=/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/R_plot   # the R that ARTP2 was built with; steps 32-34 activate it
export ARTP_REF="${ROOT}/inputs/artp_ref"                       # g1000_<panel>.chr<N>.{bed,bim,fam}
export ARTP_IN="${ROOT}/inputs/artp"                            # <ta>.summary.rds, <ta>.pathways.rds, <ta>.settings.tsv
export ARTP_OUT="${OUT}/artp"                                   # <ta>/chunk<k>.csv, then <ta>.csv
export GLGC_DATA="$LIPID_ROOT"                                  # <t>/<TRAIT>_INV_<ANC>_<SUFFIX>.meta.singlevar.results.gz
export ARTP_NCHUNK=30                                           # pathway chunks per set: array = 10 x ARTP_NCHUNK
## ---- settings Gemma decides (sARTP's own defaults shown where it has one) ------
export ARTP_NPERM=100000          # sARTP default 1e5; pathway p floor 1e-5 (Bonferroni over ~1,390 is 3.6e-5)
export ARTP_MIN_MARG_P=0          # sARTP default: drops SNPs with p < 1e-7 and everything within 500 kb; 0 turns it off
export ARTP_MAF=0.01              # sARTP default 0.05 (MAGMA's lipid run and the PascalX run used all panel SNPs / 0.01)
export ARTP_LAMBDA=1              # REQUIRED: 1 (no inflation correction, as MAGMA was given) or glgc (GLGC's own lambda, recovered from pvalue vs pvalue_GC)
export ARTP_N=median                   # REQUIRED: median or max (sARTP takes one sample size per study; GLGC N varies by SNP)
mkdir -p "$ARTP_REF" "$ARTP_IN" "$ARTP_OUT" 2>/dev/null || true
