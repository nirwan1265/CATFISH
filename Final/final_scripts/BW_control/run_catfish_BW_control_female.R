#!/usr/bin/env Rscript
# ==============================================================================
# CATFISH Pipeline — BW_control_female (Body Weight Control, Female)
# Step 1: Build gene-gene correlation pairs from .genes.raw
# Step 2: Run CATFISH omnibus pathway analysis (B=100k, MVN, GPD)
# ==============================================================================

suppressPackageStartupMessages({
  library(CATFISH)
  library(data.table)
  library(dplyr)
})

## -----------------------------------------------------------------------------
## Parameters
## -----------------------------------------------------------------------------

TRAIT      <- "BW_control_female"
CHROMS     <- c("2L", "2R", "3L", "3R", "X")

MAGMA_BASE <- "/Users/nirwantandukar/Documents/Research/results/CATFISH/MAGMA/BW_control_female"
OUT_BASE   <- "/Users/nirwantandukar/Documents/Research/results/CATFISH/DGRP/BW_control_female"
OUT_DIR    <- file.path(OUT_BASE, "CATFISH_permutation_B1000000_mvn_GPD")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

GENE_LOC_FILE     <- "/Users/nirwantandukar/Documents/Github/MAGCAT/inst/extdata/fly.genes.loc"
GENE_LENGTHS_FILE <- "/Users/nirwantandukar/Documents/Github/MAGCAT/inst/extdata/dmel.flybase.fbgn.genes.loc.tsv"
PATHWAY_FILE  <- "/Users/nirwantandukar/Documents/Github/MAGCAT/inst/extdata/pathway/Fly_Cyc.tsv"

B_PERM      <- 1000000L
PERM_MODE   <- "mvn"
OMNIBUS     <- "ACAT"
N_THREADS   <- 12
SEED        <- 123
TAU_GRID    <- c(0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.001)

TAIL_MODE          <- "hybrid_gpd"
TAIL_SWITCH_EXCEED <- 10L
TAIL_GPD_K         <- 250L
TAIL_MIN_B         <- 10000L
TAIL_MIN_TAIL      <- 50L
MIN_P              <- 1e-30

cat("=== CATFISH Pipeline:", TRAIT, "===\n")
cat("B_PERM:", B_PERM, "| Mode:", PERM_MODE, "| Tail:", TAIL_MODE, "\n")
cat("Output:", OUT_DIR, "\n\n")

## -----------------------------------------------------------------------------
## 1. Combine per-chromosome MAGMA gene results
## -----------------------------------------------------------------------------

cat("Loading MAGMA gene results...\n")

genes_list <- lapply(CHROMS, function(chr) {
  f <- file.path(MAGMA_BASE, chr, "magma_gene",
                 sprintf("%s_%s.multi_snp_wise.genes.out", TRAIT, chr))
  if (!file.exists(f)) {
    message("  Skipping missing chromosome: ", chr)
    return(NULL)
  }
  x <- read.table(f, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  message("  ", chr, ": ", nrow(x), " genes")
  x
})

genes_all <- rbindlist(Filter(Negate(is.null), genes_list))
genes_all <- genes_all[!duplicated(genes_all$GENE), ]
if (!"P" %in% names(genes_all)) names(genes_all)[9] <- "P"
genes_all <- genes_all[order(genes_all$P), ]

cat(sprintf("  Total unique genes: %d\n\n", nrow(genes_all)))

combined_file <- file.path(OUT_DIR, paste0(TRAIT, "_combined_genes.tsv"))
write.table(genes_all, combined_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved combined genes:", combined_file, "\n\n")

## -----------------------------------------------------------------------------
## 2. Adjust gene p-values
## -----------------------------------------------------------------------------

cat("Adjusting gene p-values...\n")
gene_lengths <- read.delim(GENE_LENGTHS_FILE, stringsAsFactors = FALSE)

adj <- catfish_adjust_gene_p(
  gene_results = genes_all,
  gene_lengths  = gene_lengths,
  gene_col      = "GENE",
  nsnp_col      = "NSNPS",
  p_col         = "P",
  z_col         = "ZSTAT",
  len_gene_col  = "gene_id",
  len_col       = "length"
)
genes_adj <- adj[, c("gene_id", "z_adj", "p_adj")]
colnames(genes_adj) <- c("GENE", "ZSTAT", "P")
cat(sprintf("  Adjusted %d genes\n\n", nrow(genes_adj)))

## -----------------------------------------------------------------------------
## 3. Build gene-gene correlation pairs from .genes.raw
## -----------------------------------------------------------------------------

cor_pairs_file <- file.path(OUT_DIR, paste0(TRAIT, "_gene_cor_pairs.txt"))

if (!file.exists(cor_pairs_file)) {
  cat("Building gene-gene correlation pairs...\n")

  raw_files <- sapply(CHROMS, function(chr) {
    f <- file.path(MAGMA_BASE, chr, "magma_gene",
                   sprintf("%s_%s.multi_snp_wise.genes.raw", TRAIT, chr))
    if (file.exists(f)) f else NULL
  })
  raw_files <- Filter(Negate(is.null), raw_files)
  cat(sprintf("  Found %d .genes.raw files\n", length(raw_files)))

  # Build pairs per chromosome then concatenate
  tmp_files <- c()
  for (i in seq_along(raw_files)) {
    chr_name <- names(raw_files)[i]
    tmp_out  <- file.path(OUT_DIR, sprintf("gene_cor_pairs_%s.txt", chr_name))
    tryCatch({
      magma_genesraw_to_cor_pairs_banded(
        genes_raw_file = raw_files[[i]],
        out_pairs_file = tmp_out,
        keep_abs_r_ge  = 0
      )
      tmp_files <- c(tmp_files, tmp_out)
      cat("  Done:", chr_name, "\n")
    }, error = function(e) message("  Warning for ", chr_name, ": ", e$message))
  }

  # Concatenate (header once)
  if (length(tmp_files) > 0) {
    first  <- fread(tmp_files[1])
    others <- lapply(tmp_files[-1], fread)
    all_pairs <- rbindlist(c(list(first), others))
    fwrite(all_pairs, cor_pairs_file, sep = "\t")
    cat(sprintf("  Total correlation pairs: %d\n", nrow(all_pairs)))
    file.remove(tmp_files)
  }
} else {
  cat("Gene-gene correlation pairs already exist:", cor_pairs_file, "\n")
}
cat("\n")

## -----------------------------------------------------------------------------
## 4. Load FlyCyc pathways
## -----------------------------------------------------------------------------

cat("Loading FlyCyc pathways...\n")
fly_pw <- read.delim(PATHWAY_FILE, stringsAsFactors = FALSE)
cat(sprintf("  %d pathways | %d gene-pathway entries\n\n",
            length(unique(fly_pw$pathway_id)), nrow(fly_pw)))

## -----------------------------------------------------------------------------
## 5. Run CATFISH
## -----------------------------------------------------------------------------

cat("Running CATFISH (B =", B_PERM, ")...\n\n")

omni <- catfish_omni2_pathways(
  gene_results             = genes_adj,
  pathways                 = fly_pw,
  gene_col                 = "GENE",
  p_raw_col                = "P",
  z_col                    = "ZSTAT",
  tau_grid                 = TAU_GRID,
  omnibus                  = OMNIBUS,
  B_perm                   = B_PERM,
  perm_mode                = PERM_MODE,
  magma_cor_file           = cor_pairs_file,
  mvn_marginal             = "uniform",
  mvn_calibrate_components = FALSE,
  tail_mode                = TAIL_MODE,
  tail_switch_exceed       = TAIL_SWITCH_EXCEED,
  tail_gpd_k               = TAIL_GPD_K,
  tail_min_B               = TAIL_MIN_B,
  tail_min_tail            = TAIL_MIN_TAIL,
  min_p                    = MIN_P,
  min_genes                = 2,
  n_threads                = N_THREADS,
  seed                     = SEED,
  output                   = TRUE,
  out_dir                  = OUT_DIR
)

## -----------------------------------------------------------------------------
## 6. Save and summarise
## -----------------------------------------------------------------------------

out_csv <- file.path(OUT_DIR, paste0(TRAIT, "_CATFISH_ACAT_mvn_B", B_PERM, "_GPD.csv"))
write.csv(omni, out_csv, row.names = FALSE)

omni_col <- if ("omni_p_final" %in% names(omni)) "omni_p_final" else "omni_p_mvn"
bonf     <- 0.05 / nrow(omni)
sig_bonf <- omni[omni[[omni_col]] < bonf, ]
sig_fdr  <- omni[p.adjust(omni[[omni_col]], "BH") < 0.05, ]

cat("\n=== CATFISH RESULTS:", TRAIT, "===\n")
cat("Total pathways tested:  ", nrow(omni), "\n")
cat("Bonferroni significant: ", nrow(sig_bonf), "\n")
cat("FDR < 0.05:             ", nrow(sig_fdr), "\n\n")
cat("Top 20 pathways:\n")
top20 <- omni[order(omni[[omni_col]]), ]
print(top20[1:min(20, nrow(top20)),
            c("pathway_id", "pathway_name", "n_genes", omni_col)],
      row.names = FALSE)

cat("\nResults saved to:", out_csv, "\n")
cat("=== DONE:", TRAIT, "===\n")
