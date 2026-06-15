#!/usr/bin/env Rscript
# ==============================================================================
# CATFISH Pathway Analysis — Drosophila Male Starvation Stress
# B = 100,000 | perm_mode = "mvn" | GPD tail extrapolation
# Default tau grid (DGRP moderate signal)
# ==============================================================================

suppressPackageStartupMessages({
  library(CATFISH)
  library(data.table)
  library(dplyr)
})

# ==============================================================================
# PARAMETERS
# ==============================================================================

SEX         <- "male"
TRAIT       <- "Starvation_stress_male"
B_PERM      <- 100000L
PERM_MODE   <- "both"
OMNIBUS     <- "ACAT"
N_THREADS   <- 12
SEED        <- 123

# Tau grid — default (fly DGRP is moderate signal, not strong like sorghum)
TAU_GRID <- c(0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.001)

# GPD parameters
TAIL_MODE          <- "hybrid_gpd"
TAIL_SWITCH_EXCEED <- 10L
TAIL_GPD_K         <- 250L
TAIL_MIN_B         <- 10000L
TAIL_MIN_TAIL      <- 50L
MIN_P              <- 1e-30

# Paths
MAGMA_DIR  <- "/Users/nirwantandukar/Documents/Research/results/DGRP/MAGMA/Fly_magma_genes_by_chr_male"
COR_FILE   <- "/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_multi_snp_wise_genes_by_chr_N_maize/magma_gene_cor_pairs_MLM_Fly_male.txt"
GFF3_PATH  <- "/Users/nirwantandukar/Documents/Github/MAGCAT/inst/extdata/dmel.flybase.fbgn.genes.loc.tsv"

BASE_DIR   <- "/Users/nirwantandukar/Documents/Research/results/CATFISH/DGRP/Fly_starvation_male"
OUT_SUBDIR <- paste0("CATFISH_permutation_B", B_PERM, "_mvn_GPD")
OUT_DIR    <- file.path(BASE_DIR, OUT_SUBDIR)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== CATFISH Fly Male Starvation ===\n")
cat("B_PERM:", B_PERM, "| Mode:", PERM_MODE, "| Tail:", TAIL_MODE, "\n")
cat("Output:", OUT_DIR, "\n\n")

# ==============================================================================
# 1. Load and combine MAGMA gene results
# ==============================================================================

cat("Loading MAGMA gene results...\n")

files <- list.files(
  MAGMA_DIR,
  pattern = "^Male_starvation_fly_.*\\.genes\\.out$",
  full.names = TRUE
)

if (length(files) == 0) {
  # Try alternate naming pattern
  files <- list.files(
    MAGMA_DIR,
    pattern = "\\.genes\\.out$",
    full.names = TRUE
  )
}

if (length(files) == 0) stop("No MAGMA .genes.out files found in: ", MAGMA_DIR)
cat("  Found", length(files), "chromosome files\n")

gene_list <- lapply(files, function(f) {
  x <- read.table(f, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  x
})

genes_all_raw <- do.call(rbind, gene_list)

# Standardise P column
if (!"P" %in% names(genes_all_raw)) {
  names(genes_all_raw)[9] <- "P"
}

# Keep best P per gene
genes_all_raw <- genes_all_raw[order(genes_all_raw$GENE, genes_all_raw$P), ]
genes_all     <- genes_all_raw[!duplicated(genes_all_raw$GENE), ]

if (all(c("CHR", "START") %in% names(genes_all))) {
  genes_all <- genes_all[order(genes_all$CHR, genes_all$START), ]
}

cat("  Unique genes:", nrow(genes_all), "\n\n")

# Save combined gene table for candidate gene scoring
genes_out_file <- file.path(OUT_DIR, paste0(TRAIT, "_combined_genes.tsv"))
write.table(genes_all, genes_out_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved combined genes to:", genes_out_file, "\n\n")

# ==============================================================================
# 2. Adjust gene p-values for length / NSNPS
# ==============================================================================

cat("Adjusting gene p-values for length and NSNPS...\n")

gene_lengths <- read.delim(
  GFF3_PATH,
  stringsAsFactors = FALSE
)

adj_out <- catfish_adjust_gene_p(
  gene_results = genes_all,
  gene_lengths  = gene_lengths,
  gene_col      = "GENE",
  nsnp_col      = "NSNPS",
  p_col         = "P",
  z_col         = "ZSTAT",
  len_gene_col  = "gene_id",
  len_col       = "length"
)

genes_adj <- adj_out[, c("gene_id", "z_adj", "p_adj")]
colnames(genes_adj) <- c("GENE", "ZSTAT", "P")
cat("  Adjusted", nrow(genes_adj), "genes\n\n")

# ==============================================================================
# 3. Load fly pathways (FlyCyc)
# ==============================================================================

cat("Loading FlyCyc pathways...\n")
fly_pw <- read.delim(
  "/Users/nirwantandukar/Documents/Github/MAGCAT/inst/extdata/pathway/Fly_Cyc.tsv",
  stringsAsFactors = FALSE
)
cat("  Pathways loaded:", length(unique(fly_pw$pathway_id)), "\n")
cat("  Pathway genes:", nrow(fly_pw), "\n\n")

# ==============================================================================
# 4. Run CATFISH omnibus pathway analysis
# ==============================================================================

cat("Running CATFISH omnibus (B =", B_PERM, ")...\n")
cat("This may take a while. Grab a coffee ☕\n\n")

omni_results <- catfish_omni2_pathways(
  gene_results              = genes_adj,
  species                   = NULL,
  pathways                  = fly_pw,
  gene_col                  = "GENE",
  p_raw_col                 = "P",
  z_col                     = "ZSTAT",
  tau_grid                  = TAU_GRID,
  omnibus                   = OMNIBUS,
  B_perm                    = B_PERM,
  perm_mode                 = PERM_MODE,
  magma_cor_file            = COR_FILE,
  mvn_marginal              = "uniform",
  mvn_calibrate_components  = FALSE,
  tail_mode                 = TAIL_MODE,
  tail_switch_exceed        = TAIL_SWITCH_EXCEED,
  tail_gpd_k                = TAIL_GPD_K,
  tail_min_B                = TAIL_MIN_B,
  tail_min_tail             = TAIL_MIN_TAIL,
  min_p                     = MIN_P,
  min_genes                 = 2,
  n_threads                 = N_THREADS,
  seed                      = SEED,
  output                    = TRUE,
  out_dir                   = OUT_DIR
)

# ==============================================================================
# 5. Save and summarise
# ==============================================================================

out_csv <- file.path(OUT_DIR, paste0(TRAIT, "_CATFISH_ACAT_mvn_B", B_PERM, "_GPD.csv"))
write.csv(omni_results, out_csv, row.names = FALSE)
cat("\nResults saved to:", out_csv, "\n\n")

bonf <- 0.05 / nrow(omni_results)
omni_col <- if ("omni_p_final" %in% names(omni_results)) "omni_p_final" else "omni_p_mvn"

sig_bonf <- omni_results[omni_results[[omni_col]] < bonf, ]
sig_fdr  <- omni_results[p.adjust(omni_results[[omni_col]], "BH") < 0.05, ]

cat("=== RESULTS SUMMARY ===\n")
cat("Total pathways tested:      ", nrow(omni_results), "\n")
cat("Bonferroni significant:     ", nrow(sig_bonf), "\n")
cat("FDR < 0.05 significant:     ", nrow(sig_fdr), "\n\n")

if (nrow(sig_bonf) > 0) {
  cat("Top Bonferroni-significant pathways:\n")
  top <- sig_bonf[order(sig_bonf[[omni_col]]), ]
  print(top[, c("pathway_id", "pathway_name", "n_genes", omni_col,
                "dominant_component")][seq_len(min(20, nrow(top))), ])
} else {
  cat("Top 20 pathways by p-value:\n")
  top <- omni_results[order(omni_results[[omni_col]]), ]
  print(top[, c("pathway_id", "pathway_name", "n_genes", omni_col)][1:20, ])
}

cat("\n=== DONE: Fly Male Starvation ===\n")
