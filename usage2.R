## =============================================================================
## MAGCAT Workflow: Complete Pipeline for Pathway Enrichment Analysis
## =============================================================================
##
## This script demonstrates the full MAGCAT workflow for three organisms:
##   1. Maize (Zea mays)
##   2. Drosophila (Fly)
##   3. Arabidopsis thaliana
##
## Steps:
##   1. Setup and configuration
##   2. Convert GFF3 to gene location file
##   3. Annotate SNPs to genes (MAGMA annotate)
##   4. Run MAGMA gene analysis
##   5. Combine chromosome results
##   6. (Optional) Adjust p-values for gene length/SNP count
##   7. Generate gene-gene correlation file for MVN resampling
##   8. Run pathway enrichment analysis (ACAT, Fisher, minP, etc.)
##   9. Run omnibus test combining all methods
##
## =============================================================================

## -----------------------------------------------------------------------------
## STEP 1: Setup
## -----------------------------------------------------------------------------

# Load packages
devtools::document()
devtools::load_all()
library(metapro)
library(sumFREGAT)
library(metap)
library(rtracklayer)
library(ggplot2)
library(TFisher)
library(data.table)
library(MAGCAT)
library(dplyr)

# Set MAGMA path (download from https://ctg.cncr.nl/software/magma)
options(magma.path = "/path/to/magma")  # UPDATE THIS

# Verify MAGMA is found
magma_path()

## -----------------------------------------------------------------------------
## STEP 2: Define organism-specific paths
## -----------------------------------------------------------------------------

# Choose your organism: "maize", "fly", or "arabidopsis"
ORGANISM <- "maize"  # Change this to switch organisms

# Define paths based on organism
if (ORGANISM == "maize") {
  # --- MAIZE ---
  gff_path       <- "/path/to/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3"
  gene_loc_out   <- "inst/extdata/maize.genes.loc"
  bfile          <- "/path/to/maize_genotype"  # PLINK prefix
  gwas_file      <- "/path/to/maize_gwas.txt"
  gwas_sep       <- "\t"
  n_total        <- 3107  # Sample size
  chroms         <- 1:10
  gene_regex     <- "^Zm"  # For correlation parsing
  rename_cols    <- c(CHR = "Chr", SNP = "SNP", POS = "Pos", PVALUE = "P.value")
  species_pw     <- "maize"  # Built-in pathways

} else if (ORGANISM == "fly") {
  # --- DROSOPHILA ---
  gff_path       <- "/path/to/Drosophila_melanogaster.BDGP6.gff3"
  gene_loc_out   <- "inst/extdata/fly.genes.loc"
  bfile          <- "/path/to/dgrp_genotype"
  gwas_file      <- "/path/to/fly_gwas.csv"
  gwas_sep       <- ","
  n_total        <- 166
  chroms         <- c("2L", "2R", "3L", "3R", "4", "X")
  gene_regex     <- "^FBgn"
  rename_cols    <- c(CHR = "CHR", SNP = "SNP", POS = "Positions", PVALUE = "P-Value")
  species_pw     <- "fly"  # Built-in pathways

} else if (ORGANISM == "arabidopsis") {
  # --- ARABIDOPSIS ---
  gff_path       <- "/path/to/TAIR10_GFF3_genes.gff"
  gene_loc_out   <- "inst/extdata/at.genes.loc"
  bfile          <- "/path/to/1001genomes"
  gwas_file      <- "/path/to/arabidopsis_gwas.csv"
  gwas_sep       <- ","
  n_total        <- 93
  chroms         <- 1:5
  gene_regex     <- "^AT"
  rename_cols    <- c(CHR = "CHR", SNP = "SNP", POS = "Positions", PVALUE = "P-Value")
  species_pw     <- "arabidopsis"  # Built-in pathways
}

# Output directories
out_prefix     <- paste0(ORGANISM, "_analysis")
annot_dir      <- "annot"
magma_dir      <- paste0("magma_genes_", ORGANISM)
results_dir    <- paste0("magcat_results_", ORGANISM)


## -----------------------------------------------------------------------------
## STEP 3: Convert GFF3 to MAGMA gene location file
## -----------------------------------------------------------------------------

loc <- gff3_to_geneloc(
  gff        = gff_path,
  out        = gene_loc_out,
  strict_chr = FALSE
)

cat("Gene location file written to:", gene_loc_out, "\n")


## -----------------------------------------------------------------------------
## STEP 4: Annotate SNPs to genes
## -----------------------------------------------------------------------------

ann <- magma_annotate(
  stats_file     = gwas_file,
  rename_columns = rename_cols,
  gene_loc       = gene_loc_out,
  chr_map_path   = paste0(gene_loc_out, ".chr_map.tsv"),  # If non-numeric chromosomes
  out_prefix     = out_prefix,
  out_dir        = annot_dir,
  window         = c(10, 10),  # 10kb upstream/downstream
  sep            = gwas_sep,
  nonhuman       = TRUE
)

cat("Annotation file:", ann$gene_annot, "\n")


## -----------------------------------------------------------------------------
## STEP 5: Run MAGMA gene analysis (per chromosome for correlation calculation)
## -----------------------------------------------------------------------------

magma_gene(
  bfile          = bfile,
  gene_annot     = ann$gene_annot,
  stats_file     = gwas_file,
  sep            = gwas_sep,
  n_total        = n_total,
  rename_columns = rename_cols,
  out_prefix     = out_prefix,
  out_dir        = magma_dir,
  gene_model     = "multi=snp-wise",
  chroms         = chroms,
  n_threads      = length(chroms)
)


## -----------------------------------------------------------------------------
## STEP 6: Combine chromosome results
## -----------------------------------------------------------------------------

# Find all chromosome output files
chr_pattern <- paste0("^", out_prefix, "_chr.*\\.multi_snp_wise\\.genes\\.out$")
files <- list.files(
  path       = magma_dir,
  pattern    = chr_pattern,
  full.names = TRUE
)

# Read and combine
gene_list <- lapply(files, function(f) {
  utils::read.table(f, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
})
genes_all_raw <- do.call(rbind, gene_list)

# Rename P column if needed (MAGMA uses different names based on model)
if (!"P" %in% names(genes_all_raw) && ncol(genes_all_raw) >= 9) {
  colnames(genes_all_raw)[9] <- "P"
}

# Deduplicate: keep best p-value per gene
genes_all_raw <- genes_all_raw[order(genes_all_raw$GENE, genes_all_raw$P), ]
genes_all <- genes_all_raw[!duplicated(genes_all_raw$GENE), ]

cat("Total genes:", nrow(genes_all), "\n")
head(genes_all)


## -----------------------------------------------------------------------------
## STEP 7: (Optional) Adjust p-values for gene length and SNP count
## -----------------------------------------------------------------------------

# Get gene lengths from GFF3
gene_lengths <- get_gene_lengths(
  gff3_file  = gff_path,
  output     = TRUE,
  output_dir = "inst/extdata",
  file_name  = paste0(ORGANISM, "_gene_lengths.tsv")
)

# Adjust p-values
adj_out <- magcat_adjust_gene_p(
  gene_results = genes_all,
  gene_lengths = gene_lengths,
  gene_col     = "GENE",
  nsnp_col     = "NSNPS",
  p_col        = "P",
  z_col        = "ZSTAT",
  len_gene_col = "gene_id",
  len_col      = "length"
)

# Create adjusted gene results
genes_adj <- data.frame(
  GENE  = adj_out$gene_id,
  ZSTAT = adj_out$z_adj,
  P     = adj_out$p_adj,
  stringsAsFactors = FALSE
)

head(genes_adj)


## -----------------------------------------------------------------------------
## STEP 8: Generate gene-gene correlation file (for MVN resampling)
## -----------------------------------------------------------------------------

# Find raw correlation files from MAGMA
raw_pattern <- paste0("^", out_prefix, "_chr.*\\.genes\\.raw$")
chr_raw_files <- list.files(
  path       = magma_dir,
  pattern    = raw_pattern,
  full.names = TRUE
)

# Output correlation pairs file
cor_pairs_file <- file.path(magma_dir, paste0("gene_cor_pairs_", ORGANISM, ".txt"))

# Process each chromosome and combine
if (file.exists(cor_pairs_file)) file.remove(cor_pairs_file)

first <- TRUE
for (f in chr_raw_files) {
  tmp <- tempfile(fileext = ".txt")

  magma_genesraw_to_cor_pairs_banded(
    genes_raw_file = f,
    out_pairs_file = tmp,
    gene_regex     = gene_regex,  # Organism-specific gene ID pattern
    keep_abs_r_ge  = 0,
    overwrite      = TRUE,
    verbose        = FALSE
  )

  x <- readLines(tmp, warn = FALSE)
  if (!length(x)) next
  if (!first) x <- x[-1]  # Drop header after first file

  cat(paste(x, collapse = "\n"), "\n",
      file = cor_pairs_file, append = !first)
  first <- FALSE
}

cat("Correlation file written to:", cor_pairs_file, "\n")


## -----------------------------------------------------------------------------
## STEP 9: Load pathways
## -----------------------------------------------------------------------------

# Use built-in pathways
pathways <- magcat_load_pathways(species = species_pw)
head(pathways)

# OR load custom pathways as data.frame:
# pathways <- read.delim("my_pathways.tsv")  # Must have: pathway_id, gene_id


## -----------------------------------------------------------------------------
## STEP 10: Run individual pathway tests (optional - omnibus does all)
## -----------------------------------------------------------------------------

# ACAT
acat_res <- magcat_acat_pathways(
  gene_results = genes_adj,
  species      = species_pw,
  gene_col     = "GENE",
  p_col        = "P",
  output       = TRUE,
  out_dir      = results_dir
)

# Fisher
fisher_res <- magcat_fisher_pathways(
  gene_results = genes_adj,
  species      = species_pw,
  gene_col     = "GENE",
  p_col        = "P",
  output       = TRUE,
  out_dir      = results_dir
)

# Stouffer (requires Z-scores)
stouffer_res <- magcat_stoufferZ_pathways(
  gene_results = genes_adj,
  species      = species_pw,
  gene_col     = "GENE",
  z_col        = "ZSTAT",
  alternative  = "greater",
  output       = TRUE,
  out_dir      = results_dir
)

# MinP
minp_res <- magcat_minp_pathways(
  gene_results = genes_adj,
  species      = species_pw,
  gene_col     = "GENE",
  p_col        = "P",
  output       = TRUE,
  out_dir      = results_dir
)


## -----------------------------------------------------------------------------
## STEP 11: Run Omnibus test (combines all methods with resampling)
## -----------------------------------------------------------------------------

omni_results <- magcat_omni2_pathways(
  gene_results   = genes_adj,
  species        = species_pw,
  gene_col       = "GENE",
  p_raw_col      = "P",
  z_col          = "ZSTAT",

  # Method parameters
  tau_grid       = c(0.1, 0.05, 0.02, 0.01, 0.005, 0.001),
  min_p          = 1e-15,
  do_fix         = TRUE,
  stouffer_alternative = "greater",

  # Omnibus combination
  omnibus        = "ACAT",  # or "minP"

  # Resampling for calibration
  perm_mode      = "mvn",   # "none", "global", "mvn", or "both"
  B_perm         = 10000L,
  magma_cor_file = cor_pairs_file,
  make_PD        = TRUE,
  mvn_marginal   = "uniform",

  # Output
  seed           = 123,
  output         = TRUE,
  out_dir        = results_dir
)

head(omni_results)


## -----------------------------------------------------------------------------
## STEP 12: Identify significant pathways
## -----------------------------------------------------------------------------

# Add FDR correction
omni_results$fdr <- p.adjust(omni_results$omni_p_final, method = "BH")

# Filter significant pathways
sig_pathways <- omni_results[omni_results$fdr < 0.05, ]
sig_pathways <- sig_pathways[order(sig_pathways$omni_p_final), ]

cat("\n=== Significant Pathways (FDR < 0.05) ===\n")
print(sig_pathways[, c("pathway_id", "pathway_name", "n_genes", "omni_p_final", "fdr")])


## =============================================================================
## Quick reference: Gene regex patterns by organism
## =============================================================================
##
## Organism      | Gene Regex  | Example Gene ID
## --------------|-------------|----------------
## Maize         | "^Zm"       | Zm00001eb000010
## Fly           | "^FBgn"     | FBgn0000003
## Arabidopsis   | "^AT"       | AT1G01010
## Sorghum       | "^SORBI"    | SORBI_3001G000100
## Rice          | "^Os"       | Os01g0100100
##
## =============================================================================

?MAGCAT::magcat_acat_pathways()
