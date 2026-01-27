## =============================================================================
## MAGCAT/CATFISH Workflow: Complete Pipeline for Pathway Enrichment Analysis
## =============================================================================
##
## This script demonstrates the full MAGCAT workflow for three organisms:
##   1. Maize (Zea mays)
##   2. Drosophila (Fly)
##   3. Arabidopsis thaliana
##
## Pipeline Steps:
##   1. Setup and configuration
##   2. Convert GFF3 to gene location file
##   3. Annotate SNPs to genes (MAGMA annotate)
##   4. Run MAGMA gene analysis (per chromosome)
##   5. Combine chromosome results
##   6. (Optional) Adjust p-values for gene length/SNP count
##   7. Generate gene-gene correlation file for MVN resampling
##   8. Load pathway definitions
##   9. Run omnibus pathway test
##
## Key Functions (from R/ directory):
##   - pathway_loaders.R:     magcat_load_pathways(), magcat_pmn_file()
##   - ACAT_wrappers.R:       magcat_acat_pathways(), fix_p_for_acat()
##   - Fisher_wrappers.R:     magcat_fisher_pathways()
##   - minP_wrappers.R:       magcat_minp_pathways()
##   - Stouffer_wrappers.R:   magcat_stoufferZ_pathways()
##   - adaptive_soft_TFisher_wrappers.R: magcat_soft_tfisher_adaptive_pathways()
##   - correlation_genes_wrappers.R: magma_genesraw_to_cor_pairs_banded()
##   - OMNIBUS_test.R:        magcat_omni2_pathways()
##
## =============================================================================


## -----------------------------------------------------------------------------
## STEP 1: Setup
## -----------------------------------------------------------------------------

# Install MAGCAT if needed
# devtools::install_github("nirwan1265/CATFISH")

# Load packages
library(MAGCAT)
library(dplyr)

# Optional packages (for specific features)
# library(TFisher)      # For adaptive soft TFisher
# library(ACAT)         # For ACAT combination
# library(data.table)   # For faster file reading

# Set MAGMA path (download from https://ctg.cncr.nl/software/magma)
options(magma.path = "/path/to/magma")  # UPDATE THIS

# Verify MAGMA is found
magma_path()


## -----------------------------------------------------------------------------
## STEP 2: Define organism-specific paths
## -----------------------------------------------------------------------------

# Choose your organism: "maize", "fly", or "arabidopsis"
ORGANISM <- "maize"  # <<< CHANGE THIS TO SWITCH ORGANISMS

# Define paths based on organism
if (ORGANISM == "maize") {


 # --- MAIZE (CornCyc pathways) ---
 gff_path       <- "/path/to/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3"
 gene_loc_out   <- "inst/extdata/maize.genes.loc"
 bfile          <- "/path/to/maize_genotype"      # PLINK prefix (.bed/.bim/.fam)
 gwas_file      <- "/path/to/maize_gwas.txt"
 gwas_sep       <- "\t"
 n_total        <- 3107                           # GWAS sample size
 chroms         <- 1:10
 gene_regex     <- "^Zm"                          # Gene ID pattern for correlations
 rename_cols    <- c(CHR = "Chr", SNP = "SNP", POS = "Pos", PVALUE = "P.value")
 species_pw     <- "maize"                        # Built-in pathway database

} else if (ORGANISM == "fly") {

 # --- DROSOPHILA (FlyCyc pathways) ---
 gff_path       <- "/path/to/Drosophila_melanogaster.BDGP6.gff3"
 gene_loc_out   <- "inst/extdata/fly.genes.loc"
 bfile          <- "/path/to/dgrp_genotype"
 gwas_file      <- "/path/to/fly_gwas.csv"
 gwas_sep       <- ","
 n_total        <- 166
 chroms         <- c("2L", "2R", "3L", "3R", "4", "X")
 gene_regex     <- "^FBgn"                        # FlyBase gene IDs
 rename_cols    <- c(CHR = "CHR", SNP = "SNP", POS = "Positions", PVALUE = "P-Value")
 species_pw     <- "fly"

} else if (ORGANISM == "arabidopsis") {

 # --- ARABIDOPSIS (AraCyc pathways) ---
 gff_path       <- "/path/to/TAIR10_GFF3_genes.gff"
 gene_loc_out   <- "inst/extdata/at.genes.loc"
 bfile          <- "/path/to/1001genomes"
 gwas_file      <- "/path/to/arabidopsis_gwas.csv"
 gwas_sep       <- ","
 n_total        <- 93
 chroms         <- 1:5
 gene_regex     <- "^AT"                          # TAIR gene IDs
 rename_cols    <- c(CHR = "CHR", SNP = "SNP", POS = "Positions", PVALUE = "P-Value")
 species_pw     <- "arabidopsis"

}

# Output directories
out_prefix  <- paste0(ORGANISM, "_analysis")
annot_dir   <- "annot"
magma_dir   <- paste0("magma_genes_", ORGANISM)
results_dir <- paste0("magcat_results_", ORGANISM)


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
## STEP 4: Annotate SNPs to genes (MAGMA annotate)
## -----------------------------------------------------------------------------

ann <- magma_annotate(
 stats_file     = gwas_file,
 rename_columns = rename_cols,
 gene_loc       = gene_loc_out,
 chr_map_path   = paste0(gene_loc_out, ".chr_map.tsv"),
 out_prefix     = out_prefix,
 out_dir        = annot_dir,
 window         = c(10, 10),  # 10kb upstream/downstream
 sep            = gwas_sep,
 nonhuman       = TRUE
)

cat("Annotation file:", ann$gene_annot, "\n")


## -----------------------------------------------------------------------------
## STEP 5: Run MAGMA gene analysis (per chromosome)
## -----------------------------------------------------------------------------
##
## Running per-chromosome generates .genes.raw files needed for
## gene-gene correlations (MVN resampling).

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

# Adjust p-values (removes confounding from gene size)
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
##
## This step extracts gene-gene correlations from MAGMA's .genes.raw files.
## The correlation file is required for MVN-based resampling to account for
## linkage disequilibrium between nearby genes.
##
## The function magma_genesraw_to_cor_pairs_banded() parses MAGMA's banded
## correlation format and outputs a simple 3-column file: gene1, gene2, r

# Find raw correlation files from MAGMA
raw_pattern <- paste0("^", out_prefix, "_chr.*\\.genes\\.raw$")
chr_raw_files <- list.files(
 path       = magma_dir,
 pattern    = raw_pattern,
 full.names = TRUE
)

cat("Found", length(chr_raw_files), "chromosome .raw files\n")

# Output correlation pairs file
cor_pairs_file <- file.path(magma_dir, paste0("gene_cor_pairs_", ORGANISM, ".txt"))

# Process each chromosome and combine
if (file.exists(cor_pairs_file)) file.remove(cor_pairs_file)

first <- TRUE
for (f in chr_raw_files) {
 cat("Processing:", basename(f), "\n")

 tmp <- tempfile(fileext = ".txt")

 # Extract correlations using organism-specific gene regex
 magma_genesraw_to_cor_pairs_banded(
   genes_raw_file = f,
   out_pairs_file = tmp,
   gene_regex     = gene_regex,
   keep_abs_r_ge  = 0,           # Keep all correlations (filter later if needed)
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

cat("\nCorrelation file written to:", cor_pairs_file, "\n")

# Preview the correlation file
cor_preview <- read.delim(cor_pairs_file, nrows = 5)
print(cor_preview)


## -----------------------------------------------------------------------------
## STEP 9: Load pathway definitions
## -----------------------------------------------------------------------------
##
## Built-in pathway databases (from inst/extdata/pathway/):
##   - "maize"       : CornCyc (Plant Metabolic Network)
##   - "sorghum"     : SorghumBicolorCyc (PMN)
##   - "arabidopsis" : AraCyc (PMN)
##   - "plant"       : PlantCyc general (PMN)
##   - "fly"         : FlyCyc (Drosophila)

pathways <- magcat_load_pathways(species = species_pw)

cat("Loaded", length(unique(pathways$pathway_id)), "pathways\n")
cat("Total pathway-gene mappings:", nrow(pathways), "\n")
head(pathways)

# OR load custom pathways as data.frame:
# pathways <- read.delim("my_pathways.tsv")
# Required columns: pathway_id, gene_id (optional: pathway_name)


## -----------------------------------------------------------------------------
## STEP 10: Run Omnibus Pathway Test
## -----------------------------------------------------------------------------
##
## The omnibus test combines multiple p-value combination methods:
##   - ACAT (Cauchy combination)
##   - Fisher's method
##   - Adaptive soft TFisher
##   - Minimum p-value (minP)
##   - Stouffer's Z (if Z-scores available)
##
## Results are calibrated using MVN resampling to account for gene-gene LD.

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

 # Omnibus combination method
 omnibus        = "ACAT",         # or "minP"

 # MVN resampling for LD-aware calibration
 perm_mode      = "mvn",          # "none", "global", "mvn", or "both"
 B_perm         = 10000L,         # Number of resampling iterations
 magma_cor_file = cor_pairs_file, # Gene-gene correlations from STEP 8
 make_PD        = TRUE,           # Ensure positive-definite correlation matrix
 mvn_marginal   = "uniform",      # Marginal distribution for null p-values

 # Output
 seed           = 123,
 output         = TRUE,
 out_dir        = results_dir
)

cat("\nOmnibus test complete!\n")
head(omni_results)


## -----------------------------------------------------------------------------
## STEP 11: Identify significant pathways
## -----------------------------------------------------------------------------

# The omnibus returns calibrated p-values
# Use omni_p_mvn for MVN-calibrated, or omni_p_final for best available

# Add FDR correction
omni_results$fdr <- p.adjust(omni_results$omni_p_mvn, method = "BH")

# Filter significant pathways
sig_pathways <- omni_results[!is.na(omni_results$fdr) & omni_results$fdr < 0.05, ]
sig_pathways <- sig_pathways[order(sig_pathways$omni_p_mvn), ]

cat("\n=== Significant Pathways (FDR < 0.05) ===\n")
if (nrow(sig_pathways) > 0) {
 print(sig_pathways[, c("pathway_id", "pathway_name", "n_genes", "omni_p_mvn", "fdr")])
} else {
 cat("No pathways significant at FDR < 0.05\n")
}

# Save results
write.csv(omni_results, file.path(results_dir, "omnibus_results.csv"), row.names = FALSE)
cat("\nResults saved to:", file.path(results_dir, "omnibus_results.csv"), "\n")


## =============================================================================
## Quick Reference
## =============================================================================
##
## Gene ID patterns by organism (for gene_regex parameter):
## ---------------------------------------------------------
## Organism      | Pattern   | Example
## --------------|-----------|------------------
## Maize         | "^Zm"     | Zm00001eb000010
## Fly           | "^FBgn"   | FBgn0000003
## Arabidopsis   | "^AT"     | AT1G01010
## Sorghum       | "^SORBI"  | SORBI_3001G000100
##
## Built-in pathway databases:
## ---------------------------
## Species       | Database  | Source
## --------------|-----------|------------------
## maize         | CornCyc   | Plant Metabolic Network
## sorghum       | SorghumCyc| Plant Metabolic Network
## arabidopsis   | AraCyc    | Plant Metabolic Network
## plant         | PlantCyc  | Plant Metabolic Network (general)
## fly           | FlyCyc    | FlyBase
##
## =============================================================================
