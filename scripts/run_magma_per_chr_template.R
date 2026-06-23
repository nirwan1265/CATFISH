#!/usr/bin/env Rscript

# Template: run MAGMA annotation + per-chromosome gene analysis with CATFISH.
# Edit the CONFIG section, then run:
#   Rscript run_magma_per_chr_template.R

suppressPackageStartupMessages({
  library(CATFISH)
})

## -----------------------------------------------------------------------------
## CONFIG
## -----------------------------------------------------------------------------

MAGMA_BIN <- "/full/path/to/magma"

GWAS_FILE <- "/full/path/to/gwas_results.tsv"
GWAS_SEP  <- "\t"   # "," for CSV, "\t" for TSV

BFILE_PREFIX <- "/full/path/to/plink_prefix"

GENE_LOC     <- "/full/path/to/species.genes.loc"
CHR_MAP_PATH <- paste0(GENE_LOC, ".chr_map.tsv")

OUT_PREFIX <- "trait_name"
OUT_ROOT   <- "magma_run"
ANNOT_DIR  <- file.path(OUT_ROOT, "annot")
MAGMA_DIR  <- file.path(OUT_ROOT, "magma_genes_by_chr")

N_TOTAL    <- 1000L
GENE_MODEL <- "multi=snp-wise"

# Examples:
# CHROMS <- 1:10
# CHROMS <- c("1", "2", "3", "4", "5")
# CHROMS <- c("2L", "2R", "3L", "3R", "X", "4")
CHROMS <- 1:10

N_THREADS <- min(length(CHROMS), 10L)
WINDOW_KB <- c(10, 10)

RENAME_COLUMNS <- c(
  CHR    = "CHR",
  SNP    = "SNP",
  POS    = "POS",
  PVALUE = "P"
)

## -----------------------------------------------------------------------------
## CHECKS
## -----------------------------------------------------------------------------

stopifnot(length(WINDOW_KB) == 2L)

if (!file.exists(MAGMA_BIN)) {
  stop("MAGMA binary not found: ", MAGMA_BIN, call. = FALSE)
}
if (!file.exists(GWAS_FILE)) {
  stop("GWAS file not found: ", GWAS_FILE, call. = FALSE)
}
if (!file.exists(paste0(BFILE_PREFIX, ".bed"))) {
  stop("PLINK .bed file not found: ", paste0(BFILE_PREFIX, ".bed"), call. = FALSE)
}
if (!file.exists(paste0(BFILE_PREFIX, ".bim"))) {
  stop("PLINK .bim file not found: ", paste0(BFILE_PREFIX, ".bim"), call. = FALSE)
}
if (!file.exists(paste0(BFILE_PREFIX, ".fam"))) {
  stop("PLINK .fam file not found: ", paste0(BFILE_PREFIX, ".fam"), call. = FALSE)
}
if (!file.exists(GENE_LOC)) {
  stop("Gene location file not found: ", GENE_LOC, call. = FALSE)
}

if (!file.exists(CHR_MAP_PATH)) {
  message("CHR map file not found, continuing without explicit chr_map_path: ", CHR_MAP_PATH)
  CHR_MAP_PATH <- NULL
}

dir.create(OUT_ROOT,  recursive = TRUE, showWarnings = FALSE)
dir.create(ANNOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MAGMA_DIR, recursive = TRUE, showWarnings = FALSE)

options(magma.path = MAGMA_BIN)

cat("MAGMA binary: ", CATFISH::magma_path(), "\n", sep = "")
cat("GWAS file:    ", GWAS_FILE, "\n", sep = "")
cat("BFILE prefix: ", BFILE_PREFIX, "\n", sep = "")
cat("Gene loc:     ", GENE_LOC, "\n", sep = "")
cat("Chromosomes:  ", paste(CHROMS, collapse = ", "), "\n\n", sep = "")

## -----------------------------------------------------------------------------
## STEP 1: MAGMA annotation
## -----------------------------------------------------------------------------

# magma_annotate() exists in the package but is not exported, so call it with :::.
ann <- CATFISH:::magma_annotate(
  stats_file     = GWAS_FILE,
  rename_columns = RENAME_COLUMNS,
  gene_loc       = GENE_LOC,
  chr_map_path   = CHR_MAP_PATH,
  out_prefix     = OUT_PREFIX,
  out_dir        = ANNOT_DIR,
  window         = WINDOW_KB,
  sep            = GWAS_SEP,
  nonhuman       = TRUE
)

cat("Gene annotation file:\n")
cat("  ", ann$gene_annot, "\n\n", sep = "")

## -----------------------------------------------------------------------------
## STEP 2: MAGMA gene analysis per chromosome
## -----------------------------------------------------------------------------

CATFISH::magma_gene(
  bfile          = BFILE_PREFIX,
  gene_annot     = ann$gene_annot,
  stats_file     = GWAS_FILE,
  sep            = GWAS_SEP,
  n_total        = N_TOTAL,
  rename_columns = RENAME_COLUMNS,
  out_prefix     = OUT_PREFIX,
  out_dir        = MAGMA_DIR,
  gene_model     = GENE_MODEL,
  chroms         = CHROMS,
  n_threads      = N_THREADS
)

cat("Per-chromosome MAGMA runs completed.\n\n")

## -----------------------------------------------------------------------------
## STEP 3: Combine per-chromosome .genes.out files
## -----------------------------------------------------------------------------

model_tag <- gsub("-", "_", gsub("=", "_", GENE_MODEL))

out_pattern <- paste0("^", OUT_PREFIX, "_chr.*\\.", model_tag, "\\.genes\\.out$")
out_files <- list.files(MAGMA_DIR, pattern = out_pattern, full.names = TRUE)

if (!length(out_files)) {
  out_pattern <- paste0("^", OUT_PREFIX, "_chr.*\\.genes\\.out$")
  out_files <- list.files(MAGMA_DIR, pattern = out_pattern, full.names = TRUE)
}

if (!length(out_files)) {
  stop("No per-chromosome .genes.out files were found in: ", MAGMA_DIR, call. = FALSE)
}

cat("Found ", length(out_files), " per-chromosome .genes.out files.\n", sep = "")

genes_list <- lapply(out_files, function(path) {
  utils::read.table(
    path,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
})

genes_all_raw <- do.call(rbind, genes_list)

if (!"P" %in% names(genes_all_raw) && ncol(genes_all_raw) >= 9L) {
  names(genes_all_raw)[9] <- "P"
}

if (!all(c("GENE", "P") %in% names(genes_all_raw))) {
  stop("Combined MAGMA output must contain GENE and P columns.", call. = FALSE)
}

genes_all_raw <- genes_all_raw[order(genes_all_raw$GENE, genes_all_raw$P), , drop = FALSE]
genes_all     <- genes_all_raw[!duplicated(genes_all_raw$GENE), , drop = FALSE]

combined_out <- file.path(MAGMA_DIR, paste0(OUT_PREFIX, "_ALLCHR_", model_tag, ".genes.out.tsv"))
utils::write.table(
  genes_all,
  file      = combined_out,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

cat("Combined gene-level results written to:\n")
cat("  ", combined_out, "\n", sep = "")
cat("Unique genes kept: ", nrow(genes_all), "\n\n", sep = "")

## -----------------------------------------------------------------------------
## STEP 4: Optional next step for MVN/CATFISH
## -----------------------------------------------------------------------------

raw_pattern <- paste0("^", OUT_PREFIX, "_chr.*\\.", model_tag, "\\.genes\\.raw$")
raw_files <- list.files(MAGMA_DIR, pattern = raw_pattern, full.names = TRUE)

cat("Raw MAGMA files found: ", length(raw_files), "\n", sep = "")
if (length(raw_files)) {
  cat("You can use these .genes.raw files next to build MAGMA gene correlation pairs.\n")
}
