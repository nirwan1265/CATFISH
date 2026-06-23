#!/usr/bin/env Rscript

# Run MAGMA for one chromosome at a time for Dry_tons_per_acre.
#
# Usage:
#   Rscript run_magma_single_chr_dry_tons_server.R 1
#   Rscript run_magma_single_chr_dry_tons_server.R 10
#   Rscript run_magma_single_chr_dry_tons_server.R 1 3107
#
# This script:
# 1. Picks the GWAS file for one chromosome
# 2. Builds SNP-loc from the chromosome-specific PLINK .bim
# 3. Runs MAGMA annotation
# 3. Runs MAGMA gene analysis for that chromosome's GWAS file
#
# Important:
# - This script is set up to use NMISS via: N_eff = N_TOTAL - n_miss.
# - Pass N_TOTAL as the second command-line argument, or set DEFAULT_N_TOTAL below.
# - The script auto-detects common GWAS column names for SNP / P / NMISS.

suppressPackageStartupMessages({
  library(CATFISH)
})

## -----------------------------------------------------------------------------
## USER SETTINGS
## -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
CHR_TO_RUN <- if (length(args) >= 1L) as.integer(args[[1]]) else 1L
N_TOTAL_ARG <- if (length(args) >= 2L) as.integer(args[[2]]) else NA_integer_

if (is.na(CHR_TO_RUN) || !CHR_TO_RUN %in% 1:10) {
  stop("Pass one chromosome number from 1 to 10, e.g.: Rscript script.R 1", call. = FALSE)
}

MAGMA_INPUT <- "/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/magma_v1.10"

GWAS_DIR <- "/share/maize/ntanduk/landadapt/convert/Dry_tons_per_acre"
GWAS_PATTERN <- paste0(
  "Dry_tons_per_acre_mod_sub_BAP_energy_traits_renamed_for_genotype.",
  "part1_NEW_TERRA_HA_fixed_Genotyped.recalibrated.filtered.snps_only.",
  "Chr%02d.maf01.assoc.txt"
)
# PLINK-style .assoc files are usually whitespace-delimited.
GWAS_SEP <- ""

PLINK_DIR <- "/rsstu/users/r/rrellan/DOE_CAREER/BAP/maf_01_perc/plink_bfiles"
PLINK_PATTERN <- paste0(
  "NEW_TERRA_HA_fixed_Genotyped.recalibrated.filtered.snps_only.",
  "Chr%02d.maf01"
)

GENE_LOC <- "/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/CATFISH/inst/extdata/sorghum.genes.loc"

# Optional fallback if you do not pass N_TOTAL on the command line.
# Example:
#   DEFAULT_N_TOTAL <- 3107L
DEFAULT_N_TOTAL <- NA_integer_

GENE_MODEL <- "multi=snp-wise"
WINDOW_KB  <- c(25, 25)

OUT_ROOT <- file.path(getwd(), "Dry_tons_per_acre_magma_single_chr")

## -----------------------------------------------------------------------------
## HELPERS
## -----------------------------------------------------------------------------

resolve_magma_bin <- function(path) {
  if (!file.exists(path)) {
    stop("MAGMA path does not exist: ", path, call. = FALSE)
  }

  info <- file.info(path)
  if (isTRUE(info$isdir)) {
    cand <- file.path(path, "magma")
    if (file.exists(cand)) return(cand)
    stop("MAGMA directory found, but binary was not found at: ", cand, call. = FALSE)
  }

  path
}

find_first_match <- function(choices, available) {
  hit <- choices[choices %in% available]
  if (length(hit)) hit[[1]] else NA_character_
}

detect_gwas_columns <- function(path, sep = "") {
  dat <- utils::read.table(
    path,
    header = TRUE,
    sep = sep,
    nrows = 5,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = "",
    quote = "\""
  )

  cn <- colnames(dat)

  snp_col <- find_first_match(
    c("SNP", "rs", "RS", "Marker", "marker", "ID", "id"),
    cn
  )
  p_col <- find_first_match(
    c("P", "P.value", "P_VALUE", "p_wald", "P_WALD", "pvalue", "PVAL", "Pr(>|t|)"),
    cn
  )
  nmiss_col <- find_first_match(
    c("NMISS", "n_miss", "nmiss", "N_MISS"),
    cn
  )

  list(
    colnames = cn,
    rename_columns = c(
      SNP    = snp_col,
      PVALUE = p_col,
      NMISS  = nmiss_col
    )
  )
}

write_snp_loc_from_bim <- function(bfile_prefix, out_file) {
  bim <- utils::read.table(
    paste0(bfile_prefix, ".bim"),
    header = FALSE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (ncol(bim) < 4L) {
    stop("BIM file has fewer than 4 columns: ", paste0(bfile_prefix, ".bim"), call. = FALSE)
  }

  snp_loc <- data.frame(
    SNP = bim[[2]],
    CHR = sub("^chr", "", as.character(bim[[1]]), ignore.case = TRUE),
    BP  = bim[[4]],
    stringsAsFactors = FALSE
  )

  utils::write.table(
    snp_loc,
    file = out_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )

  invisible(out_file)
}

sanitize_gene_loc_for_magma <- function(in_file, out_file) {
  x <- utils::read.table(
    in_file,
    header = TRUE,
    sep = "",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = "",
    quote = "\""
  )

  req <- c("GENE", "CHR", "START", "STOP")
  miss <- setdiff(req, names(x))
  if (length(miss)) {
    stop("Gene-loc file is missing required columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  x <- x[, req, drop = FALSE]
  x$GENE  <- as.character(x$GENE)
  x$CHR   <- sub("^chr", "", as.character(x$CHR), ignore.case = TRUE)
  x$START <- as.integer(x$START)
  x$STOP  <- as.integer(x$STOP)

  x <- x[!is.na(x$GENE) & nzchar(x$GENE), , drop = FALSE]
  n_before <- nrow(x)
  x <- x[!duplicated(x$GENE), , drop = FALSE]
  n_after <- nrow(x)

  utils::write.table(
    x,
    file = out_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  invisible(list(out_file = out_file, n_before = n_before, n_after = n_after))
}

chr_label <- sprintf("Chr%02d", CHR_TO_RUN)
chr_tag   <- sprintf("chr%02d", CHR_TO_RUN)

MAGMA_BIN <- resolve_magma_bin(MAGMA_INPUT)
GWAS_FILE <- file.path(GWAS_DIR, sprintf(GWAS_PATTERN, CHR_TO_RUN))
BFILE_PREFIX <- file.path(PLINK_DIR, sprintf(PLINK_PATTERN, CHR_TO_RUN))
N_TOTAL <- if (!is.na(N_TOTAL_ARG)) N_TOTAL_ARG else DEFAULT_N_TOTAL
GWAS_INFO <- detect_gwas_columns(GWAS_FILE, sep = GWAS_SEP)
RENAME_COLUMNS <- GWAS_INFO$rename_columns

OUT_DIR_CHR <- file.path(OUT_ROOT, chr_tag)
ANNOT_DIR   <- file.path(OUT_DIR_CHR, "annot")
MAGMA_DIR   <- file.path(OUT_DIR_CHR, "magma_gene")
OUT_PREFIX  <- paste0("Dry_tons_per_acre_", chr_tag)

## -----------------------------------------------------------------------------
## CHECKS
## -----------------------------------------------------------------------------

if (!file.exists(GWAS_FILE)) {
  stop("GWAS file not found: ", GWAS_FILE, call. = FALSE)
}
if (!file.exists(GENE_LOC)) {
  stop("Gene location file not found: ", GENE_LOC, call. = FALSE)
}
if (!dir.exists(PLINK_DIR)) {
  stop("PLINK directory not found: ", PLINK_DIR, call. = FALSE)
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
if (is.na(RENAME_COLUMNS["SNP"]) || is.na(RENAME_COLUMNS["PVALUE"])) {
  stop(
    "Could not auto-detect SNP/P columns in GWAS file.\n",
    "Columns found: ", paste(GWAS_INFO$colnames, collapse = ", "),
    call. = FALSE
  )
}
if (is.na(RENAME_COLUMNS["NMISS"])) {
  stop(
    "Could not auto-detect NMISS column in GWAS file.\n",
    "Columns found: ", paste(GWAS_INFO$colnames, collapse = ", "),
    call. = FALSE
  )
}
if (is.na(N_TOTAL) || N_TOTAL <= 0) {
  stop(
    "Using NMISS requires total sample size.\n",
    "Run like: Rscript run_magma_single_chr_dry_tons_server.R 1 3107\n",
    "or set DEFAULT_N_TOTAL in the script.",
    call. = FALSE
  )
}

dir.create(ANNOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MAGMA_DIR, recursive = TRUE, showWarnings = FALSE)

options(magma.path = MAGMA_BIN)

cat("Running MAGMA for ", chr_label, "\n", sep = "")
cat("MAGMA binary: ", CATFISH::magma_path(), "\n", sep = "")
cat("GWAS file:    ", GWAS_FILE, "\n", sep = "")
cat("PLINK prefix: ", BFILE_PREFIX, "\n", sep = "")
cat("Total N:      ", N_TOTAL, "\n", sep = "")
cat("Gene loc:     ", GENE_LOC, "\n", sep = "")
cat("GWAS columns: SNP=", RENAME_COLUMNS["SNP"], ", P=", RENAME_COLUMNS["PVALUE"],
    ", NMISS=", RENAME_COLUMNS["NMISS"], "\n", sep = "")
cat("Output dir:   ", OUT_DIR_CHR, "\n\n", sep = "")

## -----------------------------------------------------------------------------
## STEP 1: Build SNP-loc from PLINK .bim
## -----------------------------------------------------------------------------

snp_loc_file <- file.path(ANNOT_DIR, paste0(OUT_PREFIX, ".snp.loc.txt"))
write_snp_loc_from_bim(BFILE_PREFIX, snp_loc_file)

cat("SNP-loc built from BIM:\n")
cat("  ", snp_loc_file, "\n\n", sep = "")

## -----------------------------------------------------------------------------
## STEP 2: Sanitize gene-loc for MAGMA
## -----------------------------------------------------------------------------

gene_loc_clean <- file.path(ANNOT_DIR, paste0(OUT_PREFIX, ".gene_loc.clean.txt"))
gene_loc_info <- sanitize_gene_loc_for_magma(GENE_LOC, gene_loc_clean)

cat("Clean gene-loc written for MAGMA:\n")
cat("  ", gene_loc_clean, "\n", sep = "")
cat("Gene entries kept: ", gene_loc_info$n_after, " of ", gene_loc_info$n_before, "\n\n", sep = "")

## -----------------------------------------------------------------------------
## STEP 3: MAGMA annotation
## -----------------------------------------------------------------------------

annot_prefix <- file.path(ANNOT_DIR, OUT_PREFIX)
annot_args <- c(
  "--annotate",
  "nonhuman",
  paste0("window=", WINDOW_KB[1], ",", WINDOW_KB[2]),
  "--snp-loc", snp_loc_file,
  "--gene-loc", gene_loc_clean,
  "--out", annot_prefix
)

system2(MAGMA_BIN, annot_args)

gene_annot_file <- paste0(annot_prefix, ".genes.annot")
if (!file.exists(gene_annot_file)) {
  stop("MAGMA annotation did not create: ", gene_annot_file, call. = FALSE)
}

cat("Annotation completed.\n")
cat("Gene annotation file: ", gene_annot_file, "\n\n", sep = "")

## -----------------------------------------------------------------------------
## STEP 4: MAGMA gene analysis
## -----------------------------------------------------------------------------

CATFISH::magma_gene(
  bfile          = BFILE_PREFIX,
  gene_annot     = gene_annot_file,
  stats_file     = GWAS_FILE,
  sep            = GWAS_SEP,
  n_total        = N_TOTAL,
  rename_columns = RENAME_COLUMNS,
  out_prefix     = OUT_PREFIX,
  out_dir        = MAGMA_DIR,
  gene_model     = GENE_MODEL,
  n_threads      = 1
)

model_tag <- gsub("[^A-Za-z0-9]+", "_", GENE_MODEL)

cat("MAGMA gene analysis completed.\n")
cat("Expected main output:\n")
cat("  ", file.path(MAGMA_DIR, paste0(OUT_PREFIX, ".", model_tag, ".genes.out")), "\n", sep = "")
cat("Expected raw output:\n")
cat("  ", file.path(MAGMA_DIR, paste0(OUT_PREFIX, ".", model_tag, ".genes.raw")), "\n", sep = "")
