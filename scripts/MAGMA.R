#!/usr/bin/env Rscript

# Run MAGMA for one fly chromosome at a time.
#
# Usage:
#   Rscript magma_fly.R 1 172
#   Rscript magma_fly.R 5 172
#
# This script:
# 1. Picks one fly chromosome from 1..5 -> 2L, 2R, 3L, 3R, X
# 2. Builds SNP-loc from genotype.map
# 3. Runs MAGMA annotation
# 4. Runs MAGMA gene analysis for that chromosome's GWAS file
#
# Important:
# - Edit GWAS_FILE manually to First.assoc.txt / Second.assoc.txt / Third.assoc.txt
# - Pass N_TOTAL as the second command-line argument, or set DEFAULT_N_TOTAL below.
# - The script auto-detects common GWAS column names for CHR / SNP / P / NMISS.

suppressPackageStartupMessages({
  library(CATFISH)
})

## -----------------------------------------------------------------------------
## USER SETTINGS
## -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
CHR_TO_RUN <- if (length(args) >= 1L) as.integer(args[[1]]) else 1L
N_TOTAL_ARG <- if (length(args) >= 2L) as.integer(args[[2]]) else NA_integer_

CHR_LEVELS <- c("2L", "2R", "3L", "3R", "X")

if (is.na(CHR_TO_RUN) || !CHR_TO_RUN %in% seq_along(CHR_LEVELS)) {
  stop(
    "Pass one chromosome number from 1 to ", length(CHR_LEVELS),
    ", where 1=2L, 2=2R, 3=3L, 4=3R, 5=X.\n",
    "Example: Rscript magma_fly.R 1 172",
    call. = FALSE
  )
}

MAGMA_INPUT <- "/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/magma_v1.10"

# Change this manually between First / Second / Third as needed.
GWAS_FILE <- "/share/maize/ntanduk/CATFISH/coursship/First.assoc.txt"

# PLINK-style .assoc files are usually whitespace-delimited.
GWAS_SEP <- ""

PLINK_PREFIX <- "/share/maize/ntanduk/CATFISH/Genotype/genotype_DGRP"
PLINK_MAP    <- "/share/maize/ntanduk/CATFISH/Genotype/genotype.map"

GENE_LOC <- system.file("extdata", "fly.genes.loc", package = "CATFISH")

# Optional fallback if you do not pass N_TOTAL on the command line.
DEFAULT_N_TOTAL <- NA_integer_

GENE_MODEL <- "multi=snp-wise"
WINDOW_KB  <- c(25, 25)

OUT_ROOT <- file.path(getwd(), "courtship_magma_by_chr")

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

sanitize_tag <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

make_runtime_chr_map <- function() {
  # Support both fly labels and PLINK-style numeric chromosome codes.
  # DGRP PLINK files commonly use 1,2,3,4,23 for 2L,2R,3L,3R,X.
  data.frame(
    chr_original = c(
      "2L", "2R", "3L", "3R", "4", "X", "Y",
      "1", "2", "3", "4", "5", "23", "24", "25",
      "chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY",
      "chr1", "chr2", "chr3", "chr4", "chr5", "chr23", "chr24", "chr25"
    ),
    chr_magma = c(
      "1", "2", "3", "4", "5", "6", "7",
      "1", "2", "3", "4", "5", "6", "7", "8",
      "1", "2", "3", "4", "5", "6", "7",
      "1", "2", "3", "4", "5", "6", "7", "8"
    ),
    stringsAsFactors = FALSE
  )
}

write_runtime_chr_map <- function(out_file) {
  x <- make_runtime_chr_map()
  utils::write.table(
    x,
    file = out_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  invisible(out_file)
}

map_chr_to_magma <- function(chr_vec, chr_map) {
  key <- setNames(as.character(chr_map$chr_magma), as.character(chr_map$chr_original))
  out <- unname(key[as.character(chr_vec)])
  bad <- unique(as.character(chr_vec)[is.na(out)])
  if (length(bad)) {
    stop("CHR map missing entries for: ", paste(bad, collapse = ", "), call. = FALSE)
  }
  out
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

  chr_col <- find_first_match(
    c("CHR", "chr", "Chromosome", "CHROM"),
    cn
  )
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
  n_col <- find_first_match(
    c("N", "n", "NOBS", "nobs", "sample_size", "SampleSize"),
    cn
  )

  rename_columns <- c(
    CHR    = chr_col,
    SNP    = snp_col,
    PVALUE = p_col
  )

  if (!is.na(n_col) && n_col != "") {
    rename_columns["N"] <- n_col
  } else if (!is.na(nmiss_col) && nmiss_col != "") {
    rename_columns["NMISS"] <- nmiss_col
  }

  list(
    colnames = cn,
    rename_columns = rename_columns
  )
}

write_snp_loc_from_map <- function(map_file, out_file, chr_label, chr_map) {
  map_df <- utils::read.table(
    map_file,
    header = FALSE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (ncol(map_df) < 4L) {
    stop("MAP file has fewer than 4 columns: ", map_file, call. = FALSE)
  }

  map_df <- map_df[map_df[[1]] == chr_label, , drop = FALSE]
  if (!nrow(map_df)) {
    stop("No SNPs found in MAP file for chromosome: ", chr_label, call. = FALSE)
  }

  snp_loc <- data.frame(
    SNP = map_df[[2]],
    CHR = map_chr_to_magma(chr_label, chr_map),
    BP  = as.integer(map_df[[4]]),
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

sanitize_gene_loc_for_magma <- function(in_file, out_file, chr_label, chr_map) {
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

  target_chr <- map_chr_to_magma(chr_label, chr_map)

  x <- x[, req, drop = FALSE]
  x$GENE  <- as.character(x$GENE)
  x$CHR   <- as.character(x$CHR)
  x$START <- as.integer(x$START)
  x$STOP  <- as.integer(x$STOP)

  x <- x[!is.na(x$GENE) & nzchar(x$GENE), , drop = FALSE]
  x <- x[x$CHR %in% target_chr, , drop = FALSE]
  x <- x[!duplicated(x$GENE), , drop = FALSE]

  if (!nrow(x)) {
    stop("No genes left in gene-loc after filtering to chromosome: ", chr_label, call. = FALSE)
  }

  utils::write.table(
    x,
    file = out_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  invisible(list(out_file = out_file, n_after = nrow(x)))
}

## -----------------------------------------------------------------------------
## RESOLVE FILES / LABELS
## -----------------------------------------------------------------------------

chr_label <- CHR_LEVELS[[CHR_TO_RUN]]
chr_tag   <- sanitize_tag(chr_label)
trait_tag <- sanitize_tag(sub("\\.assoc\\.txt$", "", basename(GWAS_FILE)))

MAGMA_BIN <- resolve_magma_bin(MAGMA_INPUT)
N_TOTAL <- if (!is.na(N_TOTAL_ARG)) N_TOTAL_ARG else DEFAULT_N_TOTAL
GWAS_INFO <- detect_gwas_columns(GWAS_FILE, sep = GWAS_SEP)
RENAME_COLUMNS <- GWAS_INFO$rename_columns
OUT_DIR_CHR <- file.path(OUT_ROOT, trait_tag, chr_tag)
ANNOT_DIR   <- file.path(OUT_DIR_CHR, "annot")
MAGMA_DIR   <- file.path(OUT_DIR_CHR, "magma_gene")
OUT_PREFIX  <- paste0(trait_tag, "_", chr_tag)
CHR_MAP_RUNTIME <- file.path(ANNOT_DIR, paste0(OUT_PREFIX, ".chr_map.tsv"))

## -----------------------------------------------------------------------------
## CHECKS
## -----------------------------------------------------------------------------

if (!file.exists(GWAS_FILE)) {
  stop("GWAS file not found: ", GWAS_FILE, call. = FALSE)
}
if (!nzchar(GENE_LOC) || !file.exists(GENE_LOC)) {
  stop("Could not find CATFISH fly gene location file via system.file().", call. = FALSE)
}
if (!file.exists(PLINK_MAP)) {
  stop("PLINK .map file not found: ", PLINK_MAP, call. = FALSE)
}
if (!file.exists(paste0(PLINK_PREFIX, ".bed"))) {
  stop("PLINK .bed file not found: ", paste0(PLINK_PREFIX, ".bed"), call. = FALSE)
}
if (!file.exists(paste0(PLINK_PREFIX, ".bim"))) {
  stop("PLINK .bim file not found: ", paste0(PLINK_PREFIX, ".bim"), call. = FALSE)
}
if (!file.exists(paste0(PLINK_PREFIX, ".fam"))) {
  stop("PLINK .fam file not found: ", paste0(PLINK_PREFIX, ".fam"), call. = FALSE)
}
if (is.na(RENAME_COLUMNS["CHR"]) || is.na(RENAME_COLUMNS["SNP"]) || is.na(RENAME_COLUMNS["PVALUE"])) {
  stop(
    "Could not auto-detect CHR/SNP/P columns in GWAS file.\n",
    "Columns found: ", paste(GWAS_INFO$colnames, collapse = ", "),
    call. = FALSE
  )
}
if (!("N" %in% names(RENAME_COLUMNS)) && !("NMISS" %in% names(RENAME_COLUMNS)) &&
    (is.na(N_TOTAL) || N_TOTAL <= 0)) {
  stop(
    "No per-SNP N column was detected, so pass total sample size as the 2nd argument.\n",
    "Run like: Rscript magma_fly.R 1 172",
    call. = FALSE
  )
}

dir.create(ANNOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MAGMA_DIR, recursive = TRUE, showWarnings = FALSE)

write_runtime_chr_map(CHR_MAP_RUNTIME)
CHR_MAP_DF <- utils::read.delim(CHR_MAP_RUNTIME, check.names = FALSE, stringsAsFactors = FALSE)

options(magma.path = MAGMA_BIN)

cat("Running MAGMA for chromosome ", chr_label, "\n", sep = "")
cat("Trait/GWAS:    ", GWAS_FILE, "\n", sep = "")
cat("MAGMA binary:  ", CATFISH::magma_path(), "\n", sep = "")
cat("PLINK prefix:  ", PLINK_PREFIX, "\n", sep = "")
cat("PLINK map:     ", PLINK_MAP, "\n", sep = "")
cat("Total N:       ", N_TOTAL, "\n", sep = "")
cat("Gene loc:      ", GENE_LOC, "\n", sep = "")
cat("CHR map:       ", CHR_MAP_RUNTIME, "\n", sep = "")
cat(
  "GWAS columns:  CHR=", RENAME_COLUMNS["CHR"],
  ", SNP=", RENAME_COLUMNS["SNP"],
  ", P=", RENAME_COLUMNS["PVALUE"],
  sep = ""
)
if ("N" %in% names(RENAME_COLUMNS)) {
  cat(", N=", RENAME_COLUMNS["N"], sep = "")
}
if ("NMISS" %in% names(RENAME_COLUMNS)) {
  cat(", NMISS=", RENAME_COLUMNS["NMISS"], sep = "")
}
cat("\nOutput dir:    ", OUT_DIR_CHR, "\n\n", sep = "")

## -----------------------------------------------------------------------------
## STEP 1: Build SNP-loc from genotype.map
## -----------------------------------------------------------------------------

snp_loc_file <- file.path(ANNOT_DIR, paste0(OUT_PREFIX, ".snp.loc.txt"))
write_snp_loc_from_map(PLINK_MAP, snp_loc_file, chr_label, CHR_MAP_DF)

cat("SNP-loc built from MAP:\n")
cat("  ", snp_loc_file, "\n\n", sep = "")

## -----------------------------------------------------------------------------
## STEP 2: Sanitize gene-loc for MAGMA
## -----------------------------------------------------------------------------

gene_loc_clean <- file.path(ANNOT_DIR, paste0(OUT_PREFIX, ".gene_loc.clean.txt"))
gene_loc_info <- sanitize_gene_loc_for_magma(GENE_LOC, gene_loc_clean, chr_label, CHR_MAP_DF)

cat("Clean gene-loc written for MAGMA:\n")
cat("  ", gene_loc_clean, "\n", sep = "")
cat("Gene entries kept: ", gene_loc_info$n_after, "\n\n", sep = "")

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
  bfile          = PLINK_PREFIX,
  gene_annot     = gene_annot_file,
  stats_file     = GWAS_FILE,
  sep            = GWAS_SEP,
  n_total        = N_TOTAL,
  rename_columns = RENAME_COLUMNS,
  out_prefix     = OUT_PREFIX,
  out_dir        = MAGMA_DIR,
  gene_model     = GENE_MODEL,
  chr_keep       = chr_label,
  chr_map_path   = CHR_MAP_RUNTIME,
  n_threads      = 1
)

model_tag <- gsub("[^A-Za-z0-9]+", "_", GENE_MODEL)

cat("MAGMA gene analysis completed.\n")
cat("Expected main output:\n")
cat("  ", file.path(MAGMA_DIR, paste0(OUT_PREFIX, ".", model_tag, ".genes.out")), "\n", sep = "")
cat("Expected raw output:\n")
cat("  ", file.path(MAGMA_DIR, paste0(OUT_PREFIX, ".", model_tag, ".genes.raw")), "\n", sep = "")
