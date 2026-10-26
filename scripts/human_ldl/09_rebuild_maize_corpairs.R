#!/usr/bin/env Rscript
# =============================================================================
# 09_rebuild_maize_corpairs.R
#
# The existing magma_gene_cor_pairs_GLM.txt covers CHROMOSOME 1 ONLY (5,365
# genes, IDs Zm00001eb000010-065550). Chromosomes 2-10 -- 85% of maize genes --
# had no gene-gene correlations, so the Westfall-Young null treated them as
# independent.
#
# All ten per-chromosome .genes.raw files exist, so rebuild the full file the
# same way 04_corpairs_once.R does for sorghum.
#
#   Rscript 09_rebuild_maize_corpairs.R [GLM|MLM]
# =============================================================================
suppressPackageStartupMessages({ library(CATFISH) })

args  <- commandArgs(trailingOnly = TRUE)
MODEL <- if (length(args) >= 1) toupper(args[[1]]) else "GLM"

RES <- path.expand("~/Documents/Research/results/MAGMA/MAGCAT")
DIR <- if (MODEL == "GLM") {
  file.path(RES, "magma_genes_by_chr")
} else {
  file.path(RES, "magma_multi_snp_wise_genes_by_chr_N_maize")
}
OUT <- file.path(DIR, sprintf("magma_gene_cor_pairs_%s_ALLCHR.txt", MODEL))

allchr <- file.path(DIR, sprintf("N_maize_%s_ALLCHR.multi_snp_wise.genes.raw", MODEL))
raws <- if (file.exists(allchr)) {
  allchr
} else {
  file.path(DIR, sprintf("N_maize_%s_chr%d.multi_snp_wise.genes.raw", MODEL, 1:10))
}
if (file.exists(allchr)) cat("using the merged ALLCHR .genes.raw\n")
have <- file.exists(raws)
cat(sprintf("%s: %d of %d .genes.raw found\n", MODEL, sum(have), length(raws)))
if (!any(have)) stop("no .genes.raw files under ", DIR, call. = FALSE)

con <- file(OUT, "w"); header_written <- FALSE
for (f in raws[have]) {
  tmp <- tempfile(fileext = ".txt")
  cat("  ", basename(f), "... ")
  magma_genesraw_to_cor_pairs_banded(
    genes_raw_file = f,
    out_pairs_file = tmp,
    gene_regex     = "^Zm",      # maize B73 v5 gene IDs
    keep_abs_r_ge  = 0,
    overwrite      = TRUE,
    verbose        = FALSE
  )
  if (!file.exists(tmp)) { cat("no output\n"); next }
  x <- readLines(tmp, warn = FALSE)
  if (length(x) <= 1) { cat("empty\n"); next }
  if (header_written) x <- x[-1] else header_written <- TRUE
  writeLines(x, con)
  cat(length(x), "pairs\n")
  unlink(tmp)
}
close(con)

n <- length(readLines(OUT, warn = FALSE)) - 1L
cat(sprintf("\nwrote %s\n  %d pairs\n", OUT, n))
cat("Point the maize branch of 08_westfall_young.R at this file.\n")
