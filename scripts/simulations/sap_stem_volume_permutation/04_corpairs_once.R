#!/usr/bin/env Rscript

suppressMessages({ library(CATFISH) })
getcfg <- function(k){ v <- Sys.getenv(k); if (v == "") stop(paste("config missing:", k)); v }

MAGMA_OUT <- getcfg("MAGMA_OUT")
CORPAIRS <- getcfg("CORPAIRS")
GENE_MODEL <- getcfg("GENE_MODEL")
CHROMS <- strsplit(trimws(getcfg("CHROMS")), "\\s+")[[1]]
model_tag <- gsub("[^A-Za-z0-9]+", "_", GENE_MODEL)

raw_for <- function(b) file.path(
  MAGMA_OUT, sprintf("perm_%03d", b), sprintf("chr%s", CHROMS),
  sprintf("perm_%03d_chr%s.%s.genes.raw", b, CHROMS, model_tag)
)
src <- NULL
for (b in 1:100) {
  f <- raw_for(b)
  if (all(file.exists(f))) {
    src <- f
    cat(sprintf("[corpairs] using perm_%03d\n", b))
    break
  }
}
if (is.null(src)) stop("no permutation has all .genes.raw files -- check MAGMA_OUT")

if (file.exists(CORPAIRS)) file.remove(CORPAIRS)
con <- file(CORPAIRS, open = "wt")
header_written <- FALSE
for (f in src) {
  tmp <- tempfile(fileext = ".txt")
  magma_genesraw_to_cor_pairs_banded(
    genes_raw_file = f,
    out_pairs_file = tmp,
    gene_regex = "^SORBI",
    keep_abs_r_ge = 0,
    overwrite = TRUE,
    verbose = TRUE
  )
  if (!file.exists(tmp)) next
  x <- readLines(tmp, warn = FALSE)
  if (length(x) <= 1) next
  if (header_written) x <- x[-1] else header_written <- TRUE
  writeLines(x, con)
}
close(con)
n <- length(readLines(CORPAIRS, warn = FALSE))
cat(sprintf("[corpairs] wrote %d lines -> %s\n", n, CORPAIRS))
if (n <= 1) stop("CORPAIRS is empty -- extraction failed.")
