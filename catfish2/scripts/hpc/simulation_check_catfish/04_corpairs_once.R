#!/usr/bin/env Rscript
## =============================================================================
## 04_corpairs_once.R -- gene-gene correlations for CATFISH, from the first
## replicate whose merged .genes.raw exists. MAGMA's gene correlations depend on
## the panel's LD and the annotation only, not on the simulated statistics, so
## one file serves every replicate (as bap_perm/04_corpairs_once.R). Same
## settings as the lipid runs (17_prep_trait_catfish.R): Entrez IDs, all pairs.
## Writes IN/sim_cor_pairs.txt (shared by pilot and full).
## =============================================================================
suppressPackageStartupMessages(library(CATFISH))
env <- function(k) { v <- Sys.getenv(k); if (!nzchar(v)) stop(k, " is not set"); v }
RUN <- env("RUN"); OUTF <- file.path(env("IN"), "sim_cor_pairs.txt")
raws <- Sys.glob(file.path(RUN, "magma", "rep_*", "rep_*.genes.raw"))
raws <- raws[!grepl("batch", raws)]
if (!length(raws)) stop("no merged .genes.raw under ", RUN, "/magma: run step 03b")
cat("[corpairs] using", raws[1], "\n")
magma_genesraw_to_cor_pairs_banded(genes_raw_file = raws[1], out_pairs_file = OUTF,
                                   gene_regex = "^[0-9]", keep_abs_r_ge = 0,
                                   overwrite = TRUE, verbose = TRUE)
n <- length(readLines(OUTF))
cat(sprintf("[corpairs] %d lines -> %s\n", n, OUTF))
if (n <= 1) stop("no gene pairs written; do not run CATFISH")
