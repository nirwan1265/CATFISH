#!/usr/bin/env Rscript

Sys.setenv(
  CATFISH_SUBDIR = "CATFISH_permutation_B10000_mvn_GPD_paper_tau_false",
  OUT_LABEL = "B10000_GPD_paper_tau_false"
)

sys.source(
  "/Users/nirwantandukar/Documents/Github/MAGCAT/scripts/candidate_gene_analysis_dry_tons_per_acre.R",
  envir = globalenv()
)
