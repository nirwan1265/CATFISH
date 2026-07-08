#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "scripts/run_candidate_gene_analysis_dry_tons_paper_tau_false.R"
script_dir <- dirname(normalizePath(this_file, mustWork = FALSE))

Sys.setenv(
  CATFISH_SUBDIR = "CATFISH_permutation_B1000000_mvn_GPD_paper_tau_false",
  OUT_LABEL = "B1000000_GPD_paper_tau_false"
)

sys.source(
  file.path(script_dir, "candidate_gene_analysis_dry_tons_per_acre.R"),
  envir = globalenv()
)
