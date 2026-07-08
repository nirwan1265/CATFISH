#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "scripts/run_dry_tons_per_acre_arabidopsis_style.R"
script_dir <- dirname(normalizePath(this_file, mustWork = FALSE))
repo_dir <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

if (!nzchar(Sys.getenv("CATFISH_REPO"))) {
  Sys.setenv(CATFISH_REPO = repo_dir)
}

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  R_DISABLE_OPENMP = "1",
  B_PERM = "10000",
  PERM_MODE = "mvn",
  OMNIBUS = "ACAT",
  SEED = "123",
  N_THREADS = "8",
  MIN_GENES = "2",
  TAU_OPTION = "paper",
  TAU_GRID = "0.1,0.05,0.01",
  TAU_LABEL = "paper_tau_true_empirical_arabstyle",
  MVN_MARGINAL = "uniform",
  MVN_CALIBRATE_COMPONENTS = "true",
  MAKE_PD = "true",
  TAIL_MODE = "empirical",
  TAIL_SWITCH_EXCEED = "10",
  TAIL_GPD_K = "250",
  TAIL_MIN_B = "10000",
  TAIL_MIN_TAIL = "50"
)

sys.source(
  file.path(script_dir, "Usage_Dry_tons_per_acre_BAP.R"),
  envir = globalenv()
)
