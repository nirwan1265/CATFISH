repo <- "/Users/nirwantandukar/Documents/Github/MAGCAT"
run_script <- file.path(repo, "scripts/simulations/04_catfish_perm.R")
agg_script <- file.path(repo, "scripts/simulations/05_aggregate_and_plot.R")

out_root <- file.path(repo, "Final/final_main_fig_false5_gpd_seq")
per_dir <- file.path(out_root, "per_perm")

dir.create(per_dir, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(
  MAGMA_OUT = "/Users/nirwantandukar/Documents/Research/data/BAP/simulations/permutation_runs/magma",
  CORPAIRS = "/Users/nirwantandukar/Documents/Research/data/BAP/simulations/permutation_runs/gene_cor_pairs_CONSTANT.txt",
  CATFISH_OUT = per_dir,
  CATFISH_REPO = repo,
  GENE_MODEL = "multi_snp_wise",
  CHROMS = "01 02 03 04 05 06 07 08 09 10",
  MASTER_SEED = "1",
  MVN_CALIBRATE_COMPONENTS = "false",
  PATCH_EMP_NULL_LOWER = "false",
  CATFISH_THREADS = "2",
  TAU_GRID = "0.01,0.05,0.5,1",
  CATFISH_B_PERM = "2000",
  TAIL_MODE = "hybrid_gpd",
  TAIL_SWITCH_EXCEED = "10",
  TAIL_GPD_K = "250",
  TAIL_MIN_B = "1000",
  TAIL_MIN_TAIL = "50",
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  R_DISABLE_OPENMP = "1"
)

orig_commandArgs <- base::commandArgs

run_one <- function(b) {
  assign("commandArgs", function(trailingOnly = FALSE) {
    if (isTRUE(trailingOnly)) return(as.character(b))
    orig_commandArgs(trailingOnly = FALSE)
  }, envir = .GlobalEnv)
  on.exit(rm("commandArgs", envir = .GlobalEnv), add = TRUE)

  cat(sprintf("\n===== Running permutation %03d =====\n", b))
  sys.source(run_script, envir = new.env(parent = .GlobalEnv))
  cat(sprintf("===== Finished permutation %03d =====\n", b))
}

for (b in c(1L, 5L, 10L, 15L, 20L)) run_one(b)

agg_out <- system2("Rscript", c(agg_script, per_dir, out_root), stdout = TRUE, stderr = TRUE)
cat(paste(agg_out, collapse = "\n"), "\n")
