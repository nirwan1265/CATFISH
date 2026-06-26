args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("need permutation index")

repo <- "/Users/nirwantandukar/Documents/Github/MAGCAT"
out_root <- file.path(repo, "Final/final_main_fig_false5_gpd")
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

sys.source(file.path(repo, "scripts/simulations/04_catfish_perm.R"), envir = globalenv())
