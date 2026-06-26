#!/usr/bin/env Rscript

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
  TAU_GRID = "0.01,0.05,0.5,1",
  TAU_LABEL = "paper_tau_false",
  MVN_MARGINAL = "uniform",
  MVN_CALIBRATE_COMPONENTS = "false",
  MAKE_PD = "true",
  TAIL_MODE = "hybrid_gpd",
  TAIL_SWITCH_EXCEED = "10",
  TAIL_GPD_K = "250",
  TAIL_MIN_B = "10000",
  TAIL_MIN_TAIL = "50"
)

sys.source(
  "/Users/nirwantandukar/Documents/Github/MAGCAT/scripts/Usage_Dry_tons_per_acre_BAP.R",
  envir = globalenv()
)
