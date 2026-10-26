#!/usr/bin/env Rscript
## =============================================================================
## 09_component_calibration.R -- per-component calibration.
##
## Everything else in this paper reports the omnibus. This script opens it up:
## each of the five component tests is calibrated against its own null, and the
## analytic and calibrated type I error and lambda are reported per component.
##
## It answers two questions the omnibus hides:
##   1. which components are the ones LD actually breaks (the aggregation tests,
##      Fisher and Stouffer, are the prime suspects), and
##   2. whether a user who knows their pathway's signal architecture can use a
##      single calibrated component in place of the omnibus without losing
##      calibration.
##
##   Rscript 09_component_calibration.R [B_CAL] [B_EVAL]
## =============================================================================
source("00_lib.R")
args   <- commandArgs(trailingOnly = TRUE)
B_CAL  <- if (length(args) >= 1) as.integer(args[[1]]) else 500L
B_EVAL <- if (length(args) >= 2) as.integer(args[[2]]) else 500L
set.seed(20260706)

G     <- 3000L
BLOCK <- 50L
K     <- 150L
POOL  <- 600L
ALPHA <- 0.05
RHO   <- c(independent = 0.0, moderate = 0.3, strong = 0.6)
SIZE  <- c(20L, 50L)
OUT   <- "../results"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("component calibration: K=%d pathways, B_cal=%d B_eval=%d\n",
            K, B_CAL, B_EVAL))

res <- list()
for (rn in names(RHO)) {
  for (m in SIZE) {
    uni   <- make_universe(G, BLOCK, RHO[[rn]])
    paths <- make_pathways(K, m, POOL, G)
    cat(sprintf("\n  LD=%s m=%d\n", rn, m))
    cat("    calibration draws\n")
    arr_cal <- components_array(uni, paths, B_CAL)
    cat("    evaluation draws\n")
    arr_ev  <- components_array(uni, paths, B_EVAL)

    comps <- dimnames(arr_cal)[[3]]
    for (cmp in c(comps, "Omnibus")) {
      if (identical(cmp, "Omnibus")) {
        cal <- omnibus_from_components(arr_cal)
        ev  <- omnibus_from_components(arr_ev)
      } else {
        cal <- arr_cal[, , cmp]
        ev  <- arr_ev[, , cmp]
      }
      ## a component can be NA for some draws (TFisher with no gene below the
      ## smallest tau); those entries propagate as NA and are dropped by type1
      p_marg <- calibrate_cols(ev, cal)
      thr    <- wy_threshold(cal, ALPHA)
      fw     <- fwer_rates(ev, thr, K, ALPHA)

      res[[length(res) + 1L]] <- data.frame(
        ld = rn, rho = RHO[[rn]], m = m, component = cmp,
        frac_na = mean(!is.finite(ev)),
        t1_analytic = type1(exp(ev), ALPHA),
        t1_calibrated = type1(p_marg, ALPHA),
        fwer_wy = fw[["wy"]], fwer_bonferroni = fw[["bonferroni"]],
        lambda_analytic = lambda_gc(as.vector(exp(ev))),
        lambda_calibrated = lambda_gc(as.vector(p_marg)),
        stringsAsFactors = FALSE)
      cat(sprintf("    %-9s analytic=%.4f (lambda %.2f)  calibrated=%.4f  wy=%.3f\n",
                  cmp, type1(exp(ev), ALPHA), lambda_gc(as.vector(exp(ev))),
                  type1(p_marg, ALPHA), fw[["wy"]]))
    }
  }
}
out <- do.call(rbind, res)
write.csv(out, file.path(OUT, "09_component_calibration.csv"), row.names = FALSE)
cat("\nwrote results/09_component_calibration.csv\n")
print(out[out$m == 50, c("ld","component","t1_analytic","t1_calibrated",
                         "lambda_analytic","fwer_wy")], row.names = FALSE)
