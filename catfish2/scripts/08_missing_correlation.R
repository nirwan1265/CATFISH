#!/usr/bin/env Rscript
## =============================================================================
## 08_missing_correlation.R -- robustness to an INCOMPLETE gene-gene
## correlation file.
##
## WHY THIS EXISTS
##   Section 2.3 builds every correlation matrix by filling MAGMA's reported
##   gene pairs into an identity matrix, so any pair MAGMA does not report is
##   treated as independent. The Limitations section says this zero-imputation
##   is inherited identically by the marginal and the family-wise readings.
##   Until now that claim was asserted, not measured. This script measures it.
##
## DESIGN
##   Data are generated from the TRUE correlation structure. Calibration -- the
##   null bank, the marginal p-values and the Westfall-Young threshold -- is
##   built from a DEGRADED matrix in which a fraction `miss` of the off-diagonal
##   entries have been zeroed, exactly as an incomplete cor_pairs file would do.
##   The gap between the two is the error the zero-imputation causes.
##
##   miss = 0 is the matched case and serves as the baseline; miss = 1 is full
##   independence, i.e. what the analytic calibration assumes.
##
##   Rscript 08_missing_correlation.R [B_CAL] [B_EVAL]
## =============================================================================
source("00_lib.R")
args   <- commandArgs(trailingOnly = TRUE)
B_CAL  <- if (length(args) >= 1) as.integer(args[[1]]) else 400L
B_EVAL <- if (length(args) >= 2) as.integer(args[[2]]) else 400L
set.seed(20260706)

G     <- 3000L
BLOCK <- 50L
K     <- 150L
POOL  <- 600L
M     <- 50L
ALPHA <- 0.05
RHO   <- c(moderate = 0.3, strong = 0.6)
MISS  <- c(0, 0.25, 0.5, 0.75, 1)
OUT   <- "../results"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## Zero a random fraction of the off-diagonal entries, symmetrically.
degrade_R <- function(R, miss) {
  if (miss <= 0) return(R)
  idx <- which(upper.tri(R), arr.ind = TRUE)
  k <- round(miss * nrow(idx))
  if (k <= 0) return(R)
  sel <- idx[sample(nrow(idx), k), , drop = FALSE]
  R[sel] <- 0
  R[sel[, c(2, 1), drop = FALSE]] <- 0
  R
}

cat(sprintf("missing-correlation robustness: K=%d pathways of %d genes, B_cal=%d B_eval=%d\n",
            K, M, B_CAL, B_EVAL))

res <- list()
for (rn in names(RHO)) {
  rho   <- RHO[[rn]]
  Rtrue <- R_exch(BLOCK, rho)
  uni_true <- list(U = safe_chol(Rtrue), nblock = G %/% BLOCK,
                   G = G, block = BLOCK, rho = rho)
  paths <- make_pathways(K, M, POOL, G)
  cat(sprintf("\n  LD=%s : %.1f pathways per gene\n", rn,
              mean_memberships(paths, G)))

  ## evaluation draws come from the TRUE structure and are reused across all
  ## missingness levels -- only the calibration side changes.
  cat("    evaluation draws (true R)\n")
  arr_ev <- components_array(uni_true, paths, B_EVAL)
  ev_om  <- omnibus_from_components(arr_ev)

  for (miss in MISS) {
    Rd <- degrade_R(Rtrue, miss)
    uni_deg <- list(U = safe_chol(Rd), nblock = G %/% BLOCK,
                    G = G, block = BLOCK, rho = rho)
    cat(sprintf("    calibration draws (%.0f%% of pairs dropped)\n", 100 * miss))
    arr_cal <- components_array(uni_deg, paths, B_CAL, verbose = FALSE)
    cal_om  <- omnibus_from_components(arr_cal)

    p_marg <- calibrate_cols(ev_om, cal_om)
    thr    <- wy_threshold(cal_om, ALPHA)
    fw     <- fwer_rates(ev_om, thr, K, ALPHA)

    res[[length(res) + 1L]] <- data.frame(
      ld = rn, rho = rho, m = M, miss = miss,
      mean_abs_r_kept = mean(abs(Rd[upper.tri(Rd)])),
      t1_analytic = type1(exp(ev_om), ALPHA),
      t1_mvn      = type1(p_marg, ALPHA),
      fwer_wy         = fw[["wy"]],
      fwer_bonferroni = fw[["bonferroni"]],
      lambda_analytic = lambda_gc(as.vector(exp(ev_om))),
      lambda_mvn      = lambda_gc(as.vector(p_marg)),
      wy_threshold = exp(thr),
      stringsAsFactors = FALSE)
    cat(sprintf("      miss=%4.0f%%  marginal mvn=%.4f (lambda %.2f) | FWER wy=%.3f bonf=%.3f\n",
                100 * miss, type1(p_marg, ALPHA), lambda_gc(as.vector(p_marg)),
                fw[["wy"]], fw[["bonferroni"]]))
  }
}
out <- do.call(rbind, res)
write.csv(out, file.path(OUT, "08_missing_correlation.csv"), row.names = FALSE)
cat("\nwrote results/08_missing_correlation.csv\n")
print(out[, c("ld","miss","t1_mvn","fwer_wy","lambda_mvn")], row.names = FALSE)
