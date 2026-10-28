#!/usr/bin/env Rscript
## =============================================================================
## 04_wy_vs_mvn.R -- head-to-head: MVN vs Westfall-Young vs Bonferroni.
##
## Two questions:
##   (a) family-wise error under the complete null
##   (b) RANKING: with a known set of true pathways spiked in, how highly does
##       each calibration rank them? (rank is the product for polygenic traits)
##
## The spike-in is at the PATHWAY level -- a subset of pathways get genes drawn
## with a mean shift. This is NOT a simulated GWAS; no genotypes are involved.
##
## THRESHOLD SCALE. The Westfall-Young threshold is computed on the
## MVN-CALIBRATED p-values, not on the raw analytic omnibus. On the analytic
## scale the per-test null is itself badly miscalibrated under LD, so the family
## minimum is extreme because the statistic is wrong rather than because
## pathways overlap, and the Bonferroni/WY ratio is confounded: an earlier
## version of this script reported ratios up to 8e77, which measured the
## analytic statistic's heavy lower tail and nothing else. Calibrating first
## makes each test's null uniform, so the ratio isolates between-pathway
## dependence, which is what it is meant to measure. Both are reported, the
## analytic one explicitly labelled as confounded.
##
## B must exceed NP/alpha (here 300/0.05 = 6000) or the Monte Carlo floor
## 1/(B+1) sits above the Bonferroni threshold and neither can ever reject.
##
##   Rscript 04_wy_vs_mvn.R [NREP] [B]
## =============================================================================
source("00_lib.R")
args <- commandArgs(trailingOnly = TRUE)
NREP <- if (length(args) >= 1) as.integer(args[[1]]) else 200L
B    <- if (length(args) >= 2) as.integer(args[[2]]) else 1000L
B    <- max(B, 20000L)   # see the note on the Monte Carlo floor above
NP     <- 300L     # pathways per replicate
N_TRUE <- 15L      # of which this many are non-null
MU     <- 1.5      # mean shift on the latent z for true pathways
set.seed(20260706)

RHO  <- c(independent = 0.0, moderate = 0.3, strong = 0.6)
SIZE <- c(20L, 50L)
OUT <- "../results"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

res <- list()
for (rn in names(RHO)) for (m in SIZE) {
  U <- safe_chol(R_exch(m, RHO[[rn]]))
  bank <- components_from_zmat(draw_zmat(U, B))[, "Omnibus"]

  rk_mvn <- c(); rk_an <- c(); fwer_wy <- c(); fwer_bf <- c()
  minnull <- numeric(NREP); minnull_mvn <- numeric(NREP)
  for (r in seq_len(NREP)) {
    is_true <- c(rep(TRUE, N_TRUE), rep(FALSE, NP - N_TRUE))
    zz <- draw_zmat(U, NP)
    zz[is_true, ] <- zz[is_true, , drop = FALSE] + MU
    st <- components_from_zmat(zz)[, "Omnibus"]
    pm <- vapply(st, function(s) mvn_calibrate(s, bank), numeric(1))
    ## median rank of the true pathways under each ordering
    rk_an  <- c(rk_an,  stats::median(rank(st, ties.method = "average")[is_true]))
    rk_mvn <- c(rk_mvn, stats::median(rank(pm, ties.method = "average")[is_true]))
    ## family-wise error measured on the NULL pathways only
    minnull[r]     <- min(st[!is_true])
    minnull_mvn[r] <- min(pm[!is_true])
  }
  wy_t     <- stats::quantile(minnull,     0.05, names = FALSE)  # analytic scale
  wy_t_mvn <- stats::quantile(minnull_mvn, 0.05, names = FALSE)  # calibrated
  res[[length(res)+1L]] <- data.frame(
    ld = rn, rho = RHO[[rn]], m = m, n_pathways = NP, n_true = N_TRUE, mu = MU,
    median_rank_analytic = mean(rk_an),
    median_rank_mvn      = mean(rk_mvn),
    random_expectation   = NP / 2,
    B_bank               = B,
    mc_floor             = 1 / (B + 1),
    wy_threshold         = wy_t_mvn,
    wy_threshold_analytic = exp(wy_t),
    bonferroni_threshold = 0.05 / NP,
    fwer_ratio           = (0.05 / NP) / wy_t_mvn,
    fwer_ratio_analytic  = (0.05 / NP) / exp(wy_t),
    stringsAsFactors = FALSE)
  cat(sprintf("  LD=%-11s m=%-3d  median true rank: analytic=%.1f mvn=%.1f (random %.0f)\n",
              rn, m, mean(rk_an), mean(rk_mvn), NP/2))
}
out <- do.call(rbind, res)
write.csv(out, file.path(OUT, "04_wy_vs_mvn.csv"), row.names = FALSE)
cat("\nwrote results/04_wy_vs_mvn.csv\n"); print(out)
