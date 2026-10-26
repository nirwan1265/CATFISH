#!/usr/bin/env Rscript
## =============================================================================
## 02_breakage_stress.R -- component-breakage stress test, measured on BOTH the
## per-pathway (marginal) calibration AND the Westfall-Young family-wise rate.
##
## WHAT CHANGED AND WHY
##   1. Pathways are now overlapping subsets of ONE shared gene universe with
##      block LD, so a single draw moves every pathway containing a gene. The
##      earlier version drew each pathway independently, which left no
##      between-pathway dependence for Westfall-Young to control.
##   2. The headline number is the WY family-wise error rate, not only the
##      marginal type I error. Bonferroni (on the analytic p) is the competitor.
##   3. Corruption direction. For 0 < p < 1:
##        gamma < 1  =>  p^gamma > p  =>  CONSERVATIVE
##        gamma > 1  =>  p^gamma < p  =>  ANTI-CONSERVATIVE
##      The anti-conservative direction is the one that creates false positives,
##      so both are run. An earlier version used gamma = 0.7 only and called it
##      anti-conservative, which was wrong.
##
## The component log p-values are evaluated ONCE per setting and reused across
## every (gamma, k, combination) variant, since gamma and drop act on the
## components after they are computed. Recomputing them per variant would be
## roughly 25x more work.
##
##   Rscript 02_breakage_stress.R [B_CAL] [B_EVAL]
## =============================================================================
source("00_lib.R")
args   <- commandArgs(trailingOnly = TRUE)
B_CAL  <- if (length(args) >= 1) as.integer(args[[1]]) else 500L
B_EVAL <- if (length(args) >= 2) as.integer(args[[2]]) else 500L
set.seed(20260706)

G      <- 3000L     # genes in the universe
BLOCK  <- 50L       # LD block size
K      <- 150L      # pathways tested per setting
POOL   <- 600L      # pathways are drawn from the first POOL genes
ALPHA  <- 0.05
GAMMAS <- c(conservative = 0.7, anticonservative = 1.4)
MAX_COMBOS <- 3L
COMPONENTS <- c("ACAT","Fisher","Stouffer","minP","TFisher")
RHO  <- c(independent = 0.0, moderate = 0.3, strong = 0.6)
SIZE <- c(20L, 50L)
OUT  <- "../results"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("breakage stress: K=%d pathways, pool=%d of G=%d, B_cal=%d B_eval=%d\n",
            K, POOL, G, B_CAL, B_EVAL))

res <- list()
for (rn in names(RHO)) {
  for (m in SIZE) {
    uni   <- make_universe(G, BLOCK, RHO[[rn]])
    paths <- make_pathways(K, m, POOL, G)
    memb  <- mean_memberships(paths, G)
    cat(sprintf("\n  LD=%s m=%d : %.1f pathways per gene\n", rn, m, memb))

    cat("    calibration draws\n")
    arr_cal <- components_array(uni, paths, B_CAL)
    cat("    evaluation draws\n")
    arr_ev  <- components_array(uni, paths, B_EVAL)

    for (gi in seq_along(GAMMAS)) {
      gval <- GAMMAS[[gi]]; glab <- names(GAMMAS)[gi]
      for (k in 0:4) {
        combos <- if (k == 0L) list(character(0)) else
                  utils::combn(COMPONENTS, k, simplify = FALSE)
        if (length(combos) > MAX_COMBOS) combos <- combos[seq_len(MAX_COMBOS)]
        acc <- list(t1a = numeric(0), t1m = numeric(0),
                    wy = numeric(0), bonf = numeric(0),
                    lam = numeric(0), lamm = numeric(0))
        for (br in combos) {
          g <- if (length(br)) {
            stats::setNames(as.list(rep(gval, length(br))), br)
          } else {
            NULL
          }
          cal_om <- omnibus_from_components(arr_cal, gamma = g)
          ev_om  <- omnibus_from_components(arr_ev,  gamma = g)

          ## marginal: column-wise calibration
          p_marg <- calibrate_cols(ev_om, cal_om)
          ## family-wise: row-wise minimum, threshold fitted on the cal draws
          thr <- wy_threshold(cal_om, ALPHA)
          fw  <- fwer_rates(ev_om, thr, K, ALPHA)

          acc$t1a  <- c(acc$t1a,  type1(exp(ev_om), ALPHA))
          acc$t1m  <- c(acc$t1m,  type1(p_marg, ALPHA))
          acc$wy   <- c(acc$wy,   fw[["wy"]])
          acc$bonf <- c(acc$bonf, fw[["bonferroni"]])
          acc$lam  <- c(acc$lam,  lambda_gc(as.vector(exp(ev_om))))
          acc$lamm <- c(acc$lamm, lambda_gc(as.vector(p_marg)))
        }
        res[[length(res) + 1L]] <- data.frame(
          ld = rn, rho = RHO[[rn]], m = m, pathways_per_gene = memb,
          gamma = gval, direction = glab, k_broken = k, n_combos = length(combos),
          t1_analytic = mean(acc$t1a), t1_mvn = mean(acc$t1m),
          fwer_wy = mean(acc$wy), fwer_bonferroni = mean(acc$bonf),
          lambda_analytic = mean(acc$lam), lambda_mvn = mean(acc$lamm),
          stringsAsFactors = FALSE)
        cat(sprintf("    %-16s k=%d  marginal: an=%.4f mvn=%.4f | FWER: wy=%.3f bonf=%.3f\n",
                    glab, k, mean(acc$t1a), mean(acc$t1m),
                    mean(acc$wy), mean(acc$bonf)))
      }
    }
  }
}
out <- do.call(rbind, res)
write.csv(out, file.path(OUT, "02_breakage_stress.csv"), row.names = FALSE)
cat("\nwrote results/02_breakage_stress.csv\n")
print(out[out$k_broken %in% c(0L, 4L),
          c("ld","m","direction","k_broken","t1_mvn","fwer_wy","fwer_bonferroni")],
      row.names = FALSE)
