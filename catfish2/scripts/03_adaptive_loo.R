#!/usr/bin/env Rscript
## =============================================================================
## 03_adaptive_loo.R -- leave-one-out component sensitivity, measured on BOTH
## the per-pathway (marginal) calibration AND the Westfall-Young family-wise
## rate.
##
## Each component is removed in turn and the whole calibration is rebuilt from
## scratch -- the null bank, the marginal p-values and the WY threshold. The
## question is whether any single component is load-bearing, i.e. whether
## dropping it breaks family-wise control or merely costs sensitivity.
##
## Shares the overlapping-universe design of 02: K pathways drawn from a pool
## of POOL genes inside a universe of G genes with block LD, so pathways
## genuinely share genes and WY has dependence to control. Component log
## p-values are computed once per setting and reused for every drop.
##
##   Rscript 03_adaptive_loo.R [B_CAL] [B_EVAL]
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
COMPONENTS <- c("ACAT","Fisher","Stouffer","minP","TFisher")
RHO   <- c(independent = 0.0, moderate = 0.3, strong = 0.6)
SIZE  <- c(20L, 50L)
OUT   <- "../results"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("leave-one-out: K=%d pathways, pool=%d of G=%d, B_cal=%d B_eval=%d\n",
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

    for (drop in c("<none>", COMPONENTS)) {
      dr <- if (identical(drop, "<none>")) character(0) else drop
      cal_om <- omnibus_from_components(arr_cal, drop = dr)
      ev_om  <- omnibus_from_components(arr_ev,  drop = dr)

      p_marg <- calibrate_cols(ev_om, cal_om)
      thr    <- wy_threshold(cal_om, ALPHA)
      fw     <- fwer_rates(ev_om, thr, K, ALPHA)

      res[[length(res) + 1L]] <- data.frame(
        ld = rn, rho = RHO[[rn]], m = m, pathways_per_gene = memb,
        dropped = drop,
        t1_analytic = type1(exp(ev_om), ALPHA),
        t1_mvn      = type1(p_marg, ALPHA),
        fwer_wy         = fw[["wy"]],
        fwer_bonferroni = fw[["bonferroni"]],
        wy_threshold    = exp(thr),
        bonferroni_threshold = ALPHA / K,
        lambda_analytic = lambda_gc(as.vector(exp(ev_om))),
        lambda_mvn      = lambda_gc(as.vector(p_marg)),
        stringsAsFactors = FALSE)
      cat(sprintf("    drop=%-9s marginal mvn=%.4f | FWER wy=%.3f bonf=%.3f | thr %.2e vs %.2e\n",
                  drop, type1(p_marg, ALPHA), fw[["wy"]], fw[["bonferroni"]],
                  exp(thr), ALPHA / K))
    }
  }
}
out <- do.call(rbind, res)
write.csv(out, file.path(OUT, "03_adaptive_loo.csv"), row.names = FALSE)
cat("\nwrote results/03_adaptive_loo.csv\n")
print(out[, c("ld","m","dropped","t1_mvn","fwer_wy","fwer_bonferroni")], row.names = FALSE)
