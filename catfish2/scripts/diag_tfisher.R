#!/usr/bin/env Rscript
## =============================================================================
## diag_tfisher.R -- read-only diagnostic for 00_lib.R::tfisher_log.
## Changes nothing. Prints what TFisher::p.soft actually returns so we can see
## whether it underflows to 0 on strong signal, saturates at a floor, or both.
## =============================================================================
source("00_lib.R")
if (!requireNamespace("TFisher", quietly = TRUE)) stop("TFisher not installed")
set.seed(1)

one_sided <- function(z) stats::pnorm(z, lower.tail = FALSE)

cat("TAU grid:", paste(TAU_GRID_DEFAULT, collapse = ", "), "\n")
cat("packageVersion(TFisher):", as.character(packageVersion("TFisher")), "\n\n")

for (m in c(20L, 50L)) {
  scen <- list(
    null     = rep(0,   m),
    moderate = c(rep(2.0, max(1, round(0.3 * m))), rep(0, m - max(1, round(0.3 * m)))),
    strong   = c(rep(5.0, max(1, round(0.3 * m))), rep(0, m - max(1, round(0.3 * m))))
  )
  for (nm in names(scen)) {
    p <- one_sided(scen[[nm]] + stats::rnorm(m))
    cat(sprintf("== m=%d  %s  (min p = %.3e)\n", m, nm, min(p)))
    for (tau in TAU_GRID_DEFAULT) {
      st <- TFisher::stat.soft(p = p, tau1 = tau)
      pv <- TFisher::p.soft(q = st, n = m, tau1 = tau, M = NULL)
      cat(sprintf("   tau=%-6.3f  stat=%12.4f  p.soft=%-14.6g  zero=%-5s finite=%-5s\n",
                  tau, st, pv, pv == 0, is.finite(pv)))
    }
    cat(sprintf("   tfisher_log() -> %s   (exp = %s)\n\n",
                format(tfisher_log(log(p))), format(exp(tfisher_log(log(p))))))
  }
}

## How far into the tail can p.soft go at all? Push the statistic up by hand.
cat("== tail behaviour of p.soft at m=20, tau=0.05 (statistic forced upward)\n")
for (q in c(5, 10, 20, 30, 40, 60, 80, 120, 200)) {
  pv <- TFisher::p.soft(q = q, n = 20L, tau1 = 0.05, M = NULL)
  cat(sprintf("   q=%-5g  p.soft=%-14.6g zero=%s\n", q, pv, pv == 0))
}

## And what the null distribution of tfisher_log looks like over many draws,
## next to one component that is known to behave (Fisher).
cat("\n== 2000 null pathways at m=20: quantiles of log p\n")
lp_null <- replicate(2000, {
  lp <- stats::pnorm(stats::rnorm(20L), lower.tail = FALSE, log.p = TRUE)
  c(TFisher = tfisher_log(lp), Fisher = fisher_log(lp))
})
for (r in rownames(lp_null)) {
  v <- lp_null[r, ]
  cat(sprintf("   %-8s NA frac=%.3f   min=%8.3f  q01=%8.3f  q05=%8.3f  median=%8.3f\n",
              r, mean(!is.finite(v)),
              min(v, na.rm = TRUE),
              stats::quantile(v, 0.01, na.rm = TRUE, names = FALSE),
              stats::quantile(v, 0.05, na.rm = TRUE, names = FALSE),
              stats::median(v, na.rm = TRUE)))
}
