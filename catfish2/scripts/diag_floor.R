#!/usr/bin/env Rscript
## =============================================================================
## diag_floor.R -- is the family minimum in 01 tied because components saturate
## at the package min_p floor? Read-only, writes nothing. ~2 min.
## =============================================================================
suppressPackageStartupMessages(library(CATFISH))
source("00_lib.R")
set.seed(20260706)

NP <- 200L; NDATA <- 50L
for (rn in c("independent", "moderate", "strong")) {
  rho <- c(independent = 0, moderate = 0.3, strong = 0.6)[[rn]]
  for (m in c(20L, 50L)) {
    U <- safe_chol(R_exch(m, rho))
    allc <- NULL; mins <- numeric(NDATA)
    for (r in seq_len(NDATA)) {
      cc <- components_from_zmat(draw_zmat(U, NP))
      allc <- rbind(allc, cc)
      mins[r] <- min(cc[, "Omnibus"])
    }
    fl <- log(MIN_P_DEFAULT)
    at_floor <- colMeans(allc <= fl + 1e-6)          # exactly clamped
    near     <- colMeans(allc <= fl + log(10))       # within 10x of the floor
    cat(sprintf("\n== LD=%s m=%d  (%d pathways)\n", rn, m, NP * NDATA))
    for (nm in colnames(allc))
      cat(sprintf("   %-9s at floor=%.4f  within10x=%.4f  min=%.1f\n",
                  nm, at_floor[[nm]], near[[nm]], min(allc[, nm])))
    tb <- sort(table(round(mins, 9)), decreasing = TRUE)
    cat(sprintf("   family minimum: %d distinct values in %d datasets; commonest %s x%d\n",
                length(tb), NDATA, names(tb)[1], as.integer(tb[1])))
  }
}
