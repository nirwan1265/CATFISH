#!/usr/bin/env Rscript
## =============================================================================
## diag_ties.R -- what exactly is the family minimum tied AT, in 01's design?
## Read-only. ~4 min.
##
## 01 correlates ALL m genes in a pathway at rho, unlike 02/03/08/09 which use
## block LD over a shared pool. This reproduces 01's setup at the two cells that
## tie and reports the actual numbers instead of inferring them.
## =============================================================================
suppressPackageStartupMessages(library(CATFISH))
source("00_lib.R")
set.seed(20260706)

NP <- 200L; NDATA <- 120L
cat(sprintf("MIN_P_DEFAULT = %g   log = %.3f\n\n", MIN_P_DEFAULT, log(MIN_P_DEFAULT)))

for (cfg in list(list(rho = 0.3, m = 50L), list(rho = 0.6, m = 50L),
                 list(rho = 0.6, m = 20L))) {
  U <- safe_chol(R_exch(cfg$m, cfg$rho))
  mins <- numeric(NDATA); allst <- NULL
  for (r in seq_len(NDATA)) {
    cc <- components_from_zmat(draw_zmat(U, NP))
    allst <- rbind(allst, cc)
    mins[r] <- min(cc[, "Omnibus"], na.rm = TRUE)
  }
  cat(sprintf("== rho=%.1f m=%d\n", cfg$rho, cfg$m))
  cat(sprintf("   family minimum: %d distinct values in %d datasets\n",
              length(unique(round(mins, 10))), NDATA))
  tb <- sort(table(round(mins, 10)), decreasing = TRUE)
  for (i in seq_len(min(3, length(tb))))
    cat(sprintf("     log p = %-12s x%-4d  (p = %.3e)\n", names(tb)[i],
                as.integer(tb[i]), exp(as.numeric(names(tb)[i]))))
  cat(sprintf("   at MIN_P floor: %.3f of minima, %.4f of ALL pathway omnibus values\n",
              mean(mins <= log(MIN_P_DEFAULT) + 1e-6),
              mean(allst[, "Omnibus"] <= log(MIN_P_DEFAULT) + 1e-6)))
  ## where does each component bottom out?
  for (nm in colnames(allst))
    cat(sprintf("     %-9s min log p = %9.3f  (p = %.3e)   n at that value = %d\n",
                nm, min(allst[, nm], na.rm = TRUE),
                exp(min(allst[, nm], na.rm = TRUE)),
                sum(allst[, nm] <= min(allst[, nm], na.rm = TRUE) + 1e-9, na.rm = TRUE)))
  cat("\n")
}
cat("ACAT in linear space cannot resolve p below ~1e-17, because 0.5 - p\n")
cat("rounds to 0.5 in double precision; check whether the tie sits there.\n")
