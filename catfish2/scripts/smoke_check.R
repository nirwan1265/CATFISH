#!/usr/bin/env Rscript
## =============================================================================
## smoke_check.R -- checks the package-backed component layer and times it.
## Writes nothing.
##
## TFisher now contributes Wo (package column tfisher_p_min), a STATISTIC, not
## a p-value: smaller = more extreme, but not uniform under the null.
## =============================================================================
suppressPackageStartupMessages(library(CATFISH))
source("00_lib.R")
set.seed(1)

G <- 3000L; BLOCK <- 50L; RHO <- 0.3; K <- 150L; POOL <- 600L; M <- 20L
uni   <- make_universe(G, BLOCK, RHO)
paths <- make_pathways(K, M, POOL, G)
pchr  <- paths_to_chr(paths, G)

cat("[1] one null draw through the CATFISH component wrappers\n")
z  <- draw_universe_z(uni)
lp <- stats::pnorm(z, lower.tail = FALSE, log.p = TRUE)
t0 <- Sys.time()
mm <- components_from_draw(lp, z, pchr, G)
el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("    %d pathways x %d stats in %.2f s\n", nrow(mm), ncol(mm), el))
print(round(exp(mm[1:3, ]), 5))
cat(sprintf("    non-finite cells: %d\n", sum(!is.finite(mm))))

cat("\n[2] TFisher direction and saturation\n")
for (lab in c("null", "moderate", "strong")) {
  zz <- z
  shift <- switch(lab, null = 0, moderate = 2, strong = 5)
  zz[paths[[1]]] <- zz[paths[[1]]] + c(rep(shift, round(0.3 * M)),
                                       rep(0, M - round(0.3 * M)))
  ll <- stats::pnorm(zz, lower.tail = FALSE, log.p = TRUE)
  r  <- components_from_draw(ll, zz, pchr[1], G)[1, ]
  cat(sprintf("    %-9s min gene p=%.2e  TFisher Wo=%.3e  Omnibus=%.3e\n",
              lab, min(exp(ll[paths[[1]]])), exp(r[["TFisher"]]),
              exp(r[["Omnibus"]])))
}
cat("    Wo must FALL as signal rises, and must NOT sit at 1e-50\n")

cat("\n[3] how much of a null bank saturates at the floor\n")
arr <- components_array(uni, paths, 20L, verbose = FALSE)
fl  <- log(MIN_P_DEFAULT)
for (nm in dimnames(arr)[[3]])
  cat(sprintf("    %-9s at floor=%.4f   min=%8.2f\n", nm,
              mean(arr[, , nm] <= fl + 1e-6), min(arr[, , nm])))

cat("\n[4] timing\n")
t0 <- Sys.time(); invisible(components_array(uni, paths, 10L, verbose = FALSE))
per <- as.numeric(difftime(Sys.time(), t0, units = "secs")) / 10
cat(sprintf("    %.3f s per draw (K=%d, m=%d)\n", per, K, M))
for (B in c(500, 2000, 15000))
  cat(sprintf("    B=%-6d -> %.1f min per setting\n", B, per * B / 60))
