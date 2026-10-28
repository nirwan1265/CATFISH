#!/usr/bin/env Rscript
## =============================================================================
## diag_tau.R -- how often does Wo come back as exactly 1 (an empty soft-
## threshold statistic), by pathway size and tau grid? Read-only. ~3 min.
##
##   Wo = 1 happens when NO gene in the pathway beats any tau. ACAT then sees
##   1 - 1/m and is pushed toward the null, which is what wrecked lambda at m=5.
## =============================================================================
suppressPackageStartupMessages(library(CATFISH))
source("00_lib.R")
set.seed(7)

GRIDS <- list(
  current = c(0.1, 0.05, 0.01, 0.005, 0.001),
  authors = c(0.01, 0.05, 0.5),   # Zhang et al., tau=1 Fisher corner removed
  auto    = c(0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.001),
  loose   = c(0.5, 0.25, 0.1, 0.05, 0.01),
  wide    = c(0.5, 0.2, 0.05, 0.005)
)
NPATH <- 2000L

## k = number of genes below the LOOSEST tau. k>=1 is needed for the statistic
## to exist at all; k>=2 is needed for it to be COMBINING anything -- with one
## term soft TFisher is a monotone function of a single gene p, i.e. minP,
## which is already a separate component.
## lambda is a MEDIAN statistic; type I error is the tail. They can disagree,
## so both are reported. t1_05 / t1_01 are on the ANALYTIC omnibus p (no
## calibration here) -- the calibrated rates come from 01/02/09.
cat(sprintf("%-4s %-8s %7s %7s %8s %8s %8s %8s %8s\n",
            "m", "grid", "k>=1", "k>=2", "lam_Wo", "lam_omni", "medOmni",
            "t1_05", "t1_01"))
for (m in c(5L, 10L, 20L, 50L)) {
  ids <- universe_gene_ids(NPATH * m)
  pw  <- split(ids, rep(seq_len(NPATH), each = m))
  names(pw) <- sprintf("PW%06d", seq_len(NPATH))
  z   <- stats::rnorm(NPATH * m)
  gr  <- data.frame(GENE = ids, P = stats::pnorm(z, lower.tail = FALSE),
                    ZSTAT = z, stringsAsFactors = FALSE)
  for (gn in names(GRIDS)) {
    tg <- GRIDS[[gn]]
    tf <- CATFISH::catfish_soft_tfisher_adaptive_pathways(
            gene_results = gr, pathways = pw, gene_col = "GENE", p_col = "P",
            tau_grid = tg, min_p = MIN_P_DEFAULT, output = FALSE)
    Wo <- as.numeric(tf$tfisher_p_min)[match(names(pw), tf$pathway_id)]
    mm <- components_from_draw(
            stats::pnorm(z, lower.tail = FALSE, log.p = TRUE), z,
            pw, NPATH * m, tau_grid = tg)
    om <- exp(mm[, "Omnibus"])
    pmat <- matrix(gr$P, nrow = m)          # m genes x NPATH pathways
    nbelow <- colSums(pmat <= max(tg))
    cat(sprintf("%-4d %-8s %7.3f %7.3f %8.3f %8.3f %8.3f %8.4f %8.4f\n",
                m, gn, mean(nbelow >= 1), mean(nbelow >= 2), lambda_gc(Wo),
                lambda_gc(om), stats::median(om, na.rm = TRUE),
                mean(om < 0.05, na.rm = TRUE), mean(om < 0.01, na.rm = TRUE)))
  }
}
cat("\nWant: lam_omni ~1, medOmni ~0.5, t1_05 ~0.05, t1_01 ~0.01\n")
