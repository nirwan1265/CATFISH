#!/usr/bin/env Rscript
## =============================================================================
## diag_tfstat.R -- would an adaptive score built from the RAW TFisher
## statistics rank pathways the same way as Wo, while keeping resolution where
## Wo saturates? Read-only, changes nothing. ~5 min.
##
##   Wo  = min over tau of the right-tail p   (current; dies below ~1e-16)
##   S   = max over tau of the z-scored W(tau) against the null bank
##         (larger = stronger enrichment; no ceiling)
##
## If the two orderings agree where Wo is resolvable, S is the SAME test with
## working numerics rather than a different one.
## =============================================================================
suppressPackageStartupMessages(library(CATFISH))
source("00_lib.R")
set.seed(20260706)

NP <- 200L; NDATA <- 60L; TG <- TAU_GRID_DEFAULT
cat(sprintf("tau grid: %s\n\n", paste(TG, collapse = ", ")))

## W(tau) for every pathway, one package call per tau (single-tau grid)
W_by_tau <- function(gr, pw) {
  vapply(TG, function(tt) {
    o <- CATFISH::catfish_soft_tfisher_adaptive_pathways(
           gene_results = gr, pathways = pw, gene_col = "GENE", p_col = "P",
           tau_grid = tt, min_p = MIN_P_DEFAULT, output = FALSE)
    as.numeric(o$tfisher_stat_hat)[match(names(pw), o$pathway_id)]
  }, numeric(length(pw)))
}
Wo_of <- function(gr, pw) {
  o <- CATFISH::catfish_soft_tfisher_adaptive_pathways(
         gene_results = gr, pathways = pw, gene_col = "GENE", p_col = "P",
         tau_grid = TG, min_p = MIN_P_DEFAULT, output = FALSE)
  as.numeric(o$tfisher_p_min)[match(names(pw), o$pathway_id)]
}

for (cfg in list(list(rho = 0.3, m = 50L), list(rho = 0.6, m = 50L))) {
  U <- safe_chol(R_exch(cfg$m, cfg$rho))
  Wall <- NULL; Woall <- NULL
  for (r in seq_len(NDATA)) {
    z   <- draw_zmat(U, NP)
    ids <- universe_gene_ids(NP * cfg$m)
    pw  <- split(ids, rep(seq_len(NP), each = cfg$m))
    names(pw) <- sprintf("PW%06d", seq_len(NP))
    gr  <- data.frame(GENE = ids, P = stats::pnorm(as.vector(t(z)),
                      lower.tail = FALSE), stringsAsFactors = FALSE)
    Wall  <- rbind(Wall, W_by_tau(gr, pw))
    Woall <- c(Woall, Wo_of(gr, pw))
  }
  ## z-score each tau column against the pooled null, then take the max
  mu <- colMeans(Wall); sdv <- apply(Wall, 2, stats::sd)
  S  <- apply(sweep(sweep(Wall, 2, mu, "-"), 2, sdv, "/"), 1, max)

  sat <- Woall <= 1e-16
  cat(sprintf("== rho=%.1f m=%d   (%d pathways)\n", cfg$rho, cfg$m, length(Woall)))
  cat(sprintf("   Wo saturated (<=1e-16): %.4f\n", mean(sat)))
  cat(sprintf("   distinct values  Wo: %-6d   S: %d\n",
              length(unique(Woall)), length(unique(round(S, 9)))))
  ok <- !sat & is.finite(Woall) & Woall < 1
  cat(sprintf("   Spearman(S, -log Wo) where Wo resolves (n=%d): %.4f\n",
              sum(ok), stats::cor(S[ok], -log(Woall[ok]), method = "spearman")))
  if (any(sat))
    cat(sprintf("   among saturated pathways, S still takes %d distinct values\n",
                length(unique(round(S[sat], 9)))))
  cat("\n")
}
cat("Want: Spearman near 1 (same test), and S keeping many distinct values\n")
cat("where Wo has collapsed to one.\n")
