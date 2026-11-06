#!/usr/bin/env Rscript
## =============================================================================
## 24_sorghum_pcol_check.R -- which MAGMA gene p-value column to feed the
## Westfall-Young pathway test for the sorghum BAP traits?
##
## MAGMA multi=snp-wise reports three gene p-values:
##   P_SNPWISE_MEAN  mean of the SNP chi-squares       (MAGMA's default model;
##                   the one the LDL/HDL analyses used, column P)
##   P_SNPWISE_TOP1  the top SNP, corrected for gene size and LD
##   P_MULTI         the two combined
## The pathway null assumes gene p-values are uniform where there is no signal,
## so the column must be calibrated. Four checks, per trait and column:
##   1. calibration: median p, gene-level lambda, share of genes below 0.05 / 1e-3
##   2. dependence on gene size: median p by NSNPS quintile (a calibrated test
##      should not get "more significant" just because a gene has more SNPs)
##   3. the pre-specified established genes: percentile rank among all genes
##   4. QQ plot (Supplementary figure)
## LDL and HDL (mean model) are shown as a reference: both are highly
## polygenic, so some departure from 0.5 there is real signal.
##
## Usage (from catfish2/scripts):  Rscript 24_sorghum_pcol_check.R
## Env: MAGMA_ROOT (sorghum per-trait MAGMA folders), LIPID_ROOT (results/human_LDL)
## Writes ../results/sorghum_bap/24_pcol_calibration.csv
##        ../results/sorghum_bap/24_pcol_by_gene_size.csv
##        ../results/sorghum_bap/24_pcol_established_genes.csv
##        ../figures/Supp_sorghum_pcol_qq.png
## =============================================================================
MAGMA_ROOT <- Sys.getenv("MAGMA_ROOT", path.expand("~/Documents/Research/results/CATFISH/MAGMA"))
BASE       <- Sys.getenv("CATFISH_BASE", path.expand("~/Documents/Github/MAGCAT"))
LIPID_ROOT <- Sys.getenv("LIPID_ROOT", file.path(BASE, "results/human_LDL"))
OUT_DIR    <- Sys.getenv("OUT_DIR", "../results/sorghum_bap")
FIG_DIR    <- Sys.getenv("FIG_DIR", "../figures")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

TRAITS <- c("Dry_tons_per_acre", "Lignin_DM_perc")
COLS   <- c("P_SNPWISE_MEAN", "P_SNPWISE_TOP1", "P_MULTI")

## established genes, named before any pathway result was seen
ESTABLISHED <- data.frame(
  trait = c(rep("Dry_tons_per_acre", 8), rep("Lignin_DM_perc", 3)),
  gene  = c("Ma1/PRR37", "Ma3/PHYB", "Ma5/PHYC", "Ma6/GHD7", "Ma2",
            "Dw1", "Dw2", "Dw3", "bmr2/4CL", "bmr6/CAD", "bmr12/COMT"),
  id    = c("SORBI_3006G057866", "SORBI_3001G394400", "SORBI_3001G087100",
            "SORBI_3006G004400", "SORBI_3002G302700", "SORBI_3009G229800",
            "SORBI_3006G067700", "SORBI_3007G163800", "SORBI_3004G062500",
            "SORBI_3004G071000", "SORBI_3007G047300"),
  stringsAsFactors = FALSE)

read_sorghum <- function(trait) {
  f <- file.path(MAGMA_ROOT, trait,
                 sprintf("%s_chr%02d.multi_snp_wise.genes.out", trait, 1:10))
  g <- do.call(rbind, lapply(f, utils::read.table, header = TRUE,
                             stringsAsFactors = FALSE, comment.char = "#"))
  g[!duplicated(g$GENE), ]
}
read_lipid <- function(trait) {
  d <- if (trait == "ldl") "magma" else paste0("magma_", trait)
  g <- utils::read.table(file.path(LIPID_ROOT, d, sprintf("%s_eur.genes.out", trait)),
                         header = TRUE, stringsAsFactors = FALSE)
  g$P_SNPWISE_MEAN <- g$P
  g
}
lam <- function(p) stats::qchisq(stats::median(p, na.rm = TRUE), 1, lower.tail = FALSE) /
                   stats::qchisq(0.5, 1)

dat <- c(setNames(lapply(TRAITS, read_sorghum), TRAITS),
         list(LDL = read_lipid("ldl"), HDL = read_lipid("hdl")))

## ---- 1. calibration ---------------------------------------------------------
cal <- do.call(rbind, lapply(names(dat), function(t) {
  cc <- intersect(COLS, names(dat[[t]]))
  do.call(rbind, lapply(cc, function(k) {
    p <- dat[[t]][[k]]; p <- p[is.finite(p)]
    data.frame(trait = t, column = k, n_genes = length(p),
               median_p = round(stats::median(p), 3), lambda_gene = round(lam(p), 2),
               pct_below_0.05 = round(100 * mean(p < 0.05), 1),
               pct_below_1e3  = round(100 * mean(p < 1e-3), 2))
  }))
}))
utils::write.csv(cal, file.path(OUT_DIR, "24_pcol_calibration.csv"), row.names = FALSE)
cat("== 1. calibration (null: median p 0.5, lambda 1, 5% below 0.05, 0.1% below 1e-3)\n")
print(cal, row.names = FALSE)

## ---- 2. dependence on gene size ---------------------------------------------
sz <- do.call(rbind, lapply(names(dat), function(t) {
  g <- dat[[t]]
  q <- cut(g$NSNPS, unique(stats::quantile(g$NSNPS, 0:5 / 5)), include.lowest = TRUE)
  cc <- intersect(COLS, names(g))
  do.call(rbind, lapply(cc, function(k) {
    m <- tapply(g[[k]], q, stats::median, na.rm = TRUE)
    data.frame(trait = t, column = k, quintile = paste0("Q", seq_along(m)),
               nsnps_range = names(m), median_p = round(as.numeric(m), 3))
  }))
}))
utils::write.csv(sz, file.path(OUT_DIR, "24_pcol_by_gene_size.csv"), row.names = FALSE)
cat("\n== 2. median gene p by NSNPS quintile, Q1 smallest genes (flat = no size dependence)\n")
wide <- stats::reshape(sz[, c("trait", "column", "quintile", "median_p")],
                       idvar = c("trait", "column"), timevar = "quintile", direction = "wide")
names(wide) <- sub("^median_p\\.", "", names(wide))
print(wide, row.names = FALSE)

## ---- 3. established genes -------------------------------------------------
est <- do.call(rbind, lapply(seq_len(nrow(ESTABLISHED)), function(i) {
  e <- ESTABLISHED[i, ]; g <- dat[[e$trait]]
  r <- sapply(COLS, function(k) {
    j <- match(e$id, g$GENE); if (is.na(j)) return(NA_real_)
    100 * mean(g[[k]] <= g[[k]][j], na.rm = TRUE)   # top-x% of genes
  })
  p <- sapply(COLS, function(k) g[[k]][match(e$id, g$GENE)])
  data.frame(e, nsnps = g$NSNPS[match(e$id, g$GENE)],
             setNames(as.list(signif(p, 3)), paste0("p_", COLS)),
             setNames(as.list(round(r, 1)), paste0("top_pct_", COLS)))
}))
utils::write.csv(est, file.path(OUT_DIR, "24_pcol_established_genes.csv"), row.names = FALSE)
cat("\n== 3. established genes: gene p and percentile rank (top x% of all genes)\n")
print(est[, c("trait", "gene", "nsnps", grep("^(p_|top_pct_)", names(est), value = TRUE))],
      row.names = FALSE)
cat("\nmedian percentile rank of the established genes, per trait and column:\n")
print(stats::aggregate(est[, paste0("top_pct_", COLS)], list(trait = est$trait), stats::median))

## ---- 4. QQ plot -----------------------------------------------------------
col_of <- c(P_SNPWISE_MEAN = "#1b6ca8", P_SNPWISE_TOP1 = "#c0392b", P_MULTI = "#7d3c98")
png(file.path(FIG_DIR, "Supp_sorghum_pcol_qq.png"), width = 2400, height = 1150, res = 220)
op <- par(mfrow = c(1, 2), mar = c(4.2, 4.4, 2.6, 1), mgp = c(2.6, 0.7, 0))
for (t in TRAITS) {
  g <- dat[[t]]
  mx <- max(unlist(lapply(COLS, function(k) -log10(min(g[[k]], na.rm = TRUE)))))
  n <- sum(is.finite(g$P_SNPWISE_MEAN)); ex <- -log10(stats::ppoints(n))
  plot(NA, xlim = c(0, max(ex)), ylim = c(0, mx), las = 1,
       xlab = expression(Expected ~ -log[10](p)), ylab = expression(Observed ~ -log[10](p)),
       main = gsub("_", " ", t))
  abline(0, 1, col = "grey50", lty = 2)
  for (k in COLS) {
    o <- sort(-log10(g[[k]][is.finite(g[[k]])]), decreasing = TRUE)
    points(-log10(stats::ppoints(length(o))), o, pch = 16, cex = 0.35, col = col_of[k])
  }
  legend("topleft", bty = "n", pch = 16, col = col_of,
         legend = sprintf("%s (median p %.2f)", names(col_of),
                          sapply(COLS, function(k) stats::median(g[[k]], na.rm = TRUE))))
}
par(op); invisible(dev.off())
cat(sprintf("\nwrote tables to %s and %s\n", OUT_DIR, file.path(FIG_DIR, "Supp_sorghum_pcol_qq.png")))
