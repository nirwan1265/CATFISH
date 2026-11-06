#!/usr/bin/env Rscript
## =============================================================================
## 06_aggregate.R -- the permutation null, pooled over permutations.
## Writes FINAL_DIR/null_summary.csv, null_per_perm.csv and null_pathway_p.csv.
##   per-pathway calibrated p      type I error at 0.05 / 0.01, lambda
##   family-wise (Westfall-Young)  fraction of permutations with >= 1 pathway
##                                 at p_fwer_std <= 0.05 (should be <= 0.05)
##   gene level                    median gene p per permutation (about 0.5)
## =============================================================================
getcfg <- function(k) { v <- Sys.getenv(k); if (v == "") stop("config missing: ", k); v }
WY_OUT <- getcfg("WY_OUT"); FINAL_DIR <- getcfg("FINAL_DIR")
fs <- list.files(WY_OUT, pattern = "^perm_[0-9]+_wy\\.csv$", full.names = TRUE)
if (!length(fs)) stop("no per-permutation results in ", WY_OUT)
x <- do.call(rbind, lapply(fs, utils::read.csv, stringsAsFactors = FALSE))
lam <- function(p) { p <- p[is.finite(p) & p > 0 & p <= 1]; stats::qchisq(stats::median(p), 1, lower.tail = FALSE) / stats::qchisq(0.5, 1) }
ci <- function(k, n) stats::binom.test(k, n)$conf.int
per <- do.call(rbind, lapply(split(x, x$perm), function(d) data.frame(
  perm = d$perm[1], n_pathways = nrow(d), median_gene_p = d$median_gene_p[1],
  n_cal_05 = sum(d$p_calibrated < 0.05, na.rm = TRUE),
  min_p_fwer = min(d$p_fwer_std, na.rm = TRUE),
  any_fwer_05 = any(d$p_fwer_std <= 0.05, na.rm = TRUE))))
n_perm <- nrow(per); k_fw <- sum(per$any_fwer_05); pc <- x$p_calibrated[is.finite(x$p_calibrated)]
rows <- list(
  c("permutations", n_perm, NA, NA),
  c("null pathway p-values", length(pc), NA, NA),
  c("median gene p (median over permutations)", round(stats::median(per$median_gene_p), 4), NA, NA),
  c("lambda, calibrated p", round(lam(pc), 3), NA, NA),
  c("type I error at 0.05, calibrated p", round(mean(pc < 0.05), 4), round(ci(sum(pc < 0.05), length(pc)), 4)),
  c("type I error at 0.01, calibrated p", round(mean(pc < 0.01), 4), round(ci(sum(pc < 0.01), length(pc)), 4)),
  c("FWER: permutations with >=1 pathway at p_fwer_std <= 0.05",
    round(k_fw / n_perm, 4), round(ci(k_fw, n_perm), 4)))
s <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
names(s) <- c("quantity", "value", "ci_lo", "ci_hi")
dir.create(FINAL_DIR, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(s, file.path(FINAL_DIR, "null_summary.csv"), row.names = FALSE)
utils::write.csv(per, file.path(FINAL_DIR, "null_per_perm.csv"), row.names = FALSE)
utils::write.csv(x[, intersect(c("perm", "pathway_id", "n_genes", "p_calibrated", "p_fwer_std",
                                 "z_std", "p_omnibus"), names(x))],
                 file.path(FINAL_DIR, "null_pathway_p.csv"), row.names = FALSE)
print(s, row.names = FALSE)
