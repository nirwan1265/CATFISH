#!/usr/bin/env Rscript
## =============================================================================
## 09_aggregate_magma_gsa.R -- MAGMA's own gene-set analysis under the
## permutation null, summarised as 06_aggregate.R summarises CATFISH.
## Pathways: those with MIN_GENES-MAX_GENES genes present (MAGMA NGENES), the
## size rule CATFISH uses. For the competitive and the self-contained test:
##   type I error at 0.05 / 0.01 / 0.001 and lambda, over all null pathway tests
##   family-wise error: fraction of permutations with >= 1 pathway at
##     Bonferroni p <= 0.05 / K, and with >= 1 pathway at BH q <= 0.05
## Writes RUN_ROOT/final_magma/{magma_null_summary,magma_null_per_perm,magma_null_pathway_p}.csv
## =============================================================================
suppressPackageStartupMessages(library(data.table))
getcfg <- function(k) { v <- Sys.getenv(k); if (v == "") stop("config missing: ", k); v }
RUN_ROOT <- getcfg("RUN_ROOT"); MINg <- as.integer(getcfg("MIN_GENES")); MAXg <- as.integer(getcfg("MAX_GENES"))
GSA <- file.path(RUN_ROOT, "magma_gsa"); OUT <- file.path(RUN_ROOT, "final_magma")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

read_magma <- function(f) {                 # MAGMA output: '#' comment lines, then a header
  ln <- readLines(f); ln <- ln[!startsWith(ln, "#")]
  d <- fread(text = paste(ln, collapse = "\n"))
  nm <- if ("FULL_NAME" %in% names(d)) d$FULL_NAME else d[[1]]
  data.table(pathway_id = nm, ngenes = d$NGENES, p = d$P)
}
fs <- list.files(GSA, pattern = "^perm_[0-9]+\\.gsa\\.out$", full.names = TRUE)
if (!length(fs)) stop("no gsa.out in ", GSA)
x <- rbindlist(lapply(fs, function(f) {
  b <- as.integer(sub("^perm_([0-9]+)\\.gsa\\.out$", "\\1", basename(f)))
  comp <- read_magma(f)
  fself <- sub("\\.gsa\\.out$", ".gsa.self.out", f)
  self <- if (file.exists(fself)) read_magma(fself)[, .(pathway_id, p_self = p)] else data.table(pathway_id = character(), p_self = numeric())
  d <- merge(comp[, .(pathway_id, ngenes, p_comp = p)], self, by = "pathway_id", all.x = TRUE)
  d[ngenes >= MINg & ngenes <= MAXg][, perm := b]
}))
fwrite(x, file.path(OUT, "magma_null_pathway_p.csv"))

lam <- function(p) { p <- p[is.finite(p) & p > 0 & p <= 1]; qchisq(median(p), 1, lower.tail = FALSE) / qchisq(0.5, 1) }
ci  <- function(k, n) binom.test(k, n)$conf.int
per <- x[, .(K = .N,
             comp_bonf = any(p_comp <= 0.05 / .N, na.rm = TRUE),
             comp_bh   = any(p.adjust(p_comp, "BH") <= 0.05, na.rm = TRUE),
             self_bonf = any(p_self <= 0.05 / .N, na.rm = TRUE),
             self_bh   = any(p.adjust(p_self, "BH") <= 0.05, na.rm = TRUE)), by = perm]
fwrite(per, file.path(OUT, "magma_null_per_perm.csv"))
n <- nrow(per)
row <- function(test, p, fw_b, fw_q) {
  p <- p[is.finite(p)]; c5 <- ci(sum(p < 0.05), length(p)); fb <- ci(sum(fw_b), n); fq <- ci(sum(fw_q), n)
  data.table(test = test, permutations = n, null_tests = length(p), lambda = round(lam(p), 3),
             t05 = round(mean(p < 0.05), 4), t05_lo = round(c5[1], 4), t05_hi = round(c5[2], 4),
             t01 = round(mean(p < 0.01), 4), t001 = round(mean(p < 0.001), 5),
             fwer_bonf = mean(fw_b), fwer_bonf_lo = round(fb[1], 4), fwer_bonf_hi = round(fb[2], 4),
             fwer_bh = mean(fw_q), fwer_bh_lo = round(fq[1], 4), fwer_bh_hi = round(fq[2], 4))
}
s <- rbind(row("MAGMA competitive", x$p_comp, per$comp_bonf, per$comp_bh),
           row("MAGMA self-contained", x$p_self, per$self_bonf, per$self_bh))
fwrite(s, file.path(OUT, "magma_null_summary.csv"))
print(s)
cat("pathways per permutation:", paste(range(per$K), collapse = "-"), "\n")
