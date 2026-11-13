#!/usr/bin/env Rscript
## =============================================================================
## 08_aggregate.R -- power and false calls of every method, per scenario x size.
## Calls (family-wise 0.05):
##   catfish_self  p_fwer_std <= 0.05          catfish_comp  p_fwer_comp <= 0.05
##   catfish_both  max(both) <= 0.05           (Westfall-Young, B = WY_B)
##   magma_comp, magma_self, pascalx           p <= 0.05 / K (Bonferroni; K = the
##                                             pathways in SIM_SETS)
## Per replicate:
##   target_called  the target pathway is called (S2-S6)
##   n_off, any_off calls of "off-target" pathways: none of their genes within
##                  BUFFER_KB of a TARGET causal SNP. In S0 and S1 there is no
##                  target, so every call is off-target.
## Reading the off-target rate (any_off = family-wise error where it is an error):
##   S0     every method: no signal at all, any call is false
##   S1     competitive and both: background only, no pathway enriched, any call
##          is false; self-contained: background genes are real signal, so its
##          calls are expected (reported, not an error)
##   S2-S6  competitive and both: calls of pathways unrelated to the target;
##          self-contained: expected to call background pathways too
## Writes RUN/summary_by_replicate.csv and RUN/summary.csv (rates with 95% CI).
##   Rscript 08_aggregate.R
## =============================================================================
suppressPackageStartupMessages(library(data.table))
env <- function(k) { v <- Sys.getenv(k); if (!nzchar(v)) stop(k, " is not set"); v }
RUN <- env("RUN"); IN <- env("IN"); BUF <- as.numeric(env("BUFFER_KB")) * 1000
design <- fread(file.path(RUN, "design.tsv"))
causal <- fread(file.path(RUN, "causal.tsv"))
genes  <- fread(file.path(IN, "genes.tsv"), colClasses = list(character = "gene"))
sets <- list()
for (x in readLines(env("SIM_SETS"))) { f <- strsplit(x, "\t", fixed = TRUE)[[1]]; sets[[f[1]]] <- f[-1] }
K <- length(sets)
pw_genes <- rbindlist(lapply(names(sets), function(p) data.table(pathway_id = p, gene = sets[[p]])))
pw_genes <- merge(pw_genes, genes[, .(gene, chr, start, stop)], by = "gene")

## pathways with a gene within BUF of a TARGET causal SNP, per replicate
gwin <- unique(pw_genes[, .(gene, chr, lo = pmax(0, start - BUF), hi = stop + BUF)])
setkey(gwin, chr, lo, hi)
dirty <- function(r) {
  c_r <- causal[rep_id == r & grepl("target", role)]; if (!nrow(c_r)) return(character(0))
  pts <- c_r[, .(chr, lo = pos, hi = pos)]; setkey(pts, chr, lo, hi)
  g_hit <- unique(foverlaps(pts, gwin, nomatch = NULL)$gene)
  unique(pw_genes[gene %in% g_hit, pathway_id])
}
read_magma_gsa <- function(f) {
  ln <- readLines(f); ln <- ln[!startsWith(ln, "#")]
  d <- fread(text = paste(ln, collapse = "\n"), fill = TRUE)
  id <- if ("FULL_NAME" %in% names(d)) ifelse(is.na(d$FULL_NAME) | d$FULL_NAME == "", d[[1]], d$FULL_NAME) else d[[1]]
  data.table(pathway_id = id, p = d$P)
}

rows <- list()
for (i in seq_len(nrow(design))) {
  d <- design[i]; RP <- sprintf("rep_%04d", d$rep_id)
  f_cf <- file.path(RUN, "catfish", paste0(RP, ".csv"))
  f_mc <- file.path(RUN, "magma_gsa", paste0(RP, ".gsa.out"))
  f_ms <- file.path(RUN, "magma_gsa", paste0(RP, ".gsa.self.out"))
  f_px <- file.path(RUN, "pascalx", paste0(RP, ".tsv"))
  calls <- list()
  if (file.exists(f_cf)) {
    cf <- fread(f_cf)
    calls$catfish_self <- cf[p_fwer_std <= 0.05, pathway_id]
    calls$catfish_comp <- cf[p_fwer_comp <= 0.05, pathway_id]
    calls$catfish_both <- cf[pmax(p_fwer_std, p_fwer_comp) <= 0.05, pathway_id]
  }
  if (file.exists(f_mc)) { m <- read_magma_gsa(f_mc)[pathway_id %in% names(sets)]; calls$magma_comp <- m[p <= 0.05 / K, pathway_id] }
  if (file.exists(f_ms)) { m <- read_magma_gsa(f_ms)[pathway_id %in% names(sets)]; calls$magma_self <- m[p <= 0.05 / K, pathway_id] }
  if (file.exists(f_px)) {
    m <- fread(f_px, header = FALSE, col.names = c("pathway_id", "p"))[pathway_id %in% names(sets)]
    calls$pascalx <- m[p <= 0.05 / K, pathway_id]
  }
  if (!length(calls)) next
  dp <- dirty(d$rep_id)
  for (meth in names(calls)) {
    cl <- calls[[meth]]; fl <- setdiff(cl, dp)
    rows[[length(rows) + 1]] <- data.table(
      rep_id = d$rep_id, scenario = d$scenario, archetype = d$archetype, size_class = d$size_class,
      scale = d$scale, method = meth, n_called = length(cl),
      target_called = if (is.na(d$target_id) || !nzchar(d$target_id)) NA else d$target_id %in% cl,
      n_off = length(fl), any_off = length(fl) > 0)
  }
}
byrep <- rbindlist(rows)
fwrite(byrep, file.path(RUN, "summary_by_replicate.csv"))
ci <- function(k, n) if (n > 0) as.numeric(stats::binom.test(k, n)$conf.int) else rep(NA_real_, 2)
summ <- byrep[, {
  pw <- if (all(is.na(target_called))) rep(NA_real_, 3) else { k <- sum(target_called); c(k / .N, ci(k, .N)) }
  fw <- { k <- sum(any_off); c(k / .N, ci(k, .N)) }
  .(n_rep = .N, power = pw[1], power_lo = pw[2], power_hi = pw[3],
    off_rate = fw[1], off_lo = fw[2], off_hi = fw[3], mean_off = mean(n_off),
    mean_called = mean(n_called))
}, by = .(scenario, archetype, size_class, scale, method)]
setorder(summ, scenario, size_class, scale, method)
fwrite(summ, file.path(RUN, "summary.csv"))
print(summ[, .(scenario, size_class, scale, method, n_rep, power = round(power, 2),
               off_rate = round(off_rate, 2), mean_called = round(mean_called, 1))], row.names = FALSE)
