#!/usr/bin/env Rscript
## =============================================================================
## 16_ablation_multitrait.R -- does the pathway layer reorder the ranking in a
## way that helps, across traits?
##
## THE QUESTION
##   The composite score of Section 2.5 adds pathway evidence to gene-level
##   evidence. Whether that addition EARNS its place is testable: score each
##   gene with and without the pathway terms, and compare where a set of
##   pre-specified, literature-established causal genes lands under each.
##
##     genes only  :  S_GWAS + S_MAGMA
##     full        :  S_GWAS + S_MAGMA + 0.5*S_PATH + 0.3*B_hub
##
##   A single trait is an anecdote. This runs the comparison over several lipid
##   traits that share a pipeline, a reference panel and a pathway database, so
##   the only thing varying is the trait.
##
## WHAT IT CANNOT SHOW
##   The control genes are curated from the literature, not established by
##   experiment in these data. That is weaker than knockout validation and the
##   write-up should say so. It also means the controls are biased towards
##   well-studied genes, which are exactly the genes most likely to be
##   well-annotated in Reactome -- a bias that FAVOURS the pathway layer. The
##   per-trait breakdown is reported so that a reader can judge it.
##
##   Rscript 16_ablation_multitrait.R [trait ...]
## =============================================================================
suppressPackageStartupMessages({ library(data.table) })

RES <- "../results"
OUT <- file.path(RES, "16_ablation_multitrait.csv")

## ---- pre-specified positive controls ---------------------------------------
## Curated before running, from established lipid genetics: monogenic
## dyslipidaemia genes, validated drug targets, and loci with functional
## follow-up. Not derived from these data.
CONTROLS <- list(
  ldl   = c("LDLR","PCSK9","APOB","APOE","SORT1","HMGCR","NPC1L1","ABCG8",
            "LPA","CETP","LDLRAP1","MYLIP","ANGPTL3"),
  hdl   = c("CETP","LIPC","LIPG","ABCA1","APOA1","APOA5","LCAT","SCARB1",
            "PLTP","LPL","ANGPTL4","CUBN"),
  logtg = c("LPL","APOA5","APOC3","ANGPTL3","ANGPTL4","ANGPTL8","GCKR",
            "LMF1","GPIHBP1","APOE","TRIB1","MLXIPL"),
  tc    = c("LDLR","APOE","APOB","PCSK9","HMGCR","CETP","ABCG8","NPC1L1",
            "SORT1","TM6SF2","LPA","TRIB1")
)

args <- commandArgs(trailingOnly = TRUE)
traits <- if (length(args)) tolower(args) else names(CONTROLS)

rank_by <- function(d, score) {
  o <- order(-score)
  r <- integer(nrow(d)); r[o] <- seq_len(nrow(d))
  r
}

res <- list(); per_gene <- list()
for (tr in traits) {
  f <- file.path(RES, sprintf("06_candidate_genes_%s.csv", tr))
  if (!file.exists(f)) {
    cat(sprintf("skip %-6s (no %s)\n", tr, basename(f))); next
  }
  d <- data.table::fread(f, data.table = FALSE)
  if (!"symbol" %in% names(d)) { cat(sprintf("skip %-6s (no symbols)\n", tr)); next }

  s_gen  <- d$S_GWAS + d$S_MAGMA
  s_full <- d$S_GWAS + d$S_MAGMA + 0.5 * d$S_PATH + 0.3 * d$B_hub
  r_gen  <- rank_by(d, s_gen)
  r_full <- rank_by(d, s_full)

  ctrl <- CONTROLS[[tr]]
  idx  <- match(ctrl, d$symbol)
  ok   <- !is.na(idx); idx <- idx[ok]
  if (!length(idx)) { cat(sprintf("skip %-6s (no controls matched)\n", tr)); next }

  N <- nrow(d)
  per_gene[[tr]] <- data.frame(
    trait = tr, gene = ctrl[ok],
    rank_genes_only = r_gen[idx], rank_full = r_full[idx],
    change = r_gen[idx] - r_full[idx], n_genes = N,
    stringsAsFactors = FALSE)

  res[[tr]] <- data.frame(
    trait = tr, n_genes = N, n_controls = length(idx),
    median_genes_only = stats::median(r_gen[idx]),
    median_full       = stats::median(r_full[idx]),
    mean_genes_only   = mean(r_gen[idx]),
    mean_full         = mean(r_full[idx]),
    top1pct_genes_only = sum(r_gen[idx]  <= N * 0.01),
    top1pct_full       = sum(r_full[idx] <= N * 0.01),
    improved = sum(r_full[idx] < r_gen[idx]),
    worsened = sum(r_full[idx] > r_gen[idx]),
    ## median of the per-gene rank change: positive means the pathway layer
    ## moved controls up. Reported as a median because a single large move
    ## (MYLIP in the LDL analysis) dominates a mean.
    median_change = stats::median(r_gen[idx] - r_full[idx]),
    stringsAsFactors = FALSE)

  cat(sprintf("%-6s n=%d controls=%2d  median rank %5.0f -> %5.0f  (improved %d, worsened %d)\n",
              tr, N, length(idx), stats::median(r_gen[idx]),
              stats::median(r_full[idx]), res[[tr]]$improved, res[[tr]]$worsened))
}

if (!length(res)) stop("no traits had a candidate-gene table", call. = FALSE)

out <- do.call(rbind, res)
utils::write.csv(out, OUT, row.names = FALSE)
pg <- do.call(rbind, per_gene)
utils::write.csv(pg, file.path(RES, "16_ablation_per_gene.csv"), row.names = FALSE)
cat(sprintf("\nwrote %s\n       %s\n", OUT, file.path(RES, "16_ablation_per_gene.csv")))

## ---- pooled test ------------------------------------------------------------
## Paired comparison across every control gene in every trait. A one-sided
## Wilcoxon signed-rank test on the paired ranks: does adding the pathway layer
## move controls up more often than down?
cat("\npooled across traits:\n")
cat(sprintf("  control genes: %d across %d traits\n", nrow(pg), nrow(out)))
cat(sprintf("  median rank  %.0f -> %.0f\n",
            stats::median(pg$rank_genes_only), stats::median(pg$rank_full)))
cat(sprintf("  improved %d, worsened %d, unchanged %d\n",
            sum(pg$change > 0), sum(pg$change < 0), sum(pg$change == 0)))
if (sum(pg$change != 0) >= 5) {
  w <- suppressWarnings(stats::wilcox.test(pg$rank_genes_only, pg$rank_full,
                                           paired = TRUE, alternative = "greater"))
  cat(sprintf("  Wilcoxon signed-rank (paired, one-sided): V = %.0f, p = %.4g\n",
              unname(w$statistic), w$p.value))
  cat("  (tests only whether the layer helps on average; it is not a\n")
  cat("   statement about any individual gene)\n")
}
print(pg[order(-pg$change), c("trait","gene","rank_genes_only","rank_full","change")],
      row.names = FALSE)
