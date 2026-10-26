#!/usr/bin/env Rscript
## =============================================================================
## 10_pair_coverage.R -- how complete is the gene-gene correlation file?
##
## Script 08 shows that both the marginal and the Westfall-Young calibration
## degrade as reported gene-pair correlations go missing, and that WY degrades
## about twice as fast. This script answers the question that makes 08
## actionable for a real dataset: where on that axis does MY correlation file
## sit?
##
## IMPORTANT DISTINCTION. A correlation matrix built from MAGMA output is
## SUPPOSED to be sparse. MAGMA only computes correlations for gene pairs close
## enough to share LD, so a zero between two genes 40 Mb apart is correct, not
## missing. The quantity that matters is narrower: of the within-pathway gene
## pairs that are on the SAME CHROMOSOME and WITHIN `WINDOW` bp of one another
## -- the pairs that ought to have a reported correlation -- what fraction
## actually do?
##
##   Rscript 10_pair_coverage.R
## Environment overrides: COR_PAIRS, GENE_LOC, PATHWAY_SETS, WINDOW, OUT_TAG
## =============================================================================
suppressPackageStartupMessages({ library(data.table) })

BASE <- Sys.getenv("CATFISH_BASE", path.expand("~/Documents/Github/MAGCAT"))
RES  <- Sys.getenv("CATFISH_RESULTS",
                   path.expand("~/Documents/Research/results/MAGMA/MAGCAT"))
MD   <- file.path(RES, "magma_multi_snp_wise_genes_by_chr_N_maize")

COR_PAIRS <- Sys.getenv("COR_PAIRS",
  file.path(MD, "magma_gene_cor_pairs_MLM_ALLCHR.txt"))
GENE_LOC  <- Sys.getenv("GENE_LOC",
  file.path(BASE, "catfish2/data/maize.genes.loc"))
SETS      <- Sys.getenv("PATHWAY_SETS",
  file.path(BASE, "catfish2/data/maize_pathway.sets.annot.WIDE.txt"))
WINDOW    <- as.numeric(Sys.getenv("WINDOW", "1e6"))
WINDOWS   <- as.numeric(strsplit(Sys.getenv("WINDOWS",
              "1e5,5e5,1e6,5e6,1e7,Inf"), ",")[[1]])
OUT_TAG   <- Sys.getenv("OUT_TAG", "maize_N_MLM")
OUT_DIR   <- file.path(BASE, "catfish2/results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

need <- function(f, what) if (!file.exists(f))
  stop(sprintf("%s not found: %s", what, f), call. = FALSE)
need(COR_PAIRS, "correlation pairs"); need(GENE_LOC, "gene locations")
need(SETS, "pathway sets")

## ----------------------------------------------------------------- inputs ---
loc <- data.table::fread(GENE_LOC, header = FALSE, data.table = FALSE)
names(loc)[1:4] <- c("gene", "chr", "start", "end")
## MAGMA .loc files sometimes carry a header row. Left in place it becomes a
## phantom chromosome and two NA coercion warnings.
if (!grepl("^[0-9]+$", as.character(loc$start[1]))) {
  cat("  (dropping a header row from the gene location file)\n")
  loc <- loc[-1L, , drop = FALSE]
}
loc <- loc[grepl("^[0-9.]+$", as.character(loc$start)) &
           grepl("^[0-9.]+$", as.character(loc$end)), , drop = FALSE]
loc$gene <- as.character(loc$gene); loc$chr <- as.character(loc$chr)
loc$mid  <- (as.numeric(loc$start) + as.numeric(loc$end)) / 2
CHR <- stats::setNames(loc$chr, loc$gene)
MID <- stats::setNames(loc$mid, loc$gene)
cat(sprintf("gene locations: %d genes on %d chromosomes\n",
            nrow(loc), length(unique(loc$chr))))

cp <- data.table::fread(COR_PAIRS, data.table = FALSE)
names(cp)[1:3] <- c("g1", "g2", "r")
cp$g1 <- as.character(cp$g1); cp$g2 <- as.character(cp$g2)
## unordered key so lookup does not depend on the order MAGMA wrote the pair
key <- paste0(pmin(cp$g1, cp$g2), "\t", pmax(cp$g1, cp$g2))
HAVE <- new.env(hash = TRUE, size = length(key))
for (k in key) assign(k, TRUE, envir = HAVE)
cat(sprintf("correlation pairs reported: %d\n", length(key)))

lines <- readLines(SETS, warn = FALSE)
lines <- lines[nzchar(trimws(lines))]
mem <- lapply(strsplit(trimws(lines), "[ \t]+"), function(x) x)
pw  <- stats::setNames(lapply(mem, function(x) unique(x[-1L])),
                       vapply(mem, `[`, character(1), 1L))
cat(sprintf("pathways: %d\n", length(pw)))

## ------------------------------------------------------------- per pathway ---
rows <- vector("list", length(pw))
for (i in seq_along(pw)) {
  g <- pw[[i]]
  g <- g[g %in% names(CHR)]
  n <- length(g)
  if (n < 2L) {
    rows[[i]] <- data.frame(pathway = names(pw)[i], n_genes = n,
      n_pairs_total = 0L, n_pairs_expected = 0L, n_pairs_reported = 0L,
      coverage = NA_real_, stringsAsFactors = FALSE)
    next
  }
  cmb <- utils::combn(g, 2L)
  same <- CHR[cmb[1, ]] == CHR[cmb[2, ]]
  near <- same & (abs(MID[cmb[1, ]] - MID[cmb[2, ]]) <= WINDOW)
  near[is.na(near)] <- FALSE
  exp_pairs <- which(near)
  got <- 0L
  if (length(exp_pairs)) {
    kk <- paste0(pmin(cmb[1, exp_pairs], cmb[2, exp_pairs]), "\t",
                 pmax(cmb[1, exp_pairs], cmb[2, exp_pairs]))
    got <- sum(vapply(kk, function(z) exists(z, envir = HAVE, inherits = FALSE),
                      logical(1)))
  }
  rows[[i]] <- data.frame(
    pathway = names(pw)[i], n_genes = n,
    n_pairs_total = ncol(cmb),
    n_pairs_expected = length(exp_pairs),
    n_pairs_reported = got,
    coverage = if (length(exp_pairs)) got / length(exp_pairs) else NA_real_,
    stringsAsFactors = FALSE)
}
out <- do.call(rbind, rows)
f <- file.path(OUT_DIR, sprintf("10_pair_coverage_%s.csv", OUT_TAG))
utils::write.csv(out, f, row.names = FALSE)

## ---------------------------------------------------------------- summary ---
ok <- out[!is.na(out$coverage) & out$n_pairs_expected > 0, ]
overall <- sum(out$n_pairs_reported) / max(1, sum(out$n_pairs_expected))
cat(sprintf("\nwrote %s\n", f))
cat(sprintf("\nwindow for 'ought to be correlated': %.0f bp\n", WINDOW))
cat(sprintf("within-pathway pairs on the same chromosome and within window: %d\n",
            sum(out$n_pairs_expected)))
cat(sprintf("of which a correlation is reported: %d\n", sum(out$n_pairs_reported)))
cat(sprintf("\nOVERALL PAIR COVERAGE: %.1f%%\n", 100 * overall))
if (nrow(ok)) {
  q <- stats::quantile(ok$coverage, c(0, .25, .5, .75, 1))
  cat(sprintf("per-pathway coverage: min %.2f  Q1 %.2f  median %.2f  Q3 %.2f  max %.2f\n",
              q[1], q[2], q[3], q[4], q[5]))
  cat(sprintf("pathways below 50%% coverage: %d of %d\n",
              sum(ok$coverage < 0.5), nrow(ok)))
}
## ------------------------------------------------------ window sensitivity ---
## 377 qualifying pairs at 1 Mb is a thin basis for a coverage estimate, so
## report the same quantity across a range of windows. "Inf" means every
## same-chromosome within-pathway pair regardless of distance, which is the
## most pessimistic denominator: many of those zeros are correct rather than
## missing, so coverage there is a lower bound.
cat("\nwindow sensitivity:\n")
cat(sprintf("  %12s %10s %10s %9s\n", "window", "expected", "reported", "coverage"))
for (W in WINDOWS) {
  tot <- 0L; got <- 0L
  for (i in seq_along(pw)) {
    g <- pw[[i]]; g <- g[g %in% names(CHR)]
    if (length(g) < 2L) next
    cmb <- utils::combn(g, 2L)
    near <- (CHR[cmb[1, ]] == CHR[cmb[2, ]]) &
            (abs(MID[cmb[1, ]] - MID[cmb[2, ]]) <= W)
    near[is.na(near)] <- FALSE
    ii <- which(near)
    if (!length(ii)) next
    kk <- paste0(pmin(cmb[1, ii], cmb[2, ii]), "\t", pmax(cmb[1, ii], cmb[2, ii]))
    tot <- tot + length(ii)
    got <- got + sum(vapply(kk, function(z) exists(z, envir = HAVE,
                            inherits = FALSE), logical(1)))
  }
  cat(sprintf("  %12s %10d %10d %8.1f%%\n",
              if (is.infinite(W)) "any (same chr)" else format(W, scientific = TRUE),
              tot, got, if (tot) 100 * got / tot else NA_real_))
}

cat("\nRead this against results/08_missing_correlation.csv: coverage near 100%\n")
cat("places you at the left edge of that table, where Westfall-Young control is\n")
cat("close to nominal; coverage near 50% places you where it is roughly 3x\n")
cat("nominal and the family-wise p-values should be reported with that caveat.\n")
