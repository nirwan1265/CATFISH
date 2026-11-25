#!/usr/bin/env Rscript
## =============================================================================
## 07d_gsea.R <cell> -- preranked GSEA (fgsea; Korotkevich et al. 2021) on the
## MAGMA gene z-scores of every replicate of one design cell.
##   ranking    MAGMA ZSTAT of each gene (RUN/magma/rep_XXXX/rep_XXXX.genes.out),
##              largest = most associated; scoreType "pos" (enrichment at the top)
##   pathways   SIM_SETS (Reactome on chr 1-10, 10-200 genes), as every method
##   null       fgseaMultilevel: random gene sets of the same size drawn from the
##              ranked list, i.e. genes treated as exchangeable; LD between the
##              genes of a pathway is not modelled. eps = 0 for exact small p.
## Calls in 08: p <= 0.05 / K (Bonferroni), as MAGMA and PascalX. padj (BH) is
## kept in the output for reference. Competitive in kind: a pathway is compared
## with the rest of the ranked genes.
## Reads  RUN/design.tsv, RUN/magma/rep_XXXX/rep_XXXX.genes.out, SIM_SETS
## Writes RUN/gsea/rep_XXXX.tsv (pathway_id, n_genes, NES, p, padj)
## Needs  R package fgsea (Bioconductor) in FGSEA_RLIB or the active library.
##   Rscript 07d_gsea.R <cell>
## =============================================================================
suppressPackageStartupMessages(library(data.table))
env <- function(k) { v <- Sys.getenv(k); if (!nzchar(v)) stop(k, " is not set"); v }
lib <- Sys.getenv("FGSEA_RLIB"); if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(fgsea))
CELL <- as.integer(commandArgs(trailingOnly = TRUE)[1])
RUN <- env("RUN"); OUT <- file.path(RUN, "gsea"); dir.create(OUT, showWarnings = FALSE)
MINS <- as.integer(env("MIN_GENES")); MAXS <- as.integer(env("MAX_GENES"))
SEED <- as.integer(Sys.getenv("GSEA_SEED", "20261125"))

design <- fread(file.path(RUN, "design.tsv")); design <- design[design$cell == CELL]
if (!nrow(design)) stop("no replicates in cell ", CELL)
sets <- list()
for (x in readLines(env("SIM_SETS"))) { f <- strsplit(x, "\t", fixed = TRUE)[[1]]; sets[[f[1]]] <- f[-1] }
cat(sprintf("[gsea] cell %d: %d replicates, %d pathways\n", CELL, nrow(design), length(sets)))

for (r in design$rep_id) {
  RP <- sprintf("rep_%04d", r); fo <- file.path(OUT, paste0(RP, ".tsv"))
  if (file.exists(fo) && file.size(fo) > 0) { cat(RP, "already done\n"); next }
  gf <- file.path(RUN, "magma", RP, paste0(RP, ".genes.out"))
  if (!file.exists(gf)) { cat(RP, "no MAGMA genes.out, skipped\n"); next }
  g <- fread(gf, colClasses = list(character = "GENE"))
  st <- setNames(g$ZSTAT, g$GENE); st <- st[is.finite(st)]
  ## break ties so the ranking is strict and reproducible
  set.seed(SEED + r); st <- st + stats::runif(length(st), 0, 1e-9)
  res <- fgseaMultilevel(pathways = sets, stats = st, minSize = MINS, maxSize = MAXS,
                         scoreType = "pos", eps = 0, nproc = 1)
  fwrite(data.table(pathway_id = res$pathway, n_genes = res$size, NES = res$NES,
                    p = res$pval, padj = res$padj), fo, sep = "\t")
  cat(sprintf("[gsea] %s: %d pathways, min p %.3g\n", RP, nrow(res), min(res$pval, na.rm = TRUE)))
}
