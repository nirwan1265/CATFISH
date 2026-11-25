#!/usr/bin/env Rscript
## =============================================================================
## 22_gsea_lipids.R -- preranked GSEA (fgsea) on the MAGMA gene z-scores of the
## ten lipid sets, with the Reactome sets CATFISH tested (10-200 genes).
##   ranking   ZSTAT from MAGMA's merged .genes.raw (inputs/lipids/<ta>.genes.raw),
##             the gene statistics CATFISH, MAGMA gene-set and PascalX used;
##             largest z = most associated; scoreType "pos"
##   null      fgseaMultilevel (genes exchangeable; LD between genes not modelled)
## Writes out/gsea/<ta>.tsv (pathway_id, n_genes, NES, p, padj) and
##        out/gsea/summary.tsv: per set, pathways tested, Bonferroni and BH (0.05)
##        calls, and the rank of the plasma lipoprotein pathway.
## Needs  fgsea in FGSEA_RLIB (or the active R library).
##   Rscript 22_gsea_lipids.R
## =============================================================================
suppressPackageStartupMessages(library(data.table))
env <- function(k) { v <- Sys.getenv(k); if (!nzchar(v)) stop(k, " is not set"); v }
lib <- Sys.getenv("FGSEA_RLIB"); if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(fgsea))
ROOT <- env("ROOT"); OUT <- file.path(env("OUT"), "gsea"); dir.create(OUT, showWarnings = FALSE)
SEED <- as.integer(Sys.getenv("GSEA_SEED", "20261125"))
LIPO <- "REACTOME_PLASMA_LIPOPROTEIN_ASSEMBLY_REMODELING_AND_CLEARANCE"

sets <- list()
for (x in readLines(env("REACTOME_SETS"))) {
  f <- strsplit(x, "[\t ]+")[[1]]; if (length(f) > 1) sets[[f[1]]] <- f[-1]
}
## .genes.raw rows: GENE CHR START STOP NSNPS NPARAM, then one field per covariate
## listed in the "# COVAR =" header (e.g. NSAMP MAC), then ZSTAT, then correlations.
## The ZSTAT column is therefore 6 + (number of covariates) + 1.
read_raw <- function(f) {
  ln <- readLines(f)
  cv <- grep("^# COVAR", ln, value = TRUE)
  ncov <- if (length(cv)) length(strsplit(trimws(sub("^# COVAR *= *", "", cv[1])), "[\t ]+")[[1]]) else 1L
  zc <- 6L + ncov + 1L
  ln <- ln[!startsWith(ln, "#")]
  pp <- strsplit(ln, "[\t ]+")
  cat(sprintf("  %s: %d covariates, ZSTAT in field %d\n", basename(f), ncov, zc))
  data.table(gene = vapply(pp, `[`, "", 1L), z = as.numeric(vapply(pp, `[`, "", zc)))
}
summ <- list()
for (ta in strsplit(env("LIPIDS"), " ")[[1]]) {
  f <- file.path(ROOT, "inputs", "lipids", paste0(ta, ".genes.raw"))
  if (!file.exists(f)) { cat("missing", f, "\n"); next }
  g <- read_raw(f)[is.finite(z)]
  st <- setNames(g$z, g$gene)
  set.seed(SEED); st <- st + stats::runif(length(st), 0, 1e-9)   # strict, reproducible ranking
  res <- fgseaMultilevel(pathways = sets, stats = st, minSize = 10, maxSize = 200,
                         scoreType = "pos", eps = 0, nproc = 1)
  d <- data.table(pathway_id = res$pathway, n_genes = res$size, NES = res$NES, p = res$pval, padj = res$padj)
  setorder(d, p, -NES)
  fwrite(d, file.path(OUT, paste0(ta, ".tsv")), sep = "\t")
  K <- nrow(d)
  summ[[ta]] <- data.table(set = ta, n_pathways = K, n_bonf = sum(d$p <= 0.05 / K), n_bh05 = sum(d$padj <= 0.05),
                           lipo_rank = match(LIPO, d$pathway_id), lipo_p = d[pathway_id == LIPO, p])
  print(summ[[ta]])
}
fwrite(rbindlist(summ), file.path(OUT, "summary.tsv"), sep = "\t")
