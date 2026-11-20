#!/usr/bin/env Rscript
## =============================================================================
## 07c_acat.R <cell> -- two-stage ACAT gene-set test (Liu et al. 2019, ACAT;
## Liu & Xie 2020, Cauchy combination) for every replicate of one design cell.
##   stage 1  SNP -> gene:      ACAT of the SNP p-values assigned to the gene
##                              (SIM_ANNOT, the MAGMA annotation, +/- 10 kb)
##   stage 2  gene -> pathway:  ACAT of the gene p-values of the pathway (SIM_SETS)
## Equal weights; the Cauchy combination is valid under arbitrary dependence in
## the tail, so no LD correction is applied. Calls in 08: p <= 0.05 / K
## (Bonferroni), as MAGMA. Self-contained: it asks whether any gene/SNP is
## associated, like MAGMA's self-contained test.
## Reads  RUN/design.tsv, RUN/z/z_cellNN.tsv.gz, SIM_ANNOT, SIM_SETS
## Writes RUN/acat/rep_XXXX.tsv (pathway_id, n_genes, p) and
##        RUN/acat/rep_XXXX.genes.tsv.gz (gene, n_snps, p)
##   Rscript 07c_acat.R <cell>
## =============================================================================
suppressPackageStartupMessages(library(data.table))
env <- function(k) { v <- Sys.getenv(k); if (!nzchar(v)) stop(k, " is not set"); v }
CELL <- as.integer(commandArgs(trailingOnly = TRUE)[1])
RUN <- env("RUN"); OUT <- file.path(RUN, "acat"); dir.create(OUT, showWarnings = FALSE)

## Cauchy combination of p-values with equal weights; numerically safe at both ends
acat <- function(p) {
  p <- p[is.finite(p)]; if (!length(p)) return(NA_real_)
  p <- pmin(pmax(p, 1e-300), 1 - 1e-15)
  tt <- ifelse(p < 1e-15, 1 / (p * pi), tan((0.5 - p) * pi))
  T <- mean(tt)
  if (T > 1e15) 1 / (T * pi) else 0.5 - atan(T) / pi
}

design <- fread(file.path(RUN, "design.tsv")); design <- design[design$cell == CELL]
if (!nrow(design)) stop("no replicates in cell ", CELL)

## annotation: gene -> SNP row indices in the z file
z <- fread(file.path(RUN, "z", sprintf("z_cell%02d.tsv.gz", CELL)))
ln <- readLines(env("SIM_ANNOT")); ln <- ln[!startsWith(ln, "#")]
parts <- strsplit(ln, "[\t ]+")
gene_id <- vapply(parts, `[`, "", 1L)
snp_idx <- lapply(parts, function(x) { i <- match(x[-(1:2)], z$SNP); i[!is.na(i)] })
keep <- lengths(snp_idx) > 0
gene_id <- gene_id[keep]; snp_idx <- snp_idx[keep]
sets <- list()
for (x in readLines(env("SIM_SETS"))) { f <- strsplit(x, "\t", fixed = TRUE)[[1]]; sets[[f[1]]] <- f[-1] }
set_idx <- lapply(sets, function(g) { i <- match(g, gene_id); i[!is.na(i)] })
cat(sprintf("[acat] cell %d: %d replicates, %d genes, %d pathways\n",
            CELL, nrow(design), length(gene_id), length(sets)))

for (r in design$rep_id) {
  RP <- sprintf("rep_%04d", r); fo <- file.path(OUT, paste0(RP, ".tsv"))
  if (file.exists(fo) && file.size(fo) > 0) { cat(RP, "already done\n"); next }
  p_snp <- 2 * stats::pnorm(-abs(z[[RP]]))
  p_gene <- vapply(snp_idx, function(i) acat(p_snp[i]), numeric(1))
  p_set <- vapply(set_idx, function(i) acat(p_gene[i]), numeric(1))
  fwrite(data.table(gene = gene_id, n_snps = lengths(snp_idx), p = p_gene),
         file.path(OUT, paste0(RP, ".genes.tsv.gz")), sep = "\t")
  fwrite(data.table(pathway_id = names(sets), n_genes = lengths(set_idx), p = p_set), fo, sep = "\t")
  cat(sprintf("[acat] %s: min pathway p %.3g\n", RP, min(p_set, na.rm = TRUE)))
}
