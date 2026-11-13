#!/usr/bin/env Rscript
## =============================================================================
## 00_prepare_sets.R -- MAGMA annotation and Reactome sets for the simulation
## panel (run by 00_prepare_panel_LSF.sh).
##   SIM_ANNOT  MAGMA's EUR annotation (+/- 10 kb), each gene keeping only the
##              panel's SNPs; genes on chr 1-10 with at least one SNP
##   SIM_SETS   Reactome sets restricted to those genes, 10-200 genes kept
##   genes.tsv  gene, chr, start, stop (from the annotation) for step 01 and 08
## =============================================================================
suppressPackageStartupMessages(library(data.table))
env <- function(k) { v <- Sys.getenv(k); if (!nzchar(v)) stop(k, " is not set"); v }
IN <- env("IN"); CHROMS <- as.integer(strsplit(env("CHROMS"), " ")[[1]])
MIN <- as.integer(env("MIN_GENES")); MAX <- as.integer(env("MAX_GENES"))

snps <- fread(paste0(env("SIM_PANEL"), ".bim"), header = FALSE, select = 2L)[[1]]
ln <- readLines(env("ANNOT_FULL"))
hdr <- ln[startsWith(ln, "#")]; ln <- ln[!startsWith(ln, "#")]
out <- character(0); genes <- list()
for (x in ln) {
  f <- strsplit(x, "\t", fixed = TRUE)[[1]]
  loc <- strsplit(f[2], ":", fixed = TRUE)[[1]]
  if (!(as.integer(loc[1]) %in% CHROMS)) next
  s <- f[-(1:2)]; s <- s[s %in% snps]
  if (!length(s)) next
  out <- c(out, paste(c(f[1], f[2], s), collapse = "\t"))
  genes[[f[1]]] <- data.table(gene = f[1], chr = as.integer(loc[1]),
                              start = as.integer(loc[2]), stop = as.integer(loc[3]), n_snps = length(s))
}
writeLines(c(hdr, out), env("SIM_ANNOT"))
genes <- rbindlist(genes)
fwrite(genes, file.path(IN, "genes.tsv"), sep = "\t")
cat(sprintf("[sets] %d genes on chr %s with panel SNPs\n", nrow(genes), paste(range(CHROMS), collapse = "-")))

keep <- character(0); n_in <- 0L
for (x in readLines(env("REACTOME_SETS"))) {
  f <- strsplit(x, "\t", fixed = TRUE)[[1]]; if (length(f) < 2) next
  n_in <- n_in + 1L
  g <- intersect(unique(f[-1]), genes$gene)
  if (length(g) >= MIN && length(g) <= MAX) keep <- c(keep, paste(c(f[1], g), collapse = "\t"))
}
writeLines(keep, env("SIM_SETS"))
cat(sprintf("[sets] %d of %d Reactome sets have %d-%d genes on these chromosomes -> %s\n",
            length(keep), n_in, MIN, MAX, env("SIM_SETS")))
