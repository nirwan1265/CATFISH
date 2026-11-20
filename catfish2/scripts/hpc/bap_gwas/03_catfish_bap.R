#!/usr/bin/env Rscript
## =============================================================================
## 03_catfish_bap.R <trait> -- CATFISH Westfall-Young pathway analysis for one
## BAP trait, with the settings of the paper's sorghum analysis
## (scripts/hpc/gc_real/wy_gc_real.R, sorghum branch; bap_perm/05_wy_perm.R):
##   gene p-values  MAGMA snp-wise=mean, the 10 per-chromosome .genes.out of
##                  step 02 (column P), gene-level genomic control (median gene
##                  z removed), as in all the paper's real analyses
##   correlations   this trait's own .genes.raw -> gene_cor_pairs.txt (built
##                  once, banded, sorghum IDs ^SORBI), as bap_perm/04
##   pathways       SorghumBicolorCyc (.SORBI), 5-200 genes
##   test           catfish_wy_pathways(): B = 10000, five components, tau
##                  {0.01, 0.05, 0.5}, seed 20260706; self-contained and
##                  competitive (circular genomic permutation, genes in CHR/START
##                  order, B_comp = B)
## DEV_WY=<westfall_young.R> uses that file's catfish_wy_pathways() on top of the
## installed package (default: the gc_real dev copy, if present), as wy_gc_real.R.
## Writes /share/maize/ntanduk/BAP_gwas/catfish/<trait>/<trait>_wy.csv
##   Rscript 03_catfish_bap.R Stalk_height_cm
## =============================================================================
suppressPackageStartupMessages(library(CATFISH))
TRAIT <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(TRAIT)) stop("usage: Rscript 03_catfish_bap.R <trait>", call. = FALSE)
ROOT  <- Sys.getenv("BAP_ROOT", "/share/maize/ntanduk/BAP_gwas")
PWF   <- Sys.getenv("PATHWAY_FILE",
          "/share/maize/ntanduk/CATFISH/simulation/bap_perm_wy/inputs/sorghumbicolorcyc_pathways.20230103.SORBI")
DEV   <- Sys.getenv("DEV_WY", "/share/maize/ntanduk/CATFISH/gc_real/dev/westfall_young.R")
B     <- as.integer(Sys.getenv("WY_B", "10000"))
GC    <- Sys.getenv("GC", "1") == "1"
SEED  <- 20260706L; TAU <- c(0.01, 0.05, 0.5); MIN <- 5L; MAX <- 200L
CHROMS <- sprintf("%02d", 1:10)

if (nzchar(DEV) && file.exists(DEV)) {
  dev_env <- new.env(parent = asNamespace("CATFISH"))
  sys.source(DEV, envir = dev_env)
  catfish_wy_pathways <- dev_env$catfish_wy_pathways
  cat("using catfish_wy_pathways() from", DEV, "\n")
}
if (!"competitive" %in% names(formals(catfish_wy_pathways)))
  stop("this catfish_wy_pathways() has no competitive test: set DEV_WY or reinstall CATFISH", call. = FALSE)

MD  <- file.path(ROOT, "magma", TRAIT)
OUT <- file.path(ROOT, "catfish", TRAIT); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
fo  <- file.path(MD, sprintf("%s_chr%s.snp_wise_mean.genes.out", TRAIT, CHROMS))
fr  <- file.path(MD, sprintf("%s_chr%s.snp_wise_mean.genes.raw", TRAIT, CHROMS))
miss <- c(fo, fr, PWF)[!file.exists(c(fo, fr, PWF))]
if (length(miss)) stop("missing:\n  ", paste(miss, collapse = "\n  "), call. = FALSE)

## ---- gene statistics -------------------------------------------------------------
g <- do.call(rbind, lapply(fo, utils::read.table, header = TRUE, stringsAsFactors = FALSE, comment.char = "#"))
g <- g[!duplicated(g$GENE) & is.finite(g$P), ]
P_USE <- "P"; shift <- 0
if (GC) {
  z <- stats::qnorm(log(pmax(g$P, 1e-320)), lower.tail = FALSE, log.p = TRUE)
  shift <- stats::median(z)
  g$P_GC <- stats::pnorm(z - shift, lower.tail = FALSE); P_USE <- "P_GC"
}
cat(sprintf("[%s] %d genes, GC %s (median z %.3f removed), B = %d\n",
            TRAIT, nrow(g), if (GC) "on" else "off", shift, B))

## ---- gene-gene correlations from this trait's .genes.raw (once) -----------------
CORP <- file.path(OUT, "gene_cor_pairs.txt")
if (!file.exists(CORP) || file.size(CORP) == 0) {
  con <- file(CORP, open = "wt"); hdr <- FALSE
  for (f in fr) {
    tmp <- tempfile(fileext = ".txt")
    magma_genesraw_to_cor_pairs_banded(genes_raw_file = f, out_pairs_file = tmp,
      gene_regex = "^SORBI", keep_abs_r_ge = 0, overwrite = TRUE, verbose = FALSE)
    x <- if (file.exists(tmp)) readLines(tmp, warn = FALSE) else character(0)
    if (length(x) <= 1) { cat("  WARN no pairs from", basename(f), "\n"); next }
    if (hdr) x <- x[-1] else hdr <- TRUE
    writeLines(x, con)
  }
  close(con)
  if (length(readLines(CORP, warn = FALSE)) <= 1) stop("no gene-gene correlations extracted", call. = FALSE)
}
cat(sprintf("[%s] correlations: %s (%d lines)\n", TRAIT, CORP, length(readLines(CORP, warn = FALSE))))

## ---- pathways ------------------------------------------------------------------------
raw <- utils::read.delim(PWF, stringsAsFactors = FALSE, check.names = FALSE)
pw  <- unique(data.frame(pathway_id = raw[["Pathway-id"]], gene_id = raw[["Gene-name"]]))
pw  <- pw[nzchar(pw$gene_id) & pw$gene_id != "unknown", ]
pw_name <- tapply(raw[["Pathway-name"]], raw[["Pathway-id"]], function(x) x[1])

## ---- Westfall-Young: self-contained and competitive -----------------------------
t0 <- Sys.time()
res <- catfish_wy_pathways(g, pathways = pw, gene_col = "GENE", p_col = P_USE, chr_col = "CHR",
                           cor_file = CORP, B = B, use_tfisher = TRUE, tau_grid = TAU,
                           min_genes = MIN, max_genes = MAX, seed = SEED, verbose = FALSE,
                           competitive = TRUE, pos_col = "START")
res$pathway_name <- unname(pw_name[res$pathway_id])
res$gc_shift <- shift
fw <- file.path(OUT, sprintf("%s_wy.csv", TRAIT))
utils::write.csv(res, fw, row.names = FALSE)
cat(sprintf("[%s] %d pathways; FWER <= 0.05: self %d, competitive %d; top: %s (%s); %.1f min -> %s\n",
            TRAIT, nrow(res), sum(res$p_fwer_std <= 0.05), sum(res$p_fwer_comp <= 0.05),
            res$pathway_id[1], res$pathway_name[1],
            as.numeric(difftime(Sys.time(), t0, units = "mins")), fw))
