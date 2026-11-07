#!/usr/bin/env Rscript
## =============================================================================
## 05_wy_perm.R -- Westfall-Young pathway analysis for ONE permutation, with the
## settings of catfish2/scripts/23_westfall_young_sorghum.R (catfish_wy_pathways,
## B = WY_B, five components, tau grid {0.01, 0.05, 0.5}, seed 20260706,
## SorghumBicolorCyc, MIN_GENES-MAX_GENES genes).
##   Rscript 05_wy_perm.R <perm>
## Gene p-value: P_SNPWISE_MEAN when MAGMA ran multi=snp-wise, P when it ran
## snp-wise=mean (the same quantity). Writes WY_OUT/perm_<b>_wy.csv.
## =============================================================================
getcfg <- function(k) { v <- Sys.getenv(k); if (v == "") stop("config missing: ", k); v }
B_PERM   <- as.integer(commandArgs(trailingOnly = TRUE)[1]); BPAD <- sprintf("%03d", B_PERM)
suppressPackageStartupMessages(library(CATFISH))   # installed in the R these jobs use
if (!exists("catfish_wy_pathways")) stop("this CATFISH has no catfish_wy_pathways(); reinstall it from GitHub")
MAGMA_OUT <- getcfg("MAGMA_OUT"); WY_OUT <- getcfg("WY_OUT"); CORPAIRS <- getcfg("CORPAIRS")
GENE_MODEL <- getcfg("GENE_MODEL"); CHROMS <- strsplit(trimws(getcfg("CHROMS")), "\\s+")[[1]]
model_tag <- gsub("[^A-Za-z0-9]+", "_", GENE_MODEL)
WY_B <- as.integer(getcfg("WY_B")); MIN_GENES <- as.integer(getcfg("MIN_GENES"))
MAX_GENES <- as.integer(getcfg("MAX_GENES"))
SEED <- 20260706; TAU_GRID <- c(0.01, 0.05, 0.5)

files <- file.path(MAGMA_OUT, sprintf("perm_%s", BPAD), sprintf("chr%s", CHROMS),
                   sprintf("perm_%s_chr%s.%s.genes.out", BPAD, CHROMS, model_tag))
miss <- files[!file.exists(files)]
if (length(miss)) stop("missing MAGMA output:\n  ", paste(miss, collapse = "\n  "))
g <- do.call(rbind, lapply(files, function(f)
  utils::read.table(f, header = TRUE, stringsAsFactors = FALSE, comment.char = "#")))
g <- g[!duplicated(g$GENE), ]
P_COL <- if ("P_SNPWISE_MEAN" %in% names(g)) "P_SNPWISE_MEAN" else "P"
if (!file.exists(CORPAIRS)) stop("missing ", CORPAIRS, " (run step 4)")
## GC=1: gene-level genomic control -- this permutation's gene z-scores are
## centred at their median before testing (z from P on the log scale)
GC <- Sys.getenv("GC", "0") == "1"
P_USE <- P_COL; gc_shift <- 0
if (GC) {
  z <- stats::qnorm(log(pmax(g[[P_COL]], 1e-320)), lower.tail = FALSE, log.p = TRUE)
  gc_shift <- stats::median(z)
  g$P_GC <- stats::pnorm(z - gc_shift, lower.tail = FALSE)
  P_USE <- "P_GC"
  cat(sprintf("[wy] GC: median gene z %.3f removed\n", gc_shift))
}
cat(sprintf("[wy] perm %d: %d genes, column %s, median gene p %.3f\n",
            B_PERM, nrow(g), P_COL, stats::median(g[[P_COL]], na.rm = TRUE)))

pw_raw <- utils::read.delim(getcfg("PATHWAY_FILE"), stringsAsFactors = FALSE, check.names = FALSE)
pw <- unique(data.frame(pathway_id = pw_raw[["Pathway-id"]], gene_id = pw_raw[["Gene-name"]],
                        stringsAsFactors = FALSE))
pw <- pw[nzchar(pw$gene_id) & pw$gene_id != "unknown", ]

res <- catfish_wy_pathways(
  gene_results = g, pathways = pw, gene_col = "GENE", p_col = P_USE, chr_col = "CHR",
  cor_file = CORPAIRS, B = WY_B, use_tfisher = TRUE, tau_grid = TAU_GRID,
  min_genes = MIN_GENES, max_genes = MAX_GENES, seed = SEED, verbose = FALSE)
res$perm <- B_PERM
res$median_gene_p <- stats::median(g[[P_COL]], na.rm = TRUE)   # before any GC
res$gc_shift <- gc_shift
dir.create(WY_OUT, recursive = TRUE, showWarnings = FALSE)
out <- file.path(WY_OUT, sprintf("perm_%s_wy.csv", BPAD))
utils::write.csv(res, out, row.names = FALSE)
cat(sprintf("[wy] %d pathways; p_calibrated < 0.05: %d; FWER-significant: %d -> %s\n",
            nrow(res), sum(res$p_calibrated < 0.05, na.rm = TRUE),
            sum(res$p_fwer_std <= 0.05, na.rm = TRUE), out))
