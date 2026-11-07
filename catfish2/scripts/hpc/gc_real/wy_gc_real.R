#!/usr/bin/env Rscript
## =============================================================================
## wy_gc_real.R -- the paper's Westfall-Young analyses of the real data, rerun
## with gene-level genomic control (GC): each analysis's gene z-scores are
## centred at their median before testing. Everything else is as in the runs
## they replace:
##   lipids   scripts/human_ldl/18_westfall_young_trait.R (checked against the
##            package by 19): P column, genes with finite P and ZSTAT, Reactome
##            sets, 10-200 genes
##   sorghum  catfish2/scripts/23_westfall_young_sorghum.R: P_SNPWISE_MEAN of the
##            per-chromosome multi=snp-wise run, SorghumBicolorCyc (SORBI), 5-200
## all: catfish_wy_pathways(), B = 10000, 5 components, tau {0.01, 0.05, 0.5},
## seed 20260706. GC=0 reruns without centring (reproduces the original files).
##   Rscript wy_gc_real.R <analysis>    e.g. ldl_eur, hdl_afr, Dry_tons_per_acre
## Reads/writes under GC_ROOT (default: this folder's parent of inputs/, out/).
## =============================================================================
suppressPackageStartupMessages(library(CATFISH))
A    <- commandArgs(trailingOnly = TRUE)[1]
ROOT <- Sys.getenv("GC_ROOT", "/share/maize/ntanduk/CATFISH/gc_real")
GC   <- Sys.getenv("GC", "1") == "1"
B    <- as.integer(Sys.getenv("WY_B", "10000"))
IN   <- file.path(ROOT, "inputs"); OUT <- file.path(ROOT, if (GC) "out" else "out_noGC")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SEED <- 20260706L; TAU <- c(0.01, 0.05, 0.5)
sorghum <- A %in% c("Dry_tons_per_acre", "Lignin_DM_perc")

if (!sorghum) {                                   # ---- GLGC lipids
  g <- read.table(file.path(IN, "lipids", paste0(A, ".genes.out")), header = TRUE,
                  comment.char = "#", stringsAsFactors = FALSE)
  g$GENE <- as.character(g$GENE)
  g <- g[is.finite(g$P) & is.finite(g$ZSTAT), ]
  pw <- list()
  for (ln in readLines(file.path(IN, "c2.cp.reactome.v2026.2.Hs.entrez.sets.annot"), warn = FALSE)) {
    f <- strsplit(ln, "\t", fixed = TRUE)[[1]]; if (length(f) < 2) next
    gs <- unique(f[-1]); pw[[f[1]]] <- gs[nzchar(gs)]
  }
  P_COL <- "P"; cor_file <- file.path(IN, "lipids", paste0(A, ".cor_pairs.txt"))
  MIN <- 10L; pw_name <- NULL
  out_name <- paste0(A, "_reactome_5comp_westfall_young.csv")
} else {                                          # ---- sorghum BAP
  fs <- file.path(IN, "sorghum", A, sprintf("%s_chr%02d.multi_snp_wise.genes.out", A, 1:10))
  g <- do.call(rbind, lapply(fs, read.table, header = TRUE, stringsAsFactors = FALSE, comment.char = "#"))
  g <- g[!duplicated(g$GENE), ]
  raw <- utils::read.delim(file.path(IN, "sorghumbicolorcyc_pathways.20230103.SORBI"),
                           stringsAsFactors = FALSE, check.names = FALSE)
  pw <- unique(data.frame(pathway_id = raw[["Pathway-id"]], gene_id = raw[["Gene-name"]]))
  pw <- pw[nzchar(pw$gene_id) & pw$gene_id != "unknown", ]
  pw_name <- tapply(raw[["Pathway-name"]], raw[["Pathway-id"]], function(x) x[1])
  P_COL <- "P_SNPWISE_MEAN"; cor_file <- file.path(IN, "sorghum", paste0(A, "_gene_cor_pairs.txt"))
  MIN <- 5L
  out_name <- paste0(A, "_P_SNPWISE_MEAN_5comp_westfall_young.csv")
}

shift <- 0; P_USE <- P_COL
if (GC) {
  z <- stats::qnorm(log(pmax(g[[P_COL]], 1e-320)), lower.tail = FALSE, log.p = TRUE)
  shift <- stats::median(z)
  g$P_GC <- stats::pnorm(z - shift, lower.tail = FALSE)
  P_USE <- "P_GC"
}
cat(sprintf("[%s] %d genes, column %s, GC %s (median z %.3f removed), B = %d\n",
            A, nrow(g), P_COL, if (GC) "on" else "off", shift, B))
t0 <- Sys.time()
res <- catfish_wy_pathways(g, pathways = pw, gene_col = "GENE", p_col = P_USE, chr_col = "CHR",
                           z_col = if (sorghum) NULL else "ZSTAT", cor_file = cor_file, B = B,
                           use_tfisher = TRUE, tau_grid = TAU, min_genes = MIN, max_genes = 200L,
                           seed = SEED, verbose = FALSE)
if (!is.null(pw_name)) res$pathway_name <- unname(pw_name[res$pathway_id])
res$gc_shift <- shift
write.csv(res, file.path(OUT, out_name), row.names = FALSE)
cat(sprintf("[%s] %d pathways, %d at FWER <= 0.05; top: %s; %.1f min -> %s\n", A, nrow(res),
            sum(res$p_fwer_std <= 0.05), res$pathway_id[1],
            as.numeric(difftime(Sys.time(), t0, units = "mins")), file.path(OUT, out_name)))
