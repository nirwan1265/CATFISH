#!/usr/bin/env Rscript
## =============================================================================
## 23_westfall_young_sorghum.R -- Westfall-Young pathway analysis for the
## sorghum Bioenergy Association Panel (BAP) traits, using the CATFISH package.
##
## Same analysis as scripts/human_ldl/18_westfall_young_trait.R (the package
## function catfish_wy_pathways() was verified to reproduce 18 exactly), on:
##   gene results   MAGMA multi=snp-wise, run per chromosome:
##                  <MAGMA_ROOT>/<TRAIT>/<TRAIT>_chrNN.multi_snp_wise.genes.out
##   correlations   the matching per-chromosome .genes.raw files (MAGMA only
##                  computes gene-gene correlations within a chromosome, so the
##                  ten files together carry the full correlation structure)
##   pathways       SorghumBicolorCyc (PMN), SORBI gene IDs
##
## Usage (from catfish2/scripts):
##   TRAIT=Dry_tons_per_acre Rscript 23_westfall_young_sorghum.R [B] [use_tfisher] [P_COL]
##     B            null replicates             (default 10000)
##     use_tfisher  TRUE/FALSE                  (default TRUE)
##     P_COL        gene p-value column of the MAGMA output (default P_SNPWISE_MEAN)
##                  P_MULTI | P_SNPWISE_MEAN | P_SNPWISE_TOP1
##
##   P_COL matters. For both BAP traits the SNP-level GWAS is well calibrated
##   (lambda_GC 0.86 dry tons, 0.91 lignin), and so is MAGMA's SNP-wise mean
##   model (median gene p 0.447 and 0.449), but the SNP-wise top model is
##   inflated (median 0.112 and 0.263), and P_MULTI, which combines the two,
##   inherits it (0.179 and 0.322). The Westfall-Young null assumes the
##   observed gene p-values are uniform where there is no signal, so an
##   inflated gene-level p is read as pathway signal. The script prints the
##   median gene p of every column so the choice can be checked.
##
## Env overrides: TRAIT, MAGMA_ROOT, PATHWAY_FILE, MIN_GENES (5), MAX_GENES (200),
##                CHROMS (default all ten; e.g. CHROMS=8 for a quick check)
## Writes ../results/sorghum_bap/<TRAIT>_<P_COL>_<k>comp_westfall_young.csv
##        ../results/sorghum_bap/<TRAIT>_gene_cor_pairs.txt (built once, reused)
## =============================================================================
suppressPackageStartupMessages(library(CATFISH))

args <- commandArgs(trailingOnly = TRUE)
B           <- if (length(args) >= 1) as.integer(args[[1]]) else 10000L
USE_TFISHER <- if (length(args) >= 2) as.logical(args[[2]]) else TRUE
P_COL       <- if (length(args) >= 3) args[[3]] else "P_SNPWISE_MEAN"
SEED        <- 20260706                       # same seed as 18
TAU_GRID    <- c(0.01, 0.05, 0.5)             # same grid as 18

TRAIT      <- Sys.getenv("TRAIT", "Dry_tons_per_acre")
MAGMA_ROOT <- Sys.getenv("MAGMA_ROOT",
                path.expand("~/Documents/Research/results/CATFISH/MAGMA"))
BASE       <- Sys.getenv("CATFISH_BASE", path.expand("~/Documents/Github/MAGCAT"))
PATHWAY_FILE <- Sys.getenv("PATHWAY_FILE",
                  file.path(BASE, "inst/extdata/pathway/sorghumbicolorcyc_pathways.20230103.SORBI"))
## pathway size window: the same as the maize (PMN) branch of 18
MIN_GENES  <- as.integer(Sys.getenv("MIN_GENES", "5"))
MAX_GENES  <- as.integer(Sys.getenv("MAX_GENES", "200"))
OUT_DIR    <- "../results/sorghum_bap"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time()

## ------------------------------------------------------------ gene results --
tdir <- file.path(MAGMA_ROOT, TRAIT)
CHROMS <- as.integer(strsplit(Sys.getenv("CHROMS", "1,2,3,4,5,6,7,8,9,10"), ",")[[1]])  # subset only for quick checks
out_files <- file.path(tdir, sprintf("%s_chr%02d.multi_snp_wise.genes.out", TRAIT, CHROMS))
raw_files <- sub("\\.genes\\.out$", ".genes.raw", out_files)
miss <- c(out_files, raw_files)[!file.exists(c(out_files, raw_files))]
if (length(miss)) stop("missing MAGMA files:\n  ", paste(miss, collapse = "\n  "), call. = FALSE)

g <- do.call(rbind, lapply(out_files, function(f)
  utils::read.table(f, header = TRUE, stringsAsFactors = FALSE, comment.char = "#")))
if (!P_COL %in% names(g))
  stop("P_COL '", P_COL, "' not in the MAGMA output; columns: ",
       paste(names(g), collapse = ", "), call. = FALSE)
g <- g[!duplicated(g$GENE), ]
cat(sprintf("[%s] %d genes on %d chromosomes\n", TRAIT, nrow(g), length(unique(g$CHR))))
pc <- intersect(c("P_MULTI", "P_SNPWISE_MEAN", "P_SNPWISE_TOP1", "P"), names(g))
cat("  median gene p by column (about 0.5 expected where most genes are null):\n")
for (cc in pc) cat(sprintf("    %-15s %.3f%s\n", cc, stats::median(g[[cc]], na.rm = TRUE),
                           if (cc == P_COL) "   <- used" else ""))

## ------------------------------------------------- gene-gene correlations --
cor_file <- file.path(OUT_DIR, sprintf("%s_gene_cor_pairs%s.txt", TRAIT,
                    if (length(CHROMS) < 10) paste0("_chr", paste(CHROMS, collapse = "-")) else ""))
if (!file.exists(cor_file)) {
  cat("building gene correlation pairs from the .genes.raw files ...\n")
  final_file <- cor_file
  cor_file <- paste0(final_file, ".partial")   # renamed only once complete
  if (file.exists(cor_file)) unlink(cor_file)
  first <- TRUE
  for (f in raw_files) {
    tmp <- tempfile(fileext = ".txt")
    magma_genesraw_to_cor_pairs_banded(genes_raw_file = f, out_pairs_file = tmp,
                                       gene_regex = "^SORBI", keep_abs_r_ge = 0,
                                       overwrite = TRUE, verbose = FALSE)
    x <- readLines(tmp, warn = FALSE); x <- x[nzchar(x)]
    if (!first) x <- x[-1]
    cat(x, file = cor_file, sep = "\n", append = !first)
    first <- FALSE
    unlink(tmp)
  }
  file.rename(cor_file, final_file)
  cor_file <- final_file
}
cat(sprintf("  correlation pairs: %s\n", cor_file))

## ---------------------------------------------------------------- pathways --
pw_raw <- utils::read.delim(PATHWAY_FILE, stringsAsFactors = FALSE, check.names = FALSE)
pw <- unique(data.frame(pathway_id = pw_raw[["Pathway-id"]],
                        gene_id    = pw_raw[["Gene-name"]],
                        stringsAsFactors = FALSE))
pw <- pw[nzchar(pw$gene_id) & pw$gene_id != "unknown", ]
pw_name <- tapply(pw_raw[["Pathway-name"]], pw_raw[["Pathway-id"]], function(x) x[1])
cat(sprintf("  SorghumBicolorCyc: %d pathways before size filtering\n",
            length(unique(pw$pathway_id))))

## ------------------------------------------------------- Westfall-Young ----
res <- catfish_wy_pathways(
  gene_results = g, pathways = pw,
  gene_col = "GENE", p_col = P_COL, chr_col = "CHR",
  cor_file = cor_file, B = B, use_tfisher = USE_TFISHER, tau_grid = TAU_GRID,
  min_genes = MIN_GENES, max_genes = MAX_GENES, seed = SEED, verbose = TRUE)
res$pathway_name <- unname(pw_name[res$pathway_id])
res <- res[, c("pathway_id", "pathway_name", setdiff(names(res), c("pathway_id", "pathway_name")))]

tag <- sprintf("%s_%s_%dcomp", TRAIT, P_COL, if (isTRUE(USE_TFISHER)) 5L else 4L)
out <- file.path(OUT_DIR, paste0(tag, "_westfall_young.csv"))
utils::write.csv(res, out, row.names = FALSE)

cat(sprintf("\npathways tested                  : %d\n", nrow(res)))
cat(sprintf("WY 5%% FWER threshold (raw scale): %.4g\n", attr(res, "wy_threshold_05")))
cat(sprintf("significant, FWER 5%% (standardised): %d\n", sum(res$p_fwer_std <= 0.05, na.rm = TRUE)))
cat(sprintf("significant, per-pathway calibrated p < 0.05: %d\n", sum(res$p_calibrated < 0.05, na.rm = TRUE)))
cat("\ntop 15 pathways:\n")
print(utils::head(data.frame(pathway = substr(res$pathway_name, 1, 55), n = res$n_genes,
                             z_std = round(res$z_std, 1),
                             p_cal = signif(res$p_calibrated, 3),
                             p_fwer = signif(res$p_fwer_std, 3)), 15), row.names = FALSE)
cat(sprintf("\nwrote %s\ntotal runtime: %.1f min\n", out,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
