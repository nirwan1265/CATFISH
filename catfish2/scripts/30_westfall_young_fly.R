#!/usr/bin/env Rscript
## =============================================================================
## 30_westfall_young_fly.R -- Westfall-Young pathway analysis for the DGRP
## body-weight and triglyceride traits, using the CATFISH package.
##
## Same analysis and settings as 23_westfall_young_sorghum.R and
## scripts/human_ldl/18_westfall_young_trait.R (catfish_wy_pathways(),
## B = 10,000, five components, tau grid {0.01, 0.05, 0.5}, seed 20260706), on:
##   gene results   MAGMA SNP-wise mean model (column P), merged across the
##                  five chromosome arms: <MAGMA_DIR>/<TRAIT>.genes.out
##                  (GWAS: vcf2gwas/GEMMA LMM, 3 PCs; see dgrp_bw_tg/README.md)
##   correlations   the merged <TRAIT>.genes.raw (carries every arm)
##   pathways       Reactome, Drosophila melanogaster (orthology-inferred from
##                  the human curation), FlyBase FBgn IDs, all hierarchy levels
##   size window    10-200 analysed genes, as for the human Reactome run
##
## Usage (from catfish2/scripts):
##   TRAIT=BW_control_female Rscript 30_westfall_young_fly.R [B] [use_tfisher]
## Env overrides: MAGMA_DIR, PATHWAY_FILE, MIN_GENES (10), MAX_GENES (200)
## Writes ../results/dgrp/<TRAIT>_reactome_<k>comp_westfall_young.csv
##        <MAGMA_DIR>/<TRAIT>.cor_pairs.txt (built once, reused; gitignored)
## =============================================================================
suppressPackageStartupMessages(library(CATFISH))
args <- commandArgs(trailingOnly = TRUE)
B           <- if (length(args) >= 1) as.integer(args[[1]]) else 10000L
USE_TFISHER <- if (length(args) >= 2) as.logical(args[[2]]) else TRUE
SEED        <- 20260706                       # same seed as 18 and 23
TAU_GRID    <- c(0.01, 0.05, 0.5)             # same grid as 18 and 23
P_COL       <- "P"                            # SNP-wise mean: the only model run

TRAIT <- Sys.getenv("TRAIT", "BW_control_female")
BASE  <- Sys.getenv("CATFISH_BASE", path.expand("~/Documents/Github/CATFISH"))
MAGMA_DIR    <- Sys.getenv("MAGMA_DIR", file.path(BASE, "dgrp_bw_tg/magma"))
PATHWAY_FILE <- Sys.getenv("PATHWAY_FILE",
                  file.path(BASE, "inst/extdata/pathway/dme_reactome_all_levels.tsv"))
MIN_GENES <- as.integer(Sys.getenv("MIN_GENES", "10"))
MAX_GENES <- as.integer(Sys.getenv("MAX_GENES", "200"))
OUT_DIR   <- "../results/dgrp"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time()

## ------------------------------------------------------------ gene results --
go  <- file.path(MAGMA_DIR, paste0(TRAIT, ".genes.out"))
raw <- file.path(MAGMA_DIR, paste0(TRAIT, ".genes.raw"))
miss <- c(go, raw)[!file.exists(c(go, raw))]
if (length(miss)) stop("missing MAGMA files (run the merge, 04c_magma_merge_LSF.sh):\n  ",
                       paste(miss, collapse = "\n  "), call. = FALSE)
g <- utils::read.table(go, header = TRUE, stringsAsFactors = FALSE, comment.char = "#")
g <- g[!duplicated(g$GENE), ]
cat(sprintf("[%s] %d genes on %d chromosome arms; median gene p %.3f (about 0.5 expected)\n",
            TRAIT, nrow(g), length(unique(g$CHR)), stats::median(g[[P_COL]], na.rm = TRUE)))

## ------------------------------------------------- gene-gene correlations --
cor_file <- file.path(MAGMA_DIR, paste0(TRAIT, ".cor_pairs.txt"))
if (!file.exists(cor_file)) {
  cat("building gene correlation pairs from the merged .genes.raw ...\n")
  tmp <- paste0(cor_file, ".partial")          # renamed only once complete
  magma_genesraw_to_cor_pairs_banded(genes_raw_file = raw, out_pairs_file = tmp,
                                     gene_regex = "^FBgn", keep_abs_r_ge = 0,
                                     overwrite = TRUE, verbose = TRUE)
  file.rename(tmp, cor_file)
}
cat(sprintf("  correlation pairs: %s\n", cor_file))

## ---------------------------------------------------------------- pathways --
## Ensembl2Reactome_All_Levels, fly rows: gene, pathway, url, name, evidence, species
pw_raw <- utils::read.delim(PATHWAY_FILE, header = FALSE, quote = "", stringsAsFactors = FALSE,
                            col.names = c("gene_id", "pathway_id", "url", "name", "evidence", "species"))
pw <- unique(pw_raw[, c("pathway_id", "gene_id")])
pw_name <- tapply(pw_raw$name, pw_raw$pathway_id, function(x) x[1])
cat(sprintf("  Reactome (fly): %d pathways before size filtering\n", length(unique(pw$pathway_id))))

## ------------------------------------------------------- Westfall-Young ----
res <- catfish_wy_pathways(
  gene_results = g, pathways = pw,
  gene_col = "GENE", p_col = P_COL, chr_col = "CHR",
  cor_file = cor_file, B = B, use_tfisher = USE_TFISHER, tau_grid = TAU_GRID,
  min_genes = MIN_GENES, max_genes = MAX_GENES, seed = SEED, verbose = TRUE)
res$pathway_name <- unname(pw_name[res$pathway_id])
res <- res[, c("pathway_id", "pathway_name", setdiff(names(res), c("pathway_id", "pathway_name")))]
tag <- sprintf("%s_reactome_%dcomp", TRAIT, if (isTRUE(USE_TFISHER)) 5L else 4L)
out <- file.path(OUT_DIR, paste0(tag, "_westfall_young.csv"))
utils::write.csv(res, out, row.names = FALSE)

cat(sprintf("\npathways tested                       : %d\n", nrow(res)))
cat(sprintf("WY 5%% FWER threshold (raw scale)     : %.4g\n", attr(res, "wy_threshold_05")))
cat(sprintf("significant, FWER 5%% (standardised)  : %d\n", sum(res$p_fwer_std <= 0.05, na.rm = TRUE)))
cat(sprintf("per-pathway calibrated p < 0.05      : %d\n", sum(res$p_calibrated < 0.05, na.rm = TRUE)))
cat("\ntop 15 pathways:\n")
print(utils::head(data.frame(pathway = substr(res$pathway_name, 1, 55), n = res$n_genes,
                             z_std = round(res$z_std, 1),
                             p_cal = signif(res$p_calibrated, 3),
                             p_fwer = signif(res$p_fwer_std, 3)), 15), row.names = FALSE)
cat(sprintf("\nwrote %s\ntotal runtime: %.1f min\n", out,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
