#!/usr/bin/env Rscript
## =============================================================================
## 05_catfish.R <rep_id> -- CATFISH on one replicate: catfish_wy_pathways() with
## the paper's settings (five components, tau {0.01, 0.05, 0.5}, seed 20260706,
## genomic control on, the package default) and the competitive test
## (circular genomic permutation, genes in CHR/START order). B = WY_B for both
## nulls. Pathways: SIM_SETS (Reactome on chr 1-10, 10-200 genes).
## Writes RUN/catfish/rep_<r>.csv.
## =============================================================================
suppressPackageStartupMessages(library(CATFISH))
env <- function(k) { v <- Sys.getenv(k); if (!nzchar(v)) stop(k, " is not set"); v }
DEV <- Sys.getenv("DEV_WY", "")
if (nzchar(DEV)) {
  dev_env <- new.env(parent = asNamespace("CATFISH")); sys.source(DEV, envir = dev_env)
  catfish_wy_pathways <- dev_env$catfish_wy_pathways
  cat("using catfish_wy_pathways() from", DEV, "\n")
}
if (!"competitive" %in% names(formals(catfish_wy_pathways)))
  stop("this catfish_wy_pathways() has no competitive test: set DEV_WY or reinstall CATFISH")
r <- as.integer(commandArgs(trailingOnly = TRUE)[1]); RP <- sprintf("rep_%04d", r)
RUN <- env("RUN"); IN <- env("IN"); B <- as.integer(env("WY_B"))
gf <- file.path(RUN, "magma", RP, paste0(RP, ".genes.out"))
cf <- file.path(IN, "sim_cor_pairs.txt")
for (f in c(gf, cf)) if (!file.exists(f)) stop("missing ", f)
g <- utils::read.table(gf, header = TRUE, stringsAsFactors = FALSE, comment.char = "#")
g$GENE <- as.character(g$GENE)
pw <- list()
for (x in readLines(env("SIM_SETS"))) { f <- strsplit(x, "\t", fixed = TRUE)[[1]]; pw[[f[1]]] <- f[-1] }
t0 <- Sys.time()
res <- catfish_wy_pathways(g, pathways = pw, gene_col = "GENE", p_col = "P", chr_col = "CHR",
                           z_col = "ZSTAT", cor_file = cf, B = B,
                           use_tfisher = TRUE, tau_grid = c(0.01, 0.05, 0.5),
                           min_genes = as.integer(env("MIN_GENES")),
                           max_genes = as.integer(env("MAX_GENES")),
                           seed = as.integer(env("WY_SEED")),
                           competitive = TRUE, pos_col = "START", B_comp = B, verbose = FALSE)
dir.create(file.path(RUN, "catfish"), showWarnings = FALSE)
utils::write.csv(res, file.path(RUN, "catfish", paste0(RP, ".csv")), row.names = FALSE)
cat(sprintf("[catfish] %s: %d pathways, self %d, competitive %d, both %d (FWER <= 0.05); %.1f min\n",
            RP, nrow(res), sum(res$p_fwer_std <= 0.05), sum(res$p_fwer_comp <= 0.05),
            sum(pmax(res$p_fwer_std, res$p_fwer_comp) <= 0.05),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
