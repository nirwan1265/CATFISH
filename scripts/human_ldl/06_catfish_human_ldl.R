#!/usr/bin/env Rscript
# =============================================================================
# 06_catfish_human_ldl.R -- CATFISH pathway analysis for human LDL-C.
#
# Parameters mirror the sorghum stem-volume run (05_catfish_perm.R) so the two
# species are directly comparable. Three deliberate differences, all marked
# [HUMAN] below.
#
#   Rscript 06_catfish_human_ldl.R
# =============================================================================
suppressPackageStartupMessages({ library(CATFISH) })

args    <- commandArgs(trailingOnly = TRUE)
B_PERM  <- if (length(args) >= 1) as.integer(args[[1]]) else 100000L
THREADS <- if (length(args) >= 2) as.integer(args[[2]]) else 4L
cat(sprintf("B_perm = %d   threads = %d   tail_min_B = %d\n",
            B_PERM, THREADS, max(1000L, as.integer(B_PERM * 0.5))))
if (!requireNamespace("evir", quietly = TRUE))
  warning("evir is NOT installed -- GPD tail extrapolation will silently fall\n",
          "  back to empirical p-values, flooring everything at 1/(B+1).\n",
          "  install.packages(\"evir\")", call. = FALSE)

## The generic loaders (catfish_load_gmt / catfish_load_biocyc) were added to
## R/pathway_loaders_generic.R but may not be in the INSTALLED package yet.
## Source them directly if missing, so this runs without a reinstall.
if (!exists("catfish_load_gmt")) {
  .gen <- path.expand("~/Documents/Github/MAGCAT/R/pathway_loaders_generic.R")
  if (file.exists(.gen)) {
    source(.gen)
    message("sourced generic loaders from ", .gen)
  } else {
    stop("catfish_load_gmt() not found and ", .gen, " is missing.\n",
         "Run devtools::document() then devtools::install() in the repo.",
         call. = FALSE)
  }
}

MAGMA_DIR <- path.expand("~/Documents/Github/MAGCAT/results/human_LDL/magma")
DATA_DIR  <- path.expand("~/Documents/Github/MAGCAT_data")
OUT_DIR   <- path.expand("~/Documents/Github/MAGCAT/results/human_LDL/catfish")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

GENES_OUT <- file.path(MAGMA_DIR, "ldl_eur.genes.out")
GENES_RAW <- file.path(MAGMA_DIR, "ldl_eur.genes.raw")
CORPAIRS  <- file.path(OUT_DIR,   "ldl_eur.cor_pairs.txt")
GMT       <- file.path(DATA_DIR,  "c2.cp.reactome.v2026.2.Hs.entrez.sets.annot")

stopifnot(file.exists(GENES_OUT), file.exists(GENES_RAW))

## ---- 1. gene-level results -------------------------------------------------
genes <- read.table(GENES_OUT, header = TRUE, comment.char = "#",
                    stringsAsFactors = FALSE)
genes$GENE <- as.character(genes$GENE)
cat(sprintf("genes: %d   min P: %.3g\n", nrow(genes), min(genes$P, na.rm = TRUE)))

## ---- 2. gene-gene correlations from .genes.raw -----------------------------
# [HUMAN] gene_regex: the default is "^FBgn" (FlyBase) and sorghum used
# "^SORBI". Human MAGMA gene IDs are bare Entrez integers, so match digits.
if (!file.exists(CORPAIRS)) {
  magma_genesraw_to_cor_pairs_banded(
    genes_raw_file = GENES_RAW,
    out_pairs_file = CORPAIRS,
    gene_regex     = "^[0-9]",
    keep_abs_r_ge  = 0,
    overwrite      = TRUE,
    verbose        = TRUE
  )
}
cat(sprintf("correlation pairs: %d\n",
            length(readLines(CORPAIRS, warn = FALSE)) - 1L))

## ---- 3. pathways -----------------------------------------------------------
# [HUMAN] MSigDB C2:CP, size-filtered to 10-200 genes. The archetype
# simulations validated m = 5, 20, 50; unfiltered C2:CP runs to 1,497 genes,
# well outside that range.
pathways <- catfish_load_gmt(GMT, min_size = 10, max_size = 200)
cat(sprintf("pathways: %d   genes covered: %d\n",
            length(unique(pathways$pathway_id)),
            length(unique(pathways$gene_id))))

## ---- 4. omnibus ------------------------------------------------------------
omni <- catfish_omni2_pathways(
  gene_results = genes,
  pathways     = pathways,
  species      = NULL,
  gene_col     = "GENE",
  p_raw_col    = "P",
  p_adj_col    = NULL,
  z_col        = "ZSTAT",

  # [HUMAN] tau grid extended DOWNWARD. On the first run 68.6% of pathways
  # chose tau = 0.01, the smallest value in the sorghum grid -- the optimiser
  # was pressed against the boundary and wanted to truncate harder. Note the
  # MVN null re-selects tau on every draw, so widening the grid does not by
  # itself inflate; it only removes the boundary artefact.
  tau_grid     = c(0.1, 0.05, 0.01, 0.005, 0.001),
  # [HUMAN] sorghum used min_p = 1e-15, which would clamp every top LDL gene
  # (smallest gene P here is 3.1e-30). 1e-50 keeps the real tail.
  min_p        = 1e-50,
  min_genes    = 2L,
  do_fix       = TRUE,
  stouffer_alternative = "greater",

  omnibus      = "ACAT",               # ACAT-O across component tests
  # Keep MAGMA competitive OUT of the omnibus. The package default is TRUE,
  # which would fold MAGMA's gene-set p-value in as a sixth component -- that
  # makes the CATFISH-vs-MAGMA benchmark partly circular. Set explicitly so a
  # later magma_out argument cannot silently redefine the omnibus.
  magma_out             = NULL,
  include_magma_in_omni = FALSE,
  include_magma_in_perm = FALSE,
  perm_mode    = "mvn",                # LD-aware MVN null
  B_perm       = B_PERM,
  magma_cor_file = CORPAIRS,
  make_PD      = TRUE,
  mvn_marginal = "uniform",
  mvn_calibrate_components = FALSE,    # "Omnibus Alone" -- the paper default

  tail_mode          = "hybrid_gpd",   # GPD tail extrapolation
  tail_switch_exceed = 10L,
  tail_gpd_k         = 250L,
  # tail_min_B must be BELOW B_perm. The package default is 10000, and with
  # B_perm = 10000 a single NA null draw drops B to 9999, fails the
  # `B < min_B` gate, and silently skips GPD for that pathway -- which is why
  # 638 of 1348 pathways sat tied at exactly 1/(B+1).
  tail_min_B         = max(1000L, as.integer(B_PERM * 0.5)),
  tail_min_tail      = 50L,

  n_threads = THREADS,
  seed      = 20260706,
  output    = FALSE
)

saveRDS(omni, file.path(OUT_DIR, "ldl_eur_reactome_catfish_omni.rds"))
write.csv(omni, file.path(OUT_DIR, "ldl_eur_reactome_catfish_omni.csv"), row.names = FALSE)

## ---- 5. positive controls --------------------------------------------------
# omni_p_final is the MVN-calibrated result. grep("^omni_p")[1] would pick
# omni_p_analytic -- the UNCALIBRATED analytic omnibus, which saturates at the
# floor and makes every ranking meaningless.
pcol <- if ("omni_p_final" %in% names(omni)) "omni_p_final" else
        grep("^omni_p", names(omni), value = TRUE)[1]
cat(sprintf("reporting column: %s\n", pcol))
np   <- nrow(omni)
bonf <- 0.05 / np
cat(sprintf("\npathways tested: %d   Bonferroni threshold: %.3g\n", np, bonf))
cat(sprintf("significant (Bonferroni): %d\n", sum(omni[[pcol]] < bonf, na.rm = TRUE)))

cat("\n--- LDL positive controls (should rank near the top) ---\n")
ctrl <- grep("LDL_CLEARANCE|PLASMA_LIPOPROTEIN|CHOLESTEROL_BIOSYNTH|LIPOPROTEIN_REMODELING",
             omni$pathway_id, value = TRUE)
o <- omni[omni$pathway_id %in% ctrl, ]
o <- o[order(o[[pcol]]), ]
print(utils::head(data.frame(pathway = o$pathway_id,
                             n_genes = o$n_genes,
                             p = signif(o[[pcol]], 3),
                             rank = match(o$pathway_id,
                                          omni$pathway_id[order(omni[[pcol]])])), 10))

floor_p <- 1/(B_PERM + 1)
at   <- sum(abs(omni[[pcol]] - floor_p) < floor_p*1e-3, na.rm = TRUE)
below<- sum(omni[[pcol]] < floor_p*0.999, na.rm = TRUE)
cat(sprintf("\n--- Monte Carlo resolution ---\n"))
cat(sprintf("  floor 1/(B+1)          : %.3g\n", floor_p))
cat(sprintf("  tied AT the floor      : %d  (%.1f%%)  <- unrankable\n", at, 100*at/np))
cat(sprintf("  GPD-extrapolated below : %d  (%.1f%%)\n", below, 100*below/np))

ord <- order(omni[[pcol]])
lipid <- grep("LIPID|LIPOPROTEIN|CHOLESTEROL|STEROL|BILE|FATTY_ACID|LDL|SREBP|SREBF|PPARA|MEVALON",
              omni$pathway_id[ord], ignore.case = TRUE)
cat(sprintf("\n--- lipid positive controls (%d of %d) ---\n", length(lipid), np))
cat(sprintf("  median rank : %.0f / %d   (random %.0f)\n", median(lipid), np, np/2))
cat(sprintf("  in top 100  : %d   (expected %.1f)\n", sum(lipid <= 100), length(lipid)*100/np))
cat(sprintf("  in top 50   : %d   (expected %.1f)\n", sum(lipid <= 50), length(lipid)*50/np))
cat("\n--- top 12 lipid pathways by rank ---\n")
print(data.frame(rank = head(lipid, 12),
                 pathway = substr(omni$pathway_id[ord][head(lipid, 12)], 1, 56),
                 p = signif(omni[[pcol]][ord][head(lipid, 12)], 3)), row.names = FALSE)

cat("\n--- top 15 pathways overall ---\n")
t15 <- omni[order(omni[[pcol]]), ][1:15, ]
print(data.frame(pathway = t15$pathway_id, n_genes = t15$n_genes,
                 p = signif(t15[[pcol]], 3)))
