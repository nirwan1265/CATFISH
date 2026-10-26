#!/usr/bin/env Rscript
# =============================================================================
# 07_magma_competitive_compare.R
#
# The competitive arm. CATFISH's five component tests are all SELF-CONTAINED
# ("is there any signal in this set?"). At n ~ 1.3M that null is false almost
# everywhere -- 44% of genes reach p<0.05 and 11.2% pass genome-wide
# Bonferroni -- so every pathway fires and ranking collapses.
#
# MAGMA's gene-set analysis is COMPETITIVE ("are these genes MORE associated
# than other genes?"), which is the question that still discriminates at this
# power. This script runs it on identical inputs and compares.
#
#   Rscript 07_magma_competitive_compare.R
# =============================================================================
suppressPackageStartupMessages({ library(CATFISH) })

MAGMA_DIR <- path.expand("~/Documents/Github/MAGCAT/results/human_LDL/magma")
DATA_DIR  <- path.expand("~/Documents/Github/MAGCAT_data")
OUT_DIR   <- path.expand("~/Documents/Github/MAGCAT/results/human_LDL/catfish")

GENES_RAW <- file.path(MAGMA_DIR, "ldl_eur.genes.raw")
GENES_OUT <- file.path(MAGMA_DIR, "ldl_eur.genes.out")
SET_ANNOT <- file.path(DATA_DIR,  "c2.cp.reactome.v2026.2.Hs.entrez.sets.annot")
CATFISH_CSV <- file.path(OUT_DIR, "ldl_eur_reactome_catfish_omni.csv")

stopifnot(file.exists(GENES_RAW), file.exists(SET_ANNOT), file.exists(CATFISH_CSV))

genes <- read.table(GENES_OUT, header = TRUE, comment.char = "#",
                    stringsAsFactors = FALSE)
genes$GENE <- as.character(genes$GENE)

## ---- 1. MAGMA competitive gene-set analysis --------------------------------
mg <- magma_geneset_competitive(
  gene_results_raw = GENES_RAW,
  set_annot        = SET_ANNOT,
  out_prefix       = "ldl_eur_reactome_competitive",
  out_dir          = OUT_DIR,
  genes_all        = genes$GENE,
  write_tidy       = TRUE
)
saveRDS(mg, file.path(OUT_DIR, "ldl_eur_reactome_magma_competitive.rds"))

## ---- 2. side-by-side against the self-contained omnibus --------------------
cf <- read.csv(CATFISH_CSV, stringsAsFactors = FALSE)

if (!all(c("pathway_id","magma_pvalue") %in% names(mg))) {
  nm <- names(mg)
  message("columns returned: ", paste(nm, collapse = ", "))
  stop("expected pathway_id / magma_pvalue in the competitive output", call. = FALSE)
}

m <- merge(cf[, c("pathway_id","n_genes","omni_p_final")],
           mg[, c("pathway_id","magma_pvalue")],
           by = "pathway_id")
m$rank_catfish <- rank(m$omni_p_final, ties.method = "min")
m$rank_magma   <- rank(m$magma_pvalue, ties.method = "min")
n <- nrow(m)

cat(sprintf("\npathways in both: %d\n", n))
cat(sprintf("CATFISH  p<0.05: %d (%.1f%%)   Bonferroni: %d\n",
            sum(m$omni_p_final < 0.05), 100*mean(m$omni_p_final < 0.05),
            sum(m$omni_p_final < 0.05/n)))
cat(sprintf("MAGMA    p<0.05: %d (%.1f%%)   Bonferroni: %d\n",
            sum(m$magma_pvalue < 0.05, na.rm=TRUE), 100*mean(m$magma_pvalue < 0.05, na.rm=TRUE),
            sum(m$magma_pvalue < 0.05/n, na.rm=TRUE)))
cat(sprintf("Spearman rank correlation: %.3f\n",
            suppressWarnings(cor(m$rank_catfish, m$rank_magma, method = "spearman",
                                 use = "complete.obs"))))

## ---- 3. the benchmark: do the LDL positive controls rank at the top? -------
ctrl <- grep("LDL_CLEARANCE|PLASMA_LIPOPROTEIN|CHOLESTEROL_BIOSYNTH|SREBP|LIPOPROTEIN_REMODEL",
             m$pathway_id, value = TRUE)
pc <- m[m$pathway_id %in% ctrl, ]
pc <- pc[order(pc$rank_magma), ]
cat("\n--- LDL positive controls: rank under each method (of ", n, ") ---\n", sep = "")
print(data.frame(pathway   = substr(pc$pathway_id, 1, 52),
                 n_genes   = pc$n_genes,
                 magma_rank   = pc$rank_magma,
                 catfish_rank = pc$rank_catfish,
                 magma_p      = signif(pc$magma_pvalue, 3)),
      row.names = FALSE)

cat("\n--- median rank of the positive-control set ---\n")
cat(sprintf("  MAGMA competitive : %.0f / %d\n", median(pc$rank_magma), n))
cat(sprintf("  CATFISH omnibus   : %.0f / %d\n", median(pc$rank_catfish), n))
cat(sprintf("  (random expectation: %.0f)\n", n/2))

cat("\n--- MAGMA competitive top 15 ---\n")
t15 <- m[order(m$magma_pvalue), ][1:15, ]
print(data.frame(pathway = substr(t15$pathway_id, 1, 58),
                 n_genes = t15$n_genes,
                 p = signif(t15$magma_pvalue, 3)), row.names = FALSE)

write.csv(m, file.path(OUT_DIR, "ldl_eur_reactome_catfish_vs_magma.csv"), row.names = FALSE)
cat("\nwrote ldl_eur_reactome_catfish_vs_magma.csv\n")
