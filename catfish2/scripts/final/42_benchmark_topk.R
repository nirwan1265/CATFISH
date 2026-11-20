#!/usr/bin/env Rscript
## =============================================================================
## 42_benchmark_topk.R -- genome-wide gene prioritisation benchmark (the use for
## choosing genes to validate, e.g. by CRISPR): how many known genes does each
## method place in its top k genes overall, European LDL-C and HDL-C.
##
## Gene universe: genes scored by both CATFISH (06, mapped to Ensembl by unique
## symbol in the PoPS annotation) and PoPS; every method is ranked within it.
## Known-gene sets (Ensembl IDs):
##   trait      Open Targets gold-standard genes of the trait (LDL: EFO_0004611
##              and hypercholesterolaemia HP_0003124; HDL: EFO_0004612), 2019 and
##              2023 sets -- independent of this paper
##   trait+ctrl the same plus the paper's positive controls of the trait
##   any_lipid  Open Targets genes of any lipid trait plus all controls (not
##              trait-specific)
## Methods: GWAS (06 rank_gwas), MAGMA (rank_magma), CATFISH (rank), PoPS (PoPS_Score).
## Reads  catfish2/data/gold_standards/{gwas_gold_standards.191108.tsv, otg_gs_230511.json}
##        catfish2/data/pops/gene_annot_jun10.txt, catfish2/results/06_candidate_genes_{ldl,hdl}.csv,
##        catfish2/results/pops/{ldl,hdl}_eur.preds
## Writes catfish2/results/42_benchmark_topk.csv, 42_top20_lists.csv
## Usage (from catfish2/scripts/final):  Rscript 42_benchmark_topk.R   (needs jsonlite)
## =============================================================================
suppressPackageStartupMessages({ library(data.table); library(jsonlite) })
BASE <- Sys.getenv("CATFISH_BASE", path.expand("~/Documents/Github/CATFISH"))
RES  <- Sys.getenv("RES_DIR", file.path(BASE, "catfish2/results"))
GSD  <- Sys.getenv("GS_DIR", file.path(BASE, "catfish2/data/gold_standards"))
ANN  <- Sys.getenv("GENE_ANNOT", file.path(BASE, "catfish2/data/pops/gene_annot_jun10.txt"))
K    <- c(10, 20, 50, 100, 200)

ann <- fread(ANN)
sym <- ann[, .N, by = NAME][N == 1, NAME]
s2e <- setNames(ann[NAME %in% sym, ENSGID], ann[NAME %in% sym, NAME]); e2s <- setNames(ann$NAME, ann$ENSGID)
gs <- fread(file.path(GSD, "gwas_gold_standards.191108.tsv"))
J  <- lapply(readLines(file.path(GSD, "otg_gs_230511.json")), fromJSON)
ot <- function(efos) unique(c(
  gs[vapply(strsplit(trait_info.ontology, "|", fixed = TRUE), function(o) any(o %in% efos), logical(1)),
     gold_standard_info.gene_id],
  unlist(lapply(J, function(j) if (any(j$trait_info$ontology %in% efos)) j$gold_standard_info$gene_id))))
CTRL <- list(  # as 06_candidate_genes.R
  ldl = c("LDLR","PCSK9","APOB","APOE","SORT1","HMGCR","NPC1L1","ABCG8","LPA","CETP","LDLRAP1","MYLIP","ANGPTL3"),
  hdl = c("CETP","LIPC","LIPG","ABCA1","APOA1","APOA5","LCAT","SCARB1","PLTP","LPL","ANGPTL4","CUBN"))
TRAIT_EFO <- list(ldl = c("EFO_0004611", "HP_0003124"), hdl = "EFO_0004612")
LIPID_EFO <- c("EFO_0004611","EFO_0004612","HP_0003124","HMDB0000067","EFO_0004574","EFO_0004530","EFO_0004529")
any_lipid <- unique(c(ot(LIPID_EFO), s2e[unlist(CTRL)]))

out <- list(); tops <- list()
for (tr in c("ldl", "hdl")) {
  g <- fread(file.path(RES, sprintf("06_candidate_genes_%s.csv", tr)))
  g[, ENSGID := s2e[symbol]]; g <- g[!is.na(ENSGID) & !duplicated(ENSGID)]
  p <- fread(file.path(RES, "pops", sprintf("%s_eur.preds", tr)))
  u <- intersect(g$ENSGID, p$ENSGID); g <- g[ENSGID %in% u]; p <- p[ENSGID %in% u]
  sets <- list(trait = ot(TRAIT_EFO[[tr]]), `trait+ctrl` = unique(c(ot(TRAIT_EFO[[tr]]), s2e[CTRL[[tr]]])),
               any_lipid = any_lipid)
  sc <- list(GWAS = setNames(-g$rank_gwas, g$ENSGID), MAGMA = setNames(-g$rank_magma, g$ENSGID),
             CATFISH = setNames(-g$rank, g$ENSGID), PoPS = setNames(p$PoPS_Score, p$ENSGID))
  for (m in names(sc)) {
    s <- sc[[m]]; s[is.na(s)] <- -Inf
    r <- setNames(frank(-s, ties.method = "min"), names(s))
    top <- names(sort(s, decreasing = TRUE))[1:20]
    tops[[length(tops) + 1L]] <- data.table(trait = toupper(tr), method = m, rank = 1:20, gene = e2s[top],
                                            trait_gold = top %in% sets$trait, any_lipid = top %in% any_lipid)
    for (sn in names(sets)) {
      G <- intersect(sets[[sn]], u); rg <- r[G]
      out[[length(out) + 1L]] <- data.table(trait = toupper(tr), gold = sn, n_gold = length(G), method = m,
        t(setNames(vapply(K, function(k) sum(rg <= k), numeric(1)), paste0("top", K))),
        median_rank = as.numeric(median(rg)))
    }
  }
}
D <- rbindlist(out); T20 <- rbindlist(tops)
fwrite(D, file.path(RES, "42_benchmark_topk.csv")); fwrite(T20, file.path(RES, "42_top20_lists.csv"))
options(width = 180); print(D, row.names = FALSE)
cat("\nwrote", file.path(RES, "42_benchmark_topk.csv"), "and 42_top20_lists.csv\n")
