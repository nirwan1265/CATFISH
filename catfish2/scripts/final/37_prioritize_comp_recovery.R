#!/usr/bin/env Rscript
## =============================================================================
## 37_prioritize_comp_recovery.R -- does the pathway layer of the gene score
## recover the pre-specified established genes better when it is taken from the
## self-contained test (the paper), the competitive test, or the joint call?
##
## Established genes: the positive controls named before prioritisation, as in
## 06/21 (lipids: LDL-C 13, HDL-C 12 genes, the same list in every ancestry) and
## 25 (sorghum: dry biomass 8, lignin 3). Only genes that are in at least one
## pathway can be moved by the pathway layer; the others move only because
## the genes around them do.
##
## Reads  catfish2/results/06_candidate_genes_<tag>.csv           (self)
##        catfish2/results/sorghum_bap/25_<TRAIT>_gene_prioritize.csv (self)
##        catfish2/results/prioritize_comp/<comp|joint>/...         (36)
## Writes catfish2/results/prioritize_comp/37_recovery_genes.csv     one row per
##          set x established gene x pathway layer: combined and pathway rank
##        catfish2/results/prioritize_comp/37_recovery_summary.csv   one row per
##          set x pathway layer: median rank, top 10/20/50, genes moved up/down
## Usage (from catfish2/scripts/final):  Rscript 37_prioritize_comp_recovery.R
## =============================================================================
suppressPackageStartupMessages(library(data.table))

BASE <- Sys.getenv("CATFISH_BASE", path.expand("~/Documents/Github/CATFISH"))
RES  <- file.path(BASE, "catfish2/results")
OUT  <- file.path(RES, "prioritize_comp")
VARIANTS <- c(self = "self-contained", comp = "competitive", joint = "joint")

## positive controls, identical to 06_candidate_genes.R / 21_lipid_results.R
CONTROLS <- list(
  ldl = c("LDLR","PCSK9","APOB","APOE","SORT1","HMGCR","NPC1L1","ABCG8",
          "LPA","CETP","LDLRAP1","MYLIP","ANGPTL3"),
  hdl = c("CETP","LIPC","LIPG","ABCA1","APOA1","APOA5","LCAT","SCARB1",
          "PLTP","LPL","ANGPTL4","CUBN"))
## identical to 25_sorghum_gene_prioritize.R
SORGHUM_EST <- list(
  Dry_tons_per_acre = c("SORBI_3006G057866", "SORBI_3001G394400", "SORBI_3001G087100",
                        "SORBI_3006G004400", "SORBI_3002G302700", "SORBI_3009G229800",
                        "SORBI_3006G067700", "SORBI_3007G163800"),
  Lignin_DM_perc    = c("SORBI_3004G062500", "SORBI_3004G071000", "SORBI_3007G047300"))

LIPID <- c(outer(c("ldl", "hdl"), c("eur", "afr", "eas", "sas", "his"), paste, sep = "_"))

read_set <- function(s, v) {
  if (s %in% names(SORGHUM_EST)) {
    f <- if (v == "self") file.path(RES, "sorghum_bap", sprintf("25_%s_gene_prioritize.csv", s)) else
      file.path(OUT, v, "sorghum", sprintf("25_%s_gene_prioritize.csv", s))
    if (!file.exists(f)) return(NULL)
    g <- fread(f)
    data.table(id = g$gene, rank = g$rank, rank_pathway = g$pathway_rank,
               best_pathway = g$best_pathway_name, n_genes = nrow(g))[
      id %in% SORGHUM_EST[[s]]][, gene := id]
  } else {
    tr <- sub("_.*", "", s); anc <- sub(".*_", "", s)
    tag <- if (anc == "eur") tr else s
    f <- if (v == "self") file.path(RES, sprintf("06_candidate_genes_%s.csv", tag)) else
      file.path(OUT, v, sprintf("06_candidate_genes_%s.csv", tag))
    if (!file.exists(f)) return(NULL)
    g <- fread(f)
    data.table(gene = g$symbol, rank = g$rank, rank_pathway = g$rank_pathway,
               best_pathway = g$best_pathway, n_genes = nrow(g))[gene %in% CONTROLS[[tr]]]
  }
}

rows <- list()
for (s in c(LIPID, names(SORGHUM_EST))) for (v in names(VARIANTS)) {
  d <- read_set(s, v)
  if (is.null(d)) { if (v != "self") cat(sprintf("missing: %s %s (run 36 first)\n", s, v)); next }
  rows[[length(rows) + 1L]] <- d[, `:=`(set = s, layer = v)]
}
G <- rbindlist(rows, fill = TRUE)
G[, layer := factor(layer, levels = names(VARIANTS))]
G[, in_pathway := !is.na(rank_pathway)]
## change against the self-contained layer (negative = moved up)
G[, d_rank := rank - rank[layer == "self"][1], by = .(set, gene)]
setcolorder(G, c("set", "gene", "layer", "rank", "d_rank", "rank_pathway", "in_pathway",
                 "best_pathway", "n_genes"))
setorder(G, set, gene, layer)
fwrite(G, file.path(OUT, "37_recovery_genes.csv"))

S <- G[, .(n_est = .N, n_in_pathway = sum(in_pathway),
           median_rank = as.numeric(median(rank)), top10 = sum(rank <= 10), top20 = sum(rank <= 20),
           top50 = sum(rank <= 50),
           median_rank_pathway = as.numeric(median(rank_pathway, na.rm = TRUE)),
           moved_up = sum(d_rank < 0), moved_down = sum(d_rank > 0)),
       by = .(set, layer)]
S[layer == "self", c("moved_up", "moved_down") := NA_integer_]
setorder(S, set, layer)
fwrite(S, file.path(OUT, "37_recovery_summary.csv"))

## pooled over the ten lipid sets (genes repeat across ancestries, so this is a
## description, not a test)
P <- G[set %in% LIPID, .(median_rank = as.numeric(median(rank)), top10 = sum(rank <= 10),
                         top50 = sum(rank <= 50), moved_up = sum(d_rank < 0),
                         moved_down = sum(d_rank > 0), n = .N), by = layer]

options(width = 160)
cat("\nEstablished-gene recovery by pathway layer (rank among all genes, 1 = top)\n\n")
print(S, row.names = FALSE)
cat("\nAll ten lipid sets together:\n")
print(P, row.names = FALSE)
cat(sprintf("\nwrote %s\n      %s\n", file.path(OUT, "37_recovery_summary.csv"),
            file.path(OUT, "37_recovery_genes.csv")))
