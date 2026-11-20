#!/usr/bin/env Rscript
## =============================================================================
## 41_benchmark_genes.R -- locus-level gene-prioritisation benchmark on the
## Open Targets LDL-C / HDL-C gold standard (40_gold_standard_lipids.py).
##
## For every gold-standard locus (trait, causal gene, lead variant on GRCh37):
##   candidates  the protein-coding genes of the PoPS gene annotation
##               (gene_annot_jun10.txt, GRCh37) whose TSS lies within WINDOW of
##               the lead variant; loci whose causal gene is not a candidate,
##               or with a single candidate, are dropped
##   each method ranks the candidates; the locus is solved when the top-ranked
##               gene is a gold-standard gene of that trait (ties at the top
##               count 1/k); also the reciprocal rank of the gold gene
## Methods
##   nearest    distance from the lead variant to the TSS
##   gwas       best SNP p-value in the gene (06 rank_gwas)
##   magma      MAGMA gene p-value (06 rank_magma; Entrez, mapped by symbol)
##   pathway    CATFISH pathway layer alone (06 rank_pathway)
##   catfish    CATFISH score (06 rank)
##   pops       PoPS score (results/pops/<trait>_<anc>.preds), where present
## A method with no value for a candidate ranks it last.
##
## Reads  catfish2/data/gold_standards/lipid_gold_standard.tsv
##        catfish2/data/pops/gene_annot_jun10.txt
##        catfish2/results/06_candidate_genes_<tag>.csv
##        catfish2/results/pops/<trait>_<anc>.preds          (optional)
## Writes catfish2/results/41_benchmark_loci.csv, 41_benchmark_summary.csv
## Usage (from catfish2/scripts/final):  Rscript 41_benchmark_genes.R
## Env: CATFISH_BASE, WINDOW_KB (default 500), MIN_CONF (High|Medium|Low; default Low = all)
## =============================================================================
suppressPackageStartupMessages(library(data.table))
BASE  <- Sys.getenv("CATFISH_BASE", path.expand("~/Documents/Github/CATFISH"))
RES   <- Sys.getenv("RES_DIR", file.path(BASE, "catfish2/results"))
GOLD  <- Sys.getenv("GOLD", file.path(BASE, "catfish2/data/gold_standards/lipid_gold_standard.tsv"))
ANNOT <- Sys.getenv("GENE_ANNOT", file.path(BASE, "catfish2/data/pops/gene_annot_jun10.txt"))
POPS  <- Sys.getenv("POPS_DIR", file.path(RES, "pops"))
WIN   <- as.numeric(Sys.getenv("WINDOW_KB", "500")) * 1000
CONF  <- c(Low = 1, Medium = 2, High = 3)
MIN_CONF <- CONF[[Sys.getenv("MIN_CONF", "Low")]]
SETS  <- c(outer(c("ldl", "hdl"), c("eur", "afr", "eas", "sas", "his"), paste, sep = "_"))

ann  <- fread(ANNOT)
gold <- fread(GOLD, colClasses = list(character = c("chr", "gene")))[CONF[confidence] >= MIN_CONF]
ann[, CHR := as.character(CHR)]
## symbol -> Ensembl, only symbols that map to one gene
sym  <- ann[, .N, by = NAME][N == 1, NAME]
s2e  <- setNames(ann[NAME %in% sym, ENSGID], ann[NAME %in% sym, NAME])

## top-1 credit with ties, and reciprocal rank (score: larger = better)
judge <- function(score, is_gold) {
  score[is.na(score)] <- -Inf
  top <- score == max(score)
  r <- frank(-score, ties.method = "average")
  c(top1 = sum(top & is_gold) / sum(top), rr = 1 / min(r[is_gold]))
}

loci <- list()
for (A in SETS) {
  tr <- sub("_.*", "", A); anc <- sub(".*_", "", A); tag <- if (anc == "eur") tr else A
  f06 <- file.path(RES, sprintf("06_candidate_genes_%s.csv", tag)); if (!file.exists(f06)) next
  g <- fread(f06, colClasses = list(character = "gene"))
  g[, ENSGID := s2e[symbol]]; g <- g[!is.na(ENSGID) & !duplicated(ENSGID)]
  fp <- file.path(POPS, sprintf("%s_%s.preds", tr, anc))
  pops <- if (file.exists(fp)) fread(fp) else NULL
  gt <- gold[trait == tr]
  for (i in seq_len(nrow(gt))) {
    L <- gt[i]
    cand <- ann[CHR == L$chr & abs(TSS - L$pos37) <= WIN]
    if (!(L$gene %in% cand$ENSGID) || nrow(cand) < 2) next
    is_gold <- cand$ENSGID %in% gt$gene
    m <- match(cand$ENSGID, g$ENSGID)
    sc <- list(nearest = -abs(cand$TSS - L$pos37),
               gwas    = -g$rank_gwas[m],
               magma   = -g$rank_magma[m],
               pathway = -g$rank_pathway[m],
               catfish = -g$rank[m])
    if (!is.null(pops)) sc$pops <- pops$PoPS_Score[match(cand$ENSGID, pops$ENSGID)]
    for (meth in names(sc)) {
      j <- judge(sc[[meth]], is_gold)
      loci[[length(loci) + 1L]] <- data.table(set = A, trait = tr, gene = L$gene,
        symbol = cand$NAME[cand$ENSGID == L$gene][1], chr = L$chr, pos37 = L$pos37,
        confidence = L$confidence, n_cand = nrow(cand), method = meth,
        top1 = j[["top1"]], rr = j[["rr"]])
    }
  }
}
D <- rbindlist(loci)
fwrite(D, file.path(RES, "41_benchmark_loci.csv"))
S <- D[, .(n_loci = .N, top1 = mean(top1), mrr = mean(rr)), by = .(set, method)]
S[, method := factor(method, levels = c("nearest", "gwas", "magma", "pathway", "catfish", "pops"))]
setorder(S, set, method)
fwrite(S, file.path(RES, "41_benchmark_summary.csv"))
options(width = 160)
print(dcast(S, set + n_loci ~ method, value.var = "top1"), digits = 2)
cat("\nmean reciprocal rank:\n")
print(dcast(S, set ~ method, value.var = "mrr"), digits = 2)
cat("\nwrote", file.path(RES, "41_benchmark_summary.csv"), "\n")
