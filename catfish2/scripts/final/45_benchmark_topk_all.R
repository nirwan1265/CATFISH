#!/usr/bin/env Rscript
## =============================================================================
## 45_benchmark_topk_all.R -- the genome-wide gene-prioritisation benchmark of
## 42_benchmark_topk.R, extended from LDL-C and HDL-C to nine European GWAS:
## how many known genes each ranking places in its top k genes overall (the use
## of a ranking for choosing genes to validate).
##
## Traits and GWAS (EUR): LDL-C, HDL-C (GLGC 2021); coronary artery disease
## (van der Harst 2018, GCST005195), breast cancer (Michailidou 2017,
## GCST004988), Crohn's disease and ulcerative colitis (de Lange 2017,
## GCST004132/3), atrial fibrillation (Nielsen 2018, GCST006414), systemic
## lupus erythematosus (Bentham 2015, GCST003156), asthma (UK Biobank,
## GCST90014325), all through the same MAGMA (1000G EUR, +/-10 kb, SNP-wise
## mean), CATFISH (Reactome, B = 10000, gene-level GC) and PoPS v0.2 (default
## settings, Ensembl annotation) pipelines (scripts/hpc/benchmark, hpc/pops).
## Known genes: Open Targets gold-standard genes of the trait, 2019 and 2023
## releases (TRAIT_EFO; LDL-C and HDL-C exactly as in 42).
## Rankings: CATFISH (06 Score: GWAS + MAGMA + 0.5 pathway + 0.3 hub),
##   GWAS + MAGMA (06 S_GWAS + S_MAGMA: the same score with pathway and hub
##   weights 0), GWAS (best variant), MAGMA (gene p), PoPS (PoPS_Score),
##   CATFISH + PoPS (mean of the two ranks).
## Gene universes (every ranking ranked within each):
##   all_MHC   genes scored by both CATFISH (mapped to Ensembl by unique symbol)
##             and PoPS
##   no_MHC    the same outside the extended MHC (chr6 25-34 Mb, GRCh37), which
##             PoPS removes from fitting by default and whose long-range LD fills
##             the top of every list in immune traits (main analysis)
##   annotated no_MHC genes with evidence in all three CATFISH layers (in a
##             tested Reactome pathway)
## Writes catfish2/results/45_benchmark_topk_all.csv, 45_top20_lists_all.csv,
##        catfish2/tex/auto/{tab_benchmark_topk_all.tex, tab_pathway_gain_all.tex,
##        bench_all.tex (+ bench_all_stub.tex)}
## Usage (from catfish2/scripts/final):  Rscript 45_benchmark_topk_all.R  (needs jsonlite)
## =============================================================================
suppressPackageStartupMessages({ library(data.table); library(jsonlite) })
BASE <- Sys.getenv("CATFISH_BASE", path.expand("~/Documents/Github/CATFISH"))
RES  <- Sys.getenv("RES_DIR", file.path(BASE, "catfish2/results"))
GSD  <- Sys.getenv("GS_DIR", file.path(BASE, "catfish2/data/gold_standards"))
ANN  <- Sys.getenv("GENE_ANNOT", file.path(BASE, "catfish2/data/pops/gene_annot_jun10.txt"))
TEX  <- Sys.getenv("TEX_AUTO", file.path(BASE, "catfish2/tex/auto"))
K    <- c(10, 20, 50, 100, 200)
MHC  <- c(25e6, 34e6)

TRAIT_EFO <- list(ldl = c("EFO_0004611", "HP_0003124"), hdl = "EFO_0004612",
                  cad = c("EFO_0001645", "EFO_0000378"), brca = "EFO_0000305",
                  af = "EFO_0000275", cd = "EFO_0000384", uc = "EFO_0000729",
                  sle = "EFO_0002690", asthma = c("MONDO_0004979", "EFO_0000270"))
LABEL <- c(ldl = "LDL-C", hdl = "HDL-C", cad = "Coronary artery disease", brca = "Breast cancer",
           af = "Atrial fibrillation", cd = "Crohn's disease", uc = "Ulcerative colitis",
           sle = "Systemic lupus erythematosus", asthma = "Asthma")
METH <- c("CATFISH", "GWAS+MAGMA", "GWAS", "MAGMA", "PoPS", "CATFISH+PoPS")

ann <- fread(ANN)
sym <- ann[, .N, by = NAME][N == 1, NAME]
s2e <- setNames(ann[NAME %in% sym, ENSGID], ann[NAME %in% sym, NAME]); e2s <- setNames(ann$NAME, ann$ENSGID)
in_mhc <- setNames(ann$CHR == "6" & ann$START <= MHC[2] & ann$END >= MHC[1], ann$ENSGID)
gs <- fread(file.path(GSD, "gwas_gold_standards.191108.tsv"))
J  <- lapply(readLines(file.path(GSD, "otg_gs_230511.json")), fromJSON)
ot <- function(efos) unique(c(
  gs[vapply(strsplit(trait_info.ontology, "|", fixed = TRUE), function(o) any(o %in% efos), logical(1)),
     gold_standard_info.gene_id],
  unlist(lapply(J, function(j) if (any(j$trait_info$ontology %in% efos)) j$gold_standard_info$gene_id))))
cand_file <- function(tr) {
  f <- file.path(RES, sprintf("06_candidate_genes_%s.csv", tr))
  if (tr %in% c("ldl", "hdl") && file.exists(f)) f else
    file.path(RES, "benchmark", "candidates", sprintf("06_candidate_genes_%s.csv", tr))
}

out <- list(); tops <- list(); ov <- list(); GL <- list()
for (tr in names(TRAIT_EFO)) {
  g <- fread(cand_file(tr))
  g[, ENSGID := s2e[symbol]]; g <- g[!is.na(ENSGID) & !duplicated(ENSGID)]
  p <- fread(file.path(RES, "pops", sprintf("%s_eur.preds", tr)))
  for (uni in c("all_MHC", "no_MHC", "annotated")) {
    u <- intersect(g$ENSGID, p$ENSGID)
    if (uni != "all_MHC") u <- u[!in_mhc[u] %in% TRUE]
    if (uni == "annotated") u <- intersect(u, g[layers == 3, ENSGID])
    gg <- g[match(u, ENSGID)]; pp <- p[match(u, ENSGID)]
    gold <- intersect(ot(TRAIT_EFO[[tr]]), u)
    GL[[length(GL) + 1L]] <- data.table(trait = tr, universe = uni, ENSGID = gold,
      symbol = unname(e2s[gold]), entrez = gg$gene[match(gold, gg$ENSGID)])
    rk <- function(s) { s[is.na(s)] <- -Inf; frank(-s, ties.method = "min") }
    R <- list(CATFISH = rk(gg$Score), `GWAS+MAGMA` = rk(gg$S_GWAS + gg$S_MAGMA),
              GWAS = rk(-gg$rank_gwas), MAGMA = rk(-gg$rank_magma), PoPS = rk(pp$PoPS_Score))
    R[["CATFISH+PoPS"]] <- frank((R$CATFISH + R$PoPS) / 2, ties.method = "min")
    for (m in METH) {
      r <- setNames(R[[m]], u); rg <- r[gold]
      out[[length(out) + 1L]] <- data.table(trait = tr, universe = uni, n_genes = length(u),
        n_gold = length(gold), method = m,
        t(setNames(vapply(K, function(k) sum(rg <= k), numeric(1)), paste0("top", K))),
        median_rank = as.numeric(median(rg)))
      if (uni == "no_MHC" && m %in% c("CATFISH", "PoPS", "GWAS+MAGMA")) {
        top <- names(sort(r))[1:20]
        tops[[length(tops) + 1L]] <- data.table(trait = tr, method = m, rank = 1:20, gene = e2s[top],
                                                gold = top %in% gold)
      }
    }
    if (uni == "no_MHC") for (k in c(20, 50, 100)) {
      C <- u[R$CATFISH <= k]; P <- u[R$PoPS <= k]
      ov[[length(ov) + 1L]] <- data.table(trait = tr, k = k, both = sum(intersect(C, P) %in% gold),
        only_CATFISH = sum(setdiff(C, P) %in% gold), only_PoPS = sum(setdiff(P, C) %in% gold))
    }
  }
}
D <- rbindlist(out); T20 <- rbindlist(tops); OV <- rbindlist(ov)
fwrite(D, file.path(RES, "45_benchmark_topk_all.csv")); fwrite(T20, file.path(RES, "45_top20_lists_all.csv"))
fwrite(OV, file.path(RES, "45_overlap_catfish_pops.csv"))
## the gold-standard genes of each trait and universe (CONTROLS_FILE for 22, 27, 28, 29, 39)
fwrite(rbindlist(GL), file.path(RES, "45_gold_genes.csv"))

options(width = 200)
for (uni in c("no_MHC", "all_MHC", "annotated")) {
  x <- D[universe == uni]
  cat(sprintf("\n==== %s ====\n", uni))
  print(dcast(x, trait + n_gold ~ method, value.var = "top20")[order(match(trait, names(TRAIT_EFO)))], row.names = FALSE)
  tot <- x[, lapply(.SD, sum), by = method, .SDcols = paste0("top", K)]
  print(tot[match(METH, method)], row.names = FALSE)
}
print(OV[, lapply(.SD, sum), by = k, .SDcols = c("both", "only_CATFISH", "only_PoPS")])

## ---- macros ----------------------------------------------------------------------------
x  <- D[universe == "no_MHC"]; xa <- D[universe == "all_MHC"]; xl <- D[universe == "annotated"]
tot <- function(d, m, k) sum(d[method == m, get(paste0("top", k))])
M <- list()
tag <- c(CATFISH = "Cat", `GWAS+MAGMA` = "Np", GWAS = "Gw", MAGMA = "Mg", PoPS = "Pops", `CATFISH+PoPS` = "Comb")
kn  <- c(`10` = "Ten", `20` = "Twenty", `50` = "Fifty", `100` = "Hundred", `200` = "TwoHundred")
for (m in names(tag)) for (k in K) {
  M[[sprintf("nBench%s%s", tag[[m]], kn[[as.character(k)]])]] <- tot(x, m, k)
  M[[sprintf("nBenchAll%s%s", tag[[m]], kn[[as.character(k)]])]] <- tot(xa, m, k)
  M[[sprintf("nBenchAnn%s%s", tag[[m]], kn[[as.character(k)]])]] <- tot(xl, m, k)
}
M$nBenchNgold <- sum(x[method == "CATFISH", n_gold]); M$nBenchNgoldAnn <- sum(xl[method == "CATFISH", n_gold])
M$nBenchNtraits <- length(TRAIT_EFO)
M$nBenchGenesLo <- format(min(x$n_genes), big.mark = "{,}"); M$nBenchGenesHi <- format(max(x$n_genes), big.mark = "{,}")
M$nBenchAnnGenesLo <- format(min(xl$n_genes), big.mark = "{,}"); M$nBenchAnnGenesHi <- format(max(xl$n_genes), big.mark = "{,}")
for (k in c(20, 50, 100)) M[[sprintf("nBenchGain%s", kn[[as.character(k)]])]] <-
  sprintf("%.0f", 100 * (tot(x, "CATFISH", k) / tot(x, "GWAS+MAGMA", k) - 1))
o <- OV[, lapply(.SD, sum), by = k, .SDcols = c("both", "only_CATFISH", "only_PoPS")]
M$nBenchBothTwenty <- o[k == 20, both]; M$nBenchOnlyCatTwenty <- o[k == 20, only_CATFISH]; M$nBenchOnlyPopsTwenty <- o[k == 20, only_PoPS]
w <- dcast(x, trait ~ method, value.var = "top100")
M$nBenchTraitsGain <- sum(dcast(x, trait ~ method, value.var = "top50")$CATFISH >= dcast(x, trait ~ method, value.var = "top50")$`GWAS+MAGMA`)
writeLines(c("%% Placeholders for the macros written by 45_benchmark_topk_all.R",
             sprintf("\\providecommand{\\%s}{\\textbf{??}}", names(M))), file.path(TEX, "bench_all_stub.tex"))
writeLines(c("%% written by 45_benchmark_topk_all.R -- do not edit",
             sprintf("\\renewcommand{\\%s}{%s}", names(M), unlist(M))), file.path(TEX, "bench_all.tex"))

## ---- tables --------------------------------------------------------------------------------
two <- function(d, a, b, caption, label, ha, hb) {
  pair <- function(av, bv) c(if (av > bv) sprintf("\\textbf{%d}", av) else sprintf("%d", av),
                             if (bv > av) sprintf("\\textbf{%d}", bv) else sprintf("%d", bv))
  cells <- function(xx) paste(unlist(lapply(c(20, 50, 100), function(k) {
    v <- setNames(xx[, get(paste0("top", k))], xx$method); pair(v[[a]], v[[b]]) })), collapse = " & ")
  rows <- vapply(names(TRAIT_EFO), function(tr)
    sprintf("%s & %d & %s \\\\", LABEL[[tr]], d[trait == tr & method == a, n_gold], cells(d[trait == tr])), "")
  dt <- d[, lapply(.SD, sum), by = method, .SDcols = paste0("top", K)]
  c("\\begin{table}[!t]", "\\centering", sprintf("\\caption{%s}", caption), sprintf("\\label{%s}", label),
    "\\small", "\\setlength{\\tabcolsep}{4pt}", "\\begin{tabular}{lrrrrrrr}", "\\hline",
    " & & \\multicolumn{2}{c}{\\textbf{Top 20}} & \\multicolumn{2}{c}{\\textbf{Top 50}} & \\multicolumn{2}{c}{\\textbf{Top 100}} \\\\",
    sprintf("\\textbf{Trait} & $n$ & %s & %s & %s & %s & %s & %s \\\\", ha, hb, ha, hb, ha, hb), "\\hline", rows, "\\hline",
    sprintf("\\textbf{Total} & %d & %s \\\\", d[method == a, sum(n_gold)], cells(dt)), "\\hline", "\\end{tabular}", "\\end{table}")
}
common <- "Number of Open Targets gold-standard genes of each trait (2019 and 2023 releases) among the top $20$, $50$ and $100$ genes of each ranking genome-wide, within the genes scored by both CATFISH and PoPS and outside the extended MHC (chromosome~6, $25$--$34$~Mb); all rankings use the same GWAS summary statistics and MAGMA gene analysis. $n$, gold-standard genes in the gene universe. Bold, the better ranking."
writeLines(two(x, "GWAS+MAGMA", "CATFISH",
  paste("Pathway evidence and candidate genes in nine European GWAS. GWAS+MAGMA, the CATFISH gene score without its pathway layer (pathway and hub weights $0$); CATFISH, the full score.", common),
  "tab:pathway_gain_all", "\\makecell{GWAS+\\\\MAGMA}", "CATFISH"), file.path(TEX, "tab_pathway_gain_all.tex"))
writeLines(two(x, "CATFISH", "PoPS",
  paste("Known genes recovered by CATFISH and PoPS in nine European GWAS. PoPS with its default settings.", common),
  "tab:benchmark_topk_all", "CATFISH", "PoPS"), file.path(TEX, "tab_benchmark_topk_all.tex"))
cat("\nwrote", file.path(RES, "45_benchmark_topk_all.csv"), "and tex/auto/{tab_benchmark_topk_all,tab_pathway_gain_all,bench_all}.tex\n")
