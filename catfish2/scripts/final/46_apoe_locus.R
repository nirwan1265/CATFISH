#!/usr/bin/env Rscript
## =============================================================================
## 46_apoe_locus.R -- within one locus, does the gene score put the gene with
## the relevant function ahead of its linkage-disequilibrium neighbours? The
## chromosome 19q13.32 locus (APOE, APOC1, TOMM40, PVRL2/NECTIN2) in the ten
## GLGC 2021 lipid sets (LDL-C, HDL-C x EUR, AFR, EAS, SAS, HIS): each gene's
## CATFISH rank and its GWAS-layer rank (best variant).
## Reads  catfish2/results/06_candidate_genes_<trait>[_<anc>].csv
## Writes catfish2/tex/auto/tab_apoe_locus.tex, catfish2/results/46_apoe_locus.csv
##   Rscript 46_apoe_locus.R
## =============================================================================
suppressPackageStartupMessages(library(data.table))
BASE <- Sys.getenv("CATFISH_BASE", path.expand("~/Documents/Github/CATFISH"))
RES  <- file.path(BASE, "catfish2/results"); TEX <- Sys.getenv("TEX_AUTO", file.path(BASE, "catfish2/tex/auto"))
GENES <- c("APOE", "APOC1", "TOMM40", "PVRL2"); ANC <- c("eur", "afr", "eas", "sas", "his")
X <- rbindlist(lapply(c("ldl", "hdl"), function(tr) rbindlist(lapply(ANC, function(a) {
  f <- file.path(RES, sprintf("06_candidate_genes_%s%s.csv", tr, if (a == "eur") "" else paste0("_", a)))
  g <- fread(f)[symbol %in% GENES]
  data.table(trait = tr, anc = a, gene = g$symbol, rank = g$rank, rank_gwas = g$rank_gwas,
             in_pathway = !is.na(g$best_pathway) & nzchar(g$best_pathway))
}))))
fwrite(X, file.path(RES, "46_apoe_locus.csv"))
print(dcast(X, gene ~ trait + anc, value.var = "rank"), row.names = FALSE)
print(dcast(X, gene ~ trait + anc, value.var = "rank_gwas"), row.names = FALSE)
w <- dcast(X, trait + anc ~ gene, value.var = "rank")
cat(sprintf("APOE above TOMM40 in %d of %d lists; above PVRL2 in %d\n",
            sum(w$APOE < w$TOMM40), nrow(w), sum(w$APOE < w$PVRL2)))
f0 <- function(v) ifelse(abs(v - round(v)) < 1e-9, sprintf("%d", as.integer(round(v))), sprintf("%.1f", v))
row <- function(gn, col) paste(vapply(c("ldl", "hdl"), function(tr) paste(vapply(ANC, function(a)
  f0(X[trait == tr & anc == a & gene == gn, get(col)]), ""), collapse = " & "), ""), collapse = " & ")
lab <- c(APOE = "\\textit{APOE}", APOC1 = "\\textit{APOC1}", TOMM40 = "\\textit{TOMM40}", PVRL2 = "\\textit{PVRL2}")
lines <- c("\\begin{table*}[!t]", "\\centering",
  "\\caption{Ranking within one locus: the chromosome 19q13.32 apolipoprotein locus in the ten GLGC 2021 lipid sets. CATFISH, the gene's rank under the gene score among about $18{,}500$ genes; GWAS, its rank by its best variant, the GWAS layer alone. \\textit{APOE}, \\textit{APOC1} and \\textit{TOMM40} share their best variant in every set (identical GWAS ranks), so the GWAS layer cannot separate them, and \\textit{PVRL2} (\\textit{NECTIN2}) lies beside them; \\textit{APOE} and \\textit{APOC1} are members of the lipoprotein pathways that lead the pathway ranking, \\textit{TOMM40} and \\textit{PVRL2} are not.}",
  "\\label{tab:apoe_locus}", "\\small",
  "\\begin{tabular}{ll|rrrrr|rrrrr}", "\\hline",
  " & & \\multicolumn{5}{c|}{\\textbf{LDL-C}} & \\multicolumn{5}{c}{\\textbf{HDL-C}} \\\\",
  "\\textbf{Gene} & & EUR & AFR & EAS & SAS & HIS & EUR & AFR & EAS & SAS & HIS \\\\", "\\hline",
  unlist(lapply(GENES, function(gn) c(sprintf("%s & CATFISH & %s \\\\", lab[[gn]], row(gn, "rank")),
                                      sprintf(" & GWAS & %s \\\\", row(gn, "rank_gwas"))))),
  "\\hline", "\\end{tabular}", "\\end{table*}")
writeLines(lines, file.path(TEX, "tab_apoe_locus.tex"))
cat("wrote", file.path(TEX, "tab_apoe_locus.tex"), "\n")
