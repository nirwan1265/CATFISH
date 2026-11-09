#!/usr/bin/env Rscript
## =============================================================================
## 26_sorghum_results.R -- numbers and tables for the sorghum BAP section
## (Section 3.10) and its comparison with the lipid analyses.
##
## Reads only finished outputs; nothing is recomputed:
##   23  ../results/sorghum_bap/<TRAIT>_P_SNPWISE_MEAN_5comp_westfall_young.csv
##       ../results/sorghum_bap/<TRAIT>_gene_cor_pairs.txt
##   25  ../results/sorghum_bap/25_<TRAIT>_gene_prioritize.csv, 25_established_genes.csv
##   MAGMA <MAGMA_ROOT>/<TRAIT>/<TRAIT>_chrNN.multi_snp_wise.genes.out and .pval.txt
##   lipids  <BASE>/results/human_LDL/{magma,magma_hdl}/<trait>_eur.genes.out,
##           <BASE>/results/human_LDL/catfish/{ldl,hdl}_eur_reactome_5comp_westfall_young.csv,
##           <BASE>/results/human_LDL/catfish/ldl_eur.cor_pairs.txt
##
## Usage (from catfish2/scripts):  Rscript 26_sorghum_results.R
## Writes ../tex/auto/sorghum.tex            (\renewcommand macros)
##        ../tex/auto/sorghum_stub.tex       (?? placeholders, \newcommand)
##        ../tex/auto/tab_sorghum_regimes.tex
##        ../tex/auto/tab_sorghum_genes.tex
##        ../results/sorghum_bap/26_summary.csv
## =============================================================================
suppressPackageStartupMessages(library(data.table))

MAGMA_ROOT <- Sys.getenv("MAGMA_ROOT", path.expand("~/Documents/Research/results/CATFISH/MAGMA"))
BASE       <- Sys.getenv("CATFISH_BASE", path.expand("~/Documents/Github/CATFISH"))
LIPID_ROOT <- Sys.getenv("LIPID_ROOT", file.path(BASE, "results/human_LDL"))
RES_DIR    <- Sys.getenv("RES_DIR", "../results/sorghum_bap")
TEX_DIR    <- Sys.getenv("TEX_DIR", "../tex/auto")
SORGHUM_COR <- Sys.getenv("SORGHUM_COR", file.path(RES_DIR, "Lignin_DM_perc_gene_cor_pairs.txt"))
LIPID_COR   <- Sys.getenv("LIPID_COR", file.path(LIPID_ROOT, "catfish/ldl_eur.cor_pairs.txt"))
CHROMS     <- as.integer(strsplit(Sys.getenv("CHROMS", "1,2,3,4,5,6,7,8,9,10"), ",")[[1]])
PVAL_CHROMS <- as.integer(strsplit(Sys.getenv("PVAL_CHROMS", paste(CHROMS, collapse = ",")), ",")[[1]])  # quick checks only
TRAITS     <- c(dry = "Dry_tons_per_acre", lig = "Lignin_DM_perc")
P_COL      <- "P_SNPWISE_MEAN"
dir.create(TEX_DIR, showWarnings = FALSE, recursive = TRUE)

## annotations used in the text. SorghumBicolorCyc enzyme names are read from
## the pathway file; the one override is Sobic.008G029200, which carries a WRKY
## domain (InterPro IPR003657) but is listed in SorghumBicolorCyc as an
## acyl-ACP thioesterase.
PATHWAY_FILE <- Sys.getenv("PATHWAY_FILE",
                  file.path(BASE, "inst/extdata/pathway/sorghumbicolorcyc_pathways.20230103.SORBI"))
ANNOT_OVERRIDE <- c(SORBI_3008G029200 = "WRKY transcription factor$^{a}$")

`%+%` <- function(a, b) paste0(a, b)
lam   <- function(p) qchisq(median(p, na.rm = TRUE), 1, lower.tail = FALSE) / qchisq(0.5, 1)
pct   <- function(x, d = 1) formatC(100 * x, format = "f", digits = d)
num   <- function(x, d = 0) formatC(x, format = "f", digits = d, big.mark = "{,}")
f2    <- function(x, d = 2) formatC(x, format = "f", digits = d)
sci   <- function(p) { e <- floor(log10(p)); sprintf("%.1f\\times10^{%d}", p / 10^e, e) }

M <- list()                                   # macro name -> value
put <- function(name, value) M[[name]] <<- as.character(value)

## ---------------------------------------------------------------- LD summary --
ld_summary <- function(f, n_genes) {
  r <- abs(fread(f, select = 3L, showProgress = FALSE)[[1]])
  c(pairs = length(r), nb01 = 2 * sum(r >= 0.1) / n_genes,
    nb05 = 2 * sum(r >= 0.5) / n_genes, mean_abs = mean(r))
}

## --------------------------------------------------------------- per trait ----
read_genes <- function(trait) {
  f <- file.path(MAGMA_ROOT, trait, sprintf("%s_chr%02d.multi_snp_wise.genes.out", trait, CHROMS))
  g <- rbindlist(lapply(f, fread, showProgress = FALSE)); g[!duplicated(GENE)]
}
snp_lambda <- function(trait) {
  f <- file.path(MAGMA_ROOT, trait, sprintf("%s_chr%02d.pval.txt", trait, PVAL_CHROMS))
  p <- unlist(lapply(f, function(x) { d <- fread(x, showProgress = FALSE); d$P[d$N >= 50] }))
  c(n = length(p), lambda = lam(p))
}
pw_calib <- function(f) {
  w <- fread(f, showProgress = FALSE)
  out <- c(n = nrow(w), omni05 = mean(w$p_omnibus < 0.05),
           lam_omni = lam(w$p_omnibus),
           n_fwer = sum(w$p_fwer_std <= 0.05), min_fwer = min(w$p_fwer_std),
           size_rho = suppressWarnings(cor(w$n_genes, w$log10_p_omnibus, method = "spearman")))
  if ("p_calibrated" %in% names(w))
    out <- c(out, cal05 = mean(w$p_calibrated < 0.05), lam_cal = lam(w$p_calibrated),
             med_cal = median(w$p_calibrated))
  out
}

summ <- list()
genes <- list()
for (k in names(TRAITS)) {
  t <- TRAITS[[k]]
  g <- read_genes(t); genes[[k]] <- g
  s <- snp_lambda(t)
  wy <- file.path(RES_DIR, sprintf("%s_%s_5comp_westfall_young.csv", t, P_COL))
  pc <- pw_calib(wy)
  summ[[k]] <- c(snp_n = s[["n"]], snp_lambda = s[["lambda"]],
                 n_genes = nrow(g), gene_med = median(g[[P_COL]]), gene_lambda = lam(g[[P_COL]]),
                 gene_gw = mean(g[[P_COL]] < 0.05 / nrow(g)), pc)
  cat(sprintf("[%s] SNPs %s, SNP lambda %.3f | genes %d, median p %.3f, gene lambda %.2f | pathways %d: p_omni<0.05 %.1f%%, p_cal<0.05 %.1f%%, FWER sig %d (min %.3f)\n",
              t, format(s[["n"]], big.mark = ","), s[["lambda"]], nrow(g), median(g[[P_COL]]),
              lam(g[[P_COL]]), pc[["n"]], 100 * pc[["omni05"]], 100 * pc[["cal05"]],
              pc[["n_fwer"]], pc[["min_fwer"]]))
}

## lipids, for the comparison table
lip <- list()
for (k in c("ldl", "hdl")) {
  d <- if (k == "ldl") "magma" else "magma_hdl"
  g <- fread(file.path(LIPID_ROOT, d, sprintf("%s_eur.genes.out", k)), showProgress = FALSE)
  pc <- pw_calib(file.path(LIPID_ROOT, "catfish", sprintf("%s_eur_reactome_5comp_westfall_young.csv", k)))
  lip[[k]] <- c(n_genes = nrow(g), gene_med = median(g$P), gene_lambda = lam(g$P), pc)
}

ldS <- ld_summary(SORGHUM_COR, summ$lig[["n_genes"]])
ldL <- ld_summary(LIPID_COR, lip$ldl[["n_genes"]])
cat(sprintf("LD: sorghum %.0f neighbours |r|>=0.1 per gene (mean |r| %.3f); human %.0f (mean |r| %.3f)\n",
            ldS[["nb01"]], ldS[["mean_abs"]], ldL[["nb01"]], ldL[["mean_abs"]]))

## ------------------------------------------------------------ pathway ranks --
wyL <- fread(file.path(RES_DIR, sprintf("%s_%s_5comp_westfall_young.csv", TRAITS[["lig"]], P_COL)))
wyD <- fread(file.path(RES_DIR, sprintf("%s_%s_5comp_westfall_young.csv", TRAITS[["dry"]], P_COL)))
rk  <- function(w, id) match(id, w$pathway_id)
put("nSorPhenylRank", rk(wyL, "PWY-361"))
put("nSorFlavRank",   rk(wyL, "PWY1F-FLAVSYN"))
put("nSorCycloRank",  rk(wyD, "PWY-8028"))
put("nSorTopFwerDry", f2(min(wyD$p_fwer_std)))
put("nSorTopFwerLig", f2(min(wyL$p_fwer_std)))

## ------------------------------------------------------- gene prioritization --
pwf <- fread(PATHWAY_FILE, sep = "\t", quote = "", showProgress = FALSE)
pwf <- pwf[nzchar(`Protein-name`)]
enz <- pwf[, .(name = `Protein-name`[1]), by = `Gene-name`]
annot <- function(id) {
  a <- if (id %in% names(ANNOT_OVERRIDE)) ANNOT_OVERRIDE[[id]] else enz$name[match(id, enz$`Gene-name`)]
  if (is.na(a) || !nzchar(a)) "--" else gsub("&", "\\\\&", a)
}

est <- fread(file.path(RES_DIR, "25_established_genes.csv"))
cand <- list(); est_rows <- list()
for (k in names(TRAITS)) {
  t <- TRAITS[[k]]
  x <- fread(file.path(RES_DIR, sprintf("25_%s_gene_prioritize.csv", t)), showProgress = FALSE)
  ## rank with the pathway layer and hub bonus removed (GWAS + MAGMA only)
  x[, score_np := gwas_score + magma_score]
  setorder(x, -score_np, magma_p); x[, rank_np := .I]; setorder(x, rank)
  g <- genes[[k]]
  x[, mb := (g$START[match(gene, g$GENE)] + g$STOP[match(gene, g$GENE)]) / 2e6]
  x[, chr := g$CHR[match(gene, g$GENE)]]
  cand[[k]] <- x
  e <- est[trait == t]
  e[, rank_np := x$rank_np[match(id, x$gene)]]
  est_rows[[k]] <- e
  n <- nrow(x)
  better <- sum(e$combined_rank < e$gwas_rank, na.rm = TRUE)
  put(sprintf("nSor%sNgenes", if (k == "dry") "Dry" else "Lig"), num(n))
  put(sprintf("nSor%sBetter", if (k == "dry") "Dry" else "Lig"), better)
  put(sprintf("nSor%sNctrl",  if (k == "dry") "Dry" else "Lig"), nrow(e))
  put(sprintf("nSor%sCtrlMedScore", if (k == "dry") "Dry" else "Lig"), num(median(e$combined_rank)))
  put(sprintf("nSor%sCtrlMedGwas",  if (k == "dry") "Dry" else "Lig"), num(median(e$gwas_rank)))
  cat(sprintf("[%s] established genes: median rank combined %.0f, GWAS %.0f, MAGMA %.0f; better than GWAS alone %d of %d\n",
              t, median(e$combined_rank), median(e$gwas_rank), median(e$magma_rank), better, nrow(e)))
}
x <- cand$dry
put("nSorCasRank", x$rank[x$gene == "SORBI_3004G014400"])
put("nSorCasRankNp", x$rank_np[x$gene == "SORBI_3004G014400"])
put("nSorCasMagmaRank", num(x$magma_rank[x$gene == "SORBI_3004G014400"]))
put("nSorCasGwasP", sci(10^-x$gwas_neglog10p[x$gene == "SORBI_3004G014400"]))
put("nSorCasMagmaP", sci(x$magma_p[x$gene == "SORBI_3004G014400"]))
put("nSorWrkyRank", x$rank[x$gene == "SORBI_3008G029200"])
put("nSorWrkyRankNp", x$rank_np[x$gene == "SORBI_3008G029200"])
put("nSorWrkyGwasRank", num(x$gwas_rank[x$gene == "SORBI_3008G029200"]))
put("nSorWrkyLigRank", num(cand$lig$rank[cand$lig$gene == "SORBI_3008G029200"]))
put("nSorCasLigRank", num(cand$lig$rank[cand$lig$gene == "SORBI_3004G014400"]))
ln <- head(cand$lig, 6)
## genes.out coordinates include MAGMA's +/-25 kb window, so remove it once per side
put("nSorLigTopSpanKb", num((max(genes$lig$STOP[match(ln$gene, genes$lig$GENE)]) -
                             min(genes$lig$START[match(ln$gene, genes$lig$GENE)])) / 1000 - 50))
put("nSorLigTopN", nrow(ln))
put("nSorLigTopChr", ln$chr[1]); put("nSorLigTopMb", f2(median(ln$mb), 1))

## ---------------------------------------------------------- summary macros ----
for (k in names(TRAITS)) {
  K <- if (k == "dry") "Dry" else "Lig"; s <- summ[[k]]
  put(sprintf("nSor%sSnpLambda", K), f2(s[["snp_lambda"]]))
  put(sprintf("nSor%sGeneMed", K), f2(s[["gene_med"]], 3))
  put(sprintf("nSor%sGeneLambda", K), f2(s[["gene_lambda"]]))
  put(sprintf("nSor%sNpw", K), num(s[["n"]]))
  put(sprintf("nSor%sOmniFive", K), pct(s[["omni05"]]))
  put(sprintf("nSor%sCalFive", K), pct(s[["cal05"]]))
  put(sprintf("nSor%sNsig", K), s[["n_fwer"]])
  put(sprintf("nSor%sGeneGW", K), pct(s[["gene_gw"]], 3))
  put(sprintf("nSor%sSizeRho", K), f2(s[["size_rho"]], 3))
}
put("nSorNeighbours", num(ldS[["nb01"]])); put("nLipNeighbours", num(ldL[["nb01"]]))
put("nSorMeanAbsR", f2(ldS[["mean_abs"]])); put("nLipMeanAbsR", f2(ldL[["mean_abs"]]))
put("nSorStrongNb", f2(ldS[["nb05"]], 1)); put("nLipStrongNb", f2(ldL[["nb05"]], 1))
put("nSorSnps", num(summ$lig[["snp_n"]]))
put("nSorNgenes", num(summ$lig[["n_genes"]]))

## ------------------------------------------------------------- write macros --
mac_lines <- function(cmd) vapply(names(M), function(n) sprintf("\\%s{\\%s}{%s}", cmd, n, M[[n]]), "")
writeLines(c("%% written by 26_sorghum_results.R -- do not edit", mac_lines("renewcommand")),
           file.path(TEX_DIR, "sorghum.tex"))
writeLines(c("%% Placeholders for the macros written by 26_sorghum_results.R",
             vapply(names(M), function(n) sprintf("\\newcommand{\\%s}{\\textbf{??}}", n), "")),
           file.path(TEX_DIR, "sorghum_stub.tex"))

## ------------------------------------------------- table: the two regimes ----
L <- lip; S <- summ
row  <- function(lab, a, b, c, d) sprintf("%s & %s & %s & %s & %s \\\\", lab, a, b, c, d)
row2 <- function(lab, h, s) sprintf("%s & \\multicolumn{2}{c}{%s} & \\multicolumn{2}{c}{%s} \\\\", lab, h, s)
tabR <- c(
"\\begin{table*}[!t]",
"\\centering",
"\\caption{The two real-data regimes. Correlated neighbours: mean number of other genes per gene whose MAGMA gene--gene correlation is $|r|\\ge0.1$ (human values from the LDL-C run; the BAP traits share one genotype panel). $\\lambda$: genomic-control factor from the median $p$-value. Pathways at $p<0.05$: share of tested pathways below $0.05$ under the uncalibrated analytic omnibus and under the per-pathway calibration against the correlated null of Section~2.3 (not stored for the lipid runs). FWER: pathways significant at family-wise $\\alpha=0.05$.}",
"\\label{tab:sorghum_regimes}",
"\\small",
"\\begin{tabular}{lrrrr}",
"\\hline",
" & \\textbf{LDL-C} & \\textbf{HDL-C} & \\textbf{Dry tons} & \\textbf{Lignin} \\\\",
"\\hline",
row2("Samples", "$\\approx1.3$ million", "330 lines"),
row2("Correlated neighbours per gene", num(ldL[["nb01"]]), num(ldS[["nb01"]])),
row2("Mean $|r|$ between listed gene pairs", f2(ldL[["mean_abs"]]), f2(ldS[["mean_abs"]])),
row2("Pathway database", "Reactome", "SorghumBicolorCyc"),
row("Pathways tested", num(L$ldl[["n"]]), num(L$hdl[["n"]]), num(S$dry[["n"]]), num(S$lig[["n"]])),
row("SNP-level $\\lambda$", "--", "--", f2(S$dry[["snp_lambda"]]), f2(S$lig[["snp_lambda"]])),
row("Median gene $p$", f2(L$ldl[["gene_med"]], 3), f2(L$hdl[["gene_med"]], 3),
    f2(S$dry[["gene_med"]], 3), f2(S$lig[["gene_med"]], 3)),
row("Gene-level $\\lambda$", f2(L$ldl[["gene_lambda"]]), f2(L$hdl[["gene_lambda"]]),
    f2(S$dry[["gene_lambda"]]), f2(S$lig[["gene_lambda"]])),
row("Pathways at $p<0.05$, uncalibrated", pct(L$ldl[["omni05"]]) %+% "\\%", pct(L$hdl[["omni05"]]) %+% "\\%",
    pct(S$dry[["omni05"]]) %+% "\\%", pct(S$lig[["omni05"]]) %+% "\\%"),
row("Pathways at $p<0.05$, calibrated", "--", "--",
    pct(S$dry[["cal05"]]) %+% "\\%", pct(S$lig[["cal05"]]) %+% "\\%"),
row("Pathways significant, FWER", num(L$ldl[["n_fwer"]]), num(L$hdl[["n_fwer"]]),
    num(S$dry[["n_fwer"]]), num(S$lig[["n_fwer"]])),
row2("Validation", "two related traits converge", "calibration under strong LD"),
"\\hline",
"\\end{tabular}",
"\\end{table*}")
writeLines(tabR, file.path(TEX_DIR, "tab_sorghum_regimes.tex"))

## ---------------------------------------------- table: genes and candidates ----
fr <- function(v) ifelse(is.na(v), "--", num(v))
tg <- c(
"\\begin{table*}[!t]",
"\\centering",
"\\caption{Sorghum BAP gene prioritization. \\textbf{A}~Established genes, named before the analysis: rank among all " %+% M$nSorDryNgenes %+% " genes under the combined score, under the score with the pathway layer removed (GWAS + MAGMA), and under each single layer. \\textbf{B}~Five highest-scoring candidates per trait, with the gene's rank in each layer (Pathway: rank among pathway-annotated genes) and its best pathway. Annotations are SorghumBicolorCyc enzyme names. $^{a}$Carries a WRKY domain (InterPro IPR003657); SorghumBicolorCyc lists it as an acyl-ACP thioesterase, so its pathway evidence rests on an annotation that is probably wrong.}",
"\\label{tab:sorghum_genes}",
"\\small",
"\\begin{tabular}{llrrrrl}",
"\\hline",
"\\multicolumn{7}{l}{\\textbf{A}\\quad Established genes} \\\\",
"\\textbf{Gene} & \\textbf{ID} & \\textbf{Combined} & \\textbf{No pathway} & \\textbf{GWAS} & \\textbf{MAGMA} & \\textbf{In pathway} \\\\",
"\\hline")
for (k in names(TRAITS)) {
  e <- est_rows[[k]]
  tg <- c(tg, sprintf("\\multicolumn{7}{l}{\\textit{%s}} \\\\", if (k == "dry") "Dry tons per acre" else "Lignin"))
  for (i in seq_len(nrow(e)))
    tg <- c(tg, sprintf("\\textit{%s} & %s & %s & %s & %s & %s & %s \\\\",
                        e$gene[i], sub("SORBI_", "", e$id[i]), fr(e$combined_rank[i]), fr(e$rank_np[i]),
                        fr(e$gwas_rank[i]), fr(e$magma_rank[i]), if (isTRUE(e$in_pathway[i])) "yes" else "no"))
}
tg <- c(tg, "\\hline",
"\\multicolumn{7}{l}{\\textbf{B}\\quad Top candidates} \\\\",
"\\textbf{Rank} & \\textbf{Gene (locus)} & \\textbf{GWAS} & \\textbf{MAGMA} & \\textbf{Pathway} & \\multicolumn{2}{l}{\\textbf{Annotation; best pathway}} \\\\",
"\\hline")
for (k in names(TRAITS)) {
  x <- head(cand[[k]], 5)
  tg <- c(tg, sprintf("\\multicolumn{7}{l}{\\textit{%s}} \\\\", if (k == "dry") "Dry tons per acre" else "Lignin"))
  for (i in seq_len(nrow(x))) {
    bp <- if (is.na(x$best_pathway_name[i]) || !nzchar(x$best_pathway_name[i])) "" else
            paste0("; ", sub(" \\(.*\\)$", "", gsub("<[^>]+>", "", x$best_pathway_name[i])))
    tg <- c(tg, sprintf("%d & %s (%s:%s) & %s & %s & %s & \\multicolumn{2}{l}{%s%s} \\\\",
                        x$rank[i], sub("SORBI_", "", x$gene[i]), x$chr[i], f2(x$mb[i], 2),
                        fr(x$gwas_rank[i]), fr(x$magma_rank[i]), fr(x$pathway_rank[i]),
                        annot(x$gene[i]), substr(bp, 1, 42)))
  }
}
tg <- c(tg, "\\hline", "\\end{tabular}", "\\end{table*}")
writeLines(tg, file.path(TEX_DIR, "tab_sorghum_genes.tex"))

fwrite(data.table(macro = names(M), value = unlist(M)), file.path(RES_DIR, "26_summary.csv"))
cat(sprintf("wrote %s, %s, %s and %s\n", file.path(TEX_DIR, "sorghum.tex"),
            file.path(TEX_DIR, "sorghum_stub.tex"), file.path(TEX_DIR, "tab_sorghum_regimes.tex"),
            file.path(TEX_DIR, "tab_sorghum_genes.tex")))
