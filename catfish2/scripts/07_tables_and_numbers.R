#!/usr/bin/env Rscript
## 07_tables_and_numbers.R
## Turn results/01..04 CSVs into (a) LaTeX tables and (b) a numbers.tex file of
## \renewcommand macros that the manuscript prose uses.  Run AFTER 01-04.
## Usage:  Rscript 07_tables_and_numbers.R

RES <- "../results"
AUT <- "../tex/auto"
dir.create(AUT, showWarnings = FALSE, recursive = TRUE)

f2  <- function(x) formatC(x, format = "f", digits = 3)
f4  <- function(x) formatC(x, format = "f", digits = 4)
sci <- function(x) {
  e <- floor(log10(abs(x))); m <- x / 10^e
  sprintf("\\ensuremath{%.2f\\times10^{%d}}", m, e)
}
macros <- character(0)
def <- function(name, value) macros[[length(macros) + 1L]] <<-
  sprintf("\\renewcommand{\\%s}{%s}", name, value)

have <- function(f) file.exists(file.path(RES, f))

## Wrap the tabular in \resizebox so a wide table can never overflow the
## column, then write. Applied uniformly -- an earlier version patched some
## tables and not others, which left one with an unbalanced brace.
emit <- function(lines, file, star = TRUE) {
  w <- if (star) "\\textwidth" else "\\columnwidth"
  i <- grep("^\\\\begin\\{tabular\\}", lines)[1]
  j <- grep("^\\\\end\\{tabular\\}",   lines)[1]
  if (!is.na(i) && !is.na(j)) {
    lines <- append(lines, sprintf("\\resizebox{%s}{!}{%%", w), after = i - 1L)
    lines[j + 1L] <- "\\end{tabular}}"   # closes the resizebox
  }
  writeLines(lines, file)
}

## ---------------------------------------------------------------- 01 null ---
if (have("01_null_calibration.csv")) {
  d <- read.csv(file.path(RES, "01_null_calibration.csv"))
  ## Bonferroni FWER. Older runs stored type1(p, a/NP) * NP -- an expected
  ## number of false rejections per family, not a rate. 01's pathways are
  ## independent by construction, so the family-wise rate follows exactly as
  ## 1 - (1 - x/NP)^NP. Newer runs record it directly as t1_bonf_fwer.
  NPF <- if ("np" %in% names(d)) d$np else 200
  d$bonf_fwer <- if ("t1_bonf_fwer" %in% names(d)) d$t1_bonf_fwer else
                 1 - (1 - pmin(d$t1_bonferroni / NPF, 1))^NPF
  d5 <- d[d$alpha == 0.05, ]
  lines <- c(
    "\\begin{table*}[!t]", "\\centering",
    paste0("\\caption{Empirical type~I error and genomic-control inflation for the ",
           "CATFISH omnibus under the global null, by within-pathway LD level and ",
           "pathway size. Analytic: independence-based omnibus. Calibrated: per-pathway ",
           "calibration against the correlated null. Bonferroni and Westfall--Young (WY): ",
           "family-wise error rates, the probability of at least one rejection among ",
           "the pathways of one simulated family. Nominal level $\\alpha=0.05$.}"),
    "\\label{tab:null_calibration}", "\\small",
    "\\begin{tabular}{llrrrrrrr}", "\\hline",
    "\\textbf{LD} & $\\rho$ & $m$ & \\textbf{Analytic} & \\textbf{Calibrated} & \\textbf{Bonferroni} & \\textbf{WY} & $\\lambda_{\\mathrm{an}}$ & $\\lambda_{\\mathrm{cal}}$ \\\\",
    "\\hline")
  for (i in seq_len(nrow(d5))) {
    r <- d5[i, ]
    lines <- c(lines, sprintf("%s & %.1f & %d & %s & %s & %s & %s & %s & %s \\\\",
      r$ld, r$rho, r$m, f4(r$t1_analytic), f4(r$t1_mvn), f4(r$bonf_fwer),
      f4(r$t1_wy), f2(r$lambda_analytic), f2(r$lambda_mvn)))
  }
  lines <- c(lines, "\\hline", "\\end{tabular}", "\\end{table*}")
  emit(lines, file.path(AUT, "tab_null.tex"))

  g <- function(ld, m, col) d5[[col]][d5$ld == ld & d5$m == m][1]
  def("nTAnInd",   f4(g("independent", 5,  "t1_analytic")))
  def("nTAnStr",   f4(g("strong",      50, "t1_analytic")))
  def("nTMvnInd",  f4(g("independent", 5,  "t1_mvn")))
  def("nTMvnStr",  f4(g("strong",      50, "t1_mvn")))
  def("nTAnMax",   f4(max(d5$t1_analytic)))
  def("nTAnMin",   f4(min(d5$t1_analytic)))
  def("nTMvnMax",  f4(max(d5$t1_mvn)))
  def("nTMvnMin",  f4(min(d5$t1_mvn)))
  def("nTWyMax",   f4(max(d5$t1_wy)))
  def("nTWyMin",   f4(min(d5$t1_wy)))
  def("nLamAnMax", f2(max(d5$lambda_analytic)))
  def("nLamAnMin", f2(min(d5$lambda_analytic)))
  def("nLamMvnMax",f2(max(d5$lambda_mvn)))
  def("nLamMvnMin",f2(min(d5$lambda_mvn)))
  def("nTBonfInd",  f4(g("independent", 5,  "bonf_fwer")))
  def("nTBonfStr",  f4(g("strong",      50, "bonf_fwer")))
  def("nTBonfMod",  f4(g("moderate",    20, "bonf_fwer")))
  def("nLamAnInd",  f2(g("independent", 50, "lambda_analytic")))
  def("nLamAnStr",  f2(g("strong",      50, "lambda_analytic")))
  if ("wy_n_distinct" %in% names(d5))
    def("nWyDistinctMin", sprintf("%d", min(d5$wy_n_distinct)))
  d1 <- d[d$alpha == 0.01, ]
  def("nTAnMaxOne",  f4(max(d1$t1_analytic)))
  def("nTMvnMaxOne", f4(max(d1$t1_mvn)))
}

## ----------------------------------------------------------- 02 breakage ---
if (have("02_breakage_stress.csv")) {
  d <- read.csv(file.path(RES, "02_breakage_stress.csv"), stringsAsFactors = FALSE)
  lines <- c(
    "\\begin{table*}[!t]", "\\centering",
    paste0("\\caption{Component-breakage stress test on overlapping pathways. ",
           "Each of $k$ component tests is miscalibrated by raising its $p$-value ",
           "to the power $\\gamma$; $\\gamma<1$ makes a component conservative and ",
           "$\\gamma>1$ anti-conservative. Marginal columns are per-pathway type~I ",
           "error; FWER columns are family-wise error rates over the $K$ pathways ",
           "tested, measured on held-out draws. Nominal level $\\alpha=0.05$.}"),
    "\\label{tab:breakage}", "\\small",
    "\\begin{tabular}{llrlrrrrrr}", "\\hline",
    "\\textbf{LD} & $\\rho$ & $m$ & \\textbf{Direction} & $\\gamma$ & $k$ & \\multicolumn{2}{c}{\\textbf{Marginal}} & \\multicolumn{2}{c}{\\textbf{FWER}} \\\\",
    " & & & & & & Analytic & MVN & WY & Bonf. \\\\",
    "\\hline")
  d <- d[order(d$direction, d$ld, d$m, d$k_broken), ]
  for (i in seq_len(nrow(d))) {
    r <- d[i, ]
    lines <- c(lines, sprintf("%s & %.1f & %d & %s & %.1f & %d & %s & %s & %s & %s \\\\",
      r$ld, r$rho, r$m, r$direction, r$gamma, r$k_broken,
      f4(r$t1_analytic), f4(r$t1_mvn), f4(r$fwer_wy), f4(r$fwer_bonferroni)))
  }
  lines <- c(lines, "\\hline", "\\end{tabular}", "\\end{table*}")
  emit(lines, file.path(AUT, "tab_breakage.tex"))

  anti <- d[d$direction == "anticonservative", ]
  cons <- d[d$direction == "conservative", ]
  if (!nrow(anti)) anti <- cons
  def("nBrkMembership",   sprintf("%.1f", mean(d$pathways_per_gene)))
  def("nBrkAnKzero",      f4(mean(d$t1_analytic[d$k_broken == 0])))
  def("nBrkMvnKzero",     f4(mean(d$t1_mvn[d$k_broken == 0])))
  def("nBrkWyKzero",      f4(mean(d$fwer_wy[d$k_broken == 0])))
  def("nBrkBonfKzero",    f4(mean(d$fwer_bonferroni[d$k_broken == 0])))
  def("nBrkAnKfour",      f4(mean(anti$t1_analytic[anti$k_broken == 4])))
  def("nBrkMvnKfour",     f4(mean(anti$t1_mvn[anti$k_broken == 4])))
  def("nBrkWyKfour",      f4(mean(anti$fwer_wy[anti$k_broken == 4])))
  def("nBrkBonfKfour",    f4(mean(anti$fwer_bonferroni[anti$k_broken == 4])))
  def("nBrkConsAnKfour",  f4(mean(cons$t1_analytic[cons$k_broken == 4])))
  def("nBrkConsMvnKfour", f4(mean(cons$t1_mvn[cons$k_broken == 4])))
  def("nBrkConsWyKfour",  f4(mean(cons$fwer_wy[cons$k_broken == 4])))
  def("nBrkMvnMax",       f4(max(anti$t1_mvn)))
  def("nBrkAnMax",        f4(max(anti$t1_analytic)))
  def("nBrkWyRangeLo",    f4(min(d$fwer_wy)))
  def("nBrkWyRangeHi",    f4(max(d$fwer_wy)))
  def("nBrkMvnRangeLo",   f4(min(d$t1_mvn)))
  def("nBrkMvnRangeHi",   f4(max(d$t1_mvn)))
}

## ---------------------------------------------------------------- 03 LOO ---
if (have("03_adaptive_loo.csv")) {
  d <- read.csv(file.path(RES, "03_adaptive_loo.csv"), stringsAsFactors = FALSE)
  lines <- c(
    "\\begin{table*}[!t]", "\\centering",
    paste0("\\caption{Leave-one-out component sensitivity on overlapping ",
           "pathways. Each row removes a single component and rebuilds the null ",
           "bank, the marginal calibration and the Westfall--Young threshold from ",
           "scratch. Nominal level $\\alpha=0.05$.}"),
    "\\label{tab:loo}", "\\small",
    "\\begin{tabular}{llrlrrrr}", "\\hline",
    "\\textbf{LD} & $\\rho$ & $m$ & \\textbf{Dropped} & \\multicolumn{2}{c}{\\textbf{Marginal}} & \\multicolumn{2}{c}{\\textbf{FWER}} \\\\",
    " & & & & Analytic & MVN & WY & Bonf. \\\\",
    "\\hline")
  for (i in seq_len(nrow(d))) {
    r <- d[i, ]
    lines <- c(lines, sprintf("%s & %.1f & %d & %s & %s & %s & %s & %s \\\\",
      r$ld, r$rho, r$m, gsub("<|>", "", r$dropped), f4(r$t1_analytic),
      f4(r$t1_mvn), f4(r$fwer_wy), f4(r$fwer_bonferroni)))
  }
  lines <- c(lines, "\\hline", "\\end{tabular}", "\\end{table*}")
  emit(lines, file.path(AUT, "tab_loo.tex"))

  full <- d[d$dropped == "<none>", ]
  drop <- d[d$dropped != "<none>", ]
  def("nLooFullMvn", f4(mean(full$t1_mvn)))
  def("nLooFullWy",  f4(mean(full$fwer_wy)))
  def("nLooMvnMin",  f4(min(drop$t1_mvn)))
  def("nLooMvnMax",  f4(max(drop$t1_mvn)))
  def("nLooWyMin",   f4(min(drop$fwer_wy)))
  def("nLooWyMax",   f4(max(drop$fwer_wy)))
  agg <- aggregate(fwer_wy ~ dropped, data = drop, FUN = mean)
  j <- which.max(abs(agg$fwer_wy - 0.05))
  def("nLooWorst",    agg$dropped[j])
  def("nLooWorstVal", f4(agg$fwer_wy[j]))
  def("nLooThrRatio", sprintf("%.2f", mean(d$bonferroni_threshold / d$wy_threshold)))
  def("nLooBonfInd",  f4(mean(d$fwer_bonferroni[d$ld == "independent"])))
  def("nLooBonfMod",  f4(mean(d$fwer_bonferroni[d$ld == "moderate"])))
  def("nLooBonfStr",  f4(mean(d$fwer_bonferroni[d$ld == "strong"])))
  def("nLooBonfMax",  f4(max(d$fwer_bonferroni)))
  def("nLooWyInd",    f4(mean(d$fwer_wy[d$ld == "independent"])))
  def("nLooWyMod",    f4(mean(d$fwer_wy[d$ld == "moderate"])))
  def("nLooWyStr",    f4(mean(d$fwer_wy[d$ld == "strong"])))
}

## ----------------------------------------------------------- 04 ranking ---
if (have("04_wy_vs_mvn.csv")) {
  d <- read.csv(file.path(RES, "04_wy_vs_mvn.csv"))
  lines <- c(
    "\\begin{table*}[!t]", "\\centering",
    paste0("\\caption{Ranking of truly enriched pathways and family-wise significance ",
           "thresholds. Each replicate contains $N_{\\mathrm{P}}$ pathways of which ",
           "$N_{\\mathrm{true}}$ carry a latent mean shift $\\mu$. Median rank is the ",
           "median within-replicate rank of the enriched pathways (1 = best); random ",
           "expectation is $N_{\\mathrm{P}}/2$. The WY threshold is the $5$th ",
           "percentile of the minimum null omnibus $p$-value across the family; the ",
           "ratio column is Bonferroni/WY.}"),
    "\\label{tab:ranking}", "\\small",
    "\\begin{tabular}{llrrrrrr}", "\\hline",
    "\\textbf{LD} & $\\rho$ & $m$ & \\textbf{Rank (analytic)} & \\textbf{Rank (MVN)} & \\textbf{Random} & \\textbf{WY thresh.} & \\textbf{Bonf./WY} \\\\",
    "\\hline")
  for (i in seq_len(nrow(d))) {
    r <- d[i, ]
    lines <- c(lines, sprintf("%s & %.1f & %d & %.1f & %.1f & %.0f & %s & %.2f \\\\",
      r$ld, r$rho, r$m, r$median_rank_analytic, r$median_rank_mvn,
      r$random_expectation, sci(r$wy_threshold), r$fwer_ratio))
  }
  lines <- c(lines, "\\hline", "\\end{tabular}", "\\end{table*}")
  emit(lines, file.path(AUT, "tab_ranking.tex"))

  def("nRankAn",       sprintf("%.1f", mean(d$median_rank_analytic)))
  def("nRankMvn",      sprintf("%.1f", mean(d$median_rank_mvn)))
  def("nRandExp",      sprintf("%.0f", d$random_expectation[1]))
  def("nNPath",        sprintf("%d",   d$n_pathways[1]))
  def("nNTrue",        sprintf("%d",   d$n_true[1]))
  def("nMu",           sprintf("%.1f", d$mu[1]))
  def("nWyThrMin",     sci(min(d$wy_threshold)))
  def("nWyThrMax",     sci(max(d$wy_threshold)))
  def("nBonfThr",      sci(d$bonferroni_threshold[1]))
  def("nFwerRatioMin", sprintf("%.2f", min(d$fwer_ratio)))
  def("nFwerRatioMax", sprintf("%.2f", max(d$fwer_ratio)))
}

## ------------------------------------------------- 08 missing correlation ---
if (have("08_missing_correlation.csv")) {
  d <- read.csv(file.path(RES, "08_missing_correlation.csv"), stringsAsFactors = FALSE)
  lines <- c(
    "\\begin{table*}[!t]", "\\centering",
    paste0("\\caption{Robustness to an incomplete gene--gene correlation file. ",
           "Data are generated from the true correlation structure; calibration ",
           "is built from a matrix in which the stated fraction of off-diagonal ",
           "entries has been zeroed, as an unreported gene pair would be. ",
           "Marginal columns are per-pathway type~I error, FWER columns are ",
           "family-wise rates over $K$ pathways. Nominal $\\alpha=0.05$.}"),
    "\\label{tab:missing}", "\\small",
    "\\begin{tabular}{lrrrrrrr}", "\\hline",
    "\\textbf{LD} & \\textbf{Dropped} & \\multicolumn{2}{c}{\\textbf{Marginal}} & \\multicolumn{2}{c}{\\textbf{FWER}} & \\multicolumn{2}{c}{$\\lambda$} \\\\",
    " & & Analytic & MVN & WY & Bonf. & Analytic & MVN \\\\",
    "\\hline")
  for (i in seq_len(nrow(d))) {
    r <- d[i, ]
    lines <- c(lines, sprintf("%s & %.0f\\%% & %s & %s & %s & %s & %s & %s \\\\",
      r$ld, 100 * r$miss, f4(r$t1_analytic), f4(r$t1_mvn), f4(r$fwer_wy),
      f4(r$fwer_bonferroni), f2(r$lambda_analytic), f2(r$lambda_mvn)))
  }
  lines <- c(lines, "\\hline", "\\end{tabular}", "\\end{table*}")
  emit(lines, file.path(AUT, "tab_missing.tex"))

  g <- function(ld, ms, col) d[[col]][d$ld == ld & abs(d$miss - ms) < 1e-9][1]
  def("nMissMvnZero",  f4(g("moderate", 0,    "t1_mvn")))
  def("nMissWyZero",   f4(g("moderate", 0,    "fwer_wy")))
  def("nMissMvnHalf",  f4(g("moderate", 0.5,  "t1_mvn")))
  def("nMissWyHalf",   f4(g("moderate", 0.5,  "fwer_wy")))
  def("nMissMvnFull",  f4(g("moderate", 1,    "t1_mvn")))
  def("nMissWyFull",   f4(g("moderate", 1,    "fwer_wy")))
  def("nMissMvnStrFull", f4(g("strong", 1,    "t1_mvn")))
  def("nMissWyStrFull",  f4(g("strong", 1,    "fwer_wy")))
  def("nMissBonfMod",  f4(g("moderate", 0,    "fwer_bonferroni")))
  def("nMissBonfStr",  f4(g("strong",   0,    "fwer_bonferroni")))
  ## ratio to nominal, the fair comparison between a marginal and a family rate
  def("nMissRatioMvnHalf", sprintf("%.1f", g("moderate", 0.5, "t1_mvn") / 0.05))
  def("nMissRatioWyHalf",  sprintf("%.1f", g("moderate", 0.5, "fwer_wy") / 0.05))
  def("nMissLamLo", f2(min(d$lambda_mvn, na.rm = TRUE)))
  def("nMissLamHi", f2(max(d$lambda_mvn, na.rm = TRUE)))
  ## largest missingness at which WY still beats Bonferroni
  bet <- d[d$fwer_wy < d$fwer_bonferroni, ]
  def("nMissCrossover", sprintf("%.0f", 100 * max(bet$miss, na.rm = TRUE)))
}

## ------------------------------------------------ 09 component calibration ---
if (have("09_component_calibration.csv")) {
  d <- read.csv(file.path(RES, "09_component_calibration.csv"), stringsAsFactors = FALSE)
  lines <- c(
    "\\begin{table*}[!t]", "\\centering",
    paste0("\\caption{Per-component calibration. Each component test is ",
           "calibrated against its own null rather than being combined into the ",
           "omnibus, at $m=50$. Analytic columns assume gene independence; ",
           "calibrated columns use the correlated null. Nominal $\\alpha=0.05$.}"),
    "\\label{tab:components_cal}", "\\small",
    "\\begin{tabular}{llrrrrr}", "\\hline",
    "\\textbf{LD} & \\textbf{Component} & \\textbf{Analytic} & \\textbf{Calibrated} & \\textbf{WY} & $\\lambda_{\\mathrm{an}}$ & $\\lambda_{\\mathrm{cal}}$ \\\\",
    "\\hline")
  dm <- d[d$m == max(d$m), ]
  for (i in seq_len(nrow(dm))) {
    r <- dm[i, ]
    lines <- c(lines, sprintf("%s & %s & %s & %s & %s & %s & %s \\\\",
      r$ld, r$component, f4(r$t1_analytic), f4(r$t1_calibrated),
      f4(r$fwer_wy), f2(r$lambda_analytic), f2(r$lambda_calibrated)))
  }
  lines <- c(lines, "\\hline", "\\end{tabular}", "\\end{table*}")
  emit(lines, file.path(AUT, "tab_components_cal.tex"))

  gc2 <- function(ld, cmp, col) dm[[col]][dm$ld == ld & dm$component == cmp][1]
  for (cmp in c("ACAT","Fisher","Stouffer","minP","TFisher","Omnibus")) {
    def(paste0("nCmp", cmp, "Ind"), f4(gc2("independent", cmp, "t1_analytic")))
    def(paste0("nCmp", cmp, "Str"), f4(gc2("strong",      cmp, "t1_analytic")))
    def(paste0("nCmp", cmp, "Cal"), f4(gc2("strong",      cmp, "t1_calibrated")))
  }
  agg <- stats::aggregate(t1_analytic ~ component, data = dm[dm$ld != "independent", ],
                          FUN = mean)
  def("nCmpWorst",  agg$component[which.max(agg$t1_analytic)])
  def("nCmpBest",   agg$component[which.min(agg$t1_analytic)])
  def("nCmpTFisherNA", sprintf("%.1f", 100 * max(dm$frac_na, na.rm = TRUE)))
}

## ------------------------------------------------ candidate genes (human) ---
if (have("06_candidate_genes_ldl.csv")) {
  d <- read.csv(file.path(RES, "06_candidate_genes_ldl.csv"),
                stringsAsFactors = FALSE)
  d <- d[order(d$rank), ]
  lab <- function(x) if (!is.na(x$symbol) && nzchar(x$symbol)) x$symbol else x$gene
  rk  <- function(v, n) if (is.na(v)) "--" else sprintf("%.0f", v)
  shortpw <- function(p) {
    if (is.na(p) || !nzchar(p)) return("--")
    p <- sub("^REACTOME_", "", p)
    p <- gsub("_", " ", p)
    p <- tolower(p)
    if (nchar(p) > 34) p <- paste0(substr(p, 1, 33), "\\ldots")
    p
  }
  NG <- sum(!is.na(d$rank_magma)); NP2 <- sum(!is.na(d$rank_pathway))
  NW <- sum(!is.na(d$rank_gwas))

  mk <- function(dd, caption, label, star) {
    env <- if (star) "table*" else "table"
    c(sprintf("\\begin{%s}[!t]", env), "\\centering",
      paste0("\\caption{", caption, "}"), paste0("\\label{", label, "}"),
      "\\small",
      "\\begin{tabular}{rlrrrrl}", "\\hline",
      "\\textbf{Rank} & \\textbf{Gene} & \\textbf{Score} & \\textbf{GWAS} & \\textbf{MAGMA} & \\textbf{Pathway} & \\textbf{Leading pathway} \\\\",
      "\\hline",
      vapply(seq_len(nrow(dd)), function(i) {
        x <- dd[i, ]
        sprintf("%d & \\textit{%s} & %.2f & %s & %s & %s & %s \\\\",
                x$rank, lab(x), x$Score, rk(x$rank_gwas), rk(x$rank_magma),
                rk(x$rank_pathway), shortpw(x$best_pathway))
      }, character(1)),
      "\\hline", "\\end{tabular}", sprintf("\\end{%s}", env))
  }

  cap_main <- sprintf(paste0("Top five prioritized genes for LDL cholesterol. ",
    "Rank columns are the gene's position within each evidence layer: locus-level ",
    "GWAS (%d genes with an overlapping SNP), MAGMA gene-level association ",
    "(%d genes), and CATFISH pathway enrichment (%d pathway-annotated genes). ",
    "Score is the composite of Section~2.5. The full top fifty are in ",
    "Supplementary Table~S1."), NW, NG, NP2)
  emit(mk(utils::head(d, 5), cap_main, "tab:candidates", FALSE),
       file.path(AUT, "tab_candidates.tex"), star = FALSE)

  cap_supp <- sprintf(paste0("Top fifty prioritized genes for LDL cholesterol, ",
    "with per-layer ranks. Layer sizes as in Table~\\ref{tab:candidates} ",
    "(%d / %d / %d genes)."), NW, NG, NP2)
  emit(mk(utils::head(d, 50), cap_supp, "tab:candidates_supp", TRUE),
       file.path(AUT, "tab_candidates_supp.tex"), star = TRUE)

  ## positive controls named in advance
  CTRL <- c("LDLR","PCSK9","APOB","APOE","SORT1","HMGCR","NPC1L1","ABCG8",
            "LPA","CETP","LDLRAP1","MYLIP","ANGPTL3")
  idx <- match(CTRL, d$symbol); idx <- idx[!is.na(idx)]
  def("nCandN",        format(nrow(d), big.mark = ","))
  def("nCandNGwas",    format(NW,  big.mark = ","))
  def("nCandNPath",    format(NP2, big.mark = ","))
  def("nCandCtrlFound", sprintf("%d", length(idx)))
  def("nCandCtrlTotal", sprintf("%d", length(CTRL)))
  def("nCandCtrlMedian", sprintf("%.0f", stats::median(d$rank[idx])))
  def("nCandCtrlTopPct", sprintf("%d", sum(d$rank[idx] <= nrow(d) * 0.01)))
  def("nCandRandom",   format(round(nrow(d) / 2), big.mark = ","))
  j <- match("LDLR",  d$symbol); if (!is.na(j)) def("nCandLDLR",    sprintf("%d", d$rank[j]))
  j <- match("PCSK9", d$symbol); if (!is.na(j)) def("nCandPCSKnine", sprintf("%d", d$rank[j]))
  ## how many distinct values the GWAS layer resolves at its strongest point
  pg <- d$p_gwas[!is.na(d$p_gwas)]
  def("nCandGwasTies", sprintf("%d", sum(pg <= min(pg) * (1 + 1e-9))))
}

## -------------------------------------------------- 17 archetype power ---
## Verification: for each archetype and size, which component actually wins
## under the family-wise rule, and is that inside the Table 1 prediction set?
ARCHF <- c("17_archetype_power_wy_tau00100505.csv", "17_archetype_power_wy.csv")
ARCHF <- ARCHF[vapply(ARCHF, have, logical(1))]
if (length(ARCHF)) {
  d <- utils::read.csv(file.path(RES, ARCHF[1]), stringsAsFactors = FALSE)
  PRED <- list(SSS = c("ACAT", "minP"),
               MDS = c("Fisher", "Stouffer", "TFisher"),
               WDE = c("Stouffer", "Fisher"),
               MSA = c("TFisher", "Fisher", "ACAT"),
               SGP = c("minP", "ACAT"))
  ARCH <- intersect(c("SSS", "MDS", "WDE", "MSA", "SGP"), d$archetype)
  cmp  <- d[d$method != "Omnibus", ]
  rows <- list(); nin <- 0L; ntot <- 0L
  for (aa in ARCH) for (mm in sort(unique(d$m))) {
    s <- cmp[cmp$archetype == aa & cmp$m == mm, ]
    ## informative cells only: ignore effect scales where nothing or everything fires
    s <- s[ave(s$power_wy, s$effect_scale, FUN = max) > 0.05 &
           ave(s$power_wy, s$effect_scale, FUN = max) < 0.98, ]
    if (!nrow(s)) next
    tb  <- sort(table(s$method[ave(s$power_wy, s$effect_scale,
                                   FUN = function(x) x == max(x)) == 1]),
                decreasing = TRUE)
    win <- names(tb)[1]
    om  <- d[d$archetype == aa & d$m == mm & d$method == "Omnibus", ]
    bst <- tapply(s$power_wy, s$effect_scale, max)
    omp <- om$power_wy[match(names(bst), om$effect_scale)]
    ok  <- win %in% PRED[[aa]]
    nin <- nin + ok; ntot <- ntot + 1L
    rows[[length(rows) + 1L]] <- sprintf(
      "%s & %d & %s & %s & %s & %.3f \\\\",
      aa, mm, paste(PRED[[aa]], collapse = ", "), win,
      if (ok) "\\checkmark" else "--",
      max(bst - omp, na.rm = TRUE))
  }
  emit(c("\\begin{table*}[!t]", "\\centering",
         paste0("\\caption{Archetype verification under the family-wise decision ",
                "rule. For each archetype and pathway size, the component achieving ",
                "the highest family-wise power most often across effect scales, ",
                "against the set predicted in Table~\\ref{tab:archetype_summary}. ",
                "Cells in which no component or every component fires are excluded ",
                "as uninformative. Worst regret is the largest gap in a single cell ",
                "between the best component and the omnibus.}"),
         "\\label{tab:archverify}", "\\small",
         "\\begin{tabular}{llllcr}", "\\hline",
         paste("\\textbf{Archetype} & \\textbf{$m$} & \\textbf{Predicted} &",
               "\\textbf{Observed winner} & \\textbf{In set} & \\textbf{Worst regret} \\\\"),
         "\\hline", unlist(rows), "\\hline", "\\end{tabular}", "\\end{table*}"),
       file.path(AUT, "tab_archetype_verification.tex"))
  def("nArchCellsIn",  sprintf("%d", nin))
  def("nArchCellsTot", sprintf("%d", ntot))
}

## Regret summary table
REGF <- c("17_archetype_regret_wy_tau00100505.csv", "17_archetype_regret_wy.csv")
REGF <- REGF[vapply(REGF, have, logical(1))]
if (length(REGF)) {
  r <- utils::read.csv(file.path(RES, REGF[1]), stringsAsFactors = FALSE)
  r <- r[order(r$mean_regret), ]
  emit(c("\\begin{table*}[!t]", "\\centering",
         paste0("\\caption{Power regret of the ACAT-O omnibus by archetype. ",
                "Regret is the difference between the best single component and ",
                "the omnibus in the same cell, over all pathway sizes and effect ",
                "scales. The final two columns give the fraction of cells in which ",
                "the omnibus is within $0.05$ and $0.10$ of the best component.}"),
         "\\label{tab:archregret}", "\\small",
         "\\begin{tabular}{llrrrrr}", "\\hline",
         paste("\\textbf{Code} & \\textbf{Archetype} & \\textbf{Mean} &",
               "\\textbf{Median} & \\textbf{Worst} & \\textbf{$\\le0.05$} &",
               "\\textbf{$\\le0.10$} \\\\"),
         "\\hline",
         sprintf("%s & %s & %.3f & %.3f & %.3f & %.0f\\%% & %.0f\\%% \\\\",
                 r$archetype, r$archetype_label, r$mean_regret, r$median_regret,
                 r$max_regret, 100 * r$within_0.05, 100 * r$within_0.10),
         "\\hline", "\\end{tabular}", "\\end{table*}"),
       file.path(AUT, "tab_archetype_regret.tex"))
  def("nArchRegretBest",  sprintf("%.3f", min(r$mean_regret)))
  def("nArchRegretWorst", sprintf("%.3f", max(r$mean_regret)))
  def("nArchRegretBestA",  r$archetype[which.min(r$mean_regret)])
  def("nArchRegretWorstA", r$archetype[which.max(r$mean_regret)])
  def("nArchMaxRegret",   sprintf("%.3f", max(r$max_regret)))
}

## ------------------------------------- 17: component ranks across cells ---
## In each informative (m, archetype, effect) cell, rank the five components by
## family-wise power. Reports how often each finishes where; the point is that
## TFisher is never last even though it seldom wins.
ARCHF2 <- c("17_archetype_power_wy_tau00100505.csv", "17_archetype_power_wy.csv")
ARCHF2 <- ARCHF2[vapply(ARCHF2, have, logical(1))]
if (length(ARCHF2)) {
  d <- utils::read.csv(file.path(RES, ARCHF2[1]), stringsAsFactors = FALSE)
  d <- d[d$method != "Omnibus", ]
  key <- paste(d$m, d$archetype, d$effect_scale)
  mx  <- ave(d$power_wy, key, FUN = max)
  d   <- d[mx > 0.05 & mx < 0.98, ]
  key <- paste(d$m, d$archetype, d$effect_scale)
  d$rk <- ave(-d$power_wy, key, FUN = function(x) rank(x, ties.method = "min"))
  tf  <- d[d$method == "TFisher", ]
  def("nArchInfCells", sprintf("%d", length(unique(key))))
  def("nTFFirst",  sprintf("%d", sum(tf$rk == 1)))
  def("nTFSecond", sprintf("%d", sum(tf$rk == 2)))
  def("nTFThird",  sprintf("%d", sum(tf$rk == 3)))
  def("nTFWorst",  sprintf("%d", max(tf$rk)))
  fi <- d[d$method == "Fisher", ]
  def("nFisherFirst", sprintf("%d", sum(fi$rk == 1)))
  def("nFisherFourthPlus", sprintf("%d", sum(fi$rk >= 4)))
  mp <- d[d$method == "minP", ]
  def("nMinPLast", sprintf("%d", sum(mp$rk == 5)))
}

## --------------------------------- 03: WY threshold against Bonferroni ---
if (have("03_adaptive_loo.csv")) {
  d <- utils::read.csv(file.path(RES, "03_adaptive_loo.csv"), stringsAsFactors = FALSE)
  d <- d[d$dropped == "<none>", ]
  ind <- d[d$rho == 0, ]; ld <- d[d$rho > 0, ]
  ## at rho = 0 the WY threshold is MORE permissive than Bonferroni (overlap);
  ## under LD it is far stricter, because it is set on the analytic scale and
  ## the analytic p is anti-conservative there
  def("nWyGainLo", sprintf("%.2f", min(ind$wy_threshold / ind$bonferroni_threshold)))
  def("nWyGainHi", sprintf("%.2f", max(ind$wy_threshold / ind$bonferroni_threshold)))
  def("nWyStrictLo", sprintf("%.0f", min(ld$bonferroni_threshold / ld$wy_threshold)))
  def("nWyStrictHi", format(round(max(ld$bonferroni_threshold / ld$wy_threshold)),
                            big.mark = ","))
  def("nBonfThrK", sci(d$bonferroni_threshold[1]))
}

writeLines(c("%% AUTO-GENERATED by 07_tables_and_numbers.R -- do not edit by hand",
             unlist(macros)),
           file.path(AUT, "numbers.tex"))
cat(sprintf("wrote %d macros to %s/numbers.tex\n", length(macros), AUT))
cat("tables written to", AUT, "\n")
