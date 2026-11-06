#!/usr/bin/env Rscript
## =============================================================================
## 21_lipid_results.R -- figures, tables and numbers for the LDL-C / HDL-C
## section of the manuscript (Section 3.9).
##
## Run from catfish2/scripts AFTER
##   18_westfall_young_trait.R   for ldl and hdl (scripts/human_ldl)
##   06_candidate_genes.R        SPECIES=human TRAIT=ldl and TRAIT=hdl
##   20_lipid_manhattan_prep.R   (optional; without it Figure 13 shows the
##                                gene-level panels only)
##
##   Rscript 21_lipid_results.R
##
## Writes
##   ../figures/Fig13_lipid_manhattan.png     SNP and gene Manhattan, both traits
##   ../figures/Fig14_lipid_pathway_size.png  pathway statistic against size
##   ../figures/Fig15_lipid_gene_recovery.png known genes: GWAS vs MAGMA vs score
##   ../tex/auto/tab_lipid_pathways.tex       top ten pathways per trait
##   ../tex/auto/tab_lipid_candidates.tex     top ten candidate genes per trait
##   ../tex/auto/tab_lipid_recovery.tex       known-gene recovery by method
##   ../tex/auto/lipids.tex                   \renewcommand numbers for the text
## =============================================================================

BASE   <- Sys.getenv("CATFISH_BASE", path.expand("~/Documents/Github/MAGCAT"))
RES    <- "../results"
FIG    <- "../figures"
AUT    <- "../tex/auto"
WY_DIR <- Sys.getenv("LIPID_WY_DIR", file.path(BASE, "results/human_LDL/catfish"))
MAGMA_OUT <- c(
  ldl = Sys.getenv("MAGMA_GENES_OUT_LDL", file.path(BASE, "results/human_LDL/magma/ldl_eur.genes.out")),
  hdl = Sys.getenv("MAGMA_GENES_OUT_HDL", file.path(BASE, "results/human_LDL/magma_hdl/hdl_eur.genes.out")))
for (d in c(FIG, AUT)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

TRAITS <- c(ldl = "LDL-C", hdl = "HDL-C")
## positive controls, named before prioritisation (identical to 06)
CONTROLS <- list(
  ldl = c("LDLR","PCSK9","APOB","APOE","SORT1","HMGCR","NPC1L1","ABCG8",
          "LPA","CETP","LDLRAP1","MYLIP","ANGPTL3"),
  hdl = c("CETP","LIPC","LIPG","ABCA1","APOA1","APOA5","LCAT","SCARB1",
          "PLTP","LPL","ANGPTL4","CUBN"))
## Reactome "Plasma lipoprotein assembly, remodeling, and clearance" and the
## sets beneath it in the Reactome hierarchy
LIPO_RX <- "^REACTOME_(PLASMA_LIPOPROTEIN|LDL_|HDL_|VLDL_|VLDLR_|CHYLOMICRON_)"
TOPPW   <- "REACTOME_PLASMA_LIPOPROTEIN_ASSEMBLY_REMODELING_AND_CLEARANCE"
PPARA   <- "REACTOME_REGULATION_OF_LIPID_METABOLISM_BY_PPARALPHA"
SIZE_BAND <- 1.25          # size-matched = within a factor 1.25 in gene count
COL_LIPO <- "#C44E52"; COL_GWAS <- "#4C72B0"; COL_MAGMA <- "#DD8452"
COL_SCORE <- "#55A868"

need <- function(f) if (!file.exists(f)) stop("missing input: ", f, call. = FALSE)
pretty_pw <- function(id, width = 52) {
  ## readable Reactome name; truncate on plain text first, then restore
  ## capitalisation of acronyms (so a cut never lands inside TeX markup)
  x <- tolower(gsub("_", " ", sub("^REACTOME_", "", id)))
  cut <- nchar(x) > width
  x[cut] <- substr(x[cut], 1, width - 1)
  fix <- c(rrna = "rRNA", mrna = "mRNA", trna = "tRNA", er = "ER", atp = "ATP",
           pparalpha = "PPAR$\\alpha$", ldl = "LDL", hdl = "HDL", vldl = "VLDL",
           vldlr = "VLDLR", mitf = "MITF", wdr5 = "WDR5", srebp = "SREBP",
           srebf = "SREBF", nr1h2 = "NR1H2", nr1h3 = "NR1H3", tp53 = "TP53",
           gtpase = "GTPase", gtpases = "GTPases", dna = "DNA", rna = "RNA")
  x <- vapply(strsplit(x, " ", fixed = TRUE), function(w) {
    k <- match(w, names(fix)); w[!is.na(k)] <- fix[k[!is.na(k)]]
    paste(w, collapse = " ") }, character(1))
  first_fixed <- vapply(strsplit(x, " ", fixed = TRUE), function(w) w[1] %in% fix, logical(1))
  x[!first_fixed] <- sub("^([a-z])", "\\U\\1", x[!first_fixed], perl = TRUE)
  ifelse(cut, paste0(x, "\\ldots"), x)
}
tex_esc <- function(x) gsub("([&%_#$])", "\\\\\\1", x)
fmt_int <- function(x) formatC(x, format = "d", big.mark = "{,}")

macros <- character(0); macro_names <- character(0)
def <- function(name, value) {
  macro_names[[length(macro_names) + 1L]] <<- name
  macros[[length(macros) + 1L]] <<- sprintf("\\renewcommand{\\%s}{%s}", name, value)
}

## ---------------------------------------------------------------- inputs ---
W <- list(); G <- list(); M <- list()
for (t in names(TRAITS)) {
  f <- file.path(WY_DIR, sprintf("%s_eur_reactome_5comp_westfall_young.csv", t))
  need(f); w <- read.csv(f, stringsAsFactors = FALSE)
  w <- w[order(w$z_std), ]; w$rank <- seq_len(nrow(w))
  ## size-matched percentile: share of pathways of similar size (within a
  ## factor SIZE_BAND) whose standardised statistic is more extreme
  w$size_pct <- vapply(seq_len(nrow(w)), function(i) {
    pool <- w$n_genes >= w$n_genes[i] / SIZE_BAND & w$n_genes <= w$n_genes[i] * SIZE_BAND
    mean(w$z_std[pool] < w$z_std[i])
  }, numeric(1))
  W[[t]] <- w
  f <- file.path(RES, sprintf("06_candidate_genes_%s.csv", t))
  need(f); G[[t]] <- read.csv(f, stringsAsFactors = FALSE)
  need(MAGMA_OUT[[t]])
  M[[t]] <- utils::read.table(MAGMA_OUT[[t]], header = TRUE, stringsAsFactors = FALSE)
}
B   <- round(1 / min(W$ldl$p_fwer_wy)) - 1
NPW <- nrow(W$ldl)
NG  <- nrow(G$ldl)
def("nLipB", fmt_int(B)); def("nLipNpw", fmt_int(NPW)); def("nLipNgenes", fmt_int(NG))
pcol <- if (any(grepl("^p_pathway$", names(G$ldl)))) "p_pathway" else NA

## ================================================= Figure 13: Manhattan ====
chr_layout <- function(chr, pos) {
  mx <- tapply(pos, chr, max); mx <- mx[order(as.integer(names(mx)))]
  off <- c(0, cumsum(as.numeric(mx)))[seq_along(mx)]; names(off) <- names(mx)
  list(off = off, mid = off + mx / 2)
}
manhattan <- function(chr, pos, y, cap, thr, main, lab_x = NULL, lab_y = NULL,
                      lab = NULL, ylab) {
  L <- chr_layout(chr, pos)
  x <- pos + L$off[as.character(chr)]
  over <- y > cap
  col <- ifelse(chr %% 2 == 1, "#4C72B0", "#9DB4D8")
  plot(x[!over], y[!over], pch = 16, cex = 0.25, col = col[!over],
       xaxt = "n", xlab = "chromosome", ylab = ylab, ylim = c(0, cap * 1.08),
       main = main, xaxs = "i", bty = "l")
  if (any(over)) points(x[over], rep(cap, sum(over)), pch = 17, cex = 0.55,
                        col = col[over])
  axis(1, at = L$mid, labels = names(L$mid), cex.axis = 0.55, tick = FALSE)
  abline(h = thr, lty = 2, col = "grey40")
  if (length(lab)) {
    lx <- lab_x + L$off[as.character(attr(lab_x, "chr"))]
    ## alternate labels left/right of their point, in genome order, so labels
    ## at nearby positions do not overprint
    side <- integer(length(lx)); side[order(lx)] <- rep(c(2L, 4L), length.out = length(lx))
    text(lx, pmin(lab_y, cap), lab, pos = side, cex = 0.55, col = "#8B2E33",
         font = 3, offset = 0.3, xpd = NA)
  }
  if (any(over)) mtext(sprintf("capped at %g (triangles)", cap), 3, line = 0.1,
                       adj = 1, cex = 0.55, col = "grey40")
}
snp_files <- file.path(RES, sprintf("20_manhattan_snps_%s.tsv", names(TRAITS)))
have_snp <- all(file.exists(snp_files))
CAP_MAX <- 60                 # y-axis cap; stronger points drawn as triangles
cap_for <- function(y) min(CAP_MAX, ceiling(max(y, na.rm = TRUE) * 1.05) + 1)
png(file.path(FIG, "Fig13_lipid_manhattan.png"),
    width = if (have_snp) 2600 else 1500, height = 1300, res = 175)
op <- par(mfrow = c(2, if (have_snp) 2 else 1), mar = c(3.6, 4.2, 2.6, 1),
          mgp = c(2.2, 0.6, 0))
for (t in names(TRAITS)) {
  if (have_snp) {
    s <- utils::read.delim(gzfile(snp_files[match(t, names(TRAITS))]),
                           stringsAsFactors = FALSE)
    manhattan(s$CHR, s$BP, s$NEGLOG10P, cap_for(s$NEGLOG10P), -log10(5e-8),
              sprintf("%s  GWAS (variants)", TRAITS[[t]]),
              ylab = expression(-log[10]~italic(p)))
  }
  m <- M[[t]]
  y <- -log10(pmax(m$P, 1e-320))
  g <- G[[t]]
  ## label the top candidates, one per locus (skip a gene within 3 Mb of one
  ## already labelled) so neighbouring labels do not overprint
  go <- g[order(g$rank), ]
  ii <- match(as.character(go$gene), as.character(m$GENE))
  pick <- integer(0)
  for (k in seq_len(min(nrow(go), 60))) {
    j <- ii[k]; if (is.na(j)) next
    near <- any(m$CHR[ii[pick]] == m$CHR[j] &
                abs(m$START[ii[pick]] - m$START[j]) < 3e6)
    if (!near) pick <- c(pick, k)
    if (length(pick) == 6L) break
  }
  top <- go[pick, ]
  i <- ii[pick]
  lx <- (m$START[i] + m$STOP[i]) / 2; attr(lx, "chr") <- m$CHR[i]
  manhattan(m$CHR, (m$START + m$STOP) / 2, y, cap_for(y), -log10(0.05 / nrow(m)),
            sprintf("%s  MAGMA (genes); labels: top CATFISH candidate per locus", TRAITS[[t]]),
            lab_x = lx, lab_y = y[i], lab = top$symbol,
            ylab = expression(-log[10]~italic(p)[gene]))
}
par(op); dev.off(); cat("Fig13_lipid_manhattan.png\n")

## ===================================== Figure 14: pathway against size ====
png(file.path(FIG, "Fig14_lipid_pathway_size.png"), width = 2300, height = 950,
    res = 175)
op <- par(mfrow = c(1, 2), mar = c(4.4, 4.6, 2.8, 1))
for (t in names(TRAITS)) {
  w <- W[[t]]
  lipo <- grepl(LIPO_RX, w$pathway_id)
  x <- log10(w$n_genes); y <- w$z_std
  plot(x, y, pch = 16, cex = 0.45, col = "grey75", xaxt = "n",
       xlab = "pathway size (genes, log scale)",
       ylab = expression(standardised~statistic~italic(z)[std]~"(lower = stronger)"),
       main = sprintf("%s: %s Reactome pathways", TRAITS[[t]],
                      format(nrow(w), big.mark = ",")))
  axis(1, at = log10(c(10, 20, 50, 100, 200)), labels = c(10, 20, 50, 100, 200))
  lo <- stats::lowess(x, y, f = 0.3)
  lines(lo, col = "grey30", lwd = 2)
  points(x[lipo], y[lipo], pch = 16, cex = 1.0, col = COL_LIPO)
  k <- match(TOPPW, w$pathway_id)
  if (!is.na(k)) {
    points(x[k], y[k], pch = 1, cex = 2.4, lwd = 2, col = COL_LIPO)
    text(x[k], y[k], sprintf("plasma lipoprotein assembly,\nremodeling and clearance (rank %d)",
                             w$rank[k]), pos = 3, cex = 0.6, col = "#8B2E33", offset = 0.9)
  }
  legend("bottomleft", bty = "n", cex = 0.72, pch = c(16, 16, NA),
         lty = c(NA, NA, 1), lwd = c(NA, NA, 2),
         col = c("grey75", COL_LIPO, "grey30"),
         legend = c("all pathways", "plasma lipoprotein subtree",
                    "size trend (lowess)"))
}
par(op); dev.off(); cat("Fig14_lipid_pathway_size.png\n")

## ======================================= Figure 15: known-gene recovery ====
rk_of <- function(g, ctrl, col) g[[col]][match(ctrl, g$symbol)]
png(file.path(FIG, "Fig15_lipid_gene_recovery.png"), width = 2200, height = 1000,
    res = 175)
op <- par(mfrow = c(1, 2), mar = c(4.2, 4.6, 2.8, 5.5))
meth <- c(GWAS = "rank_gwas", MAGMA = "rank_magma", CATFISH = "rank")
for (t in names(TRAITS)) {
  g <- G[[t]]; ctrl <- CONTROLS[[t]]
  R <- sapply(meth, function(cc) rk_of(g, ctrl, cc))
  rownames(R) <- ctrl
  R <- R[stats::complete.cases(R), , drop = FALSE]
  yl <- c(max(R) * 1.3, 0.8)
  plot(NA, xlim = c(0.7, 3.3), ylim = yl, log = "y", xaxt = "n",
       xlab = "", ylab = "rank among all genes (log scale, 1 = top)",
       main = sprintf("%s: %d established genes", TRAITS[[t]], nrow(R)))
  axis(1, at = 1:3, labels = c("GWAS\n(best variant)", "MAGMA\n(gene test)",
                               "CATFISH\nscore"), padj = 0.5, cex.axis = 0.8)
  abline(h = c(10, 50), lty = 3, col = "grey70")
  for (i in seq_len(nrow(R)))
    lines(1:3, R[i, ], col = "grey65")
  points(rep(1, nrow(R)), R[, 1], pch = 16, col = COL_GWAS)
  points(rep(2, nrow(R)), R[, 2], pch = 16, col = COL_MAGMA)
  points(rep(3, nrow(R)), R[, 3], pch = 16, col = COL_SCORE)
  med <- apply(R, 2, stats::median)
  segments((1:3) - 0.18, med, (1:3) + 0.18, med, lwd = 3)
  text(3.15, R[, 3], rownames(R), pos = 4, cex = 0.55, font = 3, xpd = NA)
  mtext(sprintf("median rank  %s", paste(sprintf("%s %g", names(med), med),
                                          collapse = "  |  ")),
        3, line = 0.1, cex = 0.62, col = "grey30")
}
par(op); dev.off(); cat("Fig15_lipid_gene_recovery.png\n")

## =========================================== Table: top ten pathways ====
other <- c(ldl = "hdl", hdl = "ldl")
L <- c("\\begin{table*}[!t]", "\\centering",
       sprintf(paste0("\\caption{Ten highest-ranked Reactome pathways for each trait, ordered by ",
                      "the standardised Westfall--Young statistic $z_{\\mathrm{std}}$ ",
                      "(lower is stronger). All twenty have a family-wise $p$-value ",
                      "at the resolution limit $1/(B+1)=%s$ ($B=%s$). ",
                      "$-\\log_{10}p$ is the uncalibrated omnibus. Size-matched: ",
                      "percentage of pathways within a factor %.2f in size with a more ",
                      "extreme statistic. Lipoprotein-subtree sets in bold.}"),
               formatC(1 / (B + 1), format = "g", digits = 3), fmt_int(B), SIZE_BAND),
       "\\label{tab:lipid_pathways}", "\\small",
       "\\begin{tabular}{rlrrrrr}", "\\hline",
       paste0("\\textbf{Rank} & \\textbf{Pathway} & \\textbf{Genes} & ",
              "$\\boldsymbol{-\\log_{10}p}$ & $\\boldsymbol{z_{\\mathrm{std}}}$ & ",
              "\\textbf{Size-matched} & \\textbf{Rank, other trait} \\\\"),
       "\\hline")
for (t in names(TRAITS)) {
  w <- W[[t]]; wo <- W[[other[[t]]]]
  L <- c(L, sprintf("\\multicolumn{7}{l}{\\textit{%s}} \\\\", TRAITS[[t]]))
  for (i in 1:10) {
    nm <- pretty_pw(w$pathway_id[i])
    if (grepl(LIPO_RX, w$pathway_id[i])) nm <- sprintf("\\textbf{%s}", nm)
    L <- c(L, sprintf("%d & %s & %d & %.1f & %.0f & %.1f\\%% & %d \\\\",
                      i, nm, w$n_genes[i], -w$log10_p_omnibus[i], w$z_std[i],
                      100 * w$size_pct[i], wo$rank[match(w$pathway_id[i], wo$pathway_id)]))
  }
}
L <- c(L, "\\hline", "\\end{tabular}", "\\end{table*}")
writeLines(L, file.path(AUT, "tab_lipid_pathways.tex"))
cat("tab_lipid_pathways.tex\n")

## ======================================== Table: top ten candidate genes ====
L <- c("\\begin{table*}[!t]", "\\centering",
       paste0("\\caption{Ten highest-scoring candidate genes for each trait. Rank ",
              "columns give the gene's position within each evidence layer (1 is the ",
              "strongest): best variant within the MAGMA gene annotation (GWAS), MAGMA ",
              "gene-level $p$-value (MAGMA), and best pathway by $z_{\\mathrm{std}}$ ",
              "among pathway-annotated genes (Pathway). Pre-specified established genes in ",
              "bold.}"),
       "\\label{tab:lipid_candidates}", "\\small",
       "\\begin{tabular}{rlrrrrl}", "\\hline",
       paste0("\\textbf{Rank} & \\textbf{Gene} & \\textbf{Score} & \\textbf{GWAS} & ",
              "\\textbf{MAGMA} & \\textbf{Pathway} & \\textbf{Best pathway} \\\\"),
       "\\hline")
fr <- function(x) ifelse(is.na(x), "--", formatC(round(x), format = "d", big.mark = "{,}"))
for (t in names(TRAITS)) {
  g <- G[[t]][order(G[[t]]$rank), ][1:10, ]
  L <- c(L, sprintf("\\multicolumn{7}{l}{\\textit{%s}} \\\\", TRAITS[[t]]))
  for (i in seq_len(nrow(g))) {
    s <- sprintf("\\textit{%s}", tex_esc(g$symbol[i]))
    if (g$symbol[i] %in% CONTROLS[[t]]) s <- sprintf("\\textbf{%s}", s)
    bp <- if (is.na(g$best_pathway[i])) "--" else pretty_pw(g$best_pathway[i], 40)
    L <- c(L, sprintf("%d & %s & %.2f & %s & %s & %s & %s \\\\", g$rank[i], s,
                      g$Score[i], fr(g$rank_gwas[i]), fr(g$rank_magma[i]),
                      fr(g$rank_pathway[i]), bp))
  }
}
L <- c(L, "\\hline", "\\end{tabular}", "\\end{table*}")
writeLines(L, file.path(AUT, "tab_lipid_candidates.tex"))
cat("tab_lipid_candidates.tex\n")

## ============================================ Table: known-gene recovery ====
meth4 <- c(GWAS = "rank_gwas", MAGMA = "rank_magma", Pathway = "rank_pathway",
           "CATFISH score" = "rank")
L <- c("\\begin{table}[!t]", "\\centering",
       paste0("\\caption{Recovery of the pre-specified established genes by each ",
              "evidence layer alone and by the combined CATFISH score. Ranks are among ",
              "all genes in the layer (1 is the strongest). Top-$k$ columns count ",
              "established genes ranked within the first $k$.}"),
       "\\label{tab:lipid_recovery}", "\\small",
       "\\begin{tabular}{llrrrr}", "\\hline",
       paste0("\\textbf{Trait} & \\textbf{Ranking} & \\textbf{Median rank} & ",
              "\\textbf{Top 10} & \\textbf{Top 20} & \\textbf{Top 50} \\\\"),
       "\\hline")
for (t in names(TRAITS)) {
  g <- G[[t]]; ctrl <- CONTROLS[[t]][CONTROLS[[t]] %in% g$symbol]
  def(sprintf("nLipCtrl%s", toupper(t)), length(ctrl))
  for (mm in names(meth4)) {
    r <- rk_of(g, ctrl, meth4[[mm]]); r <- r[!is.na(r)]
    med <- stats::median(r)
    L <- c(L, sprintf("%s & %s & %s & %d & %d & %d \\\\",
                      if (mm == "GWAS") TRAITS[[t]] else "", mm,
                      trimws(formatC(med, format = "fg", digits = 4)),
                      sum(r <= 10), sum(r <= 20), sum(r <= 50)))
    key <- c(GWAS = "Gwas", MAGMA = "Magma", Pathway = "Path",
             "CATFISH score" = "Score")[[mm]]
    def(sprintf("nLipRec%s%sMed", toupper(t), key), trimws(formatC(med, format = "fg", digits = 4)))
    def(sprintf("nLipRec%s%sTopTen", toupper(t), key), sum(r <= 10))
    def(sprintf("nLipRec%s%sTopTwenty", toupper(t), key), sum(r <= 20))
    def(sprintf("nLipRec%s%sTopFifty", toupper(t), key), sum(r <= 50))
  }
  if (t == "ldl") L <- c(L, "\\hline")
}
L <- c(L, "\\hline", "\\end{tabular}", "\\end{table}")
writeLines(L, file.path(AUT, "tab_lipid_recovery.tex"))
cat("tab_lipid_recovery.tex\n")

## ==================================================== numbers for text ====
sp <- function(a, b) stats::cor(a, b, method = "spearman")
for (t in names(TRAITS)) {
  U <- toupper(t); w <- W[[t]]; g <- G[[t]]
  lipo <- grepl(LIPO_RX, w$pathway_id)
  def(sprintf("nLipSig%s", U), fmt_int(sum(w$p_fwer_std <= 0.05, na.rm = TRUE)))
  def(sprintf("nLipFloor%s", U), fmt_int(sum(w$p_fwer_std <= 1 / (B + 1) + 1e-12, na.rm = TRUE)))
  def(sprintf("nLipTopRank%s", U), w$rank[match(TOPPW, w$pathway_id)])
  def(sprintf("nLipPparRank%s", U), w$rank[match(PPARA, w$pathway_id)])
  def(sprintf("nLipSizeRho%s", U), sprintf("%.2f", sp(log(w$n_genes), w$z_std)))
  def(sprintf("nLipTopTwentyMedSize%s", U), stats::median(w$n_genes[1:20]))
  def(sprintf("nLipFamPctMed%s", U), sprintf("%.1f", 100 * stats::median(w$size_pct[lipo])))
  def(sprintf("nLipFamBestFive%s", U), sum(w$size_pct[lipo] < 0.05))
  def(sprintf("nLipFamTopFifty%s", U), sum(w$rank[lipo] <= 50))
  def(sprintf("nLipAllPctMedLipoFree%s", U),
      sprintf("%.1f", 100 * stats::median(w$size_pct[!lipo])))
  ## gene level
  gm <- M[[t]]
  def(sprintf("nLipGeneGW%s", U), sprintf("%.1f", 100 * mean(gm$P < 0.05 / nrow(gm))))
  def(sprintf("nLipTopTenCtrl%s", U),
      sum(g$symbol[order(g$rank)][1:10] %in% CONTROLS[[t]]))
}
def("nLipFamN", sum(grepl(LIPO_RX, W$ldl$pathway_id)))
def("nLipMedSizeAll", stats::median(W$ldl$n_genes))
## LDL vs HDL agreement
idx <- match(W$ldl$pathway_id, W$hdl$pathway_id)
def("nLipCrossRho", sprintf("%.3f", sp(W$ldl$rank, W$hdl$rank[idx])))
def("nLipCrossTopHundred", length(intersect(W$ldl$pathway_id[1:100], W$hdl$pathway_id[1:100])))
def("nLipGeneRho", sprintf("%.3f", {
  a <- M$ldl; b <- M$hdl; k <- intersect(a$GENE, b$GENE)
  sp(a$ZSTAT[match(k, a$GENE)], b$ZSTAT[match(k, b$GENE)])
}))

writeLines(c("%% AUTO-GENERATED by 21_lipid_results.R -- do not edit by hand",
             unlist(macros)), file.path(AUT, "lipids.tex"))
writeLines(c("%% Placeholders for the macros written by 21_lipid_results.R",
             sprintf("\\newcommand{\\%s}{\\textbf{??}}", unlist(macro_names))),
           file.path(AUT, "lipids_stub.tex"))
cat(sprintf("wrote %d macros to %s/lipids.tex\n", length(macros), AUT))
