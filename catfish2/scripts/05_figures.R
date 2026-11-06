#!/usr/bin/env Rscript
## =============================================================================
## 05_figures.R -- figures from the results CSVs. Run after 01-04.
##   Rscript 05_figures.R
## =============================================================================
RES <- "../results"; FIG <- "../figures"
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
has <- function(f) file.exists(file.path(RES, f))

PAL <- c(analytic = "#C44E52", mvn = "#4C72B0", wy = "#55A868", bonferroni = "#DD8452")

## --- Fig 1: type I error vs LD, by pathway size ------------------------------
if (has("01_null_calibration.csv")) {
  d <- read.csv(file.path(RES, "01_null_calibration.csv"))
  d <- d[d$alpha == 0.05, ]
  png(file.path(FIG, "Fig1_type1_error.png"), width = 1800, height = 650, res = 170)
  op <- par(mfrow = c(1, length(unique(d$m))), mar = c(4.2, 4.2, 3, 1))
  for (mm in sort(unique(d$m))) {
    s <- d[d$m == mm, ]; s <- s[order(s$rho), ]
    plot(s$rho, s$t1_analytic, type = "b", pch = 19, col = PAL["analytic"],
         ylim = c(0, max(0.12, max(s$t1_analytic, na.rm = TRUE) * 1.1)),
         xlab = expression(paste("gene-gene correlation  ", rho)),
         ylab = "type I error at alpha = 0.05",
         main = paste0("pathway size m = ", mm))
    lines(s$rho, s$t1_mvn, type = "b", pch = 17, col = PAL["mvn"])
    lines(s$rho, s$t1_wy,  type = "b", pch = 15, col = PAL["wy"])
    abline(h = 0.05, lty = 2, col = "grey40")
    legend("topleft", bty = "n", cex = 0.85,
           legend = c("analytic", "per-pathway calibrated", "Westfall-Young", "nominal"),
           col = c(PAL[["analytic"]], PAL[["mvn"]], PAL[["wy"]], "grey40"),
           pch = c(19, 17, 15, NA), lty = c(1, 1, 1, 2))
  }
  par(op); dev.off(); cat("Fig1_type1_error.png\n")

  ## --- Fig 2: genomic-control lambda ----------------------------------------
  png(file.path(FIG, "Fig2_lambda.png"), width = 1200, height = 650, res = 170)
  op <- par(mar = c(4.2, 4.2, 2, 1))
  bp <- rbind(analytic = d$lambda_analytic, mvn = d$lambda_mvn)
  colnames(bp) <- paste0(d$ld, "\nm=", d$m)
  barplot(bp, beside = TRUE, col = c(PAL[["analytic"]], PAL[["mvn"]]),
          ylab = expression(lambda[GC]), las = 2, cex.names = 0.6,
          ylim = c(0, max(2, max(bp, na.rm = TRUE) * 1.1)))
  abline(h = 1, lty = 2, col = "grey30")
  legend("topleft", c("analytic", "per-pathway calibrated"), fill = c(PAL[["analytic"]], PAL[["mvn"]]), bty = "n")
  par(op); dev.off(); cat("Fig2_lambda.png\n")
}

## --- Fig 3: component-breakage stress test -----------------------------------
## Rows = corruption direction, columns = LD level. Four methods as colours,
## two pathway sizes as line types. The earlier layout put LD in the line type
## and left 12 lines per panel in two shades of blue, which was unreadable.
if (has("02_breakage_stress.csv")) {
  d <- read.csv(file.path(RES, "02_breakage_stress.csv"), stringsAsFactors = FALSE)
  if (!"direction" %in% names(d)) d$direction <- "conservative"
  if (!"fwer_wy" %in% names(d)) { d$fwer_wy <- NA_real_; d$fwer_bonferroni <- NA_real_ }
  dirs <- sort(unique(d$direction), decreasing = TRUE)   # anticonservative first
  lds  <- unique(d$ld); ms <- sort(unique(d$m))
  SER  <- c(analytic = "t1_analytic", mvn = "t1_mvn",
            wy = "fwer_wy", bonferroni = "fwer_bonferroni")
  ymax <- max(0.1, max(unlist(d[unname(SER)]), na.rm = TRUE) * 1.08)
  png(file.path(FIG, "Fig3_breakage.png"),
      width = 640 * length(lds), height = 430 * length(dirs) + 70, res = 165)
  op <- par(mfrow = c(length(dirs), length(lds)), mar = c(4.2, 4.2, 3, 1))
  first <- TRUE
  for (dd in dirs) for (ll in lds) {
    plot(NA, xlim = range(d$k_broken), ylim = c(0, ymax),
         xlab = "corrupted components (k)", ylab = "error rate at 0.05",
         main = sprintf("%s  |  LD: %s", dd, ll))
    for (si in seq_along(SER)) {
      for (mi in seq_along(ms)) {
        ss <- d[d$direction == dd & d$ld == ll & d$m == ms[mi], ]
        ss <- ss[order(ss$k_broken), ]
        lines(ss$k_broken, ss[[SER[[si]]]], type = "b",
              pch = c(19, 17, 15, 18)[si], lty = mi, lwd = 1.6,
              col = PAL[[names(SER)[si]]])
      }
    }
    abline(h = 0.05, lty = 3, col = "grey30")
    if (first) {
      legend("topleft", bty = "n", cex = 0.78,
             c("analytic (marginal)", "MVN (marginal)",
               "WY (family-wise)", "Bonferroni (family-wise)"),
             col = PAL[names(SER)], pch = c(19, 17, 15, 18), lwd = 1.6)
      legend("bottomright", bty = "n", cex = 0.78,
             sprintf("m = %d", ms), lty = seq_along(ms))
      first <- FALSE
    }
  }
  par(op); dev.off(); cat("Fig3_breakage.png\n")

  ## Fig 3b: genomic-control lambda for the same runs. Kept separate because
  ## lambda is a property of the MARGINAL p-value distribution and has no
  ## family-wise counterpart, so it cannot share an axis with the rates above.
  if ("lambda_mvn" %in% names(d)) {
    png(file.path(FIG, "Fig3b_breakage_lambda.png"),
        width = 640 * length(lds), height = 430 * length(dirs) + 70, res = 165)
    op <- par(mfrow = c(length(dirs), length(lds)), mar = c(4.2, 4.2, 3, 1))
    ylim <- range(c(0.8, 1.2, d$lambda_analytic, d$lambda_mvn), na.rm = TRUE)
    f <- TRUE
    for (dd in dirs) for (ll in lds) {
      plot(NA, xlim = range(d$k_broken), ylim = ylim,
           xlab = "corrupted components (k)", ylab = expression(lambda),
           main = sprintf("%s  |  LD: %s", dd, ll))
      for (mi in seq_along(ms)) {
        ss <- d[d$direction == dd & d$ld == ll & d$m == ms[mi], ]
        ss <- ss[order(ss$k_broken), ]
        lines(ss$k_broken, ss$lambda_analytic, type="b", pch=19, lty=mi,
              lwd=1.6, col=PAL[["analytic"]])
        lines(ss$k_broken, ss$lambda_mvn, type="b", pch=17, lty=mi,
              lwd=1.6, col=PAL[["mvn"]])
      }
      abline(h = 1, lty = 3, col = "grey30")
      if (f) {
        legend("topright", bty="n", cex=0.78, c("analytic","per-pathway calibrated"),
               col=PAL[c("analytic","mvn")], pch=c(19,17), lwd=1.6)
        legend("bottomleft", bty="n", cex=0.78, sprintf("m = %d", ms),
               lty=seq_along(ms)); f <- FALSE
      }
    }
    par(op); dev.off(); cat("Fig3b_breakage_lambda.png\n")
  }
}

## --- Fig 4: leave-one-out ----------------------------------------------------
## Three panels: the marginal calibration, the WY family-wise rate, and the
## Bonferroni family-wise rate. The first two share a y-scale near nominal; the
## third needs its own, which is the point of the figure.
if (has("03_adaptive_loo.csv")) {
  d <- read.csv(file.path(RES, "03_adaptive_loo.csv"), stringsAsFactors = FALSE)
  if (!"fwer_wy" %in% names(d)) d$fwer_wy <- NA_real_
  s <- d[d$m == max(d$m), ]
  s <- s[order(match(s$ld, unique(d$ld)), match(s$dropped, unique(d$dropped))), ]
  lds <- unique(s$ld); drops <- unique(s$dropped)
  COLS <- c("#4C72B0", "#55A868", "#C44E52")[seq_along(lds)]
  panels <- list(
    c(col = "t1_mvn",          lab = "marginal type I error at 0.05",
      main = "per-pathway calibration"),
    c(col = "fwer_wy",         lab = "Westfall-Young FWER at 0.05",
      main = "family-wise: Westfall-Young"),
    c(col = "fwer_bonferroni", lab = "Bonferroni FWER at 0.05",
      main = "family-wise: Bonferroni"),
    c(col = "lambda_mvn",      lab = "genomic-control lambda (calibrated)",
      main = "inflation of the marginal null"))
  png(file.path(FIG, "Fig4_leave_one_out.png"), width = 3100, height = 760, res = 170)
  op <- par(mfrow = c(1, length(panels)), mar = c(6.5, 4.4, 3, 1))
  for (pn in panels) {
    if (!pn[["col"]] %in% names(s)) next
    v <- s[[pn[["col"]]]]
    lam <- grepl("lambda", pn[["col"]])
    bp <- matrix(v, nrow = length(lds), byrow = TRUE,
                 dimnames = list(lds, drops))
    barplot(bp, beside = TRUE, las = 2, col = COLS,
            ylim = if (lam) c(0, max(1.3, max(v, na.rm = TRUE) * 1.1))
                   else c(0, max(0.09, max(v, na.rm = TRUE) * 1.15)),
            ylab = pn[["lab"]], main = pn[["main"]])
    abline(h = if (lam) 1 else 0.05, lty = 2, col = "grey30")
    legend("topright", lds, fill = COLS, bty = "n", cex = 0.85)
  }
  par(op); dev.off(); cat("Fig4_leave_one_out.png\n")
}

## --- Fig 5: ranking, WY vs MVN ----------------------------------------------
if (has("04_wy_vs_mvn.csv")) {
  d <- read.csv(file.path(RES, "04_wy_vs_mvn.csv"))
  png(file.path(FIG, "Fig5_ranking.png"), width = 1300, height = 650, res = 170)
  op <- par(mar = c(6, 4.2, 2, 1))
  bp <- rbind(analytic = d$median_rank_analytic, mvn = d$median_rank_mvn)
  colnames(bp) <- paste0(d$ld, "\nm=", d$m)
  barplot(bp, beside = TRUE, col = c(PAL[["analytic"]], PAL[["mvn"]]), las = 2,
          cex.names = 0.7, ylab = "median rank of true pathways (lower is better)")
  abline(h = d$random_expectation[1], lty = 2, col = "grey30")
  legend("topright", c("analytic", "per-pathway calibrated", "random"), bty = "n",
         fill = c(PAL[["analytic"]], PAL[["mvn"]], NA), border = c(1, 1, NA))
  par(op); dev.off(); cat("Fig5_ranking.png\n")
}
cat("\nfigures written to ", normalizePath(FIG), "\n")

## --- Fig 6: robustness to an incomplete correlation file ---------------------
if (has("08_missing_correlation.csv")) {
  d <- read.csv(file.path(RES, "08_missing_correlation.csv"), stringsAsFactors = FALSE)
  lds <- unique(d$ld)
  png(file.path(FIG, "Fig6_missing_correlation.png"),
      width = 700 * (length(lds) + 1), height = 720, res = 165)
  op <- par(mfrow = c(1, length(lds) + 1), mar = c(4.4, 4.4, 3, 1))
  ymax <- max(0.1, max(c(d$t1_analytic, d$t1_mvn, d$fwer_wy, d$fwer_bonferroni),
                       na.rm = TRUE) * 1.08)
  f <- TRUE
  for (ll in lds) {
    s <- d[d$ld == ll, ]; s <- s[order(s$miss), ]
    plot(NA, xlim = c(0, 1), ylim = c(0, ymax), xaxt = "n",
         xlab = "fraction of gene-pair correlations dropped",
         ylab = "error rate at 0.05", main = sprintf("LD: %s", ll))
    axis(1, at = s$miss, labels = sprintf("%.0f%%", 100 * s$miss))
    lines(s$miss, s$t1_analytic,     type="b", pch=19, col=PAL[["analytic"]], lwd=1.6)
    lines(s$miss, s$t1_mvn,          type="b", pch=17, col=PAL[["mvn"]], lwd=1.6)
    lines(s$miss, s$fwer_wy,         type="b", pch=15, col=PAL[["wy"]], lwd=1.6)
    lines(s$miss, s$fwer_bonferroni, type="b", pch=18, col=PAL[["bonferroni"]], lwd=1.6)
    abline(h = 0.05, lty = 3, col = "grey30")
    if (f) {
      legend("topleft", bty="n", cex=0.8,
             c("analytic (marginal)","MVN (marginal)","WY (family-wise)",
               "Bonferroni (family-wise)"),
             col = PAL[c("analytic","mvn","wy","bonferroni")],
             pch = c(19,17,15,18), lwd = 1.6)
      f <- FALSE
    }
  }
  ## lambda stays near 1 throughout -- that is the point of showing it
  plot(NA, xlim = c(0, 1), ylim = c(0.5, 1.5), xaxt = "n",
       xlab = "fraction of gene-pair correlations dropped",
       ylab = expression(lambda), main = "genomic-control inflation")
  axis(1, at = unique(d$miss), labels = sprintf("%.0f%%", 100 * unique(d$miss)))
  for (i in seq_along(lds)) {
    s <- d[d$ld == lds[i], ]; s <- s[order(s$miss), ]
    lines(s$miss, s$lambda_mvn, type="b", pch=17, lty=i, col=PAL[["mvn"]], lwd=1.6)
  }
  abline(h = 1, lty = 3, col = "grey30")
  legend("topright", bty="n", cex=0.8, sprintf("LD: %s", lds), lty=seq_along(lds))
  par(op); dev.off(); cat("Fig6_missing_correlation.png\n")
}

## --- Fig 7: per-component calibration ----------------------------------------
if (has("09_component_calibration.csv")) {
  d <- read.csv(file.path(RES, "09_component_calibration.csv"), stringsAsFactors = FALSE)
  s <- d[d$m == max(d$m), ]
  lds <- unique(s$ld); comps <- unique(s$component)
  s <- s[order(match(s$ld, lds), match(s$component, comps)), ]
  COLS <- c("#4C72B0", "#55A868", "#C44E52")[seq_along(lds)]
  panels <- list(c(col="t1_analytic",   lab="type I error at 0.05",
                   main="analytic (independence assumed)"),
                 c(col="t1_calibrated", lab="type I error at 0.05",
                   main="calibrated against the correlated null"),
                 c(col="lambda_analytic", lab=NA, main="lambda, analytic"))
  png(file.path(FIG, "Fig7_component_calibration.png"),
      width = 820 * length(panels), height = 720, res = 165)
  op <- par(mfrow = c(1, length(panels)), mar = c(6.5, 4.4, 3, 1))
  for (pn in panels) {
    v <- s[[pn[["col"]]]]
    lam <- grepl("lambda", pn[["col"]])
    bp <- matrix(v, nrow = length(lds), byrow = TRUE, dimnames = list(lds, comps))
    barplot(bp, beside = TRUE, las = 2, col = COLS,
            ylim = c(0, max(if (lam) 1.3 else 0.09, max(v, na.rm = TRUE) * 1.12)),
            ylab = if (lam) expression(lambda) else pn[["lab"]],
            main = pn[["main"]])
    abline(h = if (lam) 1 else 0.05, lty = 2, col = "grey30")
    legend("topleft", lds, fill = COLS, bty = "n", cex = 0.8)
  }
  par(op); dev.off(); cat("Fig7_component_calibration.png\n")
}

## --- Fig 8: candidate-gene layer profile (human LDL) -------------------------
## For each top-ranked gene, where it sits in each evidence layer. A gene
## carried by one layer looks very different from one supported by all three,
## and the composite score alone hides that.
if (has("06_candidate_genes_ldl.csv")) {
  d <- read.csv(file.path(RES, "06_candidate_genes_ldl.csv"),
                stringsAsFactors = FALSE)
  d <- d[order(d$rank), ]
  CTRL <- c("LDLR","PCSK9","APOB","APOE","SORT1","HMGCR","NPC1L1","ABCG8",
            "LPA","CETP","LDLRAP1","MYLIP","ANGPTL3")
  N <- 20L
  s <- d[seq_len(min(N, nrow(d))), ]
  lab <- ifelse(!is.na(s$symbol) & nzchar(s$symbol), s$symbol, s$gene)
  LCOL <- c(GWAS = "#C44E52", MAGMA = "#4C72B0", pathway = "#55A868")
  m <- cbind(GWAS = s$rank_gwas, MAGMA = s$rank_magma, pathway = s$rank_pathway)
  xmax <- max(m, na.rm = TRUE)

  png(file.path(FIG, "Fig8_candidate_layers.png"),
      width = 1500, height = 1250, res = 170)
  op <- par(mar = c(4.5, 7.5, 3, 1))
  plot(NA, xlim = c(1, xmax), ylim = c(N + 0.5, 0.5), log = "x",
       xlab = "rank within layer (log scale, 1 = strongest)", ylab = "",
       yaxt = "n", main = "LDL-C: evidence layers for the top 20 genes")
  axis(2, at = seq_len(nrow(s)), labels = lab, las = 2, cex.axis = 0.8,
       font = ifelse(lab %in% CTRL, 2, 1))
  abline(h = seq_len(nrow(s)), col = "grey92")
  for (i in seq_len(nrow(s))) {
    v <- m[i, ]; ok <- !is.na(v)
    if (sum(ok) > 1) lines(range(v[ok]), c(i, i), col = "grey70")
    points(v[ok], rep(i, sum(ok)), pch = 19, cex = 1.1, col = LCOL[names(v)[ok]])
  }
  legend("bottomright", bty = "n", cex = 0.85, pch = 19, col = LCOL,
         legend = c("locus GWAS", "MAGMA gene", "CATFISH pathway"))
  mtext("bold = pre-specified positive control", side = 3, line = 0.1,
        adj = 0, cex = 0.75, col = "grey30")
  par(op); dev.off(); cat("Fig8_candidate_layers.png\n")
}

## --- Fig 9: Reactome pathway ranking for LDL-C -------------------------------
## Every pathway by rank, with lipid-related sets highlighted. Shows both that
## the correct pathway ranks first and that significance is uninformative here.
LDLWY <- file.path(RES, "..", "..", "results", "human_LDL", "catfish",
                   "ldl_eur_reactome_5comp_westfall_young.csv")
LDLWY <- Sys.getenv("LDL_WY_CSV", LDLWY)
if (file.exists(LDLWY)) {
  d <- read.csv(LDLWY, stringsAsFactors = FALSE)
  d <- d[order(d$log10_p_omnibus), ]
  d$rank <- seq_len(nrow(d))
  y <- -d$log10_p_omnibus
  lipid <- grepl("LIPO|LIPID|CHOLEST|LDL|STEROL|FATTY|SREBP|PPARALPHA",
                 d$pathway_id)
  png(file.path(FIG, "Fig9_ldl_pathway_ranking.png"),
      width = 1700, height = 950, res = 170)
  op <- par(mar = c(4.5, 4.5, 3, 1))
  plot(d$rank, y, pch = 16, cex = 0.5, col = "grey75", log = "x",
       xlab = "pathway rank (log scale)",
       ylab = expression(-log[10]~italic(p)[omnibus]),
       main = sprintf("LDL-C: %d Reactome pathways", nrow(d)))
  points(d$rank[lipid], y[lipid], pch = 16, cex = 0.9, col = "#C44E52")
  top <- which(lipid)[1:4]
  text(d$rank[top], y[top],
       labels = tolower(gsub("_", " ", sub("^REACTOME_", "", d$pathway_id[top]))),
       pos = 4, cex = 0.6, col = "#C44E52", offset = 0.4)
  if ("p_fwer_std" %in% names(d)) {
    nsig <- sum(d$p_fwer_std < 0.05, na.rm = TRUE)
    mtext(sprintf("%d of %d significant at FWER 0.05 -- significance is not the discriminating quantity here",
                  nsig, nrow(d)), side = 3, line = 0.1, adj = 0, cex = 0.72,
          col = "grey30")
  }
  legend("topright", bty = "n", pch = 16, col = c("grey75", "#C44E52"),
         legend = c("all pathways", "lipid / lipoprotein / sterol"), cex = 0.85)
  par(op); dev.off(); cat("Fig9_ldl_pathway_ranking.png\n")
}

## --- Fig 10: LDL and HDL converge on one pathway from different genes --------
## The claim this figure has to carry: the two GWAS rank genes differently, yet
## the same pathway comes first in both. Panel A shows the gene-level
## disagreement, panel B shows which genes drive the shared pathway in each
## trait, panel C shows the pathway-level agreement at the top.
DATA_D <- Sys.getenv("CATFISH_DATA", path.expand("~/Documents/Github/MAGCAT_data"))
HUM    <- Sys.getenv("HUMAN_RES", path.expand("~/Documents/Github/MAGCAT/results/human_LDL"))
f_ldl  <- file.path(HUM, "magma/ldl_eur.genes.out")
f_hdl  <- file.path(HUM, "magma_hdl/hdl_eur.genes.out")
f_wl   <- file.path(HUM, "catfish/ldl_eur_reactome_5comp_westfall_young.csv")
f_wh   <- file.path(HUM, "catfish/hdl_eur_reactome_5comp_westfall_young.csv")
f_sets <- file.path(DATA_D, "c2.cp.reactome.v2026.2.Hs.entrez.sets.annot")
f_loc  <- file.path(DATA_D, "NCBI37.3.gene.loc")
TOPPW  <- "REACTOME_PLASMA_LIPOPROTEIN_ASSEMBLY_REMODELING_AND_CLEARANCE"

if (all(file.exists(c(f_ldl, f_hdl, f_wl, f_wh, f_sets, f_loc)))) {
  rd <- function(f) {
    d <- utils::read.table(f, header = TRUE, stringsAsFactors = FALSE)
    stats::setNames(d$P, as.character(d$GENE))
  }
  L <- rd(f_ldl); H <- rd(f_hdl)
  g <- intersect(names(L), names(H))
  lx <- -log10(pmax(L[g], 1e-320)); hy <- -log10(pmax(H[g], 1e-320))

  loc <- utils::read.table(f_loc, header = FALSE, stringsAsFactors = FALSE,
                           fill = TRUE, quote = "")
  SYM <- stats::setNames(as.character(loc[[6]]), as.character(loc[[1]]))

  mem <- NULL
  for (ln in readLines(f_sets, warn = FALSE)) {
    f <- strsplit(trimws(ln), "[ \t]+")[[1]]
    if (length(f) && f[1] == TOPPW) { mem <- f[-1]; break }
  }
  mem <- intersect(mem, g)

  png(file.path(FIG, "Fig10_ldl_hdl_convergence.png"),
      width = 2500, height = 820, res = 165)
  op <- par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3.2, 1))

  ## A -- gene level
  cap <- 60
  plot(pmin(lx, cap), pmin(hy, cap), pch = 16, cex = 0.3, col = "grey80",
       xlab = expression(LDL-C~-log[10]~italic(p)~(gene)),
       ylab = expression(HDL-C~-log[10]~italic(p)~(gene)),
       main = "A  gene level: the two GWAS disagree")
  points(pmin(lx[mem], cap), pmin(hy[mem], cap), pch = 16, cex = 0.8,
         col = "#C44E52")
  lab <- c("LDLR","PCSK9","LDLRAP1","APOB","APOE","ABCA1","LIPC","LPL",
           "SCARB1","APOA1","CETP")
  for (s in lab) {
    i <- names(SYM)[match(s, SYM)]
    if (!is.na(i) && i %in% g) text(min(lx[i], cap), min(hy[i], cap), s,
                                    cex = 0.62, pos = 4, col = "#8B2E33")
  }
  rho <- stats::cor(rank(lx), rank(hy))
  legend("topright", bty = "n", cex = 0.78,
         c(sprintf("Spearman = %.3f", rho),
           sprintf("%s members", length(mem))),
         pch = c(NA, 16), col = c(NA, "#C44E52"))
  mtext(sprintf("axes capped at %d", cap), 3, line = 0.05, adj = 0,
        cex = 0.62, col = "grey40")

  ## B -- who drives the shared pathway
  nb <- 12
  ol <- mem[order(L[mem])][seq_len(nb)]
  oh <- mem[order(H[mem])][seq_len(nb)]
  vals <- c(rev(-log10(pmax(L[ol], 1e-320))), -log10(pmax(H[oh], 1e-320)))
  cols <- c(rep("#4C72B0", nb), rep("#55A868", nb))
  nm   <- c(rev(SYM[ol]), SYM[oh])
  bp <- barplot(c(rev(-log10(pmax(L[ol], 1e-320))), rep(NA, nb)),
                horiz = TRUE, col = NA, border = NA, xlim = c(-40, 40),
                xlab = expression(-log[10]~italic(p)), yaxt = "n",
                main = "B  drivers of the shared pathway")
  barplot(c(rev(-(-log10(pmax(L[ol], 1e-320)))), rep(0, nb)), horiz = TRUE,
          col = "#4C72B0", border = NA, add = TRUE, axes = FALSE, yaxt = "n")
  barplot(c(rep(0, nb), -log10(pmax(H[oh], 1e-320))), horiz = TRUE,
          col = "#55A868", border = NA, add = TRUE, axes = FALSE, yaxt = "n")
  text(rep(0, 2 * nb), bp, nm, cex = 0.6)
  legend("topleft", bty = "n", cex = 0.78, fill = c("#4C72B0", "#55A868"),
         legend = c("LDL-C drivers", "HDL-C drivers"))
  shared <- length(intersect(ol, oh))
  mtext(sprintf("%d of %d shared between traits", shared, nb), 3, line = 0.05,
        adj = 0, cex = 0.62, col = "grey40")

  ## C -- pathway level
  wl <- utils::read.csv(f_wl, stringsAsFactors = FALSE)
  wh <- utils::read.csv(f_wh, stringsAsFactors = FALSE)
  rl <- stats::setNames(rank(wl$z_std), wl$pathway_id)
  rh <- stats::setNames(rank(wh$z_std), wh$pathway_id)
  p  <- intersect(names(rl), names(rh))
  lip <- grepl("LIPO|LIPID|CHOLEST|LDL|STEROL|FATTY|SREBP|PPARALPHA|BILE", p)
  plot(rl[p], rh[p], pch = 16, cex = 0.4, col = "grey80", log = "xy",
       xlab = "LDL-C pathway rank", ylab = "HDL-C pathway rank",
       main = "C  pathway level: the same answer")
  points(rl[p][lip], rh[p][lip], pch = 16, cex = 0.9, col = "#C44E52")
  i <- match(TOPPW, p)
  if (!is.na(i)) {
    points(rl[p][i], rh[p][i], pch = 1, cex = 2.6, lwd = 2, col = "#C44E52")
    text(rl[p][i], rh[p][i], "plasma lipoprotein\nassembly / remodeling\n/ clearance",
         pos = 4, cex = 0.62, col = "#8B2E33")
  }
  abline(0, 1, lty = 3, col = "grey50")
  legend("bottomright", bty = "n", cex = 0.78, pch = 16,
         col = c("grey80", "#C44E52"),
         legend = c("all pathways", "lipid / lipoprotein / sterol"))
  par(op); dev.off(); cat("Fig10_ldl_hdl_convergence.png\n")
}

## --- Fig 11: archetype power under the Westfall-Young decision rule ----------
ARCHF <- c("17_archetype_power_wy_tau00100505.csv", "17_archetype_power_wy.csv")
ARCHF <- ARCHF[vapply(ARCHF, has, logical(1))]
if (length(ARCHF)) {
  d <- read.csv(file.path(RES, ARCHF[1]), stringsAsFactors = FALSE)
  ## Okabe-Ito palette (colour-blind safe) plus a distinct point shape per
  ## component, so the lines separate even in greyscale.
  MPAL <- c(ACAT = "#D55E00", Fisher = "#0072B2", Stouffer = "#009E73",
            minP = "#CC79A7", TFisher = "#E69F00", Omnibus = "#000000")
  MPCH <- c(ACAT = 15, Fisher = 16, Stouffer = 17, minP = 18, TFisher = 8,
            Omnibus = 1)
  MLTY <- c(ACAT = 1, Fisher = 1, Stouffer = 1, minP = 1, TFisher = 1, Omnibus = 2)
  ARCH <- c("SSS", "MDS", "WDE", "MSA", "SGP")
  ARCH <- ARCH[ARCH %in% d$archetype]
  SZ   <- sort(unique(d$m))
  meth <- intersect(names(MPAL), unique(d$method))

  png(file.path(FIG, "Fig11_archetype_power.png"),
      width = 420 * length(ARCH), height = 360 * length(SZ) + 90, res = 170)
  op <- par(mfrow = c(length(SZ), length(ARCH)), mar = c(3.6, 3.6, 2.2, 0.6),
            mgp = c(2.2, 0.7, 0), oma = c(0, 0, 2.2, 0))
  for (mm in SZ) for (aa in ARCH) {
    s <- d[d$m == mm & d$archetype == aa, ]
    plot(NA, xlim = range(s$effect_scale), ylim = c(0, 1),
         xlab = "effect scale", ylab = "power (family-wise)",
         main = sprintf("%s   m = %d", aa, mm), cex.main = 0.95)
    for (mt in meth) {
      t <- s[s$method == mt, ]; t <- t[order(t$effect_scale), ]
      lines(t$effect_scale, t$power_wy, type = "b", pch = MPCH[[mt]], cex = 0.65,
            col = MPAL[[mt]], lty = MLTY[[mt]], lwd = if (mt == "Omnibus") 2 else 1)
    }
    if (aa == ARCH[1] && mm == SZ[1])
      legend("topleft", bty = "n", cex = 0.68, legend = meth,
             col = unname(MPAL[meth]), lty = unname(MLTY[meth]), lwd = 1.4,
             pch = unname(MPCH[meth]))
  }
  mtext("Power to detect the signal pathway at family-wise alpha = 0.05",
        outer = TRUE, cex = 0.8)
  par(op); dev.off(); cat("Fig11_archetype_power.png\n")

  ## --- Fig 12: omnibus regret against the best single component -------------
  REGF <- c("17_archetype_regret_wy_tau00100505.csv", "17_archetype_regret_wy.csv")
  REGF <- REGF[vapply(REGF, has, logical(1))]
  if (length(REGF)) {
    r <- read.csv(file.path(RES, REGF[1]), stringsAsFactors = FALSE)
    r <- r[match(ARCH, r$archetype), ]
    png(file.path(FIG, "Fig12_archetype_regret.png"),
        width = 1250, height = 620, res = 170)
    op <- par(mar = c(4.0, 4.4, 2.4, 1))
    bp <- barplot(r$mean_regret, names.arg = r$archetype, col = "#4C72B0",
                  ylab = "omnibus power regret", las = 1,
                  ylim = c(0, max(0.12, max(r$max_regret, na.rm = TRUE) * 1.05)))
    arrows(bp, r$mean_regret, bp, r$max_regret, angle = 90, length = 0.04,
           col = "grey40")
    points(bp, r$max_regret, pch = 4, cex = 0.8, col = "grey30")
    abline(h = 0.05, lty = 2, col = "#C44E52")
    legend("top", bty = "n", cex = 0.8, horiz = TRUE,
           legend = c("mean regret", "worst cell", "0.05"),
           pch = c(15, 4, NA), lty = c(NA, NA, 2),
           col = c("#4C72B0", "grey30", "#C44E52"))
    mtext("power given up by the omnibus vs the best component in each cell",
          3, line = 0.4, adj = 0, cex = 0.72, col = "grey35")
    par(op); dev.off(); cat("Fig12_archetype_regret.png\n")
  }
}
