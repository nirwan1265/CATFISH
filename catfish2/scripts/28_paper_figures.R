#!/usr/bin/env Rscript
## =============================================================================
## 28_paper_figures.R -- every manuscript figure, redrawn in one style.
##
## Same data as before, same file names, so catfish_wy.tex is unchanged. No
## figure titles (journal style); panels are tagged A, B, ... The
## figures are drawn with ggplot2 and plot_theme from scripts/all_figs.R:
##   * one colour needed            -> black
##   * two colours needed           -> black and grey (black-and-white)
##   * more than two, or continuous -> viridis (rank heatmaps stay greyscale)
## Nothing is recomputed: the inputs are the result files written by 01-04,
## 06, 08, 09, 17, 18, 20, 22 and 24. (Fig0_archetypes.png is the figure from
## the first CATFISH paper and is not redrawn here; Fig17 and the replication
## supplementary figure are drawn by 27_replication_figures.R.)
##
## The old base-graphics figure code in 05, 21, 22 and 24 still writes the
## same file names; run this script last, or those files are overwritten
## with the old style.
##
## Usage (from catfish2/scripts):  Rscript 28_paper_figures.R
## Env overrides: CATFISH_BASE, CATFISH_DATA, RES_DIR, FIG_DIR, HUMAN_RES,
##                SORGHUM_MAGMA (sorghum MAGMA root, for the gene-p QQ plot)
## =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(data.table)
})

BASE    <- Sys.getenv("CATFISH_BASE", path.expand("~/Documents/Github/MAGCAT"))
DATA_D  <- Sys.getenv("CATFISH_DATA", path.expand("~/Documents/Github/MAGCAT_data"))
RES     <- Sys.getenv("RES_DIR", file.path(BASE, "catfish2/results"))
FIG     <- Sys.getenv("FIG_DIR", file.path(BASE, "catfish2/figures"))
HUM     <- Sys.getenv("HUMAN_RES", file.path(BASE, "results/human_LDL"))
SORG    <- Sys.getenv("SORGHUM_MAGMA", path.expand("~/Documents/Research/results/CATFISH/MAGMA"))
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
has  <- function(f) file.exists(file.path(RES, f))
done <- function(f) cat("  ", f, "\n")

## ---- plot_theme: verbatim from scripts/all_figs.R ---------------------------
plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.title     = element_text(
      size   = 24,
      face   = "bold",
      hjust  = 0.5,
      margin = margin(b = 10)
    ),
    axis.title.x   = element_text(
      size = 24,
      face = "bold"
    ),
    axis.title.y   = element_text(
      size = 24,
      face = "bold"
    ),
    axis.text.x    = element_text(
      size = 24,
      color = "black"
    ),
    axis.text.y    = element_text(
      size = 24,
      color = "black"
    ),
    axis.line      = element_line(color = "black"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.25),
    legend.position = "top",
    legend.title    = element_blank(),
    legend.text     = element_text(size = 24)
  )
strip_t <- theme(strip.text = element_text(size = 22, face = "bold"),
                 legend.key.width = unit(3.5, "lines"))
## no figure titles (journal style): panel identity lives in tags and the LaTeX caption
NOTITLE <- theme(plot.title = element_blank(), plot.caption = element_blank())
save <- function(p, f, w, h) {
  p <- if (inherits(p, "patchwork"))
         p + plot_annotation(tag_levels = "A") &
           theme(plot.tag = element_text(size = 30, face = "bold")) & NOTITLE
       else p + NOTITLE
  ggsave(file.path(FIG, f), p, width = w, height = h, dpi = 300, bg = "white", limitsize = FALSE)
  done(f)
}
BW   <- c("black", "grey60")                     # two-colour scale
vir  <- function(...) scale_colour_viridis_d(end = 0.85, ...)
virf <- function(...) scale_fill_viridis_d(end = 0.85, ...)
NOM  <- geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey30")
ONE  <- geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30")
LDL  <- c(independent = "Independent", moderate = "Moderate LD", strong = "Strong LD")
ld_f <- function(x) factor(LDL[x], levels = LDL)

cat("figures ->", FIG, "\n")

## ======================================================= Fig 1 and Fig 2 ===
if (has("01_null_calibration.csv")) {
  d <- fread(file.path(RES, "01_null_calibration.csv"))[alpha == 0.05]
  M3 <- c(t1_analytic = "Analytic", t1_mvn = "Per-pathway calibrated",
          t1_wy = "Westfall-Young")
  l <- melt(d, id.vars = c("rho", "ld", "m"), measure.vars = names(M3),
            variable.name = "method", value.name = "rate")
  l[, method := factor(M3[as.character(method)], levels = M3)]
  l[, m := factor(sprintf("m = %d", m), levels = sprintf("m = %d", sort(unique(d$m))))]
  p <- ggplot(l, aes(rho, rate, colour = method, shape = method)) + NOM +
    geom_line(linewidth = 1.2) + geom_point(size = 4) +
    facet_wrap(~ m, nrow = 1) + vir() +
    labs(x = expression(bold("Gene-gene correlation " * rho)),
         y = expression(bold("Type I error at " * alpha * " = 0.05"))) +
    plot_theme + strip_t
  save(p, "Fig1_type1_error.png", 20, 7.5)

  l2 <- melt(d, id.vars = c("ld", "m"), measure.vars = c("lambda_analytic", "lambda_mvn"),
             variable.name = "method", value.name = "lambda")
  l2[, method := factor(c(lambda_analytic = "Analytic", lambda_mvn = "Per-pathway calibrated")[as.character(method)],
                        levels = c("Analytic", "Per-pathway calibrated"))]
  l2[, grp := factor(sprintf("%s\nm = %d", LDL[ld], m),
                     levels = unique(sprintf("%s\nm = %d", LDL[d$ld], d$m)))]
  p <- ggplot(l2, aes(grp, lambda, fill = method)) + ONE +
    geom_col(position = position_dodge(width = 0.8), width = 0.75, colour = "black") +
    scale_fill_manual(values = c("black", "white")) +
    labs(x = NULL, y = expression(bold(lambda[GC]))) +
    plot_theme + theme(axis.text.x = element_text(size = 18, color = "black"))
  save(p, "Fig2_lambda.png", 18, 8)
}

## ====================================================== Fig 3 and Fig 3b ===
if (has("02_breakage_stress.csv")) {
  d <- fread(file.path(RES, "02_breakage_stress.csv"))
  M4 <- c(t1_analytic = "Analytic (marginal)", t1_mvn = "MVN (marginal)",
          fwer_wy = "WY (family-wise)", fwer_bonferroni = "Bonferroni (family-wise)")
  l <- melt(d, id.vars = c("ld", "m", "direction", "k_broken"), measure.vars = names(M4),
            variable.name = "method", value.name = "rate")
  l[, `:=`(method = factor(M4[as.character(method)], levels = M4), ld = ld_f(ld),
           direction = factor(tools::toTitleCase(direction), levels = c("Anticonservative", "Conservative")),
           m = factor(sprintf("m = %d", m)))]
  p <- ggplot(l, aes(k_broken, rate, colour = method, shape = method, linetype = m)) + NOM +
    geom_line(linewidth = 1) + geom_point(size = 3.2) +
    facet_grid(direction ~ ld) + vir() +
    scale_linetype_manual(values = c("solid", "dashed")) +
    labs(x = "Corrupted components (k)", y = "Error rate at 0.05") +
    plot_theme + strip_t + theme(legend.box = "vertical")
  save(p, "Fig3_breakage.png", 22, 13)

  l <- melt(d, id.vars = c("ld", "m", "direction", "k_broken"),
            measure.vars = c("lambda_analytic", "lambda_mvn"), variable.name = "method", value.name = "lambda")
  l[, `:=`(method = factor(c(lambda_analytic = "Analytic", lambda_mvn = "Per-pathway calibrated")[as.character(method)],
                           levels = c("Analytic", "Per-pathway calibrated")), ld = ld_f(ld),
           direction = factor(tools::toTitleCase(direction), levels = c("Anticonservative", "Conservative")),
           m = factor(sprintf("m = %d", m)))]
  p <- ggplot(l, aes(k_broken, lambda, colour = method, shape = method, linetype = m)) + ONE +
    geom_line(linewidth = 1) + geom_point(size = 3.2) +
    facet_grid(direction ~ ld) + scale_colour_manual(values = BW) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    labs(x = "Corrupted components (k)", y = expression(bold(lambda[GC]))) +
    plot_theme + strip_t + theme(legend.box = "vertical")
  save(p, "Fig3b_breakage_lambda.png", 22, 13)
}

## ================================================================ Fig 4 ===
if (has("03_adaptive_loo.csv")) {
  d <- fread(file.path(RES, "03_adaptive_loo.csv"))
  d <- d[m == max(m)]
  P4 <- c(t1_mvn = "Per-pathway calibration\n(type I error)",
          fwer_wy = "Westfall-Young\n(family-wise error)",
          fwer_bonferroni = "Bonferroni\n(family-wise error)",
          lambda_mvn = "Calibrated null\n(genomic-control lambda)")
  l <- melt(d, id.vars = c("ld", "dropped"), measure.vars = intersect(names(P4), names(d)),
            variable.name = "panel", value.name = "v")
  l[, `:=`(panel = factor(P4[as.character(panel)], levels = P4), ld = ld_f(ld),
           dropped = factor(sub("<none>", "None", dropped),
                            levels = sub("<none>", "None", unique(d$dropped))))]
  ref <- data.table(panel = factor(P4, levels = P4), y = c(0.05, 0.05, 0.05, 1))
  p <- ggplot(l, aes(dropped, v, fill = ld)) +
    geom_hline(data = ref, aes(yintercept = y), linetype = "dashed", colour = "grey30") +
    geom_col(position = position_dodge(width = 0.85), width = 0.8) +
    facet_wrap(~ panel, nrow = 1, scales = "free_y") + virf() +
    labs(x = "Component removed", y = NULL) +
    plot_theme + strip_t +
    theme(axis.text.x = element_text(size = 18, angle = 45, hjust = 1, color = "black"))
  save(p, "Fig4_leave_one_out.png", 26, 9)
}

## ================================================================ Fig 6 ===
if (has("08_missing_correlation.csv")) {
  d <- fread(file.path(RES, "08_missing_correlation.csv"))
  M4 <- c(t1_analytic = "Analytic (marginal)", t1_mvn = "MVN (marginal)",
          fwer_wy = "WY (family-wise)", fwer_bonferroni = "Bonferroni (family-wise)")
  l <- melt(d, id.vars = c("ld", "miss"), measure.vars = names(M4),
            variable.name = "method", value.name = "rate")
  l[, `:=`(method = factor(M4[as.character(method)], levels = M4), ld = ld_f(ld))]
  pA <- ggplot(l, aes(miss, rate, colour = method, shape = method)) + NOM +
    geom_line(linewidth = 1.2) + geom_point(size = 4) +
    facet_wrap(~ ld, nrow = 1) + vir() +
    scale_x_continuous(labels = scales::percent) +
    labs(x = "Gene-pair correlations dropped", y = "Error rate at 0.05") +
    plot_theme + strip_t
  dl <- copy(d)[, ld := ld_f(ld)]
  pB <- ggplot(dl, aes(miss, lambda_mvn, linetype = ld)) + ONE +
    geom_line(linewidth = 1.2, colour = "black") + geom_point(size = 4, colour = "black") +
    scale_x_continuous(labels = scales::percent) + coord_cartesian(ylim = c(0.5, 1.5)) +
    labs(x = "Gene-pair correlations dropped", y = expression(bold(lambda[GC]~"(calibrated)"))) +
    plot_theme
  p <- (pA | pB) + plot_layout(widths = c(length(unique(d$ld)), 1.2))
  save(p, "Fig6_missing_correlation.png", 28, 8.5)
}

## ================================================================ Fig 7 ===
if (has("09_component_calibration.csv")) {
  d <- fread(file.path(RES, "09_component_calibration.csv"))
  d <- d[m == max(m)]
  P7 <- c(t1_analytic = "Analytic\n(independence assumed)",
          t1_calibrated = "Calibrated against\nthe correlated null",
          lambda_analytic = "Genomic-control lambda,\nanalytic")
  l <- melt(d, id.vars = c("ld", "component"), measure.vars = names(P7),
            variable.name = "panel", value.name = "v")
  l[, `:=`(panel = factor(P7[as.character(panel)], levels = P7), ld = ld_f(ld),
           component = factor(component, levels = unique(d$component)))]
  ref <- data.table(panel = factor(P7, levels = P7), y = c(0.05, 0.05, 1))
  p <- ggplot(l, aes(component, v, fill = ld)) +
    geom_hline(data = ref, aes(yintercept = y), linetype = "dashed", colour = "grey30") +
    geom_col(position = position_dodge(width = 0.85), width = 0.8) +
    facet_wrap(~ panel, nrow = 1, scales = "free_y") + virf() +
    labs(x = NULL, y = NULL) + plot_theme + strip_t +
    theme(axis.text.x = element_text(size = 18, angle = 45, hjust = 1, color = "black"))
  save(p, "Fig7_component_calibration.png", 22, 9)
}

## ============================================================ Fig 11, 12 ===
ARCHF <- c("17_archetype_power_wy_tau00100505.csv", "17_archetype_power_wy.csv")
ARCHF <- ARCHF[vapply(ARCHF, has, logical(1))]
ARCH  <- c("SSS", "MDS", "WDE", "MSA", "SGP")
if (length(ARCHF)) {
  d <- fread(file.path(RES, ARCHF[1]))
  meth <- c("ACAT", "Fisher", "Stouffer", "minP", "TFisher", "Omnibus")
  d <- d[method %in% meth]
  d[, `:=`(method = factor(method, levels = meth),
           archetype = factor(archetype, levels = intersect(ARCH, unique(archetype))),
           m = factor(sprintf("m = %d", m), levels = sprintf("m = %d", sort(unique(m)))))]
  p <- ggplot(d, aes(effect_scale, power_wy, colour = method, shape = method, linetype = method)) +
    geom_line(linewidth = 1) + geom_point(size = 2.8) +
    facet_grid(m ~ archetype) +
    scale_colour_manual(values = c(viridisLite::viridis(5, end = 0.85), "black")) +
    scale_linetype_manual(values = c(rep("solid", 5), "dashed")) +
    scale_shape_manual(values = c(15, 16, 17, 18, 8, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Effect scale", y = "Power (family-wise, alpha = 0.05)") +
    plot_theme + strip_t + guides(colour = guide_legend(nrow = 1))
  save(p, "Fig11_archetype_power.png", 26, 15)

  REGF <- c("17_archetype_regret_wy_tau00100505.csv", "17_archetype_regret_wy.csv")
  REGF <- REGF[vapply(REGF, has, logical(1))]
  if (length(REGF)) {
    r <- fread(file.path(RES, REGF[1]))
    r[, archetype := factor(archetype, levels = intersect(ARCH, archetype))]
    p <- ggplot(r, aes(archetype)) +
      geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey30") +
      geom_col(aes(y = mean_regret, fill = "Mean regret"), width = 0.6) +
      geom_linerange(aes(ymin = mean_regret, ymax = max_regret), colour = "grey50", linewidth = 1) +
      geom_point(aes(y = max_regret, shape = "Worst cell"), size = 5, colour = "grey40") +
      scale_fill_manual(values = "black") + scale_shape_manual(values = 4) +
      labs(x = NULL, y = "Omnibus power regret") + plot_theme
    save(p, "Fig12_archetype_regret.png", 13, 8)
  }
}

## =============================================================== Fig 10 ===
f_ldl  <- file.path(HUM, "magma/ldl_eur.genes.out")
f_hdl  <- file.path(HUM, "magma_hdl/hdl_eur.genes.out")
f_wl   <- file.path(HUM, "catfish/ldl_eur_reactome_5comp_westfall_young.csv")
f_wh   <- file.path(HUM, "catfish/hdl_eur_reactome_5comp_westfall_young.csv")
f_sets <- file.path(DATA_D, "c2.cp.reactome.v2026.2.Hs.entrez.sets.annot")
f_loc  <- file.path(DATA_D, "NCBI37.3.gene.loc")
TOPPW  <- "REACTOME_PLASMA_LIPOPROTEIN_ASSEMBLY_REMODELING_AND_CLEARANCE"
LIPO_RX <- "^REACTOME_(PLASMA_LIPOPROTEIN|LDL_|HDL_|VLDL_|VLDLR_|CHYLOMICRON_)"
if (all(file.exists(c(f_ldl, f_hdl, f_wl, f_wh, f_sets, f_loc)))) {
  rd <- function(f) { x <- fread(f); setNames(x$P, as.character(x$GENE)) }
  L <- rd(f_ldl); H <- rd(f_hdl)
  g <- intersect(names(L), names(H))
  loc <- fread(f_loc, header = FALSE, fill = TRUE)
  SYM <- setNames(as.character(loc[[6]]), as.character(loc[[1]]))
  mem <- NULL
  for (ln in readLines(f_sets, warn = FALSE)) {
    f <- strsplit(trimws(ln), "[ \t]+")[[1]]
    if (length(f) && f[1] == TOPPW) { mem <- f[-1]; break }
  }
  mem <- intersect(mem, g)
  cap <- 60
  dA <- data.table(gene = g, x = pmin(-log10(pmax(L[g], 1e-320)), cap),
                   y = pmin(-log10(pmax(H[g], 1e-320)), cap))
  dA[, member := gene %in% mem]
  labA <- c("LDLR","PCSK9","LDLRAP1","APOB","APOE","ABCA1","LIPC","LPL","SCARB1","APOA1","CETP")
  dL <- dA[member == TRUE][, sym := SYM[gene]][sym %in% labA]
  rho <- cor(rank(dA$x), rank(dA$y))
  pA <- ggplot(dA, aes(x, y)) +
    geom_point(data = dA[member == FALSE], colour = "grey75", size = 0.8) +
    geom_point(data = dA[member == TRUE], colour = "black", size = 2.6) +
    geom_text(data = dL, aes(label = sym), hjust = -0.15, size = 6, fontface = "italic") +
    annotate("text", x = cap * 0.98, y = cap * 0.04, hjust = 1, size = 6.5,
             label = sprintf("Spearman = %.3f\n%d pathway members (black)", rho, length(mem))) +
    labs(x = expression(bold("LDL-C " * -log[10] * italic(p) * " (gene)")),
         y = expression(bold("HDL-C " * -log[10] * italic(p) * " (gene)")),
         title = "Gene level: the two GWAS disagree") + plot_theme

  nb <- 12
  ol <- mem[order(L[mem])][seq_len(nb)]; oh <- mem[order(H[mem])][seq_len(nb)]
  dB <- rbind(data.table(trait = "LDL-C", sym = SYM[ol], v = -(-log10(pmax(L[ol], 1e-320))), k = seq_len(nb)),
              data.table(trait = "HDL-C", sym = SYM[oh], v =  -log10(pmax(H[oh], 1e-320)), k = seq_len(nb)))
  dB[, trait := factor(trait, levels = c("LDL-C", "HDL-C"))]
  pB <- ggplot(dB, aes(v, -k, fill = trait)) +
    geom_col(orientation = "y", width = 0.8) +
    geom_text(aes(x = ifelse(v < 0, v - 1, v + 1), label = sym, hjust = ifelse(v < 0, 1, 0)),
              size = 5.5, fontface = "italic") +
    geom_vline(xintercept = 0, colour = "black") +
    scale_fill_manual(values = BW) +
    scale_x_continuous(labels = abs, expand = expansion(mult = 0.25)) +
    labs(x = expression(bold(-log[10] * italic(p))), y = NULL,
         title = sprintf("Drivers of the shared pathway (%d of %d shared)",
                         length(intersect(ol, oh)), nb)) +
    plot_theme + theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank(),
                       axis.line.y = element_blank())

  wl <- fread(f_wl); wh <- fread(f_wh)
  rl <- setNames(rank(wl$z_std), wl$pathway_id); rh <- setNames(rank(wh$z_std), wh$pathway_id)
  pp <- intersect(names(rl), names(rh))
  dC <- data.table(id = pp, x = rl[pp], y = rh[pp])
  dC[, lip := grepl("LIPO|LIPID|CHOLEST|LDL|STEROL|FATTY|SREBP|PPARALPHA|BILE", id)]
  pC <- ggplot(dC, aes(x, y)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey50") +
    geom_point(data = dC[lip == FALSE], colour = "grey75", size = 1.3) +
    geom_point(data = dC[lip == TRUE], colour = "black", size = 2.6) +
    geom_point(data = dC[id == TOPPW], shape = 1, size = 9, stroke = 1.6) +
    annotate("text", x = dC[id == TOPPW]$x * 1.8, y = dC[id == TOPPW]$y, hjust = 0, size = 6,
             label = "Plasma lipoprotein assembly,\nremodeling and clearance") +
    scale_x_log10() + scale_y_log10() +
    labs(x = "LDL-C pathway rank", y = "HDL-C pathway rank",
         title = "Pathway level: the same answer") + plot_theme
  p <- (pA | pB | pC)
  save(p, "Fig10_ldl_hdl_convergence.png", 30, 10)
}

## ====================================================== Fig 13, 14, 15 ===
TRAITS <- c(ldl = "LDL-C", hdl = "HDL-C")
CONTROLS <- list(
  ldl = c("LDLR","PCSK9","APOB","APOE","SORT1","HMGCR","NPC1L1","ABCG8",
          "LPA","CETP","LDLRAP1","MYLIP","ANGPTL3"),
  hdl = c("CETP","LIPC","LIPG","ABCA1","APOA1","APOA5","LCAT","SCARB1",
          "PLTP","LPL","ANGPTL4","CUBN"))
MAGMA_OUT <- c(ldl = f_ldl, hdl = f_hdl)
WYF <- c(ldl = f_wl, hdl = f_wh)
if (all(file.exists(c(MAGMA_OUT, WYF, file.path(RES, sprintf("06_candidate_genes_%s.csv", names(TRAITS))))))) {
  W <- lapply(WYF, function(f) { w <- fread(f)[order(z_std)]; w[, rank := .I]; w })
  G <- lapply(names(TRAITS), function(t) fread(file.path(RES, sprintf("06_candidate_genes_%s.csv", t))))
  names(G) <- names(TRAITS)
  M <- lapply(MAGMA_OUT, fread)

  ## ---- Fig 13: Manhattan ----
  CAP_MAX <- 60
  cap_for <- function(y) min(CAP_MAX, ceiling(max(y, na.rm = TRUE) * 1.05) + 1)
  manh <- function(chr, pos, y, thr, title, ylab, lab = NULL) {
    d <- data.table(chr = as.integer(chr), pos = as.numeric(pos), y = y)
    mx <- d[, .(mx = max(pos)), by = chr][order(chr)]
    mx[, off := c(0, cumsum(as.numeric(mx))[-.N])]
    d <- d[mx, on = "chr"][, x := pos + off]
    cap <- cap_for(d$y); d[, over := y > cap][, yy := pmin(y, cap)]
    p <- ggplot(d, aes(x, yy, colour = factor(chr %% 2))) +
      geom_point(data = d[over == FALSE], size = 0.5) +
      geom_point(data = d[over == TRUE], shape = 17, size = 2.2) +
      geom_hline(yintercept = thr, linetype = "dashed", colour = "grey30") +
      scale_colour_manual(values = BW, guide = "none") +
      scale_x_continuous(breaks = mx$off + mx$mx / 2, labels = mx$chr, expand = c(0.005, 0)) +
      labs(x = "Chromosome", y = ylab, title = title) + plot_theme +
      theme(axis.text.x = element_text(size = 13, color = "black"),
            panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
    if (!is.null(lab)) {
      lab <- merge(lab, mx, by = "chr")[, x := pos + off][, yy := pmin(y, cap)]
      setorder(lab, x)
      ## alternate labels left/right of their point so neighbours do not overprint
      lab[, hj := rep(c(1.08, -0.08), length.out = .N)]
      p <- p + geom_point(data = lab, aes(x, yy), inherit.aes = FALSE, size = 2.5) +
        geom_text(data = lab, aes(x, yy, label = sym, hjust = hj), inherit.aes = FALSE,
                  vjust = -0.3, size = 6.5, fontface = "italic") +
        coord_cartesian(clip = "off")
    }
    p
  }
  snp_files <- file.path(RES, sprintf("20_manhattan_snps_%s.tsv", names(TRAITS)))
  have_snp <- all(file.exists(snp_files))
  pl <- list()
  for (t in names(TRAITS)) {
    if (have_snp) {
      s <- fread(snp_files[match(t, names(TRAITS))])
      pl[[length(pl) + 1]] <- manh(s$CHR, s$BP, s$NEGLOG10P, -log10(5e-8),
                                   sprintf("%s: GWAS (variants)", TRAITS[[t]]),
                                   expression(bold(-log[10] * italic(p))))
    }
    m <- M[[t]]; y <- -log10(pmax(m$P, 1e-320)); go <- G[[t]][order(rank)]
    ii <- match(as.character(go$gene), as.character(m$GENE)); pick <- integer(0)
    for (k in seq_len(min(nrow(go), 60))) {
      j <- ii[k]; if (is.na(j)) next
      near <- any(m$CHR[ii[pick]] == m$CHR[j] & abs(m$START[ii[pick]] - m$START[j]) < 3e6)
      if (!near) pick <- c(pick, k)
      if (length(pick) == 6L) break
    }
    i <- ii[pick]
    lab <- data.table(chr = m$CHR[i], pos = (m$START[i] + m$STOP[i]) / 2, y = y[i], sym = go$symbol[pick])
    pl[[length(pl) + 1]] <- manh(m$CHR, (m$START + m$STOP) / 2, y, -log10(0.05 / nrow(m)),
                                 sprintf("%s: MAGMA (genes)", TRAITS[[t]]),
                                 expression(bold(-log[10] * italic(p)[gene])), lab)
  }
  p <- wrap_plots(pl, ncol = if (have_snp) 2 else 1)
  save(p, "Fig13_lipid_manhattan.png", if (have_snp) 30 else 16, 15)

  ## ---- Fig 14: pathway statistic against size ----
  pl <- lapply(names(TRAITS), function(t) {
    w <- copy(W[[t]])[, lipo := grepl(LIPO_RX, pathway_id)]
    k <- w[pathway_id == TOPPW]
    ggplot(w, aes(n_genes, z_std)) +
      geom_point(data = w[lipo == FALSE], colour = "grey75", size = 1.6) +
      geom_smooth(method = "loess", span = 0.3, se = FALSE, colour = "grey30",
                  linewidth = 1.2, formula = y ~ x) +
      geom_point(data = w[lipo == TRUE], colour = "black", size = 3.5) +
      geom_point(data = k, shape = 1, size = 10, stroke = 1.6) +
      annotate("text", x = k$n_genes, y = k$z_std, vjust = -1.6, size = 6,
               label = sprintf("Plasma lipoprotein assembly,\nremodeling and clearance (rank %d)", k$rank)) +
      scale_x_log10(breaks = c(10, 20, 50, 100, 200)) +
      labs(x = "Pathway size (genes)", y = expression(bold(italic(z)[std] * " (lower = stronger)")),
           title = sprintf("%s: %s Reactome pathways", TRAITS[[t]], format(nrow(w), big.mark = ","))) +
      plot_theme
  })
  p <- wrap_plots(pl, nrow = 1)
  save(p, "Fig14_lipid_pathway_size.png", 24, 9.5)

  ## ---- Fig 15: established-gene recovery ----
  meth <- c(rank_gwas = "GWAS\n(best variant)", rank_magma = "MAGMA\n(gene test)", rank = "CATFISH\nscore")
  pl <- lapply(names(TRAITS), function(t) {
    g <- G[[t]]; ctrl <- CONTROLS[[t]]
    R <- rbindlist(lapply(names(meth), function(cc)
      data.table(gene = ctrl, method = meth[[cc]], r = as.numeric(g[[cc]][match(ctrl, g$symbol)]))))
    R <- R[gene %in% R[, .(ok = all(!is.na(r))), by = gene][ok == TRUE]$gene]
    R[, method := factor(method, levels = meth)]
    med <- R[, .(r = median(r)), by = method]
    spread <- function(y, gap) { o <- order(y); v <- y[o]
      if (length(v) > 1) for (i in 2:length(v)) v[i] <- max(v[i], v[i - 1] + gap)
      out <- y; out[o] <- v; out }
    L3 <- R[method == meth[[3]]][, ly := 10^spread(log10(r), 0.085)]
    ggplot(R, aes(method, r)) +
      geom_line(aes(group = gene), colour = "grey75", linewidth = 0.8) +
      geom_point(aes(colour = method), size = 4.5) +
      geom_crossbar(data = med, aes(ymin = r, ymax = r), width = 0.4, linewidth = 0.9, colour = "black") +
      geom_segment(data = L3, aes(x = 3.08, xend = 3.3, y = r, yend = ly), colour = "grey55",
                   linewidth = 0.4) +
      geom_text(data = L3, aes(x = 3.33, y = ly, label = gene), hjust = 0, size = 5.5,
                fontface = "italic") +
      scale_y_log10() + vir(guide = "none") +
      scale_x_discrete(expand = expansion(add = c(0.5, 1.1))) +
      labs(x = NULL, y = "Rank among all genes (1 = top)",
           title = sprintf("%s: %d established genes", TRAITS[[t]], length(unique(R$gene)))) +
      plot_theme
  })
  p <- wrap_plots(pl, nrow = 1)
  save(p, "Fig15_lipid_gene_recovery.png", 22, 10)
}

## =============================================================== Fig 16 ===
if (has("22_weight_sensitivity_sim_summary.csv")) {
  s <- fread(file.path(RES, "22_weight_sensitivity_sim_summary.csv"))
  lip <- if (has("22_weight_sensitivity_lipids.csv")) fread(file.path(RES, "22_weight_sensitivity_lipids.csv")) else NULL
  W0 <- c(w_path = 0.5, w_hub = 0.3)
  hm <- function(d, col, title, fmt, lower_is_better = FALSE) {
    d <- copy(d)[, v := get(col)][, `:=`(xp = factor(w_path), yh = factor(w_hub))]
    d[, sh := if (lower_is_better) -v else v]
    rng <- range(d$sh); d[, dark := (sh - rng[1]) / max(diff(rng), 1e-12) > 0.55]
    def <- d[w_path == W0[["w_path"]] & w_hub == W0[["w_hub"]]]
    ggplot(d, aes(xp, yh)) +
      geom_tile(aes(fill = sh), colour = "white", linewidth = 0.6) +
      geom_text(aes(label = sprintf(fmt, v), colour = I(ifelse(dark, "white", "black"))), size = 5) +
      geom_tile(data = def, fill = NA, colour = "white", linewidth = 3.2) +
      geom_tile(data = def, fill = NA, colour = "black", linewidth = 1.3) +
      scale_fill_gradient(low = "white", high = "black", guide = "none") +
      labs(x = "Pathway weight", y = "Hub weight", title = title) + plot_theme +
      theme(axis.line = element_blank(), panel.grid = element_blank(),
            plot.title = element_text(size = 21, face = "bold", hjust = 0.5))
  }
  pl <- list(hm(s, "mean_top20_causal_fraction", "Simulation: causal genes\nin the top 20", "%.2f"),
             hm(s, "mean_rank_correlation", "Simulation: rank correlation\nwith the default weights", "%.3f"))
  if (!is.null(lip)) for (t in unique(lip$trait))
    pl[[length(pl) + 1]] <- hm(lip[trait == t], "median_rank",
      sprintf("%s: median rank of\n%d established genes", TRAITS[[t]], lip[trait == t]$n_controls[1]),
      "%g", lower_is_better = TRUE)
  p <- wrap_plots(pl, nrow = 1)
  save(p, "Fig16_weight_sensitivity.png", 7.2 * length(pl), 8.5)
}

## ======================================= Supplementary: sorghum gene-p QQ ===
STR <- c("Dry_tons_per_acre", "Lignin_DM_perc")
sf <- file.path(SORG, rep(STR, each = 10),
                sprintf("%s_chr%02d.multi_snp_wise.genes.out", rep(STR, each = 10), 1:10))
if (all(file.exists(sf))) {
  COLS <- c(P_SNPWISE_MEAN = "SNP-wise mean", P_SNPWISE_TOP1 = "SNP-wise top",
            P_MULTI = "Multi-model")
  q <- rbindlist(lapply(STR, function(t) {
    g <- rbindlist(lapply(sf[grepl(paste0("/", t, "/"), sf)], fread))
    g <- g[!duplicated(GENE)]
    rbindlist(lapply(names(COLS), function(k) {
      p <- sort(g[[k]][is.finite(g[[k]])])
      data.table(trait = gsub("_", " ", t), method = sprintf("%s (median p %.2f)", COLS[[k]], median(p)),
                 mk = COLS[[k]], exp = -log10(ppoints(length(p))), obs = -log10(p))
    }))
  }))
  q[, mk := factor(mk, levels = COLS)]
  q[, trait := factor(trait, levels = gsub("_", " ", STR))]
  lab <- unique(q[, .(trait, mk, method)])
  q <- q[, .SD[unique(c(which(obs > 2), seq(1, .N, length.out = 4000)))], by = .(trait, mk)]
  p <- ggplot(q, aes(exp, obs, colour = mk)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
    geom_point(size = 1.2) + facet_wrap(~ trait, nrow = 1, scales = "free_y") + vir() +
    labs(x = expression(bold("Expected " * -log[10] * italic(p))),
         y = expression(bold("Observed " * -log[10] * italic(p)))) +
    plot_theme + strip_t + guides(colour = guide_legend(override.aes = list(size = 4)))
  save(p, "Supp_sorghum_pcol_qq.png", 20, 9)
  print(lab)
}
cat("done\n")
