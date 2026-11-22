#!/usr/bin/env Rscript
## =============================================================================
## 48_main_figures.R -- the seven main figures and the supplementary figures of
## the condensed paper (five Results sections), as multi-panel figures.
##
## Panels that already exist are taken from the scripts that draw them: each
## script is evaluated with ggsave() replaced by a function that keeps the plot
## object (nothing is written by those calls), and the panels are pulled out of
## the patchwork. The new panels (workflow, APOE locus, nine-trait benchmark,
## overlap with PoPS, sorghum leave one gene out) are drawn here. Method names
## are set to the paper's wording: Analytic; CATFISH per-pathway p; CATFISH
## family-wise p (Westfall-Young, the same null draws).
##
##   Fig1  method and signal types       Fig5  gene score and benchmark
##   Fig2  calibration and power         Fig6  sorghum null and dry biomass
##   Fig3  stress tests                  Fig7  sorghum chr4 locus and LOGO
##   Fig4  human lipid pathways          FigS1-S7 supplementary
##
## Reads  the inputs of 28_paper_figures.R, 27_replication_figures.R,
##        33_sorghum_null_figures.R, final/35_simulation_figures.R,
##        final/43_bap_chr4_locus.R, plus
##        catfish2/results/45_benchmark_topk_all.csv, 45_overlap_catfish_pops.csv,
##        catfish2/results/46_apoe_locus.csv,
##        catfish2/results/bap_traits/logo/47_bap_logo_genes.csv,
##        catfish2/figures/Fig0_archetypes.png,
##        catfish2/figures/SuppFig_SbALDH3E1_expression.png (interim, screenshot)
## Writes catfish2/figures/main/fig1.png ... fig7.png, catfish2/figures/supp/figS1.png ...
## Env: CATFISH_BASE, CATFISH_DATA, SORGHUM_MAGMA, PATHWAY_FILE (as the panel
##      scripts), SCRIPT_DIR, OUT_DIR, DPI (default 300)
##   Rscript 48_main_figures.R
## =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(data.table) })
BASE <- Sys.getenv("CATFISH_BASE", path.expand("~/Documents/Github/CATFISH"))
SCR  <- Sys.getenv("SCRIPT_DIR", file.path(BASE, "catfish2/scripts"))
RES  <- file.path(BASE, "catfish2/results")
FIGD <- file.path(BASE, "catfish2/figures")
OUT  <- Sys.getenv("OUT_DIR", FIGD)
DPI  <- as.numeric(Sys.getenv("DPI", "300"))
EXPR <- Sys.getenv("EXPR_PNG", file.path(FIGD, "SuppFig_SbALDH3E1_expression.png"))
dir.create(file.path(OUT, "main"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "supp"), recursive = TRUE, showWarnings = FALSE)
TMP <- tempfile("fig48_"); dir.create(TMP)

## ---- plot_theme (as 28_paper_figures.R) -------------------------------------------------
plot_theme <- theme_minimal(base_size = 24) +
  theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5, margin = margin(b = 10)),
        axis.title.x = element_text(size = 24, face = "bold"),
        axis.title.y = element_text(size = 24, face = "bold"),
        axis.text.x = element_text(size = 24, color = "black"),
        axis.text.y = element_text(size = 24, color = "black"),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.25),
        legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(size = 24))
strip_t <- theme(strip.text = element_text(size = 22, face = "bold"))
NOTITLE <- theme(plot.title = element_blank(), plot.caption = element_blank())
TAG <- theme(plot.tag = element_text(size = 34, face = "bold"))

## ---- 1. existing panels: evaluate each script, keep its plots ---------------------------
CAP <- list()
grab <- function(script, vars = character()) {
  old <- setNames(Sys.getenv(names(vars), unset = NA), names(vars)); if (length(vars)) do.call(Sys.setenv, as.list(vars))
  on.exit({ for (k in names(vars)) if (is.na(old[[k]])) Sys.unsetenv(k) else do.call(Sys.setenv, setNames(list(old[[k]]), k)) })
  txt <- gsub("\\bggsave\\(", ".keep(", readLines(script))
  e <- new.env(parent = globalenv())
  e$.keep <- function(filename, plot = ggplot2::last_plot(), ...) { CAP[[basename(filename)]] <<- plot; invisible(NULL) }
  r <- try(eval(parse(text = txt), envir = e), silent = TRUE)
  if (inherits(r, "try-error")) cat("!! ", basename(script), ":", conditionMessage(attr(r, "condition")), "\n")
  else cat("panels from", basename(script), "\n")
}
grab(file.path(SCR, "28_paper_figures.R"), c(FIG_DIR = TMP))
grab(file.path(SCR, "27_replication_figures.R"), c(FIG_DIR = TMP, RES_DIR = TMP))
grab(file.path(SCR, "33_sorghum_null_figures.R"), c(FIG_DIR = TMP, TEX_AUTO = TMP))
grab(file.path(SCR, "final/35_simulation_figures.R"), c(FIG_DIR = TMP, TEX_AUTO = TMP))
grab(Sys.getenv("SCRIPT_43", file.path(SCR, "final/43_bap_chr4_locus.R")), c(FIG_DIR = TMP))

## leaves of a patchwork, in drawing order
flat <- function(p) {
  if (!inherits(p, "patchwork")) return(list(p))
  subs <- p$patches$plots; self <- p; self$patches <- NULL
  class(self) <- setdiff(class(self), "patchwork")
  out <- unlist(lapply(subs, flat), recursive = FALSE)
  if (!inherits(self, "plot_filler")) out <- c(out, list(self))
  out
}
## method names in the paper's wording (data values, facet labels, axis labels)
REN <- c("Per-pathway calibrated" = "CATFISH per-pathway p",
         "Westfall-Young" = "CATFISH family-wise p",
         "MVN (marginal)" = "CATFISH per-pathway p",
         "WY (family-wise)" = "CATFISH family-wise p",
         "Analytic (marginal)" = "Analytic",
         "Bonferroni (family-wise)" = "Bonferroni family-wise p",
         "Per-pathway calibration\n(type I error)" = "Per-pathway p\n(type I error)",
         "Westfall-Young\n(family-wise error)" = "Family-wise p\n(FWER)",
         "Bonferroni\n(family-wise error)" = "Bonferroni\n(FWER)",
         "Calibrated null\n(genomic-control lambda)" = "Per-pathway p\n(lambda)",
         "Calibrated against\nthe correlated null" = "Against the\ncorrelated null",
         "Analytic\n(independence assumed)" = "Analytic\n(independence)",
         "Genomic-control lambda,\nanalytic" = "Lambda,\nanalytic")
ren_df <- function(d) {
  if (!is.data.frame(d)) return(d)
  for (j in names(d)) {
    v <- d[[j]]
    if (is.factor(v)) { l <- levels(v); hit <- l %in% names(REN); if (any(hit)) { l[hit] <- REN[l[hit]]; levels(v) <- l } }
    else if (is.character(v)) { hit <- v %in% names(REN); v[hit] <- REN[v[hit]] }
    d[[j]] <- v
  }
  d
}
ren <- function(p) {
  p$data <- ren_df(p$data)
  for (i in seq_along(p$layers)) if (is.data.frame(p$layers[[i]]$data)) p$layers[[i]]$data <- ren_df(p$layers[[i]]$data)
  p
}
P <- function(fig, i) {                       # panel i of an existing figure
  L <- flat(CAP[[fig]]); if (i > length(L)) stop(fig, " has ", length(L), " panels", call. = FALSE)
  ren(L[[i]]) + NOTITLE
}
img <- function(f, width_px = 3000) {         # an image file as a panel
  im <- magick::image_read(f); im <- magick::image_scale(im, as.character(width_px))
  wrap_elements(full = grid::rasterGrob(as.raster(im), interpolate = TRUE))
}
fr <- if (exists("free", asNamespace("patchwork"))) patchwork::free else identity   # patchwork >= 1.2
nolegend <- theme(legend.position = "none")
save_fig <- function(p, f, w, h) {
  p <- p + plot_annotation(tag_levels = "a") & TAG
  ggsave(file.path(OUT, f), p, width = w, height = h, dpi = DPI, bg = "white", limitsize = FALSE)
  cat("  wrote", f, "\n")
}

## ---- 2. new panels -----------------------------------------------------------------------
## workflow (Fig 1a)
box <- function(x0, x1, y0, y1, lab, fill = "white", face = "plain")
  list(annotate("rect", xmin = x0, xmax = x1, ymin = y0, ymax = y1, fill = fill, colour = "black", linewidth = 0.9),
       annotate("text", x = (x0 + x1) / 2, y = (y0 + y1) / 2, label = lab, size = 7.2, lineheight = 0.95,
                fontface = face, colour = if (fill == "black") "white" else "black"))
arr <- function(x0, y0, x1, y1)
  annotate("segment", x = x0, y = y0, xend = x1, yend = y1, linewidth = 1,
           arrow = grid::arrow(length = grid::unit(0.35, "cm"), type = "closed"))
pWork <- ggplot() + theme_void() +
  box(0, 14, 18, 28, "GWAS summary\nstatistics") +
  box(18, 34, 18, 28, "MAGMA gene p-values\n(gene-level genomic\ncontrol)") +
  box(38, 55, 18, 28, "Pathway gene sets\n+ gene-gene\ncorrelation (LD)") +
  box(59, 78, 18, 28, "Five component tests\nACAT, Fisher, Stouffer,\nminP, soft TFisher") +
  box(82, 100, 18, 28, "CATFISH omnibus\n(ACAT of the five)", fill = "grey85", face = "bold") +
  box(38, 55, 2, 12, "Correlated multivariate\nnormal null: B draws,\nshared by all pathways") +
  box(59, 78, 2, 12, "Per-pathway p", fill = "black", face = "bold") +
  box(82, 100, 2, 12, "Family-wise p\n(Westfall-Young,\nsame draws)", fill = "black", face = "bold") +
  box(0, 34, 2, 12, "Gene score: GWAS + MAGMA\n+ 0.5 pathway + 0.3 hub") +
  arr(14, 23, 18, 23) + arr(34, 23, 38, 23) + arr(55, 23, 59, 23) + arr(78, 23, 82, 23) +
  arr(26, 18, 26, 12) + arr(46.5, 18, 46.5, 12) + arr(55, 7, 59, 7) + arr(78, 7, 82, 7) +
  arr(91, 18, 68.5, 12) +
  annotate("curve", x = 68.5, y = 2, xend = 17, yend = 2, curvature = -0.18, linewidth = 1,
           arrow = grid::arrow(length = grid::unit(0.35, "cm"), type = "closed")) +
  annotate("text", x = 43, y = -0.3, label = "pathway evidence", size = 6.5, fontface = "italic", vjust = 1) +
  coord_cartesian(xlim = c(0, 100), ylim = c(-3, 30), expand = FALSE, clip = "off")

## signal types (Fig 1b): one 50-gene pathway per archetype, the signal settings of the
## known-truth simulation (hpc/simulation_check_catfish/01_make_scenarios.R), gene z drawn
## at fixed normal quantiles (no noise) so the profile is the expected one
ARCH <- list(SSS = list(c(0.15, 3.00)), MDS = list(c(0.60, 1.10)), WDE = list(c(0.90, 0.55)),
             MSA = list(c(0.10, 1.80), c(0.60, 1.10)), SGP = list(c(0.05, 4.00)))
arch_lab <- c(SSS = "SSS: 15% at 3.0", MDS = "MDS: 60% at 1.1", WDE = "WDE: 90% at 0.55",
              MSA = "MSA: 10% at 1.8\n+ 60% at 1.1", SGP = "SGP: 5% at 4.0")
ad <- rbindlist(lapply(names(ARCH), function(a) {
  m <- 50; z <- numeric(0); used <- 0
  for (s in ARCH[[a]]) { k <- round(s[1] * m); z <- c(z, s[2] + qnorm(ppoints(k))); used <- used + k }
  z <- c(z, qnorm(ppoints(m - used)))
  lp <- sort(-log10(2 * pnorm(-abs(z))), decreasing = TRUE)
  data.table(arch = arch_lab[[a]], r = seq_len(m) / m, lp = lp, null = sort(-log10(ppoints(m)), decreasing = TRUE))
}))
ad[, arch := factor(arch, levels = arch_lab)]
pArch <- ggplot(ad, aes(r)) +
  geom_line(aes(y = null, linetype = "No signal"), colour = "grey55", linewidth = 1) +
  geom_line(aes(y = lp, linetype = "Pathway"), linewidth = 1.2) + geom_point(aes(y = lp), size = 1.8) +
  facet_wrap(~ arch, nrow = 1) + scale_linetype_manual(values = c(Pathway = "solid", "No signal" = "dashed")) +
  scale_x_continuous(labels = scales::percent, breaks = c(0, 0.5, 1)) +
  labs(x = "Gene rank within the pathway", y = expression(bold(-log[10] * italic(p)))) +
  plot_theme + strip_t + theme(axis.text.x = element_text(size = 18), legend.key.width = unit(3.5, "lines"))

## APOE locus (Fig 5): rank by the GWAS layer and by the gene score, 10 lipid sets
ap <- fread(file.path(RES, "46_apoe_locus.csv"))
ap[, set := paste(toupper(trait), toupper(anc))]
apl <- melt(ap, id.vars = c("set", "gene"), measure.vars = c("rank_gwas", "rank"), variable.name = "layer", value.name = "r")
apl[, layer := factor(c(rank_gwas = "GWAS\n(best variant)", rank = "CATFISH\nscore")[as.character(layer)],
                      levels = c("GWAS\n(best variant)", "CATFISH\nscore"))]
apl[, gene := factor(gene, levels = c("APOE", "APOC1", "TOMM40", "PVRL2"))]
apm <- apl[, .(r = median(r)), by = .(gene, layer)]
pApoe <- ggplot(apl, aes(layer, r, group = set)) +
  geom_line(colour = "grey65", linewidth = 0.8) + geom_point(colour = "grey40", size = 3) +
  geom_crossbar(data = apm, aes(layer, r, ymin = r, ymax = r), inherit.aes = FALSE, width = 0.45, linewidth = 1.2) +
  facet_wrap(~ gene, nrow = 1) + scale_y_log10(breaks = c(1, 2, 5, 10, 20, 50)) +
  labs(x = NULL, y = "Rank among all genes") + plot_theme + strip_t +
  theme(strip.text = element_text(size = 22, face = "bold.italic"), axis.text.x = element_text(size = 18))

## nine-trait benchmark (Fig 5): established genes in the top 20, extended MHC excluded
TRL <- c(ldl = "LDL-C", hdl = "HDL-C", cad = "CAD", brca = "Breast\ncancer", af = "AF",
         cd = "Crohn's", uc = "UC", sle = "SLE", asthma = "Asthma")
bk <- fread(file.path(RES, "45_benchmark_topk_all.csv"))[universe == "no_MHC"]
bk[, trait := factor(TRL[trait], levels = TRL)]
MB <- c("GWAS+MAGMA" = "GWAS + MAGMA", "PoPS" = "PoPS", "CATFISH" = "CATFISH", "CATFISH+PoPS" = "CATFISH + PoPS")
bk4 <- bk[method %in% names(MB)][, method := factor(MB[method], levels = MB)]
pBench <- ggplot(bk4, aes(trait, top20, fill = method)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  scale_fill_viridis_d(end = 0.85) +
  labs(x = NULL, y = "Established genes\nin the top 20") + plot_theme
## pathway layer gain per trait: CATFISH minus GWAS + MAGMA, at top 20, 50, 100
gw <- dcast(bk[method %in% c("CATFISH", "GWAS+MAGMA")], trait ~ method, value.var = c("top20", "top50", "top100"))
gl <- rbindlist(lapply(c(20, 50, 100), function(k) data.table(trait = gw$trait, k = sprintf("Top %d", k),
  gain = gw[[sprintf("top%d_CATFISH", k)]] - gw[[sprintf("top%d_GWAS+MAGMA", k)]])))
gl[, k := factor(k, levels = c("Top 20", "Top 50", "Top 100"))]
pGain <- ggplot(gl, aes(trait, gain, fill = k)) + geom_hline(yintercept = 0, colour = "grey30") +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) + scale_fill_viridis_d(end = 0.85) +
  labs(x = NULL, y = "Established genes gained\nwith the pathway layer") + plot_theme
## overlap of established genes found by CATFISH and by PoPS (top 50)
ov <- fread(file.path(RES, "45_overlap_catfish_pops.csv"))[k == 50]
ovl <- melt(ov, id.vars = "trait", measure.vars = c("both", "only_CATFISH", "only_PoPS"), variable.name = "who", value.name = "n")
ovl[, `:=`(trait = factor(TRL[trait], levels = TRL),
           who = factor(c(both = "Both", only_CATFISH = "CATFISH only", only_PoPS = "PoPS only")[as.character(who)],
                        levels = c("PoPS only", "CATFISH only", "Both")))]
pOver <- ggplot(ovl, aes(trait, n, fill = who)) + geom_col(width = 0.75) +
  scale_fill_viridis_d(end = 0.85, direction = -1, breaks = c("Both", "CATFISH only", "PoPS only")) +
  labs(x = NULL, y = "Established genes\nin the top 50") + plot_theme

## sorghum leave one gene out (Fig 7)
lg <- fread(file.path(RES, "bap_traits/logo/47_bap_logo_genes.csv"))
TRS <- c(Brix_maturity = "Brix", NDF_pct_DM = "NDF", NFC_pct_DM = "NFC")
FOC <- c(SORBI_3004G300800 = "SbALDH3E1", SORBI_3004G300900 = "TRXR")
lgl <- melt(lg, id.vars = c("trait", "gene"), measure.vars = c("rank_gwas_magma", "rank_catfish", "rank_catfish_logo"),
            variable.name = "layer", value.name = "r")
LL <- c(rank_gwas_magma = "GWAS +\nMAGMA", rank_catfish = "CATFISH", rank_catfish_logo = "CATFISH,\ngene left out")
lgl[, `:=`(layer = factor(LL[as.character(layer)], levels = LL), trait = factor(TRS[trait], levels = TRS),
           foc = ifelse(gene %in% names(FOC), FOC[gene], NA_character_))]
lab <- lgl[!is.na(foc) & layer == LL[3]]
pLogo <- ggplot(lgl[is.na(foc)], aes(layer, r, group = gene)) +
  geom_line(colour = "grey75", linewidth = 0.6) + geom_point(colour = "grey60", size = 2) +
  geom_line(data = lgl[!is.na(foc)], aes(linetype = foc), colour = "black", linewidth = 1.4) +
  geom_point(data = lgl[!is.na(foc)], aes(shape = foc), colour = "black", size = 4.5) +
  geom_text(data = lab, aes(label = r), nudge_x = 0.22, size = 7, hjust = 0) +
  geom_hline(yintercept = 20, linetype = "dashed", colour = "grey30") +
  facet_wrap(~ trait, nrow = 1) + scale_y_log10(breaks = c(1, 3, 10, 30, 100, 300)) +
  scale_linetype_manual(values = c(SbALDH3E1 = "solid", TRXR = "dashed")) +
  scale_shape_manual(values = c(SbALDH3E1 = 16, TRXR = 17)) +
  labs(x = NULL, y = "Rank among all genes") + plot_theme + strip_t +
  theme(axis.text.x = element_text(size = 17), legend.text = element_text(size = 24, face = "italic"))

## ---- 3. main figures ----------------------------------------------------------------------
cat("main figures ->", file.path(OUT, "main"), "\n")
## Fig 1: workflow, signal types, component power, omnibus regret
f1 <- wrap_plots(A = pWork, B = pArch, C = P("Fig11_archetype_power.png", 1), D = P("Fig12_archetype_regret.png", 1),
                 design = "AAAA\nBBBB\nCCCD\nCCCD") + plot_layout(heights = c(0.85, 0.75, 1, 1))
save_fig(f1, "main/fig1.png", 30, 30)

## Fig 2: type I error and lambda (analytic vs CATFISH per-pathway and family-wise), simulation
f2 <- wrap_plots(A = P("Fig1_type1_error.png", 1), B = P("Fig2_lambda.png", 1),
                 C = P("Fig22_simulation_catfish_magma.png", 1),
                 D = P("Fig22_simulation_catfish_magma.png", 2), E = P("Fig22_simulation_catfish_magma.png", 3),
                 design = "AABB\nCCCC\nDDEE")
save_fig(f2, "main/fig2.png", 34, 28)

## Fig 3: stress tests
f3 <- wrap_plots(A = P("Fig3_breakage.png", 1), B = P("Fig3b_breakage_lambda.png", 1),
                 C = P("Fig4_leave_one_out.png", 1), D = P("Fig7_component_calibration.png", 1),
                 E = P("Fig6_missing_correlation.png", 1), F = P("Fig6_missing_correlation.png", 2) +
                   labs(y = expression(bold(lambda[GC] ~ "(CATFISH per-pathway " * italic(p) * ")"))),
                 design = "AAABBB\nCCCDDD\nEEEEFF") + plot_layout(heights = c(1.25, 1, 0.9))
save_fig(f3, "main/fig3.png", 34, 34)

## Fig 4: human lipid pathways
f4 <- wrap_plots(A = fr(P("Fig13_lipid_manhattan.png", 2)), B = fr(P("Fig13_lipid_manhattan.png", 4)),
                 C = P("Fig17_replication.png", 1), D = P("Fig17_replication.png", 2),
                 E = fr(P("Fig14_lipid_pathway_size.png", 1)), F = fr(P("Fig14_lipid_pathway_size.png", 2)),
                 G = fr(P("Fig10_ldl_hdl_convergence.png", 3)),
                 design = "AAABBB\nCCCDDD\nEEFFGG") + plot_layout(heights = c(0.8, 1.15, 0.9))
save_fig(f4, "main/fig4.png", 34, 34)

## Fig 5: gene score, APOE locus, nine-trait benchmark, overlap with PoPS
f5 <- wrap_plots(A = fr(P("Fig15_lipid_gene_recovery.png", 1)), B = P("Fig15_lipid_gene_recovery.png", 2),
                 C = P("Fig17_replication.png", 3), D = P("Fig17_replication.png", 4) + nolegend,
                 E = pApoe, F = pGain, G = pBench, H = pOver,
                 design = "AABBCCDD\nEEEEFFFF\nGGGGHHHH") + plot_layout(heights = c(1.2, 0.8, 0.8))
save_fig(f5, "main/fig5.png", 38, 34)

## Fig 6: sorghum permutation null and dry biomass
f6 <- wrap_plots(A = P("Fig20_sorghum_null.png", 2), B = P("Fig20_sorghum_null.png", 3),
                 C = fr(P("Fig20_sorghum_null.png", 4) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))),
                 D = P("Fig21_sorghum_dry_tons.png", 3),
                 E = P("Fig21_sorghum_dry_tons.png", 1) + labs(x = expression(bold(-log[10] ~ "per-pathway" ~ italic(p)))), F = P("Fig21_sorghum_dry_tons.png", 2),
                 design = "AAABBCCC\nDDEEEFFF") + plot_layout(heights = c(0.8, 1.2))
save_fig(f6, "main/fig6.png", 38, 24)

## Fig 7: sorghum chromosome 4 locus and leave one gene out
f7 <- wrap_plots(A = P("Fig23_bap_chr4_locus.png", 1), B = P("Fig23_bap_chr4_locus.png", 2),
                 C = P("Fig23_bap_chr4_locus.png", 3), D = P("Fig23_bap_chr4_locus.png", 4), E = pLogo,
                 design = "AABBCC\nDDDEEE")
save_fig(f7, "main/fig7.png", 34, 22)

## ---- 4. supplementary figures --------------------------------------------------------------
cat("supplementary figures ->", file.path(OUT, "supp"), "\n")
## S1: lipid SNP-level Manhattan plots and LDL-C/HDL-C convergence of genes
s1 <- wrap_plots(A = P("Fig13_lipid_manhattan.png", 1), B = P("Fig13_lipid_manhattan.png", 3),
                 C = P("Fig10_ldl_hdl_convergence.png", 1), D = P("Fig10_ldl_hdl_convergence.png", 2),
                 design = "AB\nCD")
save_fig(s1, "supp/figS1.png", 32, 22)
## S2: established lipid genes by layer, five ancestries
s2 <- wrap_plots(P("SuppFig_replication_established.png", 1), P("SuppFig_replication_established.png", 2), ncol = 1)
save_fig(s2, "supp/figS2.png", 22, 18)
## S3: gene score weight sensitivity
s3 <- wrap_plots(lapply(1:4, function(i) P("Fig16_weight_sensitivity.png", i)), ncol = 2)
save_fig(s3, "supp/figS3.png", 26, 22)
## S4: nine-trait benchmark, every ranking and universe (top 20, 50, 100; median rank)
bka <- fread(file.path(RES, "45_benchmark_topk_all.csv"))
bka[, `:=`(trait = factor(TRL[trait], levels = TRL),
           universe = factor(c(all_MHC = "All genes", no_MHC = "Extended MHC excluded",
                               annotated = "Pathway-annotated genes")[universe],
                             levels = c("Extended MHC excluded", "All genes", "Pathway-annotated genes")),
           method = factor(method, levels = c("GWAS", "MAGMA", "GWAS+MAGMA", "PoPS", "CATFISH", "CATFISH+PoPS")))]
s4a <- ggplot(melt(bka, id.vars = c("trait", "universe", "method"), measure.vars = c("top20", "top50", "top100"))[
                , variable := factor(sub("top", "Top ", variable), levels = c("Top 20", "Top 50", "Top 100"))],
              aes(trait, value, fill = method)) +
  geom_col(position = position_dodge(width = 0.88), width = 0.85) + facet_grid(variable ~ universe, scales = "free_y") +
  scale_fill_viridis_d(end = 0.9) + guides(fill = guide_legend(nrow = 1)) +
  labs(x = NULL, y = "Established genes") + plot_theme + strip_t +
  theme(axis.text.x = element_text(size = 16, angle = 45, hjust = 1))
s4b <- ggplot(bka, aes(trait, median_rank, colour = method, shape = method)) +
  geom_point(size = 4.5, position = position_dodge(width = 0.7)) + facet_wrap(~ universe, nrow = 1) +
  scale_y_log10() + scale_colour_viridis_d(end = 0.9) + scale_shape_manual(values = c(15, 16, 17, 18, 8, 4)) +
  guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
  labs(x = NULL, y = "Median rank of\nestablished genes") + plot_theme + strip_t +
  theme(axis.text.x = element_text(size = 16, angle = 45, hjust = 1))
save_fig(s4a / s4b + plot_layout(heights = c(2.4, 1)), "supp/figS4.png", 34, 30)
## S5: sorghum null, further checks
s5 <- wrap_plots(P("Fig20_sorghum_null.png", 1), P("Fig20_sorghum_null.png", 5),
                 P("Fig21_sorghum_dry_tons.png", 4), P("Supp_sorghum_pcol_qq.png", 1), ncol = 2)
save_fig(s5, "supp/figS5.png", 26, 20)
## S6: SbALDH3E1 expression (interim: Expression Atlas screenshot until the TPM table is redrawn)
if (file.exists(EXPR)) {
  ggsave(file.path(OUT, "supp/figS6.png"), img(EXPR, 1800), width = 12, height = 10, dpi = DPI, bg = "white")
  cat("  wrote supp/figS6.png\n")
}
unlink(TMP, recursive = TRUE)
