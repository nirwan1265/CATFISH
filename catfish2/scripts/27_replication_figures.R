#!/usr/bin/env Rscript
## =============================================================================
## 27_replication_figures.R -- cross-ancestry replication of the lipid results.
##
## Five GLGC 2021 ancestry groups (EUR discovery; AFR, EAS, SAS, HIS
## replication) x two traits (LDL-C, HDL-C), analysed with the identical
## pipeline (MAGMA with the ancestry's own 1000G panel, Westfall-Young on the
## same 1,390 Reactome pathways, the same gene score and weights).
##
##   Fig17_replication.png  (main; four panels)
##     A, B  pathway ranks for LDL-C and HDL-C: the leading TOP_PW pathways of
##           every ancestry, with their rank in all five
##     C, D  gene ranks for LDL-C and HDL-C: genes in the top TOP_G in at least
##           MIN_REC of the five ancestries, with their rank in each
##   SuppFig_replication_established.png  (supplementary)
##     established genes under the combined score vs the GWAS and MAGMA layers
##
## Style: plot_theme as in scripts/all_figs.R; black-and-white for a single
## quantity (rank), viridis where more than one colour is needed.
##
## Reads (nothing is recomputed):
##   <BASE>/results/human_LDL/catfish/<trait>_<anc>_reactome_5comp_westfall_young.csv
##   <BASE>/catfish2/results/06_candidate_genes_<trait>[_<anc>].csv  (EUR has no suffix)
## Usage (from catfish2/scripts):  Rscript 27_replication_figures.R
## Writes, all inside catfish2/:
##   figures/Fig17_replication.png, figures/SuppFig_replication_established.png
##   results/27_replication_pathway_ranks.csv, 27_replication_gene_ranks.csv,
##   results/27_replication_established.csv, 27_replication_summary.csv
## =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

BASE    <- Sys.getenv("CATFISH_BASE", path.expand("~/Documents/Github/CATFISH"))
WY_DIR  <- Sys.getenv("WY_DIR",  file.path(BASE, "results/human_LDL/catfish"))
GP_DIR  <- Sys.getenv("GP_DIR",  file.path(BASE, "catfish2/results"))
FIG_DIR <- Sys.getenv("FIG_DIR", file.path(BASE, "catfish2/figures"))
RES_DIR <- Sys.getenv("RES_DIR", file.path(BASE, "catfish2/results"))
TOP_PW  <- as.integer(Sys.getenv("TOP_PW", "5"))     # leading pathways taken from each ancestry
TOP_G   <- as.integer(Sys.getenv("TOP_G", "20"))     # gene "top" cut-off
MIN_REC <- as.integer(Sys.getenv("MIN_REC", "3"))    # gene shown if in the top TOP_G in >= this many
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RES_DIR, showWarnings = FALSE, recursive = TRUE)

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
## heatmap panels: no axis lines or grid behind the tiles
heat_theme <- plot_theme + theme(
  axis.line = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.x.top = element_text(size = 20, color = "black", lineheight = 0.9),
  plot.title = element_blank())   # no figure titles (journal style)

ANC   <- c("eur", "afr", "eas", "sas", "his")
ANC_L <- c(eur = "EUR", afr = "AFR", eas = "EAS", sas = "SAS", his = "HIS")
## GLGC 2021 sample sizes per ancestry group (Graham et al. 2021)
ANC_N <- c(eur = "1.32M", afr = "99k", eas = "147k", sas = "41k", his = "48k")
TRAITS <- c(ldl = "LDL-C", hdl = "HDL-C")
TARGETS <- list(
  ldl = c("LDLR","PCSK9","APOB","APOE","SORT1","HMGCR","NPC1L1","ABCG8",
          "LPA","CETP","LDLRAP1","MYLIP","ANGPTL3"),
  hdl = c("CETP","LIPC","LIPG","ABCA1","APOA1","APOA5","LCAT","SCARB1",
          "PLTP","LPL","ANGPTL4","CUBN"))
LIPO <- "^REACTOME_(PLASMA_LIPOPROTEIN|LDL_|HDL_|VLDL_|VLDLR_|CHYLOMICRON_)"

## rank bins: one viridis scale shared by every column (rank 1 = dark purple,
## 21-50 = yellow-green); ranks beyond 50 drop to pale grey so they recede.
BINS   <- c("1", "2-3", "4-6", "7-10", "11-20", "21-50", ">50")
GREYS  <- setNames(c(viridisLite::viridis(6, end = 0.9), "#ececec"), BINS)
binr   <- function(r) factor(BINS[cut(r, c(0, 1, 3, 6, 10, 20, 50, Inf), labels = FALSE)], levels = BINS)
inkr   <- function(r) ifelse(r <= 10, "white", "black")
pretty_pw <- function(x) {
  x <- sub("^REACTOME_", "", x); x <- gsub("_", " ", x)
  x <- paste0(substr(x, 1, 1), tolower(substring(x, 2)))
  x <- gsub("(?i)\\b(ldl|hdl|vldl|vldlr|lpl|lipc|lipa|nr1h2|nr1h3|abca1|srebp|lxr|ppara|er)\\b",
            "\\U\\1", x, perl = TRUE)
  vapply(x, function(z) paste(strwrap(z, width = 34), collapse = "\n"), "", USE.NAMES = FALSE)
}

## ------------------------------------------------------------ read inputs --
wy <- list(); gp <- list()
for (t in names(TRAITS)) for (a in ANC) {
  k <- paste(t, a, sep = "_")
  w <- fread(file.path(WY_DIR, sprintf("%s_%s_reactome_5comp_westfall_young.csv", t, a)),
             showProgress = FALSE)
  w <- w[order(z_std)]; w[, rank := .I]
  wy[[k]] <- w
  f <- file.path(GP_DIR, if (a == "eur") sprintf("06_candidate_genes_%s.csv", t)
                         else sprintf("06_candidate_genes_%s_%s.csv", t, a))
  gp[[k]] <- fread(f, showProgress = FALSE)
}

## column labels: ancestry, n, and (pathway panels) number of FWER-significant pathways
xlab_anc <- function(t, with_sig) {
  lab <- sprintf("%s\nn=%s", ANC_L, ANC_N)
  if (with_sig) lab <- paste0(lab, "\n", sapply(ANC, function(a)
    sum(wy[[paste(t, a, sep = "_")]]$p_fwer_std <= 0.05)), " sig.")
  setNames(lab, ANC)
}

heat <- function(d, rows, row_face, xlabs, title, extra = NULL) {
  d[, row := factor(row, levels = rev(rows))]
  lev <- if (is.null(extra)) ANC else c(ANC, "rec")
  d[, anc := factor(anc, levels = lev)]
  p <- ggplot(d, aes(anc, row)) +
    geom_tile(aes(fill = bin), colour = "white", linewidth = 1.2) +
    geom_text(aes(label = format(rank, big.mark = ","), colour = I(ink)), size = 6.5) +
    scale_fill_manual(values = GREYS, drop = FALSE, name = "Rank",
                      guide = guide_legend(nrow = 1, title.position = "left",
                                           override.aes = list(colour = "grey60", linewidth = 0.3))) +
    scale_x_discrete(position = "top", labels = c(xlabs, rec = sprintf("Top %d\nin", TOP_G))[lev],
                     drop = FALSE) +
    labs(title = title) + heat_theme +
    theme(axis.text.y = element_text(size = 20, color = "black", face = rev(row_face), lineheight = 0.85),
          legend.title = element_text(size = 22, face = "bold"))
  if (!is.null(extra))
    p <- p + geom_text(data = extra, aes(anc, row, label = label), size = 7, inherit.aes = FALSE)
  p
}

## ------------------------------------------------- panels A, B: pathways ----
pw_tab <- list(); pw_pan <- list()
for (t in names(TRAITS)) {
  ks <- paste(t, ANC, sep = "_")
  rows <- unique(unlist(lapply(ks, function(k) head(wy[[k]]$pathway_id, TOP_PW))))
  R <- sapply(ks, function(k) wy[[k]]$rank[match(rows, wy[[k]]$pathway_id)])
  R <- matrix(R, nrow = length(rows), dimnames = list(rows, ANC))
  R <- R[order(apply(R, 1, function(v) exp(mean(log(v))))), , drop = FALSE]   # geometric-mean rank
  pw_tab[[t]] <- data.table(trait = t, pathway_id = rownames(R), R)
  d <- data.table(row = rep(pretty_pw(rownames(R)), ncol(R)), anc = rep(ANC, each = nrow(R)),
                  rank = as.vector(R))
  d[, `:=`(bin = binr(rank), ink = inkr(rank))]
  face <- ifelse(grepl(LIPO, rownames(R)), "bold", "plain")
  pw_pan[[t]] <- heat(d, pretty_pw(rownames(R)), face, xlab_anc(t, TRUE),
                      sprintf("%s: pathways", TRAITS[[t]]))
}
fwrite(rbindlist(pw_tab, fill = TRUE), file.path(RES_DIR, "27_replication_pathway_ranks.csv"))

## ----------------------------------------------------- panels C, D: genes ----
g_tab <- list(); g_pan <- list()
for (t in names(TRAITS)) {
  ks <- paste(t, ANC, sep = "_")
  top <- lapply(ks, function(k) head(gp[[k]][order(rank)]$symbol, TOP_G))
  cnt <- table(unlist(top))
  genes <- names(cnt)[cnt >= MIN_REC]
  R <- sapply(ks, function(k) gp[[k]]$rank[match(genes, gp[[k]]$symbol)])
  R <- matrix(R, nrow = length(genes), dimnames = list(genes, ANC))
  rec <- as.integer(cnt[genes])
  o <- order(-rec, apply(R, 1, function(v) exp(mean(log(v), na.rm = TRUE))))
  R <- R[o, , drop = FALSE]; rec <- rec[o]
  g_tab[[t]] <- data.table(trait = t, gene = rownames(R), n_top = rec,
                           established = rownames(R) %in% TARGETS[[t]], R)
  d <- data.table(row = rep(rownames(R), ncol(R)), anc = rep(ANC, each = nrow(R)),
                  rank = as.vector(R))
  d[, `:=`(bin = binr(rank), ink = inkr(rank))]
  extra <- data.table(row = factor(rownames(R), levels = rev(rownames(R))),
                      anc = factor("rec", levels = c(ANC, "rec")),
                      label = sprintf("%d/5", rec))
  face <- ifelse(rownames(R) %in% TARGETS[[t]], "bold.italic", "italic")
  g_pan[[t]] <- heat(d, rownames(R), face, xlab_anc(t, FALSE),
                     sprintf("%s: genes", TRAITS[[t]]), extra)
}
fwrite(rbindlist(g_tab), file.path(RES_DIR, "27_replication_gene_ranks.csv"))

## ------------------------------------------------------ assemble Fig 17 ----
n_pw <- sapply(pw_tab, nrow); n_g <- sapply(g_tab, nrow)
fig <- (pw_pan$ldl | pw_pan$hdl) / (g_pan$ldl | g_pan$hdl) +
  plot_layout(heights = c(1.25 * max(n_pw), max(n_g)), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom", plot.tag = element_text(size = 30, face = "bold"))
suppressWarnings(ggsave(file.path(FIG_DIR, "Fig17_replication.png"), fig,
                        width = 30, height = 9 + 0.62 * (1.25 * max(n_pw) + max(n_g)),
                        dpi = 250, bg = "white", limitsize = FALSE))

## ------------------------------ Supplementary: established genes vs layers ----
est <- rbindlist(lapply(names(TRAITS), function(t) rbindlist(lapply(ANC, function(a) {
  g <- gp[[paste(t, a, sep = "_")]]; x <- g[match(TARGETS[[t]], symbol)]
  data.table(trait = t, anc = a, gene = TARGETS[[t]], score = as.numeric(x$rank),
             gwas = as.numeric(x$rank_gwas), magma = as.numeric(x$rank_magma))
}))))
fwrite(est, file.path(RES_DIR, "27_replication_established.csv"))
sm <- est[, .(score_med = median(score, na.rm = TRUE), gwas_med = median(gwas, na.rm = TRUE),
              magma_med = median(magma, na.rm = TRUE),
              score_top = sum(score <= TOP_G, na.rm = TRUE), gwas_top = sum(gwas <= TOP_G, na.rm = TRUE),
              magma_top = sum(magma <= TOP_G, na.rm = TRUE), n = .N), by = .(trait, anc)]
fwrite(sm, file.path(RES_DIR, "27_replication_summary.csv"))

METH <- c(score = "Combined score", gwas = "GWAS alone", magma = "MAGMA alone")
long <- rbind(
  sm[, .(trait, anc, method = "score", med = score_med, top = score_top)],
  sm[, .(trait, anc, method = "gwas",  med = gwas_med,  top = gwas_top)],
  sm[, .(trait, anc, method = "magma", med = magma_med, top = magma_top)])
long[, `:=`(method = factor(METH[method], levels = METH),
            anc = factor(ANC_L[anc], levels = ANC_L),
            trait = factor(TRAITS[trait], levels = TRAITS))]
vir <- scale_colour_viridis_d(end = 0.8); virf <- scale_fill_viridis_d(end = 0.8)
pA <- ggplot(long, aes(anc, med, colour = method, shape = method)) +
  geom_line(aes(group = interaction(anc, trait)), colour = "grey70", linewidth = 0.8) +
  geom_point(size = 6) + scale_y_log10() + vir +
  scale_shape_manual(values = c(16, 17, 15)) +
  facet_wrap(~ trait, nrow = 1) +
  labs(x = NULL, y = "Median rank of\nestablished genes", title = "Median rank (lower is better)") +
  plot_theme + theme(strip.text = element_text(size = 24, face = "bold"))
pB <- ggplot(long, aes(anc, top, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) + virf +
  scale_y_continuous(breaks = seq(0, 12, 2)) +
  facet_wrap(~ trait, nrow = 1) +
  labs(x = NULL, y = sprintf("Established genes\nin the top %d", TOP_G),
       title = sprintf("In the top %d (of 13 LDL-C, 12 HDL-C genes)", TOP_G)) +
  plot_theme + theme(strip.text = element_text(size = 24, face = "bold"))
supp <- (pA / pB) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 30, face = "bold"), plot.title = element_blank())
ggsave(file.path(FIG_DIR, "SuppFig_replication_established.png"), supp,
       width = 20, height = 15, dpi = 300, bg = "white")

cat("established genes per set:\n"); print(sm)
cat(sprintf("\ncombined better than GWAS alone (median rank): %d of %d sets\n",
            sum(sm$score_med < sm$gwas_med), nrow(sm)))
cat(sprintf("combined >= GWAS alone (number in top %d): %d of %d sets\n",
            TOP_G, sum(sm$score_top >= sm$gwas_top), nrow(sm)))
cat("wrote", file.path(FIG_DIR, "Fig17_replication.png"), "and",
    file.path(FIG_DIR, "SuppFig_replication_established.png"), "\n")
