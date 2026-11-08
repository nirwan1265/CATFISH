#!/usr/bin/env Rscript
## =============================================================================
## 33_sorghum_null_figures.R -- sorghum BAP figures and table.
##   Fig21_sorghum_dry_tons.png  main figure for dry biomass: (A) top ten
##     SorghumBicolorCyc pathways, (B) top ten genes of the combined score and
##     their rank in each evidence layer, (C, D) the permutation null: QQ of the
##     per-pathway calibrated p over all shuffles, and the family-wise error,
##     without and with GC
##   Fig20_sorghum_null.png      diagnostics of the permutation null
##   tex/auto/tab_sorghum_null.tex
## Real-data inputs: <TRAIT>_P_SNPWISE_MEAN_5comp_westfall_young.csv (23 / the
## gc_real rerun) and 25_<TRAIT>_gene_prioritize.csv (25), both in RES_DIR.
##
## Permutation null: dry biomass yield shuffled among the 327 lines 100 times, each shuffle run
## through the full pipeline (vcf2gwas LMM per chromosome with genome-wide
## kinship and 3 PCs, MAGMA snp-wise=mean, Westfall-Young B = 10000;
## catfish2/scripts/hpc/bap_perm). Three Westfall-Young versions:
##   no GC          gene p-values as MAGMA gives them          results/sorghum_bap/null_perm/
##   GC             each shuffle's gene z centred at its median   .../null_perm_gc/
##   GC + null r    GC, and gene correlations calibrated to the     .../null_perm_gc_nullcal/
##                  shuffles' own null correlations (check only)
## Also reads null_perm/perm_gene_z.tsv.gz (gene z of every shuffle),
## null_perm/nullcal_check/rmap.tsv (MAGMA r -> null r) and
## null_perm/snp_tail_by_maf.csv (SNP-level tail, counted on the server).
## Writes ../figures/Fig21_sorghum_dry_tons.png, ../figures/Fig20_sorghum_null.png
##        and ../tex/auto/tab_sorghum_null.tex
## =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
BASE <- Sys.getenv("CATFISH_BASE", path.expand("~/Documents/Github/MAGCAT"))
RES  <- Sys.getenv("RES_DIR", file.path(BASE, "catfish2/results/sorghum_bap"))
FIG  <- Sys.getenv("FIG_DIR", file.path(BASE, "catfish2/figures"))
TEX  <- Sys.getenv("TEX_AUTO", file.path(BASE, "catfish2/tex/auto"))
TRAIT <- Sys.getenv("TRAIT", "Dry_tons_per_acre")
PATHWAY_FILE <- Sys.getenv("PATHWAY_FILE",
                  file.path(BASE, "inst/extdata/pathway/sorghumbicolorcyc_pathways.20230103.SORBI"))
N_TOP <- 10L
dir.create(FIG, showWarnings = FALSE, recursive = TRUE); dir.create(TEX, showWarnings = FALSE, recursive = TRUE)

## ---- plot_theme: verbatim from scripts/all_figs.R ---------------------------
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
NOTITLE <- theme(plot.title = element_blank(), plot.caption = element_blank())

runs <- c(nogc = "no GC", gc = "GC", nullcal = "GC + null r")
dirs <- c(nogc = "null_perm", gc = "null_perm_gc", nullcal = "null_perm_gc_nullcal")
have <- file.exists(file.path(RES, dirs, "null_pathway_p.csv")); names(have) <- names(dirs)
runs <- runs[have]; dirs <- dirs[have]
cols <- setNames(viridisLite::viridis(3, end = 0.85)[seq_along(runs)], runs)
pw  <- rbindlist(lapply(names(runs), function(k)
  fread(file.path(RES, dirs[[k]], "null_pathway_p.csv"))[, run := runs[[k]]]))
per <- rbindlist(lapply(names(runs), function(k)
  fread(file.path(RES, dirs[[k]], "null_per_perm.csv"))[, run := runs[[k]]]))
pw[, run := factor(run, levels = runs)]; per[, run := factor(run, levels = runs)]

## ---- A: QQ of per-pathway calibrated p ---------------------------------------
qq <- pw[is.finite(p_calibrated), {
  o <- sort(p_calibrated); n <- length(o)
  i <- unique(round(exp(seq(0, log(n), length.out = 400))))
  .(exp = -log10(i / (n + 1)), obs = -log10(o[i]))
}, by = run]
pA <- ggplot(qq, aes(exp, obs, colour = run)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey30") +
  geom_line(linewidth = 1.3) + scale_colour_manual(values = cols) +
  labs(x = expression(bold(Expected ~ -log[10](p))), y = expression(bold(Observed ~ -log[10](p)))) +
  plot_theme

## ---- B: the genome-wide shift of each shuffled GWAS ----------------------------
## each line: shuffle, then the MAGMA .genes.out fields (GENE CHR START STOP NSNPS NPARAM N ZSTAT P)
gz <- fread(cmd = paste("gzip -dc", shQuote(file.path(RES, "null_perm", "perm_gene_z.tsv.gz")),
                        "| awk '$2 != \"GENE\" {print $1\"\\t\"$9}'"),
            header = FALSE, sep = "\t", col.names = c("perm", "z"))
shift <- gz[, .(mean_z = mean(z, na.rm = TRUE)), by = .(perm = as.integer(perm))]
b <- merge(per[run == runs[["nogc"]]], shift, by = "perm")
b[, hit := factor(ifelse(any_fwer_05, "pathway at FWER 5%", "none"), levels = c("none", "pathway at FWER 5%"))]
pB <- ggplot(b, aes(mean_z, n_cal_05, shape = hit)) +
  geom_point(size = 3.5) + scale_shape_manual(values = c(1, 16)) +
  labs(x = "Mean gene z of the shuffle", y = "Pathways with p < 0.05") + plot_theme

## ---- C: type I error by pathway size ------------------------------------------
pw[, size := cut(n_genes, c(4, 9, 19, 49, 200), labels = c("5-9", "10-19", "20-49", "50-200"))]
c5 <- pw[, .(rate = mean(p_calibrated < 0.05)), by = .(run, size)]
pC <- ggplot(c5, aes(size, rate, fill = run)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey30") +
  scale_fill_manual(values = cols) +
  labs(x = "Genes in pathway", y = "Type I error (p < 0.05)") + plot_theme

## ---- D: family-wise error -------------------------------------------------------
fw <- per[, { k <- sum(any_fwer_05); n <- .N; ci <- binom.test(k, n)$conf.int
              .(fwer = k / n, lo = ci[1], hi = ci[2], n = n) }, by = run]
pD <- ggplot(fw, aes(run, fwer, colour = run)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey30") +
  geom_pointrange(aes(ymin = lo, ymax = hi), linewidth = 1.3, size = 1.2, show.legend = FALSE) +
  scale_colour_manual(values = cols) + labs(x = NULL, y = "Family-wise error") + plot_theme

## ---- E: MAGMA gene-gene r against the null correlation --------------------------
rm_ <- fread(file.path(RES, "null_perm", "nullcal_check", "rmap.tsv"))
pE <- ggplot(rm_, aes(magma_r, null_r)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey30") +
  geom_line(linewidth = 1.3) + geom_point(size = 3) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "MAGMA gene-gene r", y = "Null correlation") + plot_theme

## ---- F: SNP-level tail by minor allele frequency --------------------------------
st <- fread(file.path(RES, "null_perm", "snp_tail_by_maf.csv"), colClasses = c(maf_bin = "character"))
st[, ratio := n_p_below_1e4 / (n_snps * 1e-4)]
st[, maf_bin := factor(maf_bin, levels = maf_bin)]
pF <- ggplot(st, aes(maf_bin, ratio)) + geom_col(fill = "black", width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30") +
  labs(x = "Minor allele frequency", y = expression(bold(Observed / expected ~ (p < 10^-4)))) + plot_theme

## the QQ plot (pA) is panel C of the main figure below; Fig20 holds the diagnostics
p <- (pB | pC | pD) / (pE | pF | plot_spacer()) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 30, face = "bold")) & NOTITLE
ggsave(file.path(FIG, "Fig20_sorghum_null.png"), p, width = 30, height = 18, dpi = 300, bg = "white")

## =============================================================================
## Main sorghum figure (dry biomass)
## =============================================================================
clean <- function(x) { x <- gsub("<[^>]+>", "", x); x <- gsub("&alpha;", "alpha", x, fixed = TRUE)
                       x <- gsub("&beta;", "beta", x, fixed = TRUE); gsub("&[a-z]+;", "", x) }

## ---- A: top pathways ------------------------------------------------------------
wy <- fread(file.path(RES, sprintf("%s_P_SNPWISE_MEAN_5comp_westfall_young.csv", TRAIT)))
wy <- wy[order(z_std)][seq_len(N_TOP)]
wy[, lab := sprintf("%s (%d)", clean(pathway_name), n_genes)]
wy[, lab := factor(lab, levels = rev(lab))]
mA <- ggplot(wy, aes(-log10(p_calibrated), lab)) +
  geom_col(fill = "black", width = 0.7) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", colour = "grey30") +
  geom_text(aes(label = sprintf("FWER %.2f", p_fwer_std)), hjust = -0.1, size = 7) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(x = expression(bold(-log[10] ~ calibrated ~ italic(p))), y = NULL) + plot_theme +
  theme(panel.grid.major.y = element_blank())

## ---- B: top genes and their rank in each layer ----------------------------------
## annotations as in 26_sorghum_results.R: SorghumBicolorCyc enzyme names, with
## the one override for the WRKY gene that the database lists as a thioesterase
ANNOT_OVERRIDE <- c(SORBI_3008G029200 = "WRKY transcription factor")
pwf <- fread(PATHWAY_FILE, sep = "\t", quote = "", showProgress = FALSE)
pwf <- pwf[nzchar(`Protein-name`)]
enz <- pwf[, .(name = `Protein-name`[1]), by = `Gene-name`]
annot <- function(id) {
  a <- if (id %in% names(ANNOT_OVERRIDE)) ANNOT_OVERRIDE[[id]] else enz$name[match(id, enz$`Gene-name`)]
  if (is.na(a) || !nzchar(a)) "" else clean(a)
}
gp <- fread(file.path(RES, sprintf("25_%s_gene_prioritize.csv", TRAIT)))
gp <- gp[order(rank)][seq_len(N_TOP)]
gp[, ann := vapply(gene, annot, "")]
gp[, lab := sprintf("%d  %s%s", rank, sub("^SORBI_", "", gene), ifelse(nzchar(ann), paste0(" (", ann, ")"), ""))]
gp[, lab := factor(lab, levels = rev(lab))]
lay <- melt(gp[, .(lab, GWAS = gwas_rank, MAGMA = magma_rank, Pathway = pathway_rank)],
            id.vars = "lab", variable.name = "layer", value.name = "r")
lay <- lay[is.finite(r)]
mB <- ggplot(lay, aes(r, lab)) +
  geom_line(aes(group = lab), colour = "grey70", linewidth = 0.8) +
  geom_point(aes(colour = layer, shape = layer), size = 6) +
  scale_x_log10(limits = c(1, 10^ceiling(log10(max(lay$r)))), expand = expansion(mult = c(0.02, 0.05))) +
  scale_colour_viridis_d(end = 0.85) +
  labs(x = "Rank within layer (1 = strongest)", y = NULL) + plot_theme +
  theme(panel.grid.major.y = element_blank())

## ---- C: the permutation null -----------------------------------------------------
lam1 <- function(p) { p <- p[is.finite(p) & p > 0 & p <= 1]; qchisq(median(p), 1, lower.tail = FALSE) / qchisq(0.5, 1) }
two <- runs[names(runs) %in% c("nogc", "gc")]
lab_run <- vapply(two, function(r) {
  f <- per[run == r]; sprintf("%s: \u03bb %.2f, FWER %.2f", r, lam1(pw[run == r]$p_calibrated),
                              mean(f$any_fwer_05)) }, "")
qc <- qq[run %in% two][, run := factor(lab_run[match(as.character(run), two)], levels = lab_run)]
mC <- ggplot(qc, aes(exp, obs, colour = run)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey30") +
  geom_line(linewidth = 1.5) + scale_colour_manual(values = setNames(c("black", "grey60"), lab_run)) +
  guides(colour = guide_legend(ncol = 1)) +
  labs(x = expression(bold(Expected ~ -log[10](p))), y = expression(bold(Observed ~ -log[10](p)))) +
  plot_theme

fw2 <- fw[run %in% two]
mD <- ggplot(fw2, aes(run, fwer, colour = run)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey30") +
  geom_pointrange(aes(ymin = lo, ymax = hi), linewidth = 1.5, size = 1.4, show.legend = FALSE) +
  scale_colour_manual(values = setNames(c("black", "grey60"), two)) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = NULL, y = "Family-wise error") + plot_theme

m <- (mA | mB) / (free(mC) | mD) + plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 30, face = "bold")) & NOTITLE
ggsave(file.path(FIG, "Fig21_sorghum_dry_tons.png"), m, width = 32, height = 20, dpi = 300, bg = "white")

## ---- table --------------------------------------------------------------------
lam <- function(p) { p <- p[is.finite(p) & p > 0 & p <= 1]; qchisq(median(p), 1, lower.tail = FALSE) / qchisq(0.5, 1) }
tab <- pw[, .(lambda = lam(p_calibrated), t05 = mean(p_calibrated < 0.05), t01 = mean(p_calibrated < 0.01),
              t001 = mean(p_calibrated < 0.001)), by = run][fw, on = "run"]
f3 <- function(x) formatC(x, format = "f", digits = 3)
rowtex <- function(lab, v) paste0(lab, " & ", paste(v, collapse = " & "), " \\\\")
lines <- c("\\begin{tabular}{l" , paste(rep("c", nrow(tab)), collapse = ""), "}", "\\toprule",
  paste0(" & ", paste(as.character(tab$run), collapse = " & "), " \\\\"), "\\midrule",
  rowtex("Shuffles", tab$n),
  rowtex("$\\lambda$, calibrated $p$", formatC(tab$lambda, format = "f", digits = 2)),
  rowtex("Type I error, $\\alpha = 0.05$", f3(tab$t05)),
  rowtex("Type I error, $\\alpha = 0.01$", f3(tab$t01)),
  rowtex("Type I error, $\\alpha = 0.001$", formatC(tab$t001, format = "f", digits = 4)),
  rowtex("Family-wise error (95\\% CI)", sprintf("%.2f (%.2f--%.2f)", tab$fwer, tab$lo, tab$hi)),
  "\\bottomrule", "\\end{tabular}")
lines <- c(paste0(lines[1], lines[2], lines[3]), lines[-(1:3)])
writeLines(lines, file.path(TEX, "tab_sorghum_null.tex"))
print(tab)
cat("wrote", file.path(FIG, "Fig21_sorghum_dry_tons.png"), file.path(FIG, "Fig20_sorghum_null.png"),
    "and", file.path(TEX, "tab_sorghum_null.tex"), "\n")
