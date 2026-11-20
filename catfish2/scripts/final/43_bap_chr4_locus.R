#!/usr/bin/env Rscript
## =============================================================================
## 43_bap_chr4_locus.R -- the BAP chr4 63.9-64.0 Mb locus shared by stem sugar
## (Brix) and fibre (NDF, NFC, ADF): the sugar-fibre trade-off between lines and
## the allele effect on both traits, in one figure.
## Reads   catfish2/data/bap/BAP_pheno_vcf2gwas.csv
##         catfish2/results/bap_traits/chr4_locus/  (from hpc/bap_gwas/05_chr4_locus_LSF.sh)
##           leads.tsv, lead_effects.tsv, region_<trait>.tsv, leads_geno.raw,
##           leads_ld.ld, leads_ld.snplist
##         catfish2/results/bap_traits/prioritize/<trait>_gene_prioritize.csv (panel D)
## Writes  catfish2/figures/Fig23_bap_chr4_locus.png
##         catfish2/results/bap_traits/chr4_locus/43_chr4_locus_summary.csv
## Panels  A  each homozygous line's Brix (x) against NDF (y), by genotype at the
##            focal SNP, with genotype means and 95% CI: the allele shift drawn on the
##            between-line correlation
##         B  focal-SNP effect on every trait (GEMMA beta / trait SD, 95% CI),
##            for the Brix-increasing allele
##         C  every SNP in the window: GWAS z for Brix against z for NDF (one
##            shared signal puts the strong SNPs on one line through the origin)
##         D  the CATFISH gene score in the locus: each gene's rank among all
##            genes for Brix, NDF and NFC (prioritize/<trait>_gene_prioritize.csv
##            from hpc/bap_gwas/04), genes ordered by position
## Printed: r2 between the traits' lead SNPs; NDF ~ genotype with and without
## Brix in the model (how much of the fibre effect the trade-off accounts for).
##   Rscript 43_bap_chr4_locus.R
## Env: FOCAL_SNP (default the X_TRAIT lead), X_TRAIT (Brix_maturity),
##      Y_TRAIT (NDF_pct_DM), PRIO_TRAITS (Brix_maturity,NDF_pct_DM,NFC_pct_DM),
##      CATFISH_BASE, LOCUS_DIR, PRIO_DIR, PHENO_FILE, FIG_DIR,
##      LOCUS_FROM / LOCUS_TO (panel D window, bp; 63900000 / 64060000),
##      MAGMA_WINDOW (bp MAGMA added to each gene start; 25000)
## =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
BASE <- Sys.getenv("CATFISH_BASE", path.expand("~/Documents/Github/CATFISH"))
LOC  <- Sys.getenv("LOCUS_DIR", file.path(BASE, "catfish2/results/bap_traits/chr4_locus"))
PHE  <- Sys.getenv("PHENO_FILE", file.path(BASE, "catfish2/data/bap/BAP_pheno_vcf2gwas.csv"))
FIG  <- Sys.getenv("FIG_DIR", file.path(BASE, "catfish2/figures"))
PRD  <- Sys.getenv("PRIO_DIR", file.path(BASE, "catfish2/results/bap_traits/prioritize"))
PT   <- strsplit(Sys.getenv("PRIO_TRAITS", "Brix_maturity,NDF_pct_DM,NFC_pct_DM"), ",")[[1]]
LFROM <- as.numeric(Sys.getenv("LOCUS_FROM", "63900000"))
LTO   <- as.numeric(Sys.getenv("LOCUS_TO", "64060000"))
MW    <- as.numeric(Sys.getenv("MAGMA_WINDOW", "25000"))
XT   <- Sys.getenv("X_TRAIT", "Brix_maturity")
YT   <- Sys.getenv("Y_TRAIT", "NDF_pct_DM")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

LAB <- c(Anthesis_date_days = "Anthesis date", Harvest_date_days = "Harvest date",
         Total_fresh_weight_kg = "Fresh weight", Dry_tons_per_acre = "Dry biomass",
         Stalk_height_cm = "Stalk height", Brix_milk = "Brix (milk)", Brix_maturity = "Brix (maturity)",
         NFC_pct_DM = "NFC", NDF_pct_DM = "NDF", ADF_pct_DM = "ADF", Lignin_pct_DM = "Lignin")
AXIS <- c(Brix_maturity = "Brix at maturity", Brix_milk = "Brix at milk stage",
          NDF_pct_DM = "NDF (% dry matter)", ADF_pct_DM = "ADF (% dry matter)",
          NFC_pct_DM = "NFC (% dry matter)", Lignin_pct_DM = "Lignin (% dry matter)")
lab_of <- function(t, L) ifelse(t %in% names(L), L[t], t)

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

## ---- inputs --------------------------------------------------------------------
need <- file.path(LOC, c("leads.tsv", "lead_effects.tsv", "leads_geno.raw",
                         sprintf("region_%s.tsv", c(XT, YT))))
miss <- need[!file.exists(need)]
if (length(miss)) stop("missing:\n  ", paste(miss, collapse = "\n  "), call. = FALSE)

ph <- fread(PHE); setnames(ph, 1, "id")
leads <- fread(file.path(LOC, "leads.tsv"), header = FALSE, col.names = c("trait", "rs", "ps", "p_wald"))
eff   <- fread(file.path(LOC, "lead_effects.tsv"))
FOC   <- Sys.getenv("FOCAL_SNP", leads[trait == XT, rs][1])
if (is.na(FOC) || !nzchar(FOC)) stop("no lead SNP for ", XT, " in leads.tsv: set FOCAL_SNP", call. = FALSE)
options(width = 200)
cat("lead SNP per trait in the window:\n"); print(leads, row.names = FALSE)
cat("focal SNP:", FOC, "\n")

ldf <- file.path(LOC, "leads_ld.ld"); lds <- file.path(LOC, "leads_ld.snplist")
if (file.exists(ldf) && file.exists(lds)) {
  ld <- as.matrix(fread(ldf, header = FALSE)); sn <- readLines(lds)
  dimnames(ld) <- list(sn, sn)
  cat("\nr2 between the lead SNPs:\n"); print(round(ld, 2))
}

## ---- orientation: the allele that raises X_TRAIT ------------------------------
ex <- eff[rs == FOC & trait == XT]
if (nrow(ex) != 1) stop("lead_effects.tsv has no ", XT, " row for ", FOC, call. = FALSE)
UP <- if (ex$beta > 0) ex$allele1 else ex$allele0
DN <- if (ex$beta > 0) ex$allele0 else ex$allele1
cat(sprintf("allele raising %s: %s (other allele %s)\n", XT, UP, DN))

## ---- genotypes at the focal SNP -----------------------------------------------
raw  <- fread(file.path(LOC, "leads_geno.raw"))
gcol <- names(raw)[startsWith(names(raw), paste0(FOC, "_"))]
if (length(gcol) != 1) stop(FOC, " not in leads_geno.raw", call. = FALSE)
A1 <- substring(gcol, nchar(FOC) + 2)                        # allele the dosage counts
ids <- if (sum(raw$IID %in% ph$id) >= sum(raw$FID %in% ph$id)) raw$IID else raw$FID
G <- data.table(id = ids, dose = raw[[gcol]])
if (!A1 %in% c(UP, DN)) stop("plink allele ", A1, " is neither GEMMA allele (", UP, "/", DN, ")", call. = FALSE)
G[, n_up := if (A1 == UP) dose else 2 - dose]
cat(sprintf("%d of %d phenotyped lines genotyped at %s\n", sum(ph$id %in% G$id), nrow(ph), FOC))

d <- merge(ph[, c("id", XT, YT), with = FALSE], G[, .(id, n_up)], by = "id")
setnames(d, c(XT, YT), c("x", "y"))
d <- d[is.finite(x) & is.finite(y) & !is.na(n_up)]
GL <- c(paste0(DN, "/", DN), paste0(UP, "/", DN), paste0(UP, "/", UP))
d[, geno := factor(GL[n_up + 1], levels = GL)]
d[, geno := droplevels(geno)]
print(d[, .N, by = geno][order(geno)])

## ---- tests: between-line trade-off; allele effect with and without X ----------
ct  <- cor.test(d$x, d$y)
m_x <- lm(x ~ n_up, d); m_y <- lm(y ~ n_up, d); m_yx <- lm(y ~ n_up + x, d)
co  <- function(m, v) summary(m)$coefficients[v, ]
by_m <- co(m_y, "n_up"); by_c <- co(m_yx, "n_up"); bx <- co(m_x, "n_up")
explained <- 1 - by_c[1] / by_m[1]
cat(sprintf("\nlines with both traits: %d; r(%s, %s) = %.2f (p = %.2g)\n", nrow(d), XT, YT, ct$estimate, ct$p.value))
cat(sprintf("per %s allele: %s %+.3f (p = %.2g); %s %+.3f (p = %.2g)\n", UP, XT, bx[1], bx[4], YT, by_m[1], by_m[4]))
cat(sprintf("%s per allele adjusting for %s: %+.3f (p = %.2g); the trade-off accounts for %.0f%% of the %s effect\n",
            YT, XT, by_c[1], by_c[4], 100 * explained, YT))

## ---- panel A: trade-off and allele shift (homozygous lines; the few
## heterozygous calls of these inbred lines stay in the tests above) ----------------
dp <- d[n_up != 1]; dp[, geno := droplevels(geno)]
mu <- dp[, .(x = mean(x), y = mean(y), sx = sd(x) / sqrt(.N), sy = sd(y) / sqrt(.N), n = .N), by = geno]
mu <- mu[n >= 3]
COL <- setNames(c("grey60", "black")[seq_len(nlevels(dp$geno))], levels(dp$geno))
legend_lab <- sprintf("%s (n = %d)", levels(dp$geno), as.integer(table(dp$geno)))
pA <- ggplot(dp, aes(x, y)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black", linetype = "dashed", linewidth = 0.8) +
  geom_point(aes(colour = geno), size = 2.6, alpha = 0.55) +
  geom_errorbar(data = mu, aes(x = x, ymin = y - 1.96 * sy, ymax = y + 1.96 * sy, colour = geno),
                width = 0, linewidth = 1.4, inherit.aes = FALSE) +
  geom_errorbarh(data = mu, aes(y = y, xmin = x - 1.96 * sx, xmax = x + 1.96 * sx, colour = geno),
                 height = 0, linewidth = 1.4, inherit.aes = FALSE) +
  geom_point(data = mu, aes(x, y, fill = geno), shape = 21, colour = "white", size = 7, stroke = 1.5,
             inherit.aes = FALSE, show.legend = FALSE) +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.4, size = 7,
           label = sprintf("r = %.2f", ct$estimate)) +
  scale_colour_manual(values = COL, labels = legend_lab) + scale_fill_manual(values = COL) +
  labs(x = lab_of(XT, AXIS), y = lab_of(YT, AXIS)) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1, linetype = 0))) +
  plot_theme + NOTITLE

## ---- panel B: focal-SNP effect on every trait --------------------------------------
e <- eff[rs == FOC & trait %in% names(ph)]
e[, sd_t := vapply(trait, function(t) sd(ph[[t]], na.rm = TRUE), 0)]
e[, s := ifelse(allele1 == UP, 1, ifelse(allele0 == UP, -1, NA_real_))]
e[, `:=`(b = s * beta / sd_t, lo = (s * beta - 1.96 * se) / sd_t, hi = (s * beta + 1.96 * se) / sd_t)]
e[, lab := factor(lab_of(trait, LAB), levels = rev(unique(c(LAB, trait))))]
e[, ptxt := ifelse(p_wald < 1e-3, sprintf("%.0e", p_wald), sprintf("%.2f", p_wald))]
xr <- max(abs(c(e$lo, e$hi)), na.rm = TRUE)
pB <- ggplot(e, aes(b, lab)) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 1.1) +
  geom_point(size = 4.5) +
  geom_text(aes(x = xr * 1.08, label = ptxt), hjust = 0, size = 6.5) +
  scale_x_continuous(limits = c(-xr, xr * 1.6), breaks = scales::pretty_breaks(4)) +
  labs(x = sprintf("Effect of %s allele (trait SD)", UP), y = NULL) +
  plot_theme + NOTITLE

## ---- panel C: SNP-by-SNP z for X against Y across the window ---------------------
rx <- fread(file.path(LOC, sprintf("region_%s.tsv", XT)))
ry <- fread(file.path(LOC, sprintf("region_%s.tsv", YT)))
r  <- merge(rx[, .(rs, ps, a1x = allele1, a0x = allele0, zx = beta / se)],
            ry[, .(rs, a1y = allele1, zy = beta / se)], by = "rs")
r[, zy := ifelse(a1y == a1x, zy, -zy)]
r[, s := sign(zx)]                                      # orient every SNP to its X-raising allele
r[, `:=`(zx = abs(zx), zy = zy * s)]
r[, kind := ifelse(rs == FOC, "focal", ifelse(rs %in% leads$rs, "lead", "other"))]
rz <- cor(r$zx, r$zy)
cat(sprintf("\n%d SNPs in the window; cor(z %s, z %s) = %.2f\n", nrow(r), XT, YT, rz))
pC <- ggplot(r, aes(zx, zy)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_point(data = r[kind == "other"], colour = "grey60", size = 2, alpha = 0.6) +
  geom_point(data = r[kind != "other"], aes(shape = kind), colour = "black", size = 5) +
  scale_shape_manual(values = c(focal = 17, lead = 16),
                     labels = c(focal = "Focal SNP", lead = "Other traits' lead SNPs")) +
  labs(x = sprintf("|z| %s", lab_of(XT, LAB)), y = sprintf("z %s (same allele)", lab_of(YT, LAB))) +
  plot_theme + NOTITLE

## ---- panel D: the gene score across the locus -------------------------------------
GNAME <- c(SORBI_3004G300800 = "SbALDH3E1", SORBI_3004G300900 = "TRXR", SORBI_3004G301100 = "GATA")
pf <- file.path(PRD, sprintf("%s_gene_prioritize.csv", PT))
if (!all(file.exists(pf))) stop("missing:\n  ", paste(pf[!file.exists(pf)], collapse = "\n  "), call. = FALSE)
pr <- rbindlist(lapply(seq_along(PT), function(i) {
  x <- fread(pf[i], select = c("gene", "rank", "chr", "start"))
  x[, pos := start + MW][as.character(chr) == "4" & pos >= LFROM & pos <= LTO][, trait := PT[i]]
}))
pr[, lab := ifelse(gene %in% names(GNAME), sprintf("%s (%s)", sub("SORBI_3004G", "", gene), GNAME[gene]),
                   sub("SORBI_3004G", "", gene))]
pr[, lab := factor(lab, levels = unique(lab[order(pos)]))]
pr[, tl := factor(lab_of(trait, LAB), levels = lab_of(PT, LAB))]
cat("\nCATFISH rank of the genes in the locus:\n")
print(dcast(pr, lab ~ trait, value.var = "rank"), row.names = FALSE)
TCOL <- setNames(viridisLite::viridis(length(PT), end = 0.85), levels(pr$tl))
pD <- ggplot(pr, aes(lab, rank, colour = tl, shape = tl)) +
  geom_hline(yintercept = 100, linetype = "dashed", colour = "grey50") +
  geom_point(size = 5, alpha = 0.9) +
  scale_y_continuous(trans = scales::compose_trans("log10", "reverse"),
                breaks = c(1, 10, 100, 1000, 10000), labels = scales::comma) +
  scale_colour_manual(values = TCOL) + scale_shape_manual(values = c(16, 17, 15)[seq_along(PT)]) +
  labs(x = "Gene (chromosome 4, by position)", y = "CATFISH rank (all genes)") +
  plot_theme + NOTITLE + theme(axis.text.x = element_text(size = 18, angle = 60, hjust = 1, colour = "black"))

fig <- ((pA | pB) / (pC | pD)) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 28, face = "bold"))
fo <- file.path(FIG, "Fig23_bap_chr4_locus.png")
ggsave(fo, fig, width = 24, height = 20, dpi = 200, bg = "white")

sm <- data.table(focal_snp = FOC, x_trait = XT, y_trait = YT, up_allele = UP, n_lines = nrow(d),
                 r_xy = ct$estimate, p_r = ct$p.value,
                 b_x = bx[1], p_x = bx[4], b_y = by_m[1], p_y = by_m[4],
                 b_y_adj_x = by_c[1], p_y_adj_x = by_c[4], frac_y_via_x = explained,
                 n_window_snps = nrow(r), cor_z = rz)
fwrite(sm, file.path(LOC, "43_chr4_locus_summary.csv"))
cat("wrote", fo, "and", file.path(LOC, "43_chr4_locus_summary.csv"), "\n")
