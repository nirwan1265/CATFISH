#!/usr/bin/env Rscript
## =============================================================================
## 35_simulation_figures.R -- figure and table for the known-truth simulation:
## CATFISH (self-contained, competitive, both) against MAGMA's gene-set analysis
## (self-contained and competitive) on GWAS statistics simulated on real LD
## (catfish2/scripts/hpc/simulation_check_catfish; 850 replicates, 50 per cell).
## Reads   catfish2/results/simulation_check_catfish/run_full/summary.csv
##         (written by 08_aggregate.R; PascalX rows are not used)
## Writes  catfish2/figures/Fig22_simulation_catfish_magma.png
##         catfish2/tex/auto/tab_simulation_compare.tex
##         catfish2/results/simulation_check_catfish/run_full/35_simulation_compare.csv
## Error and power intervals: exact binomial 95% CI over the 50 replicates.
## =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
BASE <- Sys.getenv("CATFISH_BASE", path.expand("~/Documents/Github/CATFISH"))
RES  <- Sys.getenv("SIM_RES", file.path(BASE, "catfish2/results/simulation_check_catfish/run_full"))
FIG  <- Sys.getenv("FIG_DIR", file.path(BASE, "catfish2/figures"))
TEX  <- Sys.getenv("TEX_AUTO", file.path(BASE, "catfish2/tex/auto"))
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

METH <- c(catfish_self = "CATFISH self-contained", magma_self = "MAGMA self-contained",
          catfish_comp = "CATFISH competitive", magma_comp = "MAGMA competitive",
          catfish_both = "CATFISH both")
ARCH <- c(SSS = "SSS", MDS = "MDS", WDE = "WDE", MSA = "MSA", SGP = "SGP")
SIZE <- c(small = "Small", medium = "Medium", large = "Large")

s <- fread(file.path(RES, "summary.csv"), colClasses = list(character = "size_class"))
s <- s[method %in% names(METH)]
ci <- function(r, n) { k <- round(r * n); if (!is.finite(k)) return(c(NA_real_, NA_real_))
  as.numeric(stats::binom.test(k, n)$conf.int) }
s[, c("off_l", "off_h") := as.list(ci(off_rate, n_rep)), by = .(scenario, size_class, method)]
s[, c("pow_l", "pow_h") := as.list(ci(power, n_rep)), by = .(scenario, size_class, method)]
s[, meth := factor(METH[method], levels = METH)]
fwrite(s, file.path(RES, "35_simulation_compare.csv"))

## ---- A: family-wise error with no enriched pathway ------------------------------
##   S0 (no signal): any call is false, every test
##   S1 (background genes only, no pathway enriched): any call is false for the
##   competitive tests and "both"; the self-contained tests are right to call
##   pathways that carry background genes, so they are not shown
a <- rbind(s[scenario == "S0"][, panel := "No signal"],
           s[scenario == "S1" & method %in% c("catfish_comp", "catfish_both", "magma_comp")][
             , panel := "Background only"])
a[, panel := factor(panel, levels = c("No signal", "Background only"))]
cols5 <- setNames(viridisLite::viridis(5, end = 0.9), METH)
pA <- ggplot(a, aes(meth, off_rate, colour = meth)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey30") +
  geom_pointrange(aes(ymin = off_l, ymax = off_h), linewidth = 1.3, size = 1.2) +
  facet_grid(. ~ panel, scales = "free_x", space = "free_x") +
  scale_colour_manual(values = cols5) +
  scale_x_discrete(labels = function(x) sub(" ", "\n", x)) +
  labs(x = NULL, y = "Family-wise error") + plot_theme +
  theme(legend.position = "none", strip.text = element_text(size = 24, face = "bold"))

## ---- B, C: power to call the enriched pathway -------------------------------------
pw <- s[scenario %in% c("S2", "S3", "S4", "S5", "S6")]
pw[, arch := factor(ARCH[archetype], levels = ARCH)]
pw[, size := factor(SIZE[size_class], levels = SIZE)]
bw <- c("black", "grey65")
power_panel <- function(m) {
  d <- pw[method %in% m]; d[, meth := droplevels(meth)]
  ggplot(d, aes(arch, power, fill = meth)) +
    geom_col(position = position_dodge(0.8), width = 0.75) +
    geom_errorbar(aes(ymin = pow_l, ymax = pow_h), position = position_dodge(0.8), width = 0.25) +
    facet_grid(. ~ size) +
    scale_fill_manual(values = setNames(bw, METH[m])) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(x = "Signal architecture", y = "Power") + plot_theme +
    theme(strip.text = element_text(size = 24, face = "bold"))
}
pB <- power_panel(c("catfish_self", "magma_self"))
pC <- power_panel(c("catfish_comp", "magma_comp"))

p <- pA / pB / pC + plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 30, face = "bold")) & NOTITLE
ggsave(file.path(FIG, "Fig22_simulation_catfish_magma.png"), p, width = 26, height = 26, dpi = 300, bg = "white")

## ---- table ---------------------------------------------------------------------------
f2 <- function(x) ifelse(is.na(x), "--", formatC(x, format = "f", digits = 2))
cell <- function(sc, sz, m, what) { v <- s[scenario == sc & size_class == sz & method == m][[what]]
  if (!length(v)) NA_real_ else v }
ms <- c("catfish_self", "magma_self", "catfish_comp", "magma_comp", "catfish_both")
lines <- c("\\begin{tabular}{llccccc}", "\\toprule",
  " & & \\multicolumn{2}{c}{Self-contained} & \\multicolumn{2}{c}{Competitive} & Both \\\\",
  "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-7}",
  "Scenario & Size & CATFISH & MAGMA & CATFISH & MAGMA & CATFISH \\\\", "\\midrule",
  "\\multicolumn{7}{l}{\\emph{Family-wise error (no pathway enriched)}} \\\\",
  paste0("No signal & -- & ", paste(f2(sapply(ms, function(m) cell("S0", "", m, "off_rate"))), collapse = " & "), " \\\\"),
  paste0("Background only & -- & -- & -- & ",
         paste(f2(sapply(c("catfish_comp", "magma_comp", "catfish_both"), function(m) cell("S1", "", m, "off_rate"))), collapse = " & "), " \\\\"),
  "\\midrule", "\\multicolumn{7}{l}{\\emph{Power to call the enriched pathway}} \\\\")
for (a in names(ARCH)) for (sz in names(SIZE)) {
  sc <- s[archetype == a, scenario][1]
  lines <- c(lines, paste0(if (sz == "small") a else "", " & ", SIZE[[sz]], " & ",
    paste(f2(sapply(ms, function(m) cell(sc, sz, m, "power"))), collapse = " & "), " \\\\"))
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}")
writeLines(lines, file.path(TEX, "tab_simulation_compare.tex"))
cat("wrote", file.path(FIG, "Fig22_simulation_catfish_magma.png"), "and", file.path(TEX, "tab_simulation_compare.tex"), "\n")
