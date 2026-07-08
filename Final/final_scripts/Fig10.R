#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  source("/Users/nirwantandukar/Documents/Github/MAGCAT/scripts/all_figs.R")
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(cowplot)
  library(tidyr)
})

default_dir <- "/Users/nirwantandukar/Documents/Github/MAGCAT/Final/final_simulation_results/SAP"
out_dir <- Sys.getenv("SAP_NULL_DIR", default_dir)
raw_dir <- Sys.getenv("SAP_NULL_RAW_DIR", file.path(out_dir, "catfish"))
use_raw_panels <- tolower(Sys.getenv("SAP_NULL_USE_RAW", "false")) %in% c("true", "1", "yes")

lambda_file <- file.path(out_dir, "lambda_table.csv")
type1_file <- file.path(out_dir, "type1_error_table.csv")
hist_file <- file.path(out_dir, "hist_omnibus_null.png")
qq_file <- file.path(out_dir, "QQ_omnibus_null.png")

for (f in c(lambda_file, type1_file)) {
  if (!file.exists(f)) stop("Missing required SAP null-diagnostic file: ", f, call. = FALSE)
}

lambda_dt <- fread(lambda_file)
type1_dt <- fread(type1_file)

has_final <- "omni_p_final" %in% lambda_dt$method
has_mvn <- "omni_p_mvn" %in% lambda_dt$method
has_compcal <- "omni_p_mvn_compcal" %in% lambda_dt$method

drop_redundant_mvn <- FALSE
if (has_final && has_mvn) {
  lam_final <- lambda_dt$lambda[match("omni_p_final", lambda_dt$method)]
  lam_mvn <- lambda_dt$lambda[match("omni_p_mvn", lambda_dt$method)]
  if (isTRUE(all.equal(lam_final, lam_mvn, tolerance = 1e-15))) {
    drop_redundant_mvn <- TRUE
  }
}

pretty_method <- c(
  omni_p_final = "Omnibus Alone",
  omni_p_mvn = "Omnibus Alone",
  omni_p_analytic = "Omnibus Analytical",
  omni_p_mvn_compcal = "Omnibus Combined",
  acat_p = "ACAT",
  fisher_p = "Fisher",
  minp_p_analytic = "minP",
  stouffer_p_analytic = "Stouffer",
  tfisher_p_analytic = "TFisher"
)

method_order <- c(
  "ACAT",
  "Fisher",
  "TFisher",
  "minP",
  "Stouffer",
  "Omnibus Alone",
  "Omnibus Combined",
  "Omnibus Analytical"
)

if (drop_redundant_mvn) {
  lambda_dt <- lambda_dt %>% filter(method != "omni_p_mvn")
  type1_dt <- type1_dt %>% filter(method != "omni_p_mvn")
}

if (!has_compcal) {
  lambda_dt <- lambda_dt %>% filter(method != "omni_p_mvn_compcal")
  type1_dt <- type1_dt %>% filter(method != "omni_p_mvn_compcal")
}

lambda_dt <- lambda_dt %>%
  mutate(
    method_label = dplyr::recode(method, !!!pretty_method),
    method_label = factor(method_label, levels = method_order)
  ) %>%
  arrange(method_label)

type1_dt <- type1_dt %>%
  mutate(
    method_label = dplyr::recode(method, !!!pretty_method),
    method_label = factor(method_label, levels = method_order),
    alpha_label = factor(
      ifelse(alpha == 0.05, expression(alpha == 0.05), expression(alpha == 0.01)),
      levels = c(expression(alpha == 0.05), expression(alpha == 0.01))
    )
  ) %>%
  arrange(method_label, desc(alpha))

panel_theme <- plot_theme +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.35),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    plot.margin = margin(6, 8, 6, 8)
  )

method_colors <- c(
  ACAT = "#D55E00",
  Fisher = "#CC79A7",
  TFisher = "#7B61FF",
  minP = "#E6AB02",
  Stouffer = "#1B9E77",
  "Omnibus Alone" = "#B22222",
  "Omnibus Combined" = "#7A4E2D",
  "Omnibus Analytical" = "#8C8C8C"
)

p_lambda <- ggplot(lambda_dt, aes(method_label, lambda, fill = method_label)) +
  geom_col(width = 0.72, alpha = 0.9, show.legend = FALSE) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#2C7FB8", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f", lambda)), vjust = -0.35, size = 3.3) +
  scale_fill_manual(values = method_colors, drop = FALSE) +
  labs(x = "Method", y = expression(lambda)) +
  coord_cartesian(ylim = c(0, max(lambda_dt$lambda, na.rm = TRUE) * 1.22)) +
  panel_theme +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

p_lambda <- ggdraw(p_lambda) +
  draw_label("A", x = 0.01, y = 0.99, hjust = 0, vjust = 1, fontface = "bold", size = 20)

type1_plot <- type1_dt %>%
  mutate(alpha_f = factor(sprintf("alpha = %.2f", alpha), levels = c("alpha = 0.05", "alpha = 0.01")))

p_type1 <- ggplot(type1_plot, aes(method_label, observed, fill = method_label)) +
  geom_col(width = 0.72, alpha = 0.9, show.legend = FALSE) +
  geom_errorbar(
    aes(ymin = ci_lo, ymax = ci_hi),
    width = 0.16,
    linewidth = 0.45
  ) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "#D7301F", linewidth = 0.7) +
  geom_hline(yintercept = 0.01, linetype = "dotted", color = "#252525", linewidth = 0.7) +
  scale_fill_manual(values = method_colors, drop = FALSE) +
  labs(x = "Method", y = "Type I error") +
  facet_wrap(~alpha_f, ncol = 2) +
  coord_cartesian(ylim = c(0, max(type1_plot$ci_hi, na.rm = TRUE) * 1.18)) +
  panel_theme +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.text = element_text(size = 11, face = "bold")
  )

p_type1 <- ggdraw(p_type1) +
  draw_label("B", x = 0.01, y = 0.99, hjust = 0, vjust = 1, fontface = "bold", size = 20)

perm_files <- if (use_raw_panels && dir.exists(raw_dir)) {
  list.files(raw_dir, pattern = "^pathway_pvals_perm_.*\\.csv$", full.names = TRUE)
} else {
  character(0)
}

if (length(perm_files) > 0) {
  dat <- rbindlist(lapply(perm_files, fread), fill = TRUE)
  head_p <- if ("omni_p_final" %in% names(dat)) "omni_p_final" else names(dat)[1]
  pvals <- dat[[head_p]]
  pvals <- pvals[is.finite(pvals) & pvals > 0 & pvals <= 1]

  qd <- data.table(
    expected = -log10(ppoints(length(pvals))),
    observed = -log10(sort(pvals))
  )
  k <- seq_len(nrow(qd))
  ci_band <- data.table(
    expected = -log10(ppoints(nrow(qd))),
    lo = -log10(qbeta(0.975, k, nrow(qd) - k + 1)),
    hi = -log10(qbeta(0.025, k, nrow(qd) - k + 1))
  )

  p_hist_inner <- ggplot(data.table(p = pvals), aes(p)) +
    geom_histogram(
      breaks = seq(0, 1, 0.05),
      fill = "#4C78A8",
      colour = "white",
      linewidth = 0.25
    ) +
    geom_hline(yintercept = length(pvals) / 20, colour = "red", linetype = "dashed", linewidth = 0.7) +
    labs(
      title = "Null p-values should be ~Uniform(0,1)",
      x = "Omnibus p-value (null)",
      y = "Count"
    ) +
    panel_theme

  p_hist <- ggdraw(p_hist_inner) +
    draw_label("C", x = 0.01, y = 0.99, hjust = 0, vjust = 1, fontface = "bold", size = 20)

  p_qq_inner <- ggplot() +
    geom_ribbon(data = ci_band, aes(expected, ymin = lo, ymax = hi), fill = "grey80", alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed", linewidth = 0.7) +
    geom_point(data = qd, aes(expected, observed), size = 0.7, alpha = 0.6) +
    labs(
      title = "CATFISH omnibus under SAP stem-volume pheno",
      subtitle = sprintf(
        "%d permutations, %d pooled pathway p-values; lambda = %.2f",
        length(perm_files), nrow(qd),
        lambda_dt$lambda[match("omni_p_final", lambda_dt$method)]
      ),
      x = expression(Expected ~ -log[10](p)),
      y = expression(Observed ~ -log[10](p))
    ) +
    panel_theme +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5)
    )

  p_qq <- ggdraw(p_qq_inner) +
    draw_label("D", x = 0.01, y = 0.99, hjust = 0, vjust = 1, fontface = "bold", size = 20)
} else {
  if (!file.exists(hist_file) || !file.exists(qq_file)) {
    stop(
      "Pre-rendered SAP QQ/hist images are missing.",
      call. = FALSE
    )
  }

  p_hist <- ggdraw() +
    draw_image(hist_file, scale = 0.96) +
    draw_label("C", x = 0.01, y = 0.99, hjust = 0, vjust = 1, fontface = "bold", size = 20)

  p_qq <- ggdraw() +
    draw_image(qq_file, scale = 0.96) +
    draw_label("D", x = 0.01, y = 0.99, hjust = 0, vjust = 1, fontface = "bold", size = 20)
}

full_fig <- plot_grid(
  p_lambda, p_type1,
  p_hist, p_qq,
  ncol = 2,
  rel_widths = c(1, 1.12),
  rel_heights = c(0.95, 1.05)
)

png_file <- file.path(out_dir, "Figure_SAP_Stem_volume_Null_Diagnostics.png")
pdf_file <- file.path(out_dir, "Figure_SAP_Stem_volume_Null_Diagnostics.pdf")

ggsave(png_file, full_fig, width = 15, height = 11, dpi = 300, bg = "white")
ggsave(pdf_file, full_fig, width = 15, height = 11, bg = "white")

main_png <- "/Users/nirwantandukar/Documents/Github/MAGCAT/Final/final_main_fig/SAP_Stem_volume_Null_Diagnostics.png"
main_pdf <- "/Users/nirwantandukar/Documents/Github/MAGCAT/Final/final_main_fig/SAP_Stem_volume_Null_Diagnostics.pdf"
file.copy(png_file, main_png, overwrite = TRUE)
file.copy(pdf_file, main_pdf, overwrite = TRUE)

cat("Saved figure:\n")
cat("  ", png_file, "\n", sep = "")
cat("  ", pdf_file, "\n", sep = "")
cat("Copied to:\n")
cat("  ", main_png, "\n", sep = "")
cat("  ", main_pdf, "\n", sep = "")
