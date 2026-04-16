# QQ Plot for CATFISH B=1M GPD Results - Sorghum Stem Volume
# For "Comparison of CATFISH component tests" figure

library(ggplot2)
library(dplyr)
library(patchwork)

# Load results
res <- read.csv("sorghum_catfish_results/sorghum_stem_vol_CATFISH_B1000000_GPD.csv",
                stringsAsFactors = FALSE)

n_pathways <- nrow(res)
bonf_thresh <- 0.05 / n_pathways

# Calculate expected and observed -log10(p)
res <- res %>%
  arrange(omni_p_final) %>%
  mutate(
    rank = row_number(),
    expected_p = rank / (n_pathways + 1),
    expected_neglog10 = -log10(expected_p),
    observed_neglog10 = -log10(omni_p_final),
    significance = ifelse(omni_p_final < bonf_thresh,
                          "Bonferroni significant",
                          "Non-significant")
  )

res$significance <- factor(res$significance,
                           levels = c("Bonferroni significant", "Non-significant"))

n_sig <- sum(res$omni_p_final < bonf_thresh)
max_obs <- max(res$observed_neglog10)

cat("Summary:\n")
cat("  Total pathways:", n_pathways, "\n")
cat("  Bonferroni threshold:", format(bonf_thresh, digits = 3), "\n")
cat("  Bonferroni significant:", n_sig, "\n")
cat("  Max observed -log10(p):", round(max_obs, 1), "\n")

# =============================================================================
# QQ Plot with inset for extreme values
# =============================================================================

# Main panel - lower range (0-6)
panel_main <- ggplot(res, aes(x = expected_neglog10, y = observed_neglog10)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(color = significance), size = 2.5, alpha = 0.8) +
  scale_color_manual(
    values = c("Bonferroni significant" = "#D7191C", "Non-significant" = "gray60"),
    name = NULL
  ) +
  geom_hline(yintercept = -log10(bonf_thresh), linetype = "dotted", color = "#D7191C", linewidth = 0.8) +
  annotate("text", x = 0.1, y = -log10(bonf_thresh) + 0.3,
           label = sprintf("Bonferroni (P < %.2e)", bonf_thresh),
           hjust = 0, size = 3, color = "#D7191C") +
  labs(
    x = expression(Expected ~ -log[10](P)),
    y = expression(Observed ~ -log[10](P))
  ) +
  coord_cartesian(xlim = c(0, 3), ylim = c(0, 6)) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  )

# Inset panel - extreme values
panel_inset <- ggplot(res %>% filter(observed_neglog10 > 5),
                      aes(x = expected_neglog10, y = observed_neglog10)) +
  geom_point(aes(color = significance), size = 2, alpha = 0.8) +
  scale_color_manual(
    values = c("Bonferroni significant" = "#D7191C", "Non-significant" = "gray60"),
    guide = "none"
  ) +
  scale_y_continuous(breaks = seq(10, 35, by = 10)) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(xlim = c(2, 2.7), ylim = c(5, 36)) +
  theme_bw(base_size = 9) +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", color = "gray40", linewidth = 0.5),
    plot.margin = margin(3, 3, 3, 3),
    axis.text = element_text(size = 7)
  )

# Combine with inset
qq_inset <- panel_main +
  inset_element(panel_inset, left = 0.4, bottom = 0.45, right = 0.98, top = 0.98) +
  # Add annotation pointing to inset
 plot_annotation(
    caption = sprintf("n = %d pathways | %d Bonferroni significant | max -log10(P) = %.1f",
                     n_pathways, n_sig, max_obs)
  ) &
  theme(plot.caption = element_text(size = 9, hjust = 0.5, color = "gray40"))

# =============================================================================
# Alternative: Broken axis version
# =============================================================================

# Transform y-axis: compress values > 5
res_broken <- res %>%
  mutate(
    y_display = case_when(
      observed_neglog10 <= 5 ~ observed_neglog10,
      observed_neglog10 > 5 ~ 5 + (observed_neglog10 - 5) * 0.12
    )
  )

max_y_display <- max(res_broken$y_display)

# Custom breaks
y_breaks <- c(0, 1, 2, 3, 4, 5, 5 + (10-5)*0.12, 5 + (20-5)*0.12, 5 + (30-5)*0.12)
y_labels <- c("0", "1", "2", "3", "4", "5", "10", "20", "30")

qq_broken <- ggplot(res_broken, aes(x = expected_neglog10, y = y_display)) +
  # Diagonal in lower region
  annotate("segment", x = 0, xend = 3, y = 0, yend = 3,
           linetype = "dashed", color = "gray50") +
  # Break indicator
  annotate("rect", xmin = -0.15, xmax = 3.05, ymin = 4.95, ymax = 5.15,
           fill = "white", color = NA) +
  annotate("segment", x = -0.05, xend = 0.1, y = 4.9, yend = 5.2, color = "gray30", linewidth = 0.5) +
  annotate("segment", x = 0.1, xend = 0.25, y = 4.9, yend = 5.2, color = "gray30", linewidth = 0.5) +
  # Points
  geom_point(aes(color = significance), size = 2.5, alpha = 0.8) +
  scale_color_manual(
    values = c("Bonferroni significant" = "#D7191C", "Non-significant" = "gray60"),
    name = NULL
  ) +
  # Bonferroni line
  geom_hline(yintercept = -log10(bonf_thresh), linetype = "dotted", color = "#D7191C", linewidth = 0.8) +
  scale_y_continuous(breaks = y_breaks, labels = y_labels) +
  labs(
    x = expression(Expected ~ -log[10](P)),
    y = expression(Observed ~ -log[10](P)),
    caption = sprintf("n = %d pathways | %d Bonferroni significant (P < %.2e) | max -log10(P) = %.1f",
                     n_pathways, n_sig, bonf_thresh, max_obs)
  ) +
  coord_cartesian(xlim = c(0, 3), ylim = c(0, max_y_display + 0.2), clip = "off") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = c(0.75, 0.15),
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.text = element_text(size = 9),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(size = 9, hjust = 0.5, color = "gray40")
  )

# =============================================================================
# Save figures
# =============================================================================

fig_dir <- "sorghum_catfish_results/figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# Also save to main Figures folder for the component test figure
main_fig_dir <- "Figures/Fig_Sorghum_Component_Tests"
if (!dir.exists(main_fig_dir)) dir.create(main_fig_dir, recursive = TRUE)

# Save inset version
ggsave(file.path(fig_dir, "qq_plot_GPD_B1M_pathways.png"), qq_inset, width = 6, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "qq_plot_GPD_B1M_pathways.pdf"), qq_inset, width = 6, height = 6)

# Save broken axis version
ggsave(file.path(fig_dir, "qq_plot_GPD_B1M_broken.png"), qq_broken, width = 6, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "qq_plot_GPD_B1M_broken.pdf"), qq_broken, width = 6, height = 6)

# Copy to main figures folder
ggsave(file.path(main_fig_dir, "qq_omnibus_pathways_sorghum.png"), qq_inset, width = 6, height = 6, dpi = 300)
ggsave(file.path(main_fig_dir, "qq_omnibus_pathways_sorghum.pdf"), qq_inset, width = 6, height = 6)

cat("\nFigures saved:\n")
cat("  ", fig_dir, "/qq_plot_GPD_B1M_pathways.png/pdf (inset version)\n")
cat("  ", fig_dir, "/qq_plot_GPD_B1M_broken.png/pdf (broken axis version)\n")
cat("  ", main_fig_dir, "/qq_omnibus_pathways_sorghum.png/pdf\n")
