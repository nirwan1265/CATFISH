suppressPackageStartupMessages({
  source("/Users/nirwantandukar/Documents/Github/MAGCAT/all_figs.R")
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# Final Figure 5 script
# Combined component-breakage stress test and adaptive/leave-one-out diagnostics
# Panels:
#   A. Component breakage lambda
#   B. Component breakage type I error
#   C. Adaptive/leave-one-out lambda
#   D. Adaptive/leave-one-out type I error
# ------------------------------------------------------------------------------

results_dir <- "/Users/nirwantandukar/Documents/Github/MAGCAT/simulation_results"
final_fig_dir <- "/Users/nirwantandukar/Documents/Github/MAGCAT/Final/final_main_fig"
dir.create(final_fig_dir, recursive = TRUE, showWarnings = FALSE)

block_d <- readRDS(file.path(results_dir, "block_d.rds"))
block_e_adaptive <- readRDS(file.path(results_dir, "block_e_adaptive.rds"))
block_e_leave1out <- readRDS(file.path(results_dir, "block_e_leave1out.rds"))

p_break_lambda <- strip_plot_headers(plot_block_d_bar_metric(block_d, metric = "lambda")) +
  labs(y = expression(lambda), x = NULL) +
  theme(
    legend.position = "top",
    plot.margin = margin(5, 5, 5, 5)
  )

p_break_type1 <- strip_plot_headers(plot_block_d_bar_metric(block_d, metric = "type1_05")) +
  labs(y = "Type I error", x = "Number of broken components") +
  theme(
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  )

adaptive_df <- block_e_adaptive %>%
  mutate(
    source = "Adaptive",
    method_label = case_when(
      method == "omnibus_analytical" ~ "Analytic",
      method == "omnibus_mvn" ~ "MVN",
      method == "omnibus_adaptive" ~ "Adaptive+Analytic",
      method == "omnibus_adaptive_mvn" ~ "Adaptive+MVN",
      TRUE ~ method
    )
  ) %>%
  group_by(source, method_label) %>%
  summarise(
    n_obs = n(),
    lambda = mean(lambda, na.rm = TRUE),
    lambda_se = sd(lambda, na.rm = TRUE) / sqrt(n_obs),
    type1_05 = mean(type1_05, na.rm = TRUE),
    type1_05_se = sd(type1_05, na.rm = TRUE) / sqrt(n_obs),
    .groups = "drop"
  ) %>%
  dplyr::select(-n_obs)

leave1_df <- block_e_leave1out %>%
  mutate(
    source = "Leave-one-out",
    method_label = case_when(
      method == "omnibus_minus_acat" ~ "-ACAT",
      method == "omnibus_minus_fisher" ~ "-Fisher",
      method == "omnibus_minus_tfisher" ~ "-TFisher",
      method == "omnibus_minus_minp" ~ "-minP",
      method == "omnibus_minus_stouffer" ~ "-Stouffer",
      method == "omnibus_all" ~ "All",
      TRUE ~ method
    )
  ) %>%
  group_by(source, method_label) %>%
  summarise(
    n_obs = n(),
    lambda = mean(lambda, na.rm = TRUE),
    lambda_se = sd(lambda, na.rm = TRUE) / sqrt(n_obs),
    type1_05 = mean(type1_05, na.rm = TRUE),
    type1_05_se = sd(type1_05, na.rm = TRUE) / sqrt(n_obs),
    .groups = "drop"
  ) %>%
  dplyr::select(-n_obs)

block_e_df <- bind_rows(adaptive_df, leave1_df) %>%
  mutate(
    source = factor(source, levels = c("Adaptive", "Leave-one-out")),
    lambda_se = coalesce(lambda_se, 0),
    type1_05_se = coalesce(type1_05_se, 0)
  )

p_adapt_lambda <-
  ggplot(block_e_df, aes(x = method_label, y = lambda, fill = source)) +
  geom_col(alpha = 0.82, width = 0.72) +
  geom_errorbar(
    aes(ymin = pmax(lambda - lambda_se, 0), ymax = lambda + lambda_se),
    width = 0.14, linewidth = 0.7, color = "black"
  ) +
  facet_wrap(~ source, scales = "free_x", ncol = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#1F77B4", linewidth = 0.7) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 4),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(x = NULL, y = expression(lambda), fill = "Family") +
  block_a_theme +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    strip.text = element_text(size = 13, face = "bold"),
    legend.position = "top",
    legend.text = element_text(size = 13),
    plot.margin = margin(5, 5, 5, 5)
  )

p_adapt_type1 <-
  ggplot(block_e_df, aes(x = method_label, y = type1_05, fill = source)) +
  geom_col(alpha = 0.82, width = 0.72) +
  geom_errorbar(
    aes(ymin = pmax(type1_05 - type1_05_se, 0), ymax = type1_05 + type1_05_se),
    width = 0.14, linewidth = 0.7, color = "black"
  ) +
  facet_wrap(~ source, scales = "free_x", ncol = 1) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.6) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 4),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(x = "Variant", y = "Type I error", fill = "Family") +
  block_a_theme +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    strip.text = element_text(size = 13, face = "bold"),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  )

p_adapt_lambda <- strip_plot_headers(p_adapt_lambda)
p_adapt_type1 <- strip_plot_headers(p_adapt_type1)

combined_plot <-
  (p_break_lambda / p_break_type1 / p_adapt_lambda / p_adapt_type1) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 18))

ggsave(
  filename = file.path(final_fig_dir, "Fig5.png"),
  plot = combined_plot,
  width = 16,
  height = 24,
  dpi = 300,
  bg = "white"
)

message("Fig5 output written to Final/final_main_fig/Fig5.png")
