suppressPackageStartupMessages({
  library(ggplot2)
})

# Focused helper layer for the manuscript simulation figures that depend on the
# broader simulation backend in scripts/all_figs.R.

.repo_dir <- "/Users/nirwantandukar/Documents/Github/MAGCAT"
.final_dir <- file.path(.repo_dir, "Final")
.sim_results_dir <- file.path(.repo_dir, "simulation_results")

.ensure_sim_backend <- local({
  loaded <- FALSE

  function() {
    if (!loaded || !exists("run_block_a", inherits = TRUE) || !exists("run_block_h", inherits = TRUE)) {
      source(file.path(.repo_dir, "scripts", "all_figs.R"))
      loaded <<- TRUE
    }
    invisible(TRUE)
  }
})

refresh_fig3_outputs <- function() {
  .ensure_sim_backend()

  write_block_a_by_m_compact(
    summary_csv = file.path(.sim_results_dir, "block_a_archetype_summary_by_m.csv"),
    results_dir = file.path(.final_dir, "final_main_fig"),
    supp_dir = file.path(.final_dir, "final_sup")
  )

  file.copy(
    from = file.path(.final_dir, "final_main_fig", "block_a_archetype_recovery_by_m_compact.png"),
    to = file.path(.final_dir, "final_main_fig", "Fig3.png"),
    overwrite = TRUE
  )

  write_block_a_by_m_tables(
    summary_csv = file.path(.sim_results_dir, "block_a_archetype_summary_by_m.csv"),
    tables_dir = file.path(.final_dir, "final_sup", "tables")
  )

  message("Fig3 outputs written to Final/final_main_fig and Final/final_sup/tables")
}

refresh_suppfig4_output <- function() {
  dir.create(file.path(.final_dir, "final_sup", "supp_figs"), recursive = TRUE, showWarnings = FALSE)

  file.copy(
    from = file.path(.sim_results_dir, "block_h_power.png"),
    to = file.path(.final_dir, "final_sup", "supp_figs", "SuppFig4_block_h_power.png"),
    overwrite = TRUE
  )

  message("SuppFig4 output written to Final/final_sup/supp_figs/SuppFig4_block_h_power.png")
}

refresh_suppfig5_output <- function() {
  dir.create(file.path(.final_dir, "final_sup", "supp_figs"), recursive = TRUE, showWarnings = FALSE)

  file.copy(
    from = file.path(.sim_results_dir, "block_h_omnibus_regret.png"),
    to = file.path(.final_dir, "final_sup", "supp_figs", "SuppFig5_block_h_omnibus_regret.png"),
    overwrite = TRUE
  )

  message("SuppFig5 output written to Final/final_sup/supp_figs/SuppFig5_block_h_omnibus_regret.png")
}

refresh_suppfig6_omnibus_compare_output <- function() {
  .ensure_sim_backend()
  suppressPackageStartupMessages(library(patchwork))

  block_b <- readRDS(file.path(.sim_results_dir, "block_b.rds"))
  block_c <- readRDS(file.path(.sim_results_dir, "block_c.rds"))

  variant_colors <- c("Combined" = "#377EB8", "Alone" = "#111111")
  strategy_colors <- c(
    "MVN True Cor" = "#377EB8",
    "MVN Imputed (0)" = "#4DAF4A"
  )

  b_df <- block_b %>%
    dplyr::filter(
      method %in% c("omnibus_combined", "omnibus_alone"),
      calibration %in% c("mvn", "mvn_alone")
    ) %>%
    dplyr::mutate(
      cor_structure = factor(cor_structure, levels = c("LD_moderate", "LD_strong", "LD_independent")),
      variant = factor(dplyr::case_when(
        method == "omnibus_combined" ~ "Combined",
        method == "omnibus_alone" ~ "Alone",
        TRUE ~ NA_character_
      ), levels = c("Combined", "Alone")),
      lambda_dist = abs(lambda - 1),
      type1_dist = abs(type1_05 - 0.05)
    )

  c_df <- block_c %>%
    dplyr::filter(
      method %in% c("omnibus_combined", "omnibus_alone"),
      strategy %in% c("mvn_true_cor", "mvn_imputed_cor")
    ) %>%
    dplyr::mutate(
      pathway_size = as.integer(pathway_size),
      facet_label = factor(paste0("m=", pathway_size),
                           levels = paste0("m=", sort(unique(pathway_size)))),
      variant = factor(dplyr::case_when(
        method == "omnibus_combined" ~ "Combined",
        method == "omnibus_alone" ~ "Alone",
        TRUE ~ NA_character_
      ), levels = c("Combined", "Alone")),
      strategy_label = factor(dplyr::case_when(
        strategy == "mvn_true_cor" ~ "MVN True Cor",
        strategy == "mvn_imputed_cor" ~ "MVN Imputed (0)",
        TRUE ~ strategy
      ), levels = c("MVN True Cor", "MVN Imputed (0)")),
      lambda_dist = abs(lambda - 1),
      type1_dist = abs(type1_05 - 0.05)
    )

  pA <- ggplot(b_df, aes(x = variant, y = lambda_dist, fill = variant)) +
    geom_col(width = 0.72, alpha = 0.9) +
    geom_text(aes(label = sprintf("%.3f", lambda_dist)), vjust = -0.25, size = 2.8) +
    facet_grid(pathway_size ~ cor_structure,
               labeller = labeller(pathway_size = function(x) paste0("m=", x))) +
    scale_fill_manual(values = variant_colors) +
    labs(
      title = "A. Null Calibration: Distance to lambda = 1",
      subtitle = "Smaller is better",
      x = "Omnibus variant",
      y = expression("|" * lambda - 1 * "|")
    ) +
    block_a_theme +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 10),
      strip.text = element_text(size = 11, face = "bold")
    )

  pB <- ggplot(b_df, aes(x = variant, y = type1_dist, fill = variant)) +
    geom_col(width = 0.72, alpha = 0.9) +
    geom_text(aes(label = sprintf("%.3f", type1_dist)), vjust = -0.25, size = 2.8) +
    facet_grid(pathway_size ~ cor_structure,
               labeller = labeller(pathway_size = function(x) paste0("m=", x))) +
    scale_fill_manual(values = variant_colors) +
    labs(
      title = "B. Null Calibration: Distance to Type I = 0.05",
      subtitle = "Smaller is better",
      x = "Omnibus variant",
      y = expression("|" * "Type I error" - 0.05 * "|")
    ) +
    block_a_theme +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 10),
      strip.text = element_text(size = 11, face = "bold")
    )

  pC <- ggplot(c_df, aes(x = missing_frac, y = lambda_dist,
                         color = strategy_label, linetype = variant, shape = variant)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    facet_wrap(~ facet_label, ncol = 3) +
    scale_color_manual(values = strategy_colors) +
    scale_linetype_manual(values = c("Combined" = "solid", "Alone" = "dashed")) +
    scale_shape_manual(values = c("Combined" = 16, "Alone" = 17)) +
    labs(
      title = "C. Missing Correlation: Distance to lambda = 1",
      subtitle = "Smaller is better",
      x = "Fraction missing",
      y = expression("|" * lambda - 1 * "|"),
      color = "Strategy",
      linetype = "Omnibus",
      shape = "Omnibus"
    ) +
    block_a_theme +
    theme(
      legend.position = "top",
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 10),
      strip.text = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 9)
    )

  pD <- ggplot(c_df, aes(x = missing_frac, y = type1_dist,
                         color = strategy_label, linetype = variant, shape = variant)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    facet_wrap(~ facet_label, ncol = 3) +
    scale_color_manual(values = strategy_colors) +
    scale_linetype_manual(values = c("Combined" = "solid", "Alone" = "dashed")) +
    scale_shape_manual(values = c("Combined" = 16, "Alone" = 17)) +
    labs(
      title = "D. Missing Correlation: Distance to Type I = 0.05",
      subtitle = "Smaller is better",
      x = "Fraction missing",
      y = expression("|" * "Type I error" - 0.05 * "|"),
      color = "Strategy",
      linetype = "Omnibus",
      shape = "Omnibus"
    ) +
    block_a_theme +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 10),
      strip.text = element_text(size = 11, face = "bold")
    )

  combined_plot <- (pA / pB) | (pC / pD)
  combined_plot <- combined_plot + patchwork::plot_annotation(
    title = "Omnibus Combined vs Omnibus Alone",
    subtitle = "Comparison by distance to ideal calibration targets",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    )
  )

  out_dir <- file.path(.final_dir, "final_sup", "supp_figs")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(out_dir, "SuppFig6_omnibus_combined_vs_alone.png")
  ggsave(out_path, combined_plot, width = 20, height = 12, dpi = 300, bg = "white")

  message("SuppFig6 output written to Final/final_sup/supp_figs/SuppFig6_omnibus_combined_vs_alone.png")
}

refresh_suppfig7_stress_diag_compare_output <- function() {
  .ensure_sim_backend()
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(patchwork)
  })

  block_d <- readRDS(file.path(.sim_results_dir, "block_d.rds"))
  block_e_adaptive <- readRDS(file.path(.sim_results_dir, "block_e_adaptive.rds"))
  block_e_leave1out <- readRDS(file.path(.sim_results_dir, "block_e_leave1out.rds"))

  stress_df <- block_d %>%
    dplyr::filter(method %in% c("omnibus_mvn_combined", "omnibus_mvn_alone")) %>%
    dplyr::mutate(
      cor_structure = factor(
        dplyr::case_when(
          rho == 0 ~ "LD_independent",
          rho == 0.3 ~ "LD_moderate",
          rho == 0.7 ~ "LD_strong",
          TRUE ~ paste0("rho=", rho)
        ),
        levels = c("LD_independent", "LD_moderate", "LD_strong")
      ),
      lambda_dist = abs(lambda - 1),
      type1_dist = abs(type1_05 - 0.05)
    ) %>%
    dplyr::select(pathway_size, cor_structure, n_broken, method, lambda_dist, type1_dist) %>%
    tidyr::pivot_wider(names_from = method, values_from = c(lambda_dist, type1_dist)) %>%
    dplyr::mutate(
      delta_lambda = lambda_dist_omnibus_mvn_combined - lambda_dist_omnibus_mvn_alone,
      delta_type1 = type1_dist_omnibus_mvn_combined - type1_dist_omnibus_mvn_alone
    )

  adaptive_df <- block_e_adaptive %>%
    dplyr::filter(method %in% c(
      "omnibus_mvn_combined", "omnibus_mvn_alone",
      "omnibus_adaptive_mvn_combined", "omnibus_adaptive_mvn_alone"
    )) %>%
    dplyr::mutate(
      comparison = dplyr::case_when(
        method %in% c("omnibus_mvn_combined", "omnibus_mvn_alone") ~ "Full omnibus",
        method %in% c("omnibus_adaptive_mvn_combined", "omnibus_adaptive_mvn_alone") ~ "Adaptive omnibus",
        TRUE ~ method
      ),
      mode = dplyr::case_when(
        grepl("_combined$", method) ~ "combined",
        grepl("_alone$", method) ~ "alone",
        TRUE ~ NA_character_
      ),
      setting = factor(
        paste0(
          "m=", pathway_size, "\n",
          dplyr::case_when(
            rho == 0 ~ "LD_independent",
            rho == 0.3 ~ "LD_moderate",
            rho == 0.7 ~ "LD_strong",
            TRUE ~ paste0("rho=", rho)
          )
        ),
        levels = c(
          "m=5\nLD_independent", "m=5\nLD_moderate", "m=5\nLD_strong",
          "m=25\nLD_independent", "m=25\nLD_moderate", "m=25\nLD_strong",
          "m=50\nLD_independent", "m=50\nLD_moderate", "m=50\nLD_strong"
        )
      ),
      lambda_dist = abs(lambda - 1),
      type1_dist = abs(type1_05 - 0.05)
    ) %>%
    dplyr::select(comparison, setting, mode, lambda_dist, type1_dist) %>%
    tidyr::pivot_wider(names_from = mode, values_from = c(lambda_dist, type1_dist)) %>%
    dplyr::mutate(
      delta_lambda = lambda_dist_combined - lambda_dist_alone,
      delta_type1 = type1_dist_combined - type1_dist_alone
    ) %>%
    tidyr::pivot_longer(
      cols = c(delta_lambda, delta_type1),
      names_to = "metric",
      values_to = "delta"
    ) %>%
    dplyr::mutate(
      metric = factor(metric, levels = c("delta_lambda", "delta_type1"),
                      labels = c(expression("|" * lambda - 1 * "|"),
                                 expression("|" * "Type I" - 0.05 * "|"))),
      comparison = factor(comparison, levels = c("Full omnibus", "Adaptive omnibus"))
    )

  leave1_df <- block_e_leave1out %>%
    dplyr::mutate(
      setting = factor(
        paste0("m=", pathway_size, "\n", cor_structure),
        levels = c(
          "m=5\nLD_independent", "m=5\nLD_moderate", "m=5\nLD_strong",
          "m=25\nLD_independent", "m=25\nLD_moderate", "m=25\nLD_strong",
          "m=50\nLD_independent", "m=50\nLD_moderate", "m=50\nLD_strong"
        )
      ),
      variant = dplyr::case_when(
        method == "omnibus_all" ~ "All",
        method == "omnibus_minus_acat" ~ "-ACAT",
        method == "omnibus_minus_fisher" ~ "-Fisher",
        method == "omnibus_minus_tfisher" ~ "-TFisher",
        method == "omnibus_minus_minp" ~ "-minP",
        method == "omnibus_minus_stouffer" ~ "-Stouffer",
        TRUE ~ method
      ),
      lambda_dist = abs(lambda - 1),
      type1_dist = abs(type1_05 - 0.05)
    ) %>%
    dplyr::select(setting, variant, calibration, lambda_dist, type1_dist) %>%
    tidyr::pivot_wider(names_from = calibration, values_from = c(lambda_dist, type1_dist)) %>%
    dplyr::mutate(
      delta_lambda = lambda_dist_combined - lambda_dist_alone,
      delta_type1 = type1_dist_combined - type1_dist_alone
    ) %>%
    tidyr::pivot_longer(
      cols = c(delta_lambda, delta_type1),
      names_to = "metric",
      values_to = "delta"
    ) %>%
    dplyr::mutate(
      metric = factor(metric, levels = c("delta_lambda", "delta_type1"),
                      labels = c(expression("|" * lambda - 1 * "|"),
                                 expression("|" * "Type I" - 0.05 * "|"))),
      variant = factor(variant, levels = c("All", "-ACAT", "-Fisher", "-TFisher", "-minP", "-Stouffer"))
    )

  all_delta_vals <- c(
    stress_df$delta_lambda, stress_df$delta_type1,
    adaptive_df$delta, leave1_df$delta
  )
  fill_lim <- max(abs(all_delta_vals), na.rm = TRUE)

  pA <- ggplot(
    stress_df,
    aes(x = factor(n_broken), y = factor(pathway_size), fill = delta_lambda)
  ) +
    geom_tile(color = "white", linewidth = 0.4) +
    facet_wrap(~ cor_structure, nrow = 1) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-fill_lim, fill_lim)
    ) +
    labs(
      title = "A. Stress Test: Delta in |lambda - 1|",
      subtitle = "Combined - Alone; negative means Combined is closer to ideal",
      x = "Number of broken components",
      y = "Pathway size",
      fill = "Delta"
    ) +
    block_a_theme +
    theme(
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 10),
      strip.text = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 9)
    )

  pB <- ggplot(
    stress_df,
    aes(x = factor(n_broken), y = factor(pathway_size), fill = delta_type1)
  ) +
    geom_tile(color = "white", linewidth = 0.4) +
    facet_wrap(~ cor_structure, nrow = 1) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-fill_lim, fill_lim)
    ) +
    labs(
      title = "B. Stress Test: Delta in |Type I - 0.05|",
      subtitle = "Combined - Alone; negative means Combined is closer to ideal",
      x = "Number of broken components",
      y = "Pathway size",
      fill = "Delta"
    ) +
    block_a_theme +
    theme(
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 10),
      strip.text = element_text(size = 11, face = "bold"),
      legend.position = "none"
    )

  pC <- ggplot(
    adaptive_df,
    aes(x = setting, y = comparison, fill = delta)
  ) +
    geom_tile(color = "white", linewidth = 0.4) +
    facet_wrap(~ metric, ncol = 1, scales = "free_y", labeller = label_parsed) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-fill_lim, fill_lim)
    ) +
    labs(
      title = "C. Adaptive Diagnostics: Combined - Alone",
      subtitle = "Negative means Combined is closer to ideal",
      x = "Setting",
      y = NULL,
      fill = "Delta"
    ) +
    block_a_theme +
    theme(
      axis.text.x = element_text(size = 8, angle = 35, hjust = 1),
      axis.text.y = element_text(size = 10),
      strip.text = element_text(size = 11, face = "bold"),
      legend.position = "none"
    )

  pD <- ggplot(
    leave1_df,
    aes(x = setting, y = variant, fill = delta)
  ) +
    geom_tile(color = "white", linewidth = 0.35) +
    facet_wrap(~ metric, ncol = 1, scales = "free_y", labeller = label_parsed) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-fill_lim, fill_lim)
    ) +
    labs(
      title = "D. Leave-one-out Diagnostics: Combined - Alone",
      subtitle = "Negative means Combined is closer to ideal",
      x = "Setting",
      y = NULL,
      fill = "Delta"
    ) +
    block_a_theme +
    theme(
      axis.text.x = element_text(size = 8, angle = 35, hjust = 1),
      axis.text.y = element_text(size = 9),
      strip.text = element_text(size = 11, face = "bold"),
      legend.position = "none"
    )

  combined_plot <- (pA / pB) | (pC / pD)
  combined_plot <- combined_plot + patchwork::plot_annotation(
    title = "Where Does Omnibus Combined Beat Omnibus Alone?",
    subtitle = "Signed distance difference: Combined - Alone. Blue favors Combined; red favors Alone.",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    )
  )

  out_dir <- file.path(.final_dir, "final_sup", "supp_figs")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(out_dir, "SuppFig7_stress_diag_combined_vs_alone.png")
  ggsave(out_path, combined_plot, width = 20, height = 14, dpi = 300, bg = "white")

  message("SuppFig7 output written to Final/final_sup/supp_figs/SuppFig7_stress_diag_combined_vs_alone.png")
}

refresh_simulation_figure_sources <- function(run_block_a_flag = TRUE,
                                              run_block_h_flag = TRUE) {
  .ensure_sim_backend()

  if (isTRUE(run_block_a_flag)) {
    message("Refreshing Block A with pathway-size scaling and current oTFisher code")
    run_block_a()
  }

  if (isTRUE(run_block_h_flag)) {
    message("Refreshing Block H with pathway-size scaling and current oTFisher code")
    run_block_h()
  }
}

refresh_final_simulation_figures <- function(run_block_a_flag = TRUE,
                                             run_block_h_flag = TRUE) {
  refresh_simulation_figure_sources(
    run_block_a_flag = run_block_a_flag,
    run_block_h_flag = run_block_h_flag
  )

  message("Refreshing Final figure copies")
  refresh_fig3_outputs()
  refresh_suppfig4_output()
  refresh_suppfig5_output()
  message("Done refreshing simulation-derived Final figures")
}
