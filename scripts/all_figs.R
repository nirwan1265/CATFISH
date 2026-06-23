suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(MASS)
})

# ==============================================================================
# ALL FIGURES
#   Block A: Archetype recovery simulation at the pathway level
# ==============================================================================

if (requireNamespace("TFisher", quietly = TRUE)) {
  use_tfisher_pkg <- TRUE
} else {
  use_tfisher_pkg <- FALSE
  warning(
    "TFisher package not installed. Using a fallback approximation for TFisher."
  )
}

BLOCK_A_CONFIG <- list(
  seed = 20260312L,
  m = 20L,
  rho = 0.20,
  n_reps = 500L,
  n_null_cal = 2000L,
  alpha = 0.05,
  tau_grid = c(0.10, 0.05, 0.01),
  min_p = 1e-15,
  output_dir = file.path("Figures", "Fig2.Simulation_validation", "Ind_figs"),
  results_dir = "simulation_results",
  save_to_ind_figs = TRUE
)

BLOCK_H_CONFIG <- list(
  seed = 20260317L,
  m_values = c(5L, 20L, 50L),
  rho = 0.20,
  effect_scales = c(0.5, 0.8, 1.0, 1.2, 1.5, 1.8),
  n_reps = 250L,
  n_null_cal = 1500L,
  alpha = 0.05,
  tau_grid = c(0.10, 0.05, 0.01),
  results_dir = "simulation_results"
)

# Use the same plotting theme style used in the diagnostics scripts.
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

block_a_theme <- plot_theme

block_a_method_order <- c("ACAT", "Fisher", "TFisher", "minP", "Stouffer", "Omnibus")
block_a_archetype_order <- c(
  "Archetype I - SDA",
  "Archetype II - CME",
  "Archetype III - DPS",
  "Archetype IV - HDS",
  "Archetype V - SGP"
)

make_cor_exchangeable <- function(m, rho = 0.20) {
  sigma <- matrix(rho, nrow = m, ncol = m)
  diag(sigma) <- 1
  sigma
}

ensure_pd <- function(sigma, tol = 1e-8) {
  eig <- eigen(sigma, symmetric = TRUE)
  vals <- pmax(eig$values, tol)
  sigma_pd <- eig$vectors %*% diag(vals, nrow = length(vals)) %*% t(eig$vectors)
  sigma_pd <- (sigma_pd + t(sigma_pd)) / 2

  d <- sqrt(diag(sigma_pd))
  sigma_pd <- sweep(sigma_pd, 1, d, "/")
  sigma_pd <- sweep(sigma_pd, 2, d, "/")
  diag(sigma_pd) <- 1
  sigma_pd
}

sanitize_p <- function(p_vals, min_p = BLOCK_A_CONFIG$min_p) {
  p_vals <- p_vals[is.finite(p_vals)]
  pmax(pmin(p_vals, 1 - min_p), min_p)
}

z_to_p_onesided <- function(z_vals, min_p = BLOCK_A_CONFIG$min_p) {
  sanitize_p(stats::pnorm(z_vals, lower.tail = FALSE), min_p = min_p)
}

compute_acat <- function(p_vals, min_p = BLOCK_A_CONFIG$min_p) {
  p_use <- sanitize_p(p_vals, min_p = min_p)
  if (!length(p_use)) {
    return(1)
  }

  if (requireNamespace("ACAT", quietly = TRUE)) {
    p_acat <- ACAT::ACAT(Pvals = p_use)
  } else {
    tan_vals <- tan((0.5 - p_use) * pi)
    p_acat <- 0.5 - atan(mean(tan_vals)) / pi
  }

  max(min(p_acat, 1 - min_p), min_p)
}

compute_fisher <- function(p_vals, min_p = BLOCK_A_CONFIG$min_p) {
  p_use <- sanitize_p(p_vals, min_p = min_p)
  stat <- -2 * sum(log(p_use))
  p_fisher <- stats::pchisq(stat, df = 2 * length(p_use), lower.tail = FALSE)
  max(min(p_fisher, 1 - min_p), min_p)
}

compute_stouffer <- function(z_vals, alternative = "greater",
                             min_p = BLOCK_A_CONFIG$min_p) {
  z_stouffer <- sum(z_vals) / sqrt(length(z_vals))

  if (alternative == "greater") {
    p_val <- stats::pnorm(z_stouffer, lower.tail = FALSE)
  } else if (alternative == "less") {
    p_val <- stats::pnorm(z_stouffer, lower.tail = TRUE)
  } else {
    p_val <- 2 * stats::pnorm(-abs(z_stouffer))
  }

  max(min(p_val, 1 - min_p), min_p)
}

compute_minp <- function(p_vals, min_p = BLOCK_A_CONFIG$min_p) {
  p_use <- sanitize_p(p_vals, min_p = min_p)
  p_min <- min(p_use)
  p_minp <- 1 - (1 - p_min)^length(p_use)
  max(min(p_minp, 1 - min_p), min_p)
}

compute_tfisher <- function(p_vals, tau_grid = BLOCK_A_CONFIG$tau_grid,
                            min_p = BLOCK_A_CONFIG$min_p) {
  p_use <- sanitize_p(p_vals, min_p = min_p)
  tau_grid <- sort(unique(tau_grid[tau_grid > 0 & tau_grid < 1]), decreasing = FALSE)

  if (!length(p_use) || !length(tau_grid)) {
    return(1)
  }

  if (use_tfisher_pkg) {
    tfisher_fit <- tryCatch(
      {
        stat_out <- TFisher::stat.soft.omni(p = p_use, TAU1 = tau_grid, M = NULL)
        stat_val <- as.numeric(stat_out$omni)
        as.numeric(TFisher::p.soft.omni(q = stat_val, n = length(p_use),
                                        TAU1 = tau_grid, M = NULL))
      },
      error = function(e) NA_real_
    )

    if (is.finite(tfisher_fit) && !is.na(tfisher_fit)) {
      return(max(min(tfisher_fit, 1 - min_p), min_p))
    }
  }

  tfisher_ps <- vapply(tau_grid, function(tau) {
    keep <- p_use[p_use < tau]
    if (!length(keep)) {
      return(1)
    }
    stat <- -2 * sum(log(keep))
    stats::pchisq(stat, df = 2 * length(keep), lower.tail = FALSE)
  }, numeric(1))

  p_tfisher <- min(tfisher_ps) * length(tau_grid)
  max(min(p_tfisher, 1 - min_p), min_p)
}

compute_components_analytic <- function(z_vals, tau_grid = BLOCK_A_CONFIG$tau_grid,
                                        stouffer_alt = "greater",
                                        min_p = BLOCK_A_CONFIG$min_p) {
  p_vals <- z_to_p_onesided(z_vals, min_p = min_p)

  component_p <- c(
    ACAT = compute_acat(p_vals, min_p = min_p),
    Fisher = compute_fisher(p_vals, min_p = min_p),
    TFisher = compute_tfisher(p_vals, tau_grid = tau_grid, min_p = min_p),
    minP = compute_minp(p_vals, min_p = min_p),
    Stouffer = compute_stouffer(z_vals, alternative = stouffer_alt, min_p = min_p)
  )

  c(component_p, Omnibus = compute_acat(component_p, min_p = min_p))
}

precompute_null_distribution <- function(m, sigma, n_null = BLOCK_A_CONFIG$n_null_cal,
                                         seed = BLOCK_A_CONFIG$seed,
                                         tau_grid = BLOCK_A_CONFIG$tau_grid) {
  set.seed(seed)
  sigma_pd <- ensure_pd(sigma)

  null_stats <- matrix(NA_real_, nrow = n_null, ncol = length(block_a_method_order))
  colnames(null_stats) <- block_a_method_order

  for (i in seq_len(n_null)) {
    z_null <- as.numeric(MASS::mvrnorm(1, mu = rep(0, m), Sigma = sigma_pd))
    null_stats[i, ] <- compute_components_analytic(z_null, tau_grid = tau_grid)
  }

  null_stats
}

calibrate_from_null <- function(obs_stats, null_stats) {
  calibrated <- vapply(names(obs_stats), function(method) {
    obs_val <- obs_stats[[method]]
    null_vals <- null_stats[, method]
    (1 + sum(null_vals <= obs_val, na.rm = TRUE)) / (sum(!is.na(null_vals)) + 1)
  }, numeric(1))

  calibrated[names(obs_stats)]
}

make_mu_from_layers <- function(m, layers) {
  mu <- rep(0, m)
  available <- sample.int(m, size = m, replace = FALSE)
  start <- 1L

  for (layer in layers) {
    n_layer <- min(as.integer(layer$n), m - start + 1L)
    if (n_layer <= 0) {
      next
    }

    idx <- available[start:(start + n_layer - 1L)]
    mu[idx] <- layer$delta
    start <- start + n_layer
  }

  mu
}

block_a_archetypes <- list(
  list(
    code = "SDA",
    label = "Archetype I - SDA",
    layers = list(list(n = 3L, delta = 3.0)),
    description = "Few strong genes; sparse driver pattern."
  ),
  list(
    code = "CME",
    label = "Archetype II - CME",
    layers = list(list(n = 12L, delta = 1.1)),
    description = "Many moderate genes with shared direction."
  ),
  list(
    code = "DPS",
    label = "Archetype III - DPS",
    layers = list(list(n = 18L, delta = 0.55)),
    description = "Diffuse weak positive shift across most genes."
  ),
  list(
    code = "HDS",
    label = "Archetype IV - HDS",
    layers = list(list(n = 2L, delta = 1.8), list(n = 12L, delta = 1.10)),
    description = "Two strong drivers plus a moderate support layer."
  ),
  list(
    code = "SGP",
    label = "Archetype V - SGP",
    layers = list(list(n = 1L, delta = 4.0)),
    description = "Single-gene proxy boundary case."
  )
)

simulate_archetype_once <- function(archetype, sigma, null_stats,
                                    tau_grid = BLOCK_A_CONFIG$tau_grid) {
  mu <- make_mu_from_layers(m = nrow(sigma), layers = archetype$layers)
  z_alt <- as.numeric(MASS::mvrnorm(1, mu = mu, Sigma = sigma))
  obs_stats <- compute_components_analytic(z_alt, tau_grid = tau_grid)
  p_cal <- calibrate_from_null(obs_stats, null_stats)
  method_rank <- rank(p_cal, ties.method = "average")
  top_method <- names(p_cal)[which.min(p_cal)]

  tibble(
    archetype = archetype$label,
    method = names(p_cal),
    p_value = as.numeric(p_cal),
    rank = as.numeric(method_rank[names(p_cal)]),
    is_top = names(p_cal) == top_method
  )
}

scale_archetype_layers <- function(archetype, effect_scale) {
  archetype_scaled <- archetype
  archetype_scaled$layers <- lapply(archetype$layers, function(layer) {
    layer$delta <- as.numeric(layer$delta) * effect_scale
    layer
  })
  archetype_scaled
}

simulate_archetype_once_scaled <- function(archetype, sigma, null_stats,
                                           effect_scale,
                                           tau_grid = BLOCK_A_CONFIG$tau_grid) {
  arch_scaled <- scale_archetype_layers(archetype, effect_scale = effect_scale)
  mu <- make_mu_from_layers(m = nrow(sigma), layers = arch_scaled$layers)
  z_alt <- as.numeric(MASS::mvrnorm(1, mu = mu, Sigma = sigma))
  obs_stats <- compute_components_analytic(z_alt, tau_grid = tau_grid)
  p_cal <- calibrate_from_null(obs_stats, null_stats)
  method_rank <- rank(p_cal, ties.method = "average")
  top_method <- names(p_cal)[which.min(p_cal)]

  tibble(
    archetype = archetype$label,
    effect_scale = effect_scale,
    method = names(p_cal),
    p_value = as.numeric(p_cal),
    rank = as.numeric(method_rank[names(p_cal)]),
    is_top = names(p_cal) == top_method
  )
}

summarize_archetype_results <- function(results_long,
                                        alpha = BLOCK_A_CONFIG$alpha) {
  results_long %>%
    group_by(archetype, method) %>%
    summarize(
      power = mean(p_value < alpha),
      mean_rank = mean(rank),
      top_probability = mean(is_top),
      median_p = stats::median(p_value),
      .groups = "drop"
    ) %>%
    mutate(
      archetype = factor(archetype, levels = block_a_archetype_order),
      method = factor(method, levels = block_a_method_order)
    ) %>%
    arrange(archetype, method)
}

plot_metric_heatmap <- function(summary_df, value_col, title, fill_limits,
                                low_color, high_color, reverse = FALSE,
                                digits = 2,
                                legend_title = NULL) {
  values <- summary_df[[value_col]]
  labels <- sprintf(paste0("%.", digits, "f"), values)
  if (is.null(legend_title)) {
    legend_title <- title
  }

  p <- ggplot(summary_df, aes(x = archetype, y = method, fill = .data[[value_col]])) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(label = labels, size = 4.2, fontface = "bold") +
    labs(x = NULL, y = NULL, title = title, fill = legend_title) +
    block_a_theme +
    theme(
      axis.text.x = element_text(angle = 20, hjust = 1),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(size = 12)
    )

  if (reverse) {
    p + scale_fill_gradient(
      limits = fill_limits,
      low = high_color,
      high = low_color
    )
  } else {
    p + scale_fill_gradient(
      limits = fill_limits,
      low = low_color,
      high = high_color
    )
  }
}

block_a_method_colors <- c(
  "ACAT" = "#E41A1C",
  "Fisher" = "#377EB8",
  "TFisher" = "#4DAF4A",
  "minP" = "#984EA3",
  "Stouffer" = "#FF7F00",
  "Omnibus" = "#222222"
)

compact_archetype_labels <- c(
  "Archetype I - SDA" = "Type I",
  "Archetype II - CME" = "Type II",
  "Archetype III - DPS" = "Type III",
  "Archetype IV - HDS" = "Type IV",
  "Archetype V - SGP" = "Type V"
)

build_block_a_by_m_compact_summary <- function(summary_by_m) {
  power_best <- summary_by_m %>%
    group_by(pathway_size, archetype) %>%
    arrange(desc(power), method) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(
      pathway_size,
      archetype,
      metric = "Power",
      best_method = method,
      value = power,
      label = sprintf("%.2f\n%s", power, method)
    )

  rank_best <- summary_by_m %>%
    group_by(pathway_size, archetype) %>%
    arrange(mean_rank, method) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(
      pathway_size,
      archetype,
      metric = "Mean rank",
      best_method = method,
      value = mean_rank,
      label = sprintf("%.2f\n%s", mean_rank, method)
    )

  top_best <- summary_by_m %>%
    group_by(pathway_size, archetype) %>%
    arrange(desc(top_probability), method) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(
      pathway_size,
      archetype,
      metric = "Top method",
      best_method = method,
      value = top_probability,
      label = sprintf("%.2f\n%s", top_probability, method)
    )

  bind_rows(power_best, rank_best, top_best) %>%
    mutate(
      pathway_size = factor(pathway_size, levels = c(50, 20, 5), labels = c("G=50", "G=20", "G=5")),
      archetype_short = factor(
        compact_archetype_labels[as.character(archetype)],
        levels = c("Type I", "Type II", "Type III", "Type IV", "Type V")
      ),
      metric = factor(metric, levels = c("Power", "Mean rank", "Top method")),
      x_group = interaction(archetype_short, metric, lex.order = TRUE),
      best_method = factor(best_method, levels = names(block_a_method_colors))
    )
}

plot_block_a_by_m_compact <- function(summary_by_m) {
  compact_df <- build_block_a_by_m_compact_summary(summary_by_m)

  ggplot(compact_df, aes(x = x_group, y = pathway_size, fill = best_method)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = label), size = 4.2, fontface = "bold", lineheight = 0.95) +
    facet_grid(~ archetype_short, scales = "free_x", space = "free_x", switch = "x") +
    scale_fill_manual(values = block_a_method_colors, drop = FALSE) +
    labs(
      x = NULL,
      y = NULL,
      fill = "Winning method",
      title = "Pathway-size sensitivity of archetype verification",
      subtitle = "Each cell shows the best summary for that G x archetype setting:\nPower = highest power, Mean rank = lowest mean rank, Top method = highest top-method probability"
    ) +
    scale_x_discrete(labels = rep(c("Power", "Mean rank", "Top method"), 5)) +
    block_a_theme +
    theme(
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.x = element_text(face = "bold", size = 18),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
      axis.text.y = element_text(face = "bold", size = 16),
      panel.spacing.x = unit(0.8, "lines"),
      legend.position = "top",
      legend.title = element_text(face = "bold", size = 16),
      legend.text = element_text(size = 14)
    )
}

write_block_a_by_m_compact <- function(
    summary_csv = file.path("simulation_results", "block_a_archetype_summary_by_m.csv"),
    results_dir = "simulation_results",
    supp_dir = file.path("fig", "supp")) {
  if (!file.exists(summary_csv)) {
    stop("Could not find summary file: ", summary_csv, call. = FALSE)
  }

  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

  summary_by_m <- read.csv(summary_csv, stringsAsFactors = FALSE)
  compact_plot <- plot_block_a_by_m_compact(summary_by_m)

  results_path <- file.path(results_dir, "block_a_archetype_recovery_by_m_compact.png")
  supp_path <- file.path(supp_dir, "SuppFig2_block_a_archetype_recovery_by_m_compact.png")

  ggsave(results_path, compact_plot, width = 18, height = 7.5, dpi = 300, bg = "white")
  ggsave(supp_path, compact_plot, width = 18, height = 7.5, dpi = 300, bg = "white")

  invisible(list(
    plot = compact_plot,
    results_path = results_path,
    supp_path = supp_path
  ))
}

write_block_a_by_m_tables <- function(
    summary_csv = file.path("simulation_results", "block_a_archetype_summary_by_m.csv"),
    tables_dir = file.path("fig", "supp", "tables")) {
  if (!file.exists(summary_csv)) {
    stop("Could not find summary file: ", summary_csv, call. = FALSE)
  }

  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

  summary_by_m <- read.csv(summary_csv, stringsAsFactors = FALSE)
  compact_df <- build_block_a_by_m_compact_summary(summary_by_m) %>%
    dplyr::select(pathway_size, archetype_short, metric, best_method, value) %>%
    dplyr::mutate(
      value = dplyr::case_when(
        metric == "Mean rank" ~ sprintf("%.2f", value),
        TRUE ~ sprintf("%.2f", value)
      )
    )

  compact_wide <- compact_df %>%
    tidyr::pivot_wider(
      names_from = metric,
      values_from = c(best_method, value),
      names_glue = "{metric}_{.value}"
    ) %>%
    dplyr::rename(
      G = pathway_size,
      Archetype = archetype_short
    ) %>%
    dplyr::arrange(factor(G, levels = c("G=5", "G=20", "G=50")), Archetype)

  compact_csv <- file.path(tables_dir, "TableS_Archetype_Verification_ByG_Compact.csv")
  full_csv <- file.path(tables_dir, "TableS_Archetype_Verification_ByG_Full.csv")
  latex_path <- file.path(tables_dir, "TableS_Archetype_Verification_ByG_Compact.tex")

  utils::write.csv(compact_wide, compact_csv, row.names = FALSE)

  full_export <- summary_by_m %>%
    dplyr::mutate(
      G = paste0("G=", pathway_size),
      Archetype = compact_archetype_labels[archetype]
    ) %>%
    dplyr::select(G, Archetype, method, power, mean_rank, top_probability, median_p)
  utils::write.csv(full_export, full_csv, row.names = FALSE)

  latex_lines <- c(
    "\\begin{table*}[!t]",
    "\\centering",
    "\\caption{Compact summary of fixed-parameter archetype verification across pathway sizes. For each pathway size (G) and archetype, entries report the best-performing method for power, mean rank, and top-method probability, together with the corresponding value.}",
    "\\label{tab:supp_archetype_verification_byG}",
    "\\small",
    "\\begin{tabular}{llcccccc}",
    "\\hline",
    "G & Archetype & Power method & Power & Rank method & Mean rank & Top method & Top prob.\\\\",
    "\\hline"
  )

  for (i in seq_len(nrow(compact_wide))) {
    row <- compact_wide[i, ]
    latex_lines <- c(
      latex_lines,
      sprintf(
        "%s & %s & %s & %s & %s & %s & %s & %s\\\\",
        row$G,
        row$Archetype,
        row$Power_best_method,
        row$Power_value,
        row$`Mean rank_best_method`,
        row$`Mean rank_value`,
        row$`Top method_best_method`,
        row$`Top method_value`
      )
    )
  }

  latex_lines <- c(
    latex_lines,
    "\\hline",
    "\\end{tabular}",
    "\\end{table*}"
  )

  writeLines(latex_lines, latex_path)

  invisible(list(
    compact_csv = compact_csv,
    full_csv = full_csv,
    latex = latex_path
  ))
}

run_block_a <- function(config = BLOCK_A_CONFIG) {
  dir.create(config$results_dir, recursive = TRUE, showWarnings = FALSE)
  if (isTRUE(config$save_to_ind_figs)) {
    dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  sigma <- ensure_pd(make_cor_exchangeable(config$m, config$rho))

  message(
    sprintf(
      "Block A: precomputing null calibration (m=%d, rho=%.2f, B=%d)",
      config$m, config$rho, config$n_null_cal
    )
  )
  null_stats <- precompute_null_distribution(
    m = config$m,
    sigma = sigma,
    n_null = config$n_null_cal,
    seed = config$seed,
    tau_grid = config$tau_grid
  )

  message(sprintf("Block A: running %d replicates per archetype", config$n_reps))
  results_long <- bind_rows(lapply(block_a_archetypes, function(archetype) {
    bind_rows(lapply(seq_len(config$n_reps), function(rep_idx) {
      simulate_archetype_once(
        archetype = archetype,
        sigma = sigma,
        null_stats = null_stats,
        tau_grid = config$tau_grid
      ) %>%
        mutate(replicate = rep_idx)
    }))
  })) %>%
    mutate(
      archetype = factor(archetype, levels = block_a_archetype_order),
      method = factor(method, levels = block_a_method_order)
    )

  summary_df <- summarize_archetype_results(results_long, alpha = config$alpha)

  power_plot <- plot_metric_heatmap(
    summary_df = summary_df,
    value_col = "power",
    title = sprintf("Power at alpha = %.2f", config$alpha),
    fill_limits = c(0, 1),
    low_color = "#F7FBFF",
    high_color = "#08519C",
    legend_title = "Power"
  )

  rank_plot <- plot_metric_heatmap(
    summary_df = summary_df,
    value_col = "mean_rank",
    title = "Mean rank (lower is better)",
    fill_limits = c(1, length(block_a_method_order)),
    low_color = "#F7FCF5",
    high_color = "#00441B",
    reverse = TRUE,
    legend_title = "Mean rank"
  )

  top_plot <- plot_metric_heatmap(
    summary_df = summary_df,
    value_col = "top_probability",
    title = "Top-method probability",
    fill_limits = c(0, 1),
    low_color = "#FFF5EB",
    high_color = "#A63603",
    legend_title = "Top-method probability"
  )

  combined_plot <- (power_plot / rank_plot / top_plot) +
    plot_annotation(tag_levels = "A")

  summary_path <- file.path(config$results_dir, "block_a_archetype_summary.csv")
  raw_path <- file.path(config$results_dir, "block_a_archetype_results.rds")
  figure_path_results <- file.path(config$results_dir, "block_a_archetype_recovery.png")

  write.csv(summary_df, summary_path, row.names = FALSE)
  saveRDS(
    list(
      config = config,
      sigma = sigma,
      null_stats = null_stats,
      replicate_results = results_long,
      summary = summary_df
    ),
    raw_path
  )

  ggsave(
    filename = figure_path_results,
    plot = combined_plot,
    width = 16,
    height = 18,
    dpi = 300
  )

  if (isTRUE(config$save_to_ind_figs)) {
    figure_path_ind_figs <- file.path(config$output_dir, "block_a_archetype_recovery.png")
    ggsave(
      filename = figure_path_ind_figs,
      plot = combined_plot,
      width = 16,
      height = 18,
      dpi = 300
    )
  }

  message("Block A outputs written:")
  message("  - ", summary_path)
  message("  - ", raw_path)
  message("  - ", figure_path_results)
  if (isTRUE(config$save_to_ind_figs)) {
    message("  - ", file.path(config$output_dir, "block_a_archetype_recovery.png"))
  }

  invisible(
    list(
      config = config,
      sigma = sigma,
      null_stats = null_stats,
      replicate_results = results_long,
      summary = summary_df,
      plot = combined_plot
    )
  )
}

# ------------------------------------------------------------------------------
# Block B-G wrappers (delegated to simulation_null_diagnostics.R)
# Block mapping:
#   Block B <- legacy Block A (null calibration across LD / pathway sizes)
#   Block C <- legacy Block C (missing-correlation omnibus robustness)
#   Block D <- legacy Block B profile (component breakage stress test)
#   Block E <- legacy Block D + Block E combined (adaptive + leave-one-out)
#   Block F <- legacy Block A power curve (effect-size power)
#   Block G <- legacy Block F (random two-component dropout; optional)
# ------------------------------------------------------------------------------

.legacy_sim_env <- NULL

get_legacy_sim_env <- function(script_path = "simulation_null_diagnostics.R") {
  if (!is.null(.legacy_sim_env) && exists("run_all_simulations", envir = .legacy_sim_env)) {
    return(.legacy_sim_env)
  }

  env <- new.env(parent = globalenv())
  sys.source(script_path, envir = env)
  env$sim_theme <- plot_theme
  .legacy_sim_env <<- env
  env
}

run_legacy_blocks <- function(run_block_b_flag = TRUE,
                              run_block_c_flag = TRUE,
                              run_block_d_flag = TRUE,
                              run_block_e_flag = TRUE,
                              run_block_f_flag = TRUE,
                              run_block_g_flag = TRUE,
                              output_dir = "simulation_results",
                              reduced = FALSE) {
  out <- list()
  if (isTRUE(run_block_b_flag)) out$block_b <- run_block_b(output_dir = output_dir, reduced = reduced)
  if (isTRUE(run_block_c_flag)) out$block_c <- run_block_c(output_dir = output_dir, reduced = reduced)
  if (isTRUE(run_block_d_flag)) out$block_d <- run_block_d(output_dir = output_dir, reduced = reduced)
  if (isTRUE(run_block_e_flag)) out$block_e <- run_block_e(output_dir = output_dir, reduced = reduced)
  if (isTRUE(run_block_f_flag)) out$block_f <- run_block_f(output_dir = output_dir, reduced = reduced)
  if (isTRUE(run_block_g_flag)) out$block_g <- run_block_g(output_dir = output_dir, reduced = reduced)
  invisible(out)
}

legacy_run_params <- function(env, reduced = FALSE) {
  if (isTRUE(reduced)) {
    list(
      n_null = 50L,
      b_sim = 100L,
      pathway_sizes_run = c(20L, 50L),
      pathway_sizes_cde = c(10L, 30L, 50L)
    )
  } else {
    list(
      n_null = env$n_sims_null,
      b_sim = env$b_perm,
      pathway_sizes_run = env$pathway_sizes,
      pathway_sizes_cde = env$pathway_sizes
    )
  }
}

retitle_block_plot <- function(plot_obj, old_block_label, new_block_label) {
  if (inherits(plot_obj, "ggplot") &&
      !is.null(plot_obj$labels$title) &&
      is.character(plot_obj$labels$title)) {
    plot_obj$labels$title <- gsub(old_block_label, new_block_label, plot_obj$labels$title)
  }
  plot_obj
}

strip_geom_text_layers <- function(plot_obj) {
  if (!inherits(plot_obj, "ggplot")) {
    return(plot_obj)
  }
  keep_layer <- vapply(plot_obj$layers, function(layer) {
    !inherits(layer$geom, "GeomText")
  }, logical(1))
  plot_obj$layers <- plot_obj$layers[keep_layer]
  plot_obj
}

strip_plot_headers <- function(plot_obj) {
  if (!inherits(plot_obj, "ggplot")) {
    return(plot_obj)
  }
  plot_obj$labels$title <- NULL
  plot_obj$labels$subtitle <- NULL
  plot_obj$labels$caption <- NULL
  plot_obj
}

make_panel_ab <- function(plot_a, plot_b, title = NULL) {
  out <- plot_a / plot_b +
    patchwork::plot_annotation(tag_levels = "A")
  if (!is.null(title)) {
    out <- out + patchwork::plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(face = "bold", hjust = 0.5))
    )
  }
  out
}

plot_linebar_block_b <- function(results) {
  method_order <- c("ACAT", "FISHER", "TFISHER", "MINP", "STOUFFER", "OMNIBUS")
  df <- results %>%
    dplyr::mutate(
      method_label = factor(toupper(method), levels = method_order),
      cor_structure = factor(cor_structure, levels = c("LD_moderate", "LD_strong", "LD_independent")),
      pathway_size = as.integer(pathway_size)
    )

  scale_fac <- max(df$lambda, na.rm = TRUE) / max(df$type1_05, na.rm = TRUE)

  ggplot(df, aes(x = method_label, group = calibration)) +
    geom_col(aes(y = type1_05, fill = calibration),
             position = position_dodge(width = 0.75), alpha = 0.75, width = 0.65) +
    geom_line(aes(y = lambda / scale_fac, color = calibration),
              position = position_dodge(width = 0.75), linewidth = 0.85) +
    geom_point(aes(y = lambda / scale_fac, color = calibration),
               position = position_dodge(width = 0.75), size = 1.8) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.6) +
    geom_hline(yintercept = 1 / scale_fac, linetype = "dashed", color = "#1F77B4", linewidth = 0.6) +
    facet_grid(pathway_size ~ cor_structure,
               labeller = labeller(pathway_size = function(x) paste0("m=", x))) +
    scale_y_continuous(
      name = "Type I error",
      sec.axis = sec_axis(~ . * scale_fac, name = expression(lambda))
    ) +
    scale_fill_manual(values = c("analytic" = "#E41A1C", "mvn" = "#377EB8")) +
    scale_color_manual(values = c("analytic" = "#B22222", "mvn" = "#1F5AA6")) +
    labs(
      title = "Block B: Null calibration (line+bar)",
      subtitle = "Bars = Type I error, Line = genomic control lambda.\nRed dashed: Type I=0.05; Blue dashed: lambda=1",
      x = "Method", fill = "Calibration", color = "Calibration"
    ) +
    block_a_theme +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))
}

plot_linebar_block_c <- function(results) {
  df <- results %>%
    dplyr::filter(method == "omnibus") %>%
    dplyr::mutate(
      pathway_size = as.integer(pathway_size),
      facet_label = factor(paste0("m=", pathway_size),
                           levels = paste0("m=", sort(unique(pathway_size))))
    )

  scale_fac <- max(df$lambda, na.rm = TRUE) / max(df$type1_05, na.rm = TRUE)

  ggplot(df, aes(x = factor(missing_frac), group = strategy)) +
    geom_col(aes(y = type1_05, fill = strategy),
             position = position_dodge(width = 0.75), alpha = 0.75, width = 0.65) +
    geom_line(aes(y = lambda / scale_fac, color = strategy),
              position = position_dodge(width = 0.75), linewidth = 0.85) +
    geom_point(aes(y = lambda / scale_fac, color = strategy),
               position = position_dodge(width = 0.75), size = 1.8) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.6) +
    geom_hline(yintercept = 1 / scale_fac, linetype = "dashed", color = "#1F77B4", linewidth = 0.6) +
    facet_wrap(~ facet_label, ncol = 3) +
    scale_y_continuous(
      name = "Type I error (omnibus)",
      sec.axis = sec_axis(~ . * scale_fac, name = expression(lambda))
    ) +
    scale_fill_manual(values = c(
      "analytic_fallback" = "#E41A1C",
      "mvn_true_cor" = "#377EB8",
      "mvn_imputed_cor" = "#4DAF4A"
    )) +
    scale_color_manual(values = c(
      "analytic_fallback" = "#B22222",
      "mvn_true_cor" = "#1F5AA6",
      "mvn_imputed_cor" = "#2E8B57"
    )) +
    labs(
      title = "Block C: Missing-correlation robustness (line+bar)",
      subtitle = "Bars = Type I error, Line = genomic control lambda.\nRed dashed: Type I=0.05; Blue dashed: lambda=1",
      x = "Missing correlation fraction", fill = "Strategy", color = "Strategy"
    ) +
    block_a_theme +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
}

prepare_block_d_plot_df <- function(results) {
  results %>%
    dplyr::filter(method %in% c("acat", "fisher", "tfisher", "minp", "stouffer", "omnibus_mvn")) %>%
    dplyr::mutate(
      cor_structure = dplyr::case_when(
        rho == 0 ~ "LD_independent",
        rho == 0.3 ~ "LD_moderate",
        rho == 0.7 ~ "LD_strong",
        TRUE ~ paste0("rho=", rho)
      ),
      cor_structure = factor(cor_structure, levels = c("LD_independent", "LD_moderate", "LD_strong")),
      method_label = factor(dplyr::case_when(
        method == "acat" ~ "ACAT",
        method == "fisher" ~ "Fisher",
        method == "tfisher" ~ "TFisher",
        method == "minp" ~ "minP",
        method == "stouffer" ~ "Stouffer",
        method == "omnibus_mvn" ~ "Omnibus (MVN)",
        TRUE ~ method
      ), levels = c("ACAT", "Fisher", "TFisher", "minP", "Stouffer", "Omnibus (MVN)"))
    )
}

plot_block_d_bar_metric <- function(results, metric = c("lambda", "type1_05")) {
  metric <- match.arg(metric)
  df <- prepare_block_d_plot_df(results) %>%
    dplyr::mutate(
      y = .data[[metric]],
      y_se = dplyr::case_when(
        metric == "lambda" ~ lambda_se,
        metric == "type1_05" ~ type1_05_se,
        TRUE ~ NA_real_
      ),
      ymin = pmax(y - y_se, 0),
      ymax = y + y_se
    )

  dodge <- position_dodge(width = 0.85)
  y_lab <- if (metric == "lambda") expression(lambda) else "Type I error"

  p <- ggplot(df, aes(x = factor(n_broken), y = y, fill = method_label)) +
    geom_col(position = dodge, width = 0.78, alpha = 0.82) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), position = dodge, width = 0.14, linewidth = 0.42) +
    facet_grid(pathway_size ~ cor_structure,
               scales = "free_y",
               labeller = labeller(pathway_size = function(x) paste0("m=", x))) +
    scale_fill_manual(values = c(
      "ACAT" = "#E41A1C", "Fisher" = "#377EB8", "TFisher" = "#4DAF4A",
      "minP" = "#984EA3", "Stouffer" = "#FF7F00", "Omnibus (MVN)" = "#000000"
    )) +
    labs(x = "Number of broken components", y = y_lab, fill = "Method") +
    block_a_theme +
    theme(legend.position = "top")

  if (metric == "lambda") {
    p <- p + geom_hline(yintercept = 1, linetype = "dashed", color = "#1F77B4", linewidth = 0.6)
  } else {
    p <- p + geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.6)
  }

  p
}

plot_linebar_block_d <- function(results) {
  p_lambda <- plot_block_d_bar_metric(results, metric = "lambda")
  p_type1 <- plot_block_d_bar_metric(results, metric = "type1_05")
  make_panel_ab(p_lambda, p_type1, title = "Block D: Component breakage stress test")
}

run_block_b <- function(output_dir = "simulation_results", reduced = FALSE) {
  env <- get_legacy_sim_env()
  prm <- legacy_run_params(env, reduced = reduced)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  results <- env$run_block_a_null(
    n_sims = prm$n_null,
    b = prm$b_sim,
    seed = env$seed_base + 1000L,
    pathway_sizes_arg = prm$pathway_sizes_run
  )

  saveRDS(results, file.path(output_dir, "block_b.rds"))
  plots <- env$plot_block_a_null(results)
  plots$lambda <- retitle_block_plot(plots$lambda, "Block A", "Block B")
  plots$type1 <- retitle_block_plot(plots$type1, "Block A", "Block B")
  # Remove embedded numeric labels for cleaner readability.
  plots$lambda <- strip_geom_text_layers(plots$lambda)
  plots$type1 <- strip_geom_text_layers(plots$type1)
  # Remove panel headers and keep the legend only in the upper panel.
  plots$lambda <- strip_plot_headers(plots$lambda)
  plots$type1 <- strip_plot_headers(plots$type1) + theme(legend.position = "none")
  # Improve readability: allow each facet to use its own y-axis range.
  plots$lambda <- plots$lambda +
    facet_grid(pathway_size ~ cor_structure,
               scales = "free_y",
               labeller = labeller(pathway_size = function(x) paste0("m=", x)))
  plots$type1 <- plots$type1 +
    facet_grid(pathway_size ~ cor_structure,
               scales = "free_y",
               labeller = labeller(pathway_size = function(x) paste0("m=", x)))
  ggsave(file.path(output_dir, "block_b_lambda.png"),
         plots$lambda, width = 14, height = 10, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, "block_b_type1.png"),
         plots$type1, width = 14, height = 10, dpi = 300, bg = "white")

  panel_ab <- make_panel_ab(plots$lambda, plots$type1)
  ggsave(file.path(output_dir, "block_b_panel.png"),
         panel_ab, width = 14, height = 16, dpi = 300, bg = "white")

  # Keep the "linebar" filename for backward compatibility, but save the
  # same A/B panel style used in the main output for readability.
  linebar <- panel_ab
  ggsave(file.path(output_dir, "block_b_linebar.png"),
         linebar, width = 16, height = 12, dpi = 300, bg = "white")

  invisible(results)
}

run_block_c <- function(output_dir = "simulation_results", reduced = FALSE) {
  env <- get_legacy_sim_env()
  prm <- legacy_run_params(env, reduced = reduced)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  results <- env$run_block_c(
    n_sims = prm$n_null,
    b = prm$b_sim,
    seed = env$seed_base + 3000L,
    pathway_sizes_arg = prm$pathway_sizes_cde,
    missing_fracs = c(0, 0.1, 0.3, 0.5, 0.7)
  )

  saveRDS(results, file.path(output_dir, "block_c.rds"))
  plots <- env$plot_block_c(results)
  plots$lambda <- retitle_block_plot(plots$lambda, "Block C", "Block C")
  plots$type1 <- retitle_block_plot(plots$type1, "Block C", "Block C")
  plots$lambda <- strip_plot_headers(plots$lambda)
  plots$type1 <- strip_plot_headers(plots$type1) + theme(legend.position = "none")
  ggsave(file.path(output_dir, "block_c_lambda.png"),
         plots$lambda, width = 10, height = 8, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, "block_c_type1.png"),
         plots$type1, width = 10, height = 8, dpi = 300, bg = "white")

  panel_ab <- make_panel_ab(plots$lambda, plots$type1)
  ggsave(file.path(output_dir, "block_c_panel.png"),
         panel_ab, width = 12, height = 14, dpi = 300, bg = "white")

  linebar <- plot_linebar_block_c(results)
  ggsave(file.path(output_dir, "block_c_linebar.png"),
         linebar, width = 12, height = 10, dpi = 300, bg = "white")

  invisible(results)
}

run_block_d <- function(output_dir = "simulation_results", reduced = FALSE) {
  env <- get_legacy_sim_env()
  prm <- legacy_run_params(env, reduced = reduced)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  results <- env$run_block_b_profile(
    n_sims = prm$n_null,
    b = prm$b_sim,
    seed = env$seed_base + 5200L,
    pathway_sizes_arg = if (isTRUE(reduced)) c(5L, 25L) else c(5L, 25L, 50L),
    rho_values = c(0, 0.3, 0.7),
    n_broken_options = 1:4,
    broken_power = 0.7
  )

  saveRDS(results, file.path(output_dir, "block_d.rds"))
  p_lambda <- strip_plot_headers(plot_block_d_bar_metric(results, metric = "lambda"))
  p_type1 <- strip_plot_headers(plot_block_d_bar_metric(results, metric = "type1_05")) +
    theme(legend.position = "none")
  ggsave(file.path(output_dir, "block_d_lambda.png"),
         p_lambda, width = 16, height = 10, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, "block_d_type1.png"),
         p_type1, width = 16, height = 10, dpi = 300, bg = "white")

  panel_ab <- make_panel_ab(p_lambda, p_type1)
  ggsave(file.path(output_dir, "block_d_panel.png"),
         panel_ab, width = 16, height = 16, dpi = 300, bg = "white")

  # Keep legacy filename but use the same two-panel style.
  linebar <- panel_ab
  ggsave(file.path(output_dir, "block_d_linebar.png"),
         linebar, width = 16, height = 12, dpi = 300, bg = "white")

  invisible(results)
}

plot_block_e_linebar_panel <- function(adaptive_results, leave1out_results) {
  adaptive_df <- adaptive_results %>%
    dplyr::filter(method %in% c("omnibus_analytical", "omnibus_mvn",
                                "omnibus_adaptive", "omnibus_adaptive_mvn")) %>%
    dplyr::mutate(
      source = "Adaptive omnibus",
      method_label = dplyr::case_when(
        method == "omnibus_analytical" ~ "Analytic",
        method == "omnibus_mvn" ~ "MVN",
        method == "omnibus_adaptive" ~ "Adaptive+Analytic",
        method == "omnibus_adaptive_mvn" ~ "Adaptive+MVN",
        TRUE ~ method
      )
    ) %>%
    dplyr::group_by(source, method_label) %>%
    dplyr::summarise(
      lambda = mean(lambda, na.rm = TRUE),
      type1_05 = mean(type1_05, na.rm = TRUE),
      .groups = "drop"
    )

  leave1out_df <- leave1out_results %>%
    dplyr::filter(method %in% c("omnibus_minus_acat", "omnibus_minus_fisher",
                                "omnibus_minus_tfisher", "omnibus_minus_minp",
                                "omnibus_minus_stouffer", "omnibus_all")) %>%
    dplyr::mutate(
      source = "Leave-one-out omnibus",
      method_label = dplyr::case_when(
        method == "omnibus_minus_acat" ~ "-ACAT",
        method == "omnibus_minus_fisher" ~ "-Fisher",
        method == "omnibus_minus_tfisher" ~ "-TFisher",
        method == "omnibus_minus_minp" ~ "-minP",
        method == "omnibus_minus_stouffer" ~ "-Stouffer",
        method == "omnibus_all" ~ "All",
        TRUE ~ method
      )
    ) %>%
    dplyr::group_by(source, method_label) %>%
    dplyr::summarise(
      lambda = mean(lambda, na.rm = TRUE),
      type1_05 = mean(type1_05, na.rm = TRUE),
      .groups = "drop"
    )

  panel_df <- dplyr::bind_rows(adaptive_df, leave1out_df)
  panel_df <- panel_df %>%
    dplyr::mutate(
      source = factor(source, levels = c("Adaptive omnibus", "Leave-one-out omnibus"))
    )

  max_type1 <- max(panel_df$type1_05, na.rm = TRUE)
  max_lambda <- max(panel_df$lambda, na.rm = TRUE)
  lambda_scale <- ifelse(max_type1 > 0, max_lambda / max_type1, 1)
  panel_df$lambda_scaled <- panel_df$lambda / lambda_scale

  ggplot(panel_df, aes(x = method_label)) +
    geom_col(aes(y = type1_05), fill = "#377EB8", alpha = 0.75, width = 0.7) +
    geom_line(aes(y = lambda_scaled, group = 1), color = "#E41A1C", linewidth = 0.9) +
    geom_point(aes(y = lambda_scaled), color = "#E41A1C", size = 2) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.6) +
    geom_hline(yintercept = 1 / lambda_scale, linetype = "dashed", color = "#1F77B4", linewidth = 0.6) +
    facet_wrap(~ source, scales = "free_x", ncol = 1) +
    scale_y_continuous(
      name = "Type I error (alpha=0.05)",
      sec.axis = sec_axis(~ . * lambda_scale, name = expression(lambda))
    ) +
    labs(
      title = "Block E: Adaptive and Leave-one-out Omnibus (line+bar)",
      subtitle = "Bars: Type I error; Line: genomic control lambda.\nRed dashed: Type I=0.05; Blue dashed: lambda=1",
      x = "Omnibus variant"
    ) +
    block_a_theme +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      panel.grid.major.x = element_blank()
    )
}

plot_block_e_ab_panel <- function(adaptive_results, leave1out_results) {
  adaptive_df <- adaptive_results %>%
    dplyr::mutate(
      source = "Adaptive",
      method_label = dplyr::case_when(
        method == "omnibus_analytical" ~ "Analytic",
        method == "omnibus_mvn" ~ "MVN",
        method == "omnibus_adaptive" ~ "Adaptive+Analytic",
        method == "omnibus_adaptive_mvn" ~ "Adaptive+MVN",
        TRUE ~ method
      )
    ) %>%
    dplyr::group_by(source, method_label) %>%
    dplyr::summarise(
      n_obs = dplyr::n(),
      lambda = mean(lambda, na.rm = TRUE),
      lambda_se = stats::sd(lambda, na.rm = TRUE) / sqrt(n_obs),
      type1_05 = mean(type1_05, na.rm = TRUE),
      type1_05_se = stats::sd(type1_05, na.rm = TRUE) / sqrt(n_obs),
      .groups = "drop"
    ) %>%
    dplyr::select(-n_obs)

  leave1_df <- leave1out_results %>%
    dplyr::mutate(
      source = "Leave-one-out",
      method_label = dplyr::case_when(
        method == "omnibus_minus_acat" ~ "-ACAT",
        method == "omnibus_minus_fisher" ~ "-Fisher",
        method == "omnibus_minus_tfisher" ~ "-TFisher",
        method == "omnibus_minus_minp" ~ "-minP",
        method == "omnibus_minus_stouffer" ~ "-Stouffer",
        method == "omnibus_all" ~ "All",
        TRUE ~ method
      )
    ) %>%
    dplyr::group_by(source, method_label) %>%
    dplyr::summarise(
      n_obs = dplyr::n(),
      lambda = mean(lambda, na.rm = TRUE),
      lambda_se = stats::sd(lambda, na.rm = TRUE) / sqrt(n_obs),
      type1_05 = mean(type1_05, na.rm = TRUE),
      type1_05_se = stats::sd(type1_05, na.rm = TRUE) / sqrt(n_obs),
      .groups = "drop"
    ) %>%
    dplyr::select(-n_obs)

  df <- dplyr::bind_rows(adaptive_df, leave1_df) %>%
    dplyr::mutate(
      source = factor(source, levels = c("Adaptive", "Leave-one-out")),
      lambda_se = dplyr::coalesce(lambda_se, 0),
      type1_05_se = dplyr::coalesce(type1_05_se, 0)
    )

  p_lambda <- ggplot(df, aes(x = method_label, y = lambda, fill = source)) +
    geom_col(alpha = 0.82, width = 0.72) +
    geom_errorbar(aes(ymin = pmax(lambda - lambda_se, 0), ymax = lambda + lambda_se),
                  width = 0.14, linewidth = 0.7, color = "black") +
    facet_wrap(~ source, scales = "free_x", ncol = 1) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "#1F77B4", linewidth = 0.7) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4),
                       expand = expansion(mult = c(0, 0.08))) +
    labs(x = "Variant", y = expression(lambda), fill = "Family") +
    block_a_theme +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      axis.title = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(angle = 35, hjust = 1, size = 13),
      axis.text.y = element_text(size = 13),
      strip.text = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 14)
    )

  p_type1 <- ggplot(df, aes(x = method_label, y = type1_05, fill = source)) +
    geom_col(alpha = 0.82, width = 0.72) +
    geom_errorbar(aes(ymin = pmax(type1_05 - type1_05_se, 0), ymax = type1_05 + type1_05_se),
                  width = 0.14, linewidth = 0.7, color = "black") +
    facet_wrap(~ source, scales = "free_x", ncol = 1) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.6) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4),
                       expand = expansion(mult = c(0, 0.08))) +
    labs(x = "Variant", y = "Type I error", fill = "Family") +
    block_a_theme +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      axis.title = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(angle = 35, hjust = 1, size = 13),
      axis.text.y = element_text(size = 13),
      strip.text = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 14),
      legend.position = "none"
    )

  p_lambda <- strip_plot_headers(p_lambda)
  p_type1 <- strip_plot_headers(p_type1)
  make_panel_ab(p_lambda, p_type1)
}

run_block_e <- function(output_dir = "simulation_results", reduced = FALSE) {
  env <- get_legacy_sim_env()
  prm <- legacy_run_params(env, reduced = reduced)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  results_adaptive <- env$run_block_d(
    n_train = if (isTRUE(reduced)) 50L else 200L,
    n_test = if (isTRUE(reduced)) 100L else 300L,
    b = prm$b_sim,
    seed = env$seed_base + 4000L,
    pathway_sizes_arg = prm$pathway_sizes_cde,
    rho_values = c(0, 0.3, 0.7)
  )

  results_leave1out <- env$run_block_e(
    n_sims = prm$n_null,
    b = prm$b_sim,
    seed = env$seed_base + 5000L,
    pathway_sizes_arg = if (isTRUE(reduced)) c(20L, 50L) else prm$pathway_sizes_run
  )

  combined_results <- list(
    adaptive = results_adaptive,
    leave1out = results_leave1out
  )

  saveRDS(combined_results, file.path(output_dir, "block_e.rds"))
  saveRDS(results_adaptive, file.path(output_dir, "block_e_adaptive.rds"))
  saveRDS(results_leave1out, file.path(output_dir, "block_e_leave1out.rds"))

  # Use the same two-panel style as other calibration blocks.
  panel_plot <- plot_block_e_ab_panel(results_adaptive, results_leave1out)
  ggsave(file.path(output_dir, "block_e_panel.png"),
         panel_plot, width = 14, height = 12, dpi = 300, bg = "white")

  panel_ab <- panel_plot
  ggsave(file.path(output_dir, "block_e_panel_ab.png"),
         panel_ab, width = 14, height = 14, dpi = 300, bg = "white")

  ggsave(file.path(output_dir, "block_e_linebar.png"),
         panel_plot, width = 14, height = 12, dpi = 300, bg = "white")

  invisible(combined_results)
}

run_block_f <- function(output_dir = "simulation_results", reduced = FALSE) {
  env <- get_legacy_sim_env()
  prm <- legacy_run_params(env, reduced = reduced)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  results <- env$run_block_a(
    n_sims = if (isTRUE(reduced)) 100L else 200L,
    b = prm$b_sim,
    seed = env$seed_base + 1100L,
    pathway_sizes_arg = if (isTRUE(reduced)) c(20L, 50L) else prm$pathway_sizes_run
  )

  saveRDS(results, file.path(output_dir, "block_f_power.rds"))
  power_plot <- env$plot_block_a_power(results)
  power_plot <- retitle_block_plot(power_plot, "Block A", "Block F")
  ggsave(file.path(output_dir, "block_f_power.png"),
         power_plot, width = 16, height = 12, dpi = 300, bg = "white")

  invisible(results)
}

run_block_g <- function(output_dir = "simulation_results", reduced = FALSE) {
  env <- get_legacy_sim_env()
  prm <- legacy_run_params(env, reduced = reduced)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  results_full <- env$run_block_a_null(
    n_sims = prm$n_null,
    b = prm$b_sim,
    seed = env$seed_base + 6000L,
    pathway_sizes_arg = if (isTRUE(reduced)) c(20L, 50L) else prm$pathway_sizes_run
  )

  results <- results_full %>%
    dplyr::filter(method %in% c("acat", "fisher", "tfisher", "minp", "stouffer"))

  saveRDS(results, file.path(output_dir, "block_g.rds"))
  saveRDS(results_full, file.path(output_dir, "block_g_raw.rds"))
  plots <- env$plot_block_a_null(results)
  plots$lambda <- retitle_block_plot(plots$lambda, "Block A", "Block G")
  plots$type1 <- retitle_block_plot(plots$type1, "Block A", "Block G")
  plots$lambda <- strip_plot_headers(strip_geom_text_layers(plots$lambda)) +
    facet_grid(pathway_size ~ cor_structure,
               scales = "free_y",
               labeller = labeller(pathway_size = function(x) paste0("m=", x))) +
    theme(legend.position = "top") +
    guides(fill = guide_legend(title = NULL))
  plots$type1 <- strip_plot_headers(strip_geom_text_layers(plots$type1)) +
    facet_grid(pathway_size ~ cor_structure,
               scales = "free_y",
               labeller = labeller(pathway_size = function(x) paste0("m=", x))) +
    theme(legend.position = "none")
  ggsave(file.path(output_dir, "block_g_lambda.png"),
         plots$lambda, width = 14, height = 10, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, "block_g_type1.png"),
         plots$type1, width = 14, height = 10, dpi = 300, bg = "white")

  panel_ab <- make_panel_ab(plots$lambda, plots$type1)
  ggsave(file.path(output_dir, "block_g_panel.png"),
         panel_ab, width = 14, height = 16, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, "block_g_linebar.png"),
         panel_ab, width = 14, height = 12, dpi = 300, bg = "white")

  invisible(results)
}

run_block_h <- function(config = BLOCK_H_CONFIG, reduced = FALSE) {
  cfg <- config
  if (isTRUE(reduced)) {
    cfg$m_values <- intersect(cfg$m_values, c(20L, 50L))
    cfg$effect_scales <- c(0.6, 1.0, 1.4, 1.8)
    cfg$n_reps <- 120L
    cfg$n_null_cal <- 800L
  }

  dir.create(cfg$results_dir, recursive = TRUE, showWarnings = FALSE)

  method_colors <- c(
    "ACAT" = "#E41A1C",
    "Fisher" = "#377EB8",
    "TFisher" = "#4DAF4A",
    "minP" = "#984EA3",
    "Stouffer" = "#FF7F00",
    "Omnibus" = "#000000"
  )

  set.seed(cfg$seed)
  all_results <- list()
  iter <- 1L
  total <- length(cfg$m_values) * length(block_a_archetypes) * length(cfg$effect_scales)

  for (m in cfg$m_values) {
    sigma <- ensure_pd(make_cor_exchangeable(m = m, rho = cfg$rho))
    message(sprintf("Block H: precomputing null stats for m=%d", m))
    null_stats <- precompute_null_distribution(
      m = m,
      sigma = sigma,
      n_null = cfg$n_null_cal,
      seed = cfg$seed + m,
      tau_grid = cfg$tau_grid
    )

    for (arch in block_a_archetypes) {
      for (eff in cfg$effect_scales) {
        message(sprintf(
          "Block H: %d/%d | m=%d | %s | effect_scale=%.2f",
          iter, total, m, arch$code, eff
        ))
        res <- bind_rows(lapply(seq_len(cfg$n_reps), function(rep_idx) {
          simulate_archetype_once_scaled(
            archetype = arch,
            sigma = sigma,
            null_stats = null_stats,
            effect_scale = eff,
            tau_grid = cfg$tau_grid
          ) %>%
            mutate(replicate = rep_idx, pathway_size = m)
        }))
        all_results[[length(all_results) + 1L]] <- res
        iter <- iter + 1L
      }
    }
  }

  results_long <- bind_rows(all_results) %>%
    mutate(
      pathway_size = factor(pathway_size, levels = sort(unique(pathway_size))),
      archetype = factor(archetype, levels = block_a_archetype_order),
      method = factor(method, levels = block_a_method_order)
    )

  summary_df <- results_long %>%
    group_by(pathway_size, archetype, effect_scale, method) %>%
    summarize(
      power = mean(p_value < cfg$alpha),
      mean_rank = mean(rank),
      top_probability = mean(is_top),
      .groups = "drop"
    )

  regret_df <- summary_df %>%
    group_by(pathway_size, archetype, effect_scale) %>%
    mutate(
      best_power = max(power, na.rm = TRUE),
      power_regret = best_power - power
    ) %>%
    ungroup()

  omnibus_regret <- regret_df %>%
    filter(method == "Omnibus") %>%
    group_by(archetype) %>%
    summarize(
      mean_regret = mean(power_regret, na.rm = TRUE),
      median_regret = median(power_regret, na.rm = TRUE),
      max_regret = max(power_regret, na.rm = TRUE),
      within_0.05 = mean(power_regret <= 0.05, na.rm = TRUE),
      within_0.10 = mean(power_regret <= 0.10, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(archetype)

  p_power <- ggplot(
    summary_df,
    aes(x = effect_scale, y = power, color = method, linetype = method)
  ) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.6) +
    facet_grid(pathway_size ~ archetype) +
    scale_color_manual(values = method_colors) +
    scale_linetype_manual(values = c(
      "ACAT" = "solid",
      "Fisher" = "dashed",
      "TFisher" = "dotdash",
      "minP" = "longdash",
      "Stouffer" = "twodash",
      "Omnibus" = "solid"
    )) +
    labs(
      x = "Effect-size scale multiplier",
      y = "Power (alpha = 0.05)",
      color = "Method",
      linetype = "Method"
    ) +
    block_a_theme +
    theme(
      legend.position = "top",
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      strip.text = element_text(size = 13, face = "bold")
    )

  p_regret <- ggplot(
    regret_df %>% filter(method == "Omnibus"),
    aes(x = effect_scale, y = power_regret, group = pathway_size, color = pathway_size)
  ) +
    geom_hline(yintercept = 0.05, color = "grey50", linetype = "dashed", linewidth = 0.4) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.8) +
    facet_wrap(~ archetype, ncol = 3) +
    scale_color_manual(values = c("5" = "#1b9e77", "20" = "#7570b3", "50" = "#d95f02")) +
    labs(
      x = "Effect-size scale multiplier",
      y = "Omnibus power regret",
      color = "Pathway size (m)"
    ) +
    block_a_theme +
    theme(
      legend.position = "top",
      axis.text.x = element_text(size = 13),
      axis.text.y = element_text(size = 13),
      strip.text = element_text(size = 13, face = "bold")
    )

  out_summary <- file.path(cfg$results_dir, "block_h_archetype_power_summary.csv")
  out_regret <- file.path(cfg$results_dir, "block_h_omnibus_regret_by_archetype.csv")
  out_rds <- file.path(cfg$results_dir, "block_h_archetype_power.rds")
  out_plot_power <- file.path(cfg$results_dir, "block_h_power.png")
  out_plot_regret <- file.path(cfg$results_dir, "block_h_omnibus_regret.png")

  write.csv(summary_df, out_summary, row.names = FALSE)
  write.csv(omnibus_regret, out_regret, row.names = FALSE)
  saveRDS(
    list(config = cfg, results_long = results_long, summary = summary_df, regret = regret_df,
         omnibus_regret = omnibus_regret),
    out_rds
  )

  ggsave(out_plot_power, p_power, width = 20, height = 12, dpi = 300, bg = "white")
  ggsave(out_plot_regret, p_regret, width = 16, height = 9, dpi = 300, bg = "white")

  message("Block H outputs written:")
  message("  - ", out_summary)
  message("  - ", out_regret)
  message("  - ", out_rds)
  message("  - ", out_plot_power)
  message("  - ", out_plot_regret)

  invisible(list(
    config = cfg,
    summary = summary_df,
    regret = regret_df,
    omnibus_regret = omnibus_regret,
    power_plot = p_power,
    regret_plot = p_regret
  ))
}

run_all_figs <- function(run_block_a_flag = TRUE,
                         run_block_b_flag = FALSE,
                         run_block_c_flag = FALSE,
                         run_block_d_flag = FALSE,
                         run_block_e_flag = FALSE,
                         run_block_f_flag = FALSE,
                         run_block_g_flag = FALSE,
                         run_block_h_flag = FALSE,
                         run_block_i_flag = FALSE,
                         reduced = FALSE,
                         output_dir = "simulation_results") {
  out <- list()
  if (isTRUE(run_block_a_flag)) {
    out$block_a <- run_block_a()
  }
  if (any(c(run_block_b_flag, run_block_c_flag, run_block_d_flag, run_block_e_flag, run_block_f_flag, run_block_g_flag))) {
    out$legacy <- run_legacy_blocks(
      run_block_b_flag = run_block_b_flag,
      run_block_c_flag = run_block_c_flag,
      run_block_d_flag = run_block_d_flag,
      run_block_e_flag = run_block_e_flag,
      run_block_f_flag = run_block_f_flag,
      run_block_g_flag = run_block_g_flag,
      output_dir = output_dir,
      reduced = reduced
    )
  }
  if (isTRUE(run_block_h_flag)) {
    out$block_h <- run_block_h(reduced = reduced)
  }
  if (isTRUE(run_block_i_flag)) {
    out$block_i <- run_block_i(reduced = reduced)
  }
  invisible(out)
}

# ==============================================================================
# Block I: Computational Benchmark (Time and Memory vs B and Pathway Size)
# ==============================================================================

BLOCK_I_CONFIG <- list(
  seed = 20260318L,
  pathway_sizes = c(10L, 20L, 50L),
  n_pathways_values = c(1L, 10L, 100L),
  n_pathways_exact_max = 10L,
  B_values = c(50000L, 100000L, 250000L, 500000L, 1000000L),
  thread_values = c(1L, 4L),
  rho = 0.20,
  tau_grid = c(0.10, 0.05, 0.01),
  min_p = 1e-15,
  stouffer_alternative = "greater",
  n_reps = 1L,
  results_dir = "simulation_results",
  output_dir = file.path("Figures", "Fig_Benchmark")
)

# Lightweight memory probe for block I benchmarks.
.block_i_mem_used_mb <- local({
  has_lobstr <- requireNamespace("lobstr", quietly = TRUE)
  function() {
    if (has_lobstr) {
      as.numeric(lobstr::mem_used()) / (1024^2)
    } else {
      # Fallback: approximate from GC summary.
      sum(gc()[, 2])
    }
  }
})

#' Run a single benchmark iteration
#' @param m Pathway size (number of genes)
#' @param B Number of permutations
#' @param sigma Correlation matrix
#' @param n_pathways Number of pathways processed in one run
#' @param threads Number of threads for inner null calculations
#' @param tau_grid Adaptive TFisher tau grid
#' @param min_p Lower p clipping bound for numerical stability
#' @param stouffer_alternative Stouffer alternative ("greater", "less", "two.sided")
#' @param seed Random seed
#' @return List with time and memory metrics
run_benchmark_single <- function(m, B, sigma, n_pathways = 1L, threads = 1L,
                                 tau_grid = c(0.10, 0.05, 0.01), min_p = 1e-15,
                                 stouffer_alternative = "greater", seed = 123L) {
  set.seed(seed)
  n_pathways <- suppressWarnings(as.integer(n_pathways))
  if (!is.finite(n_pathways) || is.na(n_pathways) || n_pathways < 1L) n_pathways <- 1L
  threads <- suppressWarnings(as.integer(threads))
  if (!is.finite(threads) || is.na(threads) || threads < 1L) threads <- 1L
  if (.Platform$OS.type == "windows") threads <- 1L
  tau_grid <- sort(unique(as.numeric(tau_grid)))
  tau_grid <- tau_grid[is.finite(tau_grid) & !is.na(tau_grid) & tau_grid > 0 & tau_grid < 1]
  if (!length(tau_grid)) tau_grid <- c(0.10, 0.05, 0.01)
  stouffer_alternative <- match.arg(stouffer_alternative, c("greater", "less", "two.sided"))

  # Track memory baseline (used for approximate peak delta over main steps).
  gc(reset = TRUE)
  mem_base <- .block_i_mem_used_mb()
  mem_trace <- c(mem_base)

  # Time the computation
  t_start <- proc.time()[["elapsed"]]

  # Process n_pathways sequentially, matching CATFISH's per-pathway MVN workload.
  checksum <- 0
  for (p_idx in seq_len(n_pathways)) {
    Z_null <- MASS::mvrnorm(n = B, mu = rep(0, m), Sigma = sigma)
    mem_trace <- c(mem_trace, .block_i_mem_used_mb())

    P_null <- pnorm(Z_null, lower.tail = FALSE)
    P_null <- pmax(pmin(P_null, 1 - min_p), min_p)
    mem_trace <- c(mem_trace, .block_i_mem_used_mb())

    # Full five-component workload per null draw:
    # ACAT, Fisher, adaptive TFisher, minP, Stouffer, then omnibus across components.
    acat_null <- {
      tan_vals <- tan((0.5 - P_null) * pi)
      tstat <- rowMeans(tan_vals)
      p <- 0.5 - atan(tstat) / pi
      pmax(pmin(p, 1 - min_p), min_p)
    }

    fisher_null <- {
      stat <- -2 * rowSums(log(P_null))
      p <- stats::pchisq(stat, df = 2 * ncol(P_null), lower.tail = FALSE)
      pmax(pmin(p, 1 - min_p), min_p)
    }

    tfisher_row <- function(b) {
      pb <- P_null[b, ]
      pb <- pmax(pmin(pb, 1 - min_p), min_p)
      if (length(pb) < 2L) return(NA_real_)

      # Exact adaptive soft-TFisher when available; fallback approximation otherwise.
      if (use_tfisher_pkg) {
        best <- Inf
        for (tau in tau_grid) {
          st <- TFisher::stat.soft(p = pb, tau1 = tau)
          pv <- 1 - as.numeric(TFisher::p.soft(q = st, n = length(pb), tau1 = tau, M = NULL))
          if (is.finite(pv) && pv < best) best <- pv
        }
        if (is.finite(best)) {
          return(pmax(pmin(best, 1 - min_p), min_p))
        }
      }

      tfisher_ps <- vapply(tau_grid, function(tau) {
        keep <- pb[pb < tau]
        if (!length(keep)) return(1)
        stat <- -2 * sum(log(keep))
        stats::pchisq(stat, df = 2 * length(keep), lower.tail = FALSE)
      }, numeric(1))
      p_tf <- min(tfisher_ps) * length(tau_grid)
      pmax(pmin(p_tf, 1 - min_p), min_p)
    }

    if (threads > 1L && B >= 1000L && .Platform$OS.type != "windows") {
      tfisher_null <- unlist(parallel::mclapply(seq_len(B), tfisher_row, mc.cores = threads))
    } else {
      tfisher_null <- vapply(seq_len(B), tfisher_row, numeric(1))
    }

    minp_null <- {
      pmin_gene <- apply(P_null, 1, min)
      p <- 1 - (1 - pmin_gene)^ncol(P_null)
      pmax(pmin(p, 1 - min_p), min_p)
    }

    stouffer_null <- {
      z_st <- rowSums(Z_null) / sqrt(ncol(Z_null))
      if (stouffer_alternative == "greater") {
        p <- stats::pnorm(z_st, lower.tail = FALSE)
      } else if (stouffer_alternative == "less") {
        p <- stats::pnorm(z_st, lower.tail = TRUE)
      } else {
        p <- 2 * stats::pnorm(-abs(z_st))
      }
      pmax(pmin(p, 1 - min_p), min_p)
    }

    method_mat <- cbind(acat_null, fisher_null, tfisher_null, minp_null, stouffer_null)
    method_mat <- pmax(pmin(method_mat, 1 - min_p), min_p)
    tan_method <- tan((0.5 - method_mat) * pi)
    valid_counts <- rowSums(is.finite(tan_method) & !is.na(tan_method))
    tan_sum <- rowSums(ifelse(is.finite(tan_method) & !is.na(tan_method), tan_method, 0))
    omni_null <- rep(NA_real_, B)
    ok <- valid_counts > 0
    omni_null[ok] <- 0.5 - atan(tan_sum[ok] / valid_counts[ok]) / pi
    omni_null <- pmax(pmin(omni_null, 1 - min_p), min_p)

    mem_trace <- c(mem_trace, .block_i_mem_used_mb())

    checksum <- checksum +
      sum(acat_null, na.rm = TRUE) +
      sum(fisher_null, na.rm = TRUE) +
      sum(tfisher_null, na.rm = TRUE) +
      sum(minp_null, na.rm = TRUE) +
      sum(stouffer_null, na.rm = TRUE) +
      sum(omni_null, na.rm = TRUE)
  }

  elapsed_sec <- as.numeric(proc.time()[["elapsed"]] - t_start)

  # Approximate peak memory consumed during benchmarked steps.
  mem_peak <- max(mem_trace, na.rm = TRUE) - mem_base
  mem_peak <- max(as.numeric(mem_peak), 0)

  # Keep outputs alive through timing scope and prevent accidental elision.
  if (!is.finite(checksum) ||
      length(acat_null) != B ||
      length(fisher_null) != B ||
      length(tfisher_null) != B ||
      length(minp_null) != B ||
      length(stouffer_null) != B ||
      length(omni_null) != B) {
    stop("Block I benchmark internal size mismatch.", call. = FALSE)
  }

  list(
    m = m,
    n_pathways = n_pathways,
    B = B,
    threads = threads,
    time_sec = elapsed_sec,
    mem_mb = mem_peak
  )
}

#' Run Block I: Computational Benchmark
#' @param config Configuration list
#' @param reduced If TRUE, use reduced parameter set for quick testing
#' @return Invisible list with results and plots
run_block_i <- function(config = BLOCK_I_CONFIG, reduced = FALSE) {
  cfg <- config
  if (isTRUE(reduced)) {
    cfg$pathway_sizes <- c(10L, 50L)
    cfg$n_pathways_values <- c(1L, 10L, 100L)
    cfg$n_pathways_exact_max <- 10L
    cfg$B_values <- c(50000L, 100000L)
    cfg$thread_values <- c(1L, 4L)
    cfg$n_reps <- 1L
  }

  dir.create(cfg$results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

  n_cores <- suppressWarnings(as.integer(parallel::detectCores(logical = FALSE)))
  if (!is.finite(n_cores) || is.na(n_cores) || n_cores < 1L) n_cores <- 1L

  cfg$thread_values <- suppressWarnings(as.integer(cfg$thread_values))
  cfg$thread_values <- sort(unique(cfg$thread_values[is.finite(cfg$thread_values) & !is.na(cfg$thread_values) & cfg$thread_values >= 1L]))
  if (!length(cfg$thread_values)) cfg$thread_values <- 1L
  cfg$thread_values <- cfg$thread_values[cfg$thread_values <= n_cores]
  if (!length(cfg$thread_values)) cfg$thread_values <- 1L
  if (.Platform$OS.type == "windows") cfg$thread_values <- 1L
  cfg$n_pathways_values <- suppressWarnings(as.integer(cfg$n_pathways_values))
  cfg$n_pathways_values <- sort(unique(cfg$n_pathways_values[is.finite(cfg$n_pathways_values) & !is.na(cfg$n_pathways_values) & cfg$n_pathways_values >= 1L]))
  if (!length(cfg$n_pathways_values)) cfg$n_pathways_values <- 1L
  cfg$n_pathways_exact_max <- suppressWarnings(as.integer(cfg$n_pathways_exact_max))
  if (!is.finite(cfg$n_pathways_exact_max) || is.na(cfg$n_pathways_exact_max) || cfg$n_pathways_exact_max < 1L) {
    cfg$n_pathways_exact_max <- max(cfg$n_pathways_values)
  }
  n_pathways_bench <- cfg$n_pathways_values[cfg$n_pathways_values <= cfg$n_pathways_exact_max]
  if (!length(n_pathways_bench)) n_pathways_bench <- min(cfg$n_pathways_values)

  set.seed(cfg$seed)

  # Run benchmarks
  results_list <- list()
  total_iters <- length(cfg$pathway_sizes) * length(n_pathways_bench) * length(cfg$B_values) * length(cfg$thread_values) * cfg$n_reps
  iter <- 1L

  for (m in cfg$pathway_sizes) {
    # Create correlation matrix for this pathway size
    sigma <- ensure_pd(make_cor_exchangeable(m = m, rho = cfg$rho))

    for (n_pathways in n_pathways_bench) {
      for (B in cfg$B_values) {
        for (threads in cfg$thread_values) {
          for (rep_idx in seq_len(cfg$n_reps)) {
            message(sprintf(
              "Block I: %d/%d | m=%d | pathways=%d | B=%s | threads=%d | rep=%d",
              iter, total_iters, m, n_pathways, format(B, big.mark = ","), threads, rep_idx
            ))

            res <- run_benchmark_single(
              m = m,
              n_pathways = n_pathways,
              B = B,
              sigma = sigma,
              threads = threads,
              tau_grid = cfg$tau_grid,
              min_p = cfg$min_p,
              stouffer_alternative = cfg$stouffer_alternative,
              seed = cfg$seed + iter
            )
            res$replicate <- rep_idx
            results_list[[length(results_list) + 1L]] <- res
            iter <- iter + 1L
          }
        }
      }
    }
  }

  # Combine results
  results_df <- do.call(rbind, lapply(results_list, as.data.frame))

  # Summarize across replicates
  summary_df <- results_df %>%
    group_by(m, n_pathways, B, threads) %>%
    summarize(
      time_sec_mean = mean(time_sec),
      time_sec_median = stats::median(time_sec),
      time_sec_sd = sd(time_sec, na.rm = TRUE),
      mem_mb_mean = mean(mem_mb),
      mem_mb_median = stats::median(mem_mb),
      mem_mb_sd = sd(mem_mb, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      estimate_type = "observed",
      time_sec_sd = ifelse(is.na(time_sec_sd), 0, time_sec_sd),
      mem_mb_sd = ifelse(is.na(mem_mb_sd), 0, mem_mb_sd)
    )

  # Extrapolate larger pathway-count settings (e.g., 100) from observed settings (e.g., 1 and 10).
  missing_np <- setdiff(cfg$n_pathways_values, sort(unique(summary_df$n_pathways)))
  if (length(missing_np)) {
    pred_rows <- list()
    idx <- 1L
    key_df <- unique(summary_df[, c("m", "B", "threads"), drop = FALSE])
    for (k in seq_len(nrow(key_df))) {
      sub <- summary_df[
        summary_df$m == key_df$m[k] &
          summary_df$B == key_df$B[k] &
          summary_df$threads == key_df$threads[k],
        , drop = FALSE
      ]
      if (nrow(sub) < 2L || length(unique(sub$n_pathways)) < 2L) next

      fit_time_mean <- stats::lm(time_sec_mean ~ n_pathways, data = sub)
      fit_time_med  <- stats::lm(time_sec_median ~ n_pathways, data = sub)
      fit_mem_mean  <- stats::lm(mem_mb_mean ~ n_pathways, data = sub)
      fit_mem_med   <- stats::lm(mem_mb_median ~ n_pathways, data = sub)

      for (np in missing_np) {
        pred_rows[[idx]] <- data.frame(
          m = key_df$m[k],
          n_pathways = as.integer(np),
          B = key_df$B[k],
          threads = key_df$threads[k],
          time_sec_mean = pmax(as.numeric(stats::predict(fit_time_mean, newdata = data.frame(n_pathways = np))), 0),
          time_sec_median = pmax(as.numeric(stats::predict(fit_time_med, newdata = data.frame(n_pathways = np))), 0),
          time_sec_sd = NA_real_,
          mem_mb_mean = pmax(as.numeric(stats::predict(fit_mem_mean, newdata = data.frame(n_pathways = np))), 0),
          mem_mb_median = pmax(as.numeric(stats::predict(fit_mem_med, newdata = data.frame(n_pathways = np))), 0),
          mem_mb_sd = NA_real_,
          estimate_type = "projected",
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
    if (length(pred_rows)) {
      summary_df <- dplyr::bind_rows(summary_df, do.call(rbind, pred_rows))
    }
  }

  summary_df <- summary_df %>%
    group_by(m, n_pathways, B) %>%
    mutate(
      time_baseline_1thread = ifelse(any(threads == 1L), time_sec_median[threads == 1L][1], min(time_sec_median, na.rm = TRUE)),
      speedup_vs_1thread = time_baseline_1thread / time_sec_median
    ) %>%
    ungroup() %>%
    mutate(
      pathway_size = factor(paste0("m = ", m), levels = paste0("m = ", sort(unique(m)))),
      n_pathways_label = factor(
        paste0("P = ", n_pathways),
        levels = paste0("P = ", sort(unique(cfg$n_pathways_values)))
      ),
      B_label = format(B, big.mark = ",", scientific = FALSE),
      total_null_draws = B * n_pathways,
      thread_label = factor(
        ifelse(threads == 1L, "1 thread", paste0(threads, " threads")),
        levels = unique(ifelse(sort(unique(threads)) == 1L, "1 thread", paste0(sort(unique(threads)), " threads")))
      ),
      time_plot = pmax(time_sec_median, 1e-8)
    )

  # Paper-ready benchmark table.
  table_df <- summary_df %>%
    transmute(
      pathway_size = m,
      n_pathways = n_pathways,
      B = B,
      total_null_draws = B * n_pathways,
      threads = threads,
      estimate_type = estimate_type,
      time_sec_mean = round(time_sec_mean, 4),
      time_sec_median = round(time_sec_median, 4),
      time_sec_sd = round(time_sec_sd, 4),
      speedup_vs_1thread = round(speedup_vs_1thread, 3),
      mem_mb_mean = round(mem_mb_mean, 2),
      mem_mb_median = round(mem_mb_median, 2),
      mem_mb_sd = round(mem_mb_sd, 2)
    ) %>%
    arrange(pathway_size, n_pathways, B, threads)

  # Linear model for runtime prediction: time_sec_mean = intercept + slope * B
  time_model_df <- summary_df %>%
    group_by(m, n_pathways, threads, pathway_size, n_pathways_label, thread_label) %>%
    summarize(
      intercept = {
        fit <- stats::lm(time_sec_median ~ B)
        as.numeric(stats::coef(fit)[1])
      },
      slope_per_B = {
        fit <- stats::lm(time_sec_median ~ B)
        as.numeric(stats::coef(fit)[2])
      },
      r2 = {
        fit <- stats::lm(time_sec_median ~ B)
        as.numeric(summary(fit)$r.squared)
      },
      .groups = "drop"
    ) %>%
    mutate(
      sec_per_10k_B = slope_per_B * 10000,
      sec_per_100k_B = slope_per_B * 100000
    )

  # Prediction table at practical B values.
  pred_B_targets <- sort(unique(c(cfg$B_values, 500000L, 1000000L)))
  prediction_df <- merge(
    time_model_df,
    data.frame(B_pred = pred_B_targets, stringsAsFactors = FALSE),
    by = NULL
  ) %>%
    mutate(
      pred_time_sec = pmax(intercept + slope_per_B * B_pred, 0),
      pred_time_min = pred_time_sec / 60
    ) %>%
    dplyr::select(
      m, n_pathways, threads, B_pred, pred_time_sec, pred_time_min,
      intercept, slope_per_B, r2
    ) %>%
    arrange(m, n_pathways, threads, B_pred)

  # Smooth line data for linear-fit overlay.
  pred_line_df <- merge(
    time_model_df,
    data.frame(B = seq(min(cfg$B_values), max(cfg$B_values), length.out = 100), stringsAsFactors = FALSE),
    by = NULL
  ) %>%
    mutate(pred_time_sec = pmax(intercept + slope_per_B * B, 0))

  # Linear model for memory prediction: mem_mb_mean = intercept + slope * B
  memory_model_df <- summary_df %>%
    group_by(m, n_pathways, threads, pathway_size, n_pathways_label, thread_label) %>%
    summarize(
      intercept_mem = {
        fit <- stats::lm(mem_mb_mean ~ B)
        as.numeric(stats::coef(fit)[1])
      },
      slope_mem_per_B = {
        fit <- stats::lm(mem_mb_mean ~ B)
        as.numeric(stats::coef(fit)[2])
      },
      r2_mem = {
        fit <- stats::lm(mem_mb_mean ~ B)
        as.numeric(summary(fit)$r.squared)
      },
      .groups = "drop"
    )

  mem_pred_line_df <- merge(
    memory_model_df,
    data.frame(B = seq(min(cfg$B_values), max(cfg$B_values), length.out = 100), stringsAsFactors = FALSE),
    by = NULL
  ) %>%
    mutate(pred_mem_mb = pmax(intercept_mem + slope_mem_per_B * B, 0))

  # Create time plot (faceted by pathway size and number of pathways; color = threads).
  p_time <- ggplot(summary_df, aes(x = B, y = time_plot, color = thread_label, group = thread_label)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    geom_errorbar(
      aes(ymin = pmax(time_plot - time_sec_sd, 1e-8), ymax = time_plot + time_sec_sd),
      width = 0.05 * max(summary_df$B),
      linewidth = 0.6
    ) +
    facet_grid(rows = vars(n_pathways_label), cols = vars(pathway_size), scales = "free_y", switch = "y") +
    scale_x_log10(
      breaks = cfg$B_values,
      labels = function(x) format(x, big.mark = ",", scientific = FALSE)
    ) +
    scale_y_log10() +
    scale_color_brewer(palette = "Dark2") +
    labs(
      title = "Block I: Runtime vs Resampling Budget",
      subtitle = "Full five-component omnibus workload (ACAT, Fisher, TFisher, minP, Stouffer), faceted by pathway size and number of pathways",
      x = "Number of null draws per pathway (B)",
      y = "Time (seconds)",
      color = "Threads"
    ) +
    block_a_theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
      legend.position = "top",
      strip.placement = "outside",
      strip.text.x = element_text(size = 18, face = "bold"),
      strip.text.y.left = element_text(size = 16, face = "bold", angle = 0),
      strip.background = element_rect(fill = "grey95", color = "grey80")
    )

  # Create memory plot (same facets for direct read-across).
  p_memory <- ggplot(summary_df, aes(x = B, y = mem_mb_mean, color = thread_label, group = thread_label)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    geom_errorbar(
      aes(ymin = pmax(mem_mb_mean - mem_mb_sd, 0), ymax = mem_mb_mean + mem_mb_sd),
      width = 0.05 * max(summary_df$B),
      linewidth = 0.6
    ) +
    facet_grid(rows = vars(n_pathways_label), cols = vars(pathway_size), scales = "free_y", switch = "y") +
    scale_x_log10(
      breaks = cfg$B_values,
      labels = function(x) format(x, big.mark = ",", scientific = FALSE)
    ) +
    scale_color_brewer(palette = "Dark2") +
    labs(
      title = "Block I: Memory vs Resampling Budget",
      subtitle = "Approximate peak resident memory during full five-component omnibus workload",
      x = "Number of null draws per pathway (B)",
      y = "Memory (MB)",
      color = "Threads"
    ) +
    block_a_theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
      legend.position = "top"
    )

  # Linear-fit runtime plot for prediction use.
  axis_breaks_linear <- c(50000L, 100000L, 250000L, 500000L, 1000000L)
  axis_labels_linear <- c("50,000", "100,000", "250,000", "500,000", "1,000,000")

  p_time_linear <- ggplot(summary_df, aes(x = B, y = time_sec_mean, color = thread_label, group = thread_label)) +
    geom_point(size = 2.8) +
    geom_line(linewidth = 0.8, alpha = 0.5) +
    geom_line(
      data = pred_line_df,
      aes(x = B, y = pred_time_sec, color = thread_label, group = thread_label),
      linewidth = 1.1,
      linetype = "dashed"
    ) +
    facet_grid(n_pathways_label ~ pathway_size, scales = "free_y", switch = "y") +
    scale_x_continuous(
      breaks = axis_breaks_linear,
      labels = axis_labels_linear
    ) +
    scale_color_brewer(palette = "Dark2") +
    labs(
      x = "Number of null draws per pathway (B)",
      y = "Time (seconds)",
      color = "Threads"
    ) +
    block_a_theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
      legend.position = "top",
      strip.placement = "outside",
      strip.text.x = element_text(size = 18, face = "bold"),
      strip.text.y.left = element_text(size = 16, face = "bold", angle = 90),
      strip.background = element_rect(fill = "grey95", color = "grey80")
    )

  p_memory_linear <- ggplot(summary_df, aes(x = B, y = mem_mb_mean, color = thread_label, group = thread_label)) +
    geom_point(size = 2.8) +
    geom_line(linewidth = 0.8, alpha = 0.5) +
    geom_line(
      data = mem_pred_line_df,
      aes(x = B, y = pred_mem_mb, color = thread_label, group = thread_label),
      linewidth = 1.1,
      linetype = "dashed"
    ) +
    facet_grid(n_pathways_label ~ pathway_size, scales = "free_y", switch = "y") +
    scale_x_continuous(
      breaks = axis_breaks_linear,
      labels = axis_labels_linear
    ) +
    scale_color_brewer(palette = "Dark2") +
    labs(
      x = "Number of null draws per pathway (B)",
      y = "Memory (MB)",
      color = "Threads"
    ) +
    block_a_theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
      legend.position = "none",
      strip.placement = "outside",
      strip.text.x = element_text(size = 18, face = "bold"),
      strip.text.y.left = element_text(size = 16, face = "bold", angle = 90),
      strip.background = element_rect(fill = "grey95", color = "grey80")
    )

  # Speedup plot: makes threading benefit explicit.
  p_speedup <- ggplot(
    summary_df %>% dplyr::filter(threads > 1L),
    aes(x = B, y = speedup_vs_1thread, color = thread_label, group = thread_label)
  ) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.8) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    facet_grid(n_pathways_label ~ pathway_size, scales = "free_y") +
    scale_x_log10(
      breaks = cfg$B_values,
      labels = function(x) format(x, big.mark = ",", scientific = FALSE)
    ) +
    scale_color_brewer(palette = "Dark2") +
    labs(
      title = "Block I: Thread Speedup vs B",
      subtitle = "Speedup > 1 means faster than 1-thread baseline",
      x = "Number of null draws per pathway (B)",
      y = "Speedup (vs 1 thread)",
      color = "Threads"
    ) +
    block_a_theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
      legend.position = "top"
    )

  # Create combined panel
  combined_plot <- (p_time / p_memory) +
    plot_annotation(tag_levels = "A")

  # Combined linear panel (A=time, B=memory) for manuscript.
  linear_panel_plot <- (p_time_linear / p_memory_linear) +
    plot_annotation(tag_levels = "A")

  # Save outputs
  out_csv <- file.path(cfg$results_dir, "block_i_benchmark_summary.csv")
  out_table_csv <- file.path(cfg$results_dir, "block_i_benchmark_table.csv")
  out_model_csv <- file.path(cfg$results_dir, "block_i_time_linear_model.csv")
  out_pred_csv <- file.path(cfg$results_dir, "block_i_time_predictions.csv")
  out_rds <- file.path(cfg$results_dir, "block_i_benchmark.rds")
  out_time_plot <- file.path(cfg$output_dir, "block_i_time.png")
  out_time_linear_plot <- file.path(cfg$output_dir, "block_i_time_linear_fit.png")
  out_memory_linear_plot <- file.path(cfg$output_dir, "block_i_memory_linear_fit.png")
  out_linear_panel_plot <- file.path(cfg$output_dir, "block_i_linear_panel_AB.png")
  out_speedup_plot <- file.path(cfg$output_dir, "block_i_speedup.png")
  out_memory_plot <- file.path(cfg$output_dir, "block_i_memory.png")
  out_combined_plot <- file.path(cfg$output_dir, "block_i_benchmark_panel.png")

  write.csv(summary_df, out_csv, row.names = FALSE)
  write.csv(table_df, out_table_csv, row.names = FALSE)
  write.csv(time_model_df, out_model_csv, row.names = FALSE)
  write.csv(prediction_df, out_pred_csv, row.names = FALSE)
  saveRDS(
    list(
      config = cfg,
      results_raw = results_df,
      summary = summary_df,
      table = table_df,
      time_model = time_model_df,
      time_predictions = prediction_df,
      memory_model = memory_model_df
    ),
    out_rds
  )

  ggsave(out_time_plot, p_time, width = 12, height = 8, dpi = 300, bg = "white")
  ggsave(out_time_linear_plot, p_time_linear, width = 12, height = 8, dpi = 300, bg = "white")
  ggsave(out_memory_linear_plot, p_memory_linear, width = 12, height = 8, dpi = 300, bg = "white")
  ggsave(out_linear_panel_plot, linear_panel_plot, width = 12, height = 14, dpi = 300, bg = "white")
  ggsave(out_speedup_plot, p_speedup, width = 12, height = 8, dpi = 300, bg = "white")
  ggsave(out_memory_plot, p_memory, width = 12, height = 8, dpi = 300, bg = "white")
  ggsave(out_combined_plot, combined_plot, width = 12, height = 14, dpi = 300, bg = "white")

  message("Block I outputs written:")
  message("  - ", out_csv)
  message("  - ", out_table_csv)
  message("  - ", out_model_csv)
  message("  - ", out_pred_csv)
  message("  - ", out_rds)
  message("  - ", out_time_plot)
  message("  - ", out_time_linear_plot)
  message("  - ", out_memory_linear_plot)
  message("  - ", out_linear_panel_plot)
  message("  - ", out_speedup_plot)
  message("  - ", out_memory_plot)
  message("  - ", out_combined_plot)

  invisible(list(
    config = cfg,
    results_raw = results_df,
    summary = summary_df,
    table = table_df,
    time_model = time_model_df,
    time_predictions = prediction_df,
    memory_model = memory_model_df,
    time_plot = p_time,
    time_linear_plot = p_time_linear,
    memory_linear_plot = p_memory_linear,
    linear_panel_plot = linear_panel_plot,
    speedup_plot = p_speedup,
    memory_plot = p_memory,
    combined_plot = combined_plot
  ))
}

if (sys.nframe() == 0L) {
  run_all_figs()
}
