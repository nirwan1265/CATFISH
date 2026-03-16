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

run_all_figs <- function(run_block_a_flag = TRUE,
                         run_block_b_flag = FALSE,
                         run_block_c_flag = FALSE,
                         run_block_d_flag = FALSE,
                         run_block_e_flag = FALSE,
                         run_block_f_flag = FALSE,
                         run_block_g_flag = FALSE,
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
  invisible(out)
}

if (sys.nframe() == 0L) {
  run_all_figs()
}
