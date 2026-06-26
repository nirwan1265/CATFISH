# ==============================================================================
# SIMULATION FRAMEWORK FOR CATFISH PATHWAY TESTS (v2)
# ==============================================================================
# This script implements comprehensive simulations to evaluate:
#   Block A: Power across different signal shapes (sparse vs dense)
#   Block B: "Broken component" stress test
#   Block C: Missing correlation / incomplete Sigma scenarios
#   Block D: Adaptive omnibus with train/test split
#   Block E: Leave-one-out omnibus sensitivity (MVN)
#   Block F: Random two-component dropout (MVN)
#
# Key fixes from v1:
#   - Uses actual TFisher package for correct null distribution
#   - Precomputes null distributions per (m, Sigma) for massive speedup
#   - Block B stress-tests robustness to broken component tests
#   - Block D implements adaptive drop/downweight with proper train/test split
#   - Uses two-sided p-values for mixed direction scenario
# ==============================================================================

devtools::load_all(".")
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(Matrix)
library(MASS)
library(data.table)

# Check for TFisher package
if (!requireNamespace("TFisher", quietly = TRUE)) {
  warning("TFisher package not installed. Using approximation. ",
          "Install with: install.packages('TFisher')")
  use_tfisher_pkg <- FALSE
} else {
  library(TFisher)
  use_tfisher_pkg <- TRUE
}

set.seed(42)

# ==============================================================================
# GLOBAL PARAMETERS
# ==============================================================================

n_sims_null <- 500L
n_sims_alt <- 500L
b_perm <- 1000L
#tau_grid <- c(0.000000005)
tau_grid <- c(0.05, 0.01, 0.1)
min_p <- 1e-15
seed_base <- 12345L
lambda_boot <- 200L

# Pathway sizes to test
pathway_sizes <- c(5L, 25L, 50L)

# Effect sizes (delta) for alternative hypotheses
effect_sizes <- c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0)

# Alpha levels for testing
alpha_levels <- c(0.05, 0.01)

# Cap for plotting lambda (keeps figures readable)
lambda_cap <- 2.5

# ==============================================================================
# PLOTTING THEME
# ==============================================================================

sim_theme <- theme_minimal(base_size = 24) +
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

    ## <-- CHANGED (was element_blank())
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.25),

    legend.position = "top",
    legend.title    = element_blank(),
    legend.text     = element_text(size = 24)
  )
# ==============================================================================
# HELPER FUNCTIONS: CORRELATION MATRIX GENERATION
# ==============================================================================

make_cor_independent <- function(m) {
  diag(m)
}

make_cor_block <- function(m, block_size = 20L, rho = 0.3) {
  sigma <- diag(m)
  n_blocks <- ceiling(m / block_size)

  for (b in seq_len(n_blocks)) {
    start_idx <- (b - 1) * block_size + 1
    end_idx <- min(b * block_size, m)
    block_genes <- start_idx:end_idx

    for (i in block_genes) {
      for (j in block_genes) {
        if (i != j) {
          sigma[i, j] <- rho
        }
      }
    }
  }
  sigma
}

make_cor_ar1 <- function(m, rho = 0.5) {
  sigma <- matrix(0, m, m)
  for (i in seq_len(m)) {
    for (j in seq_len(m)) {
      sigma[i, j] <- rho^abs(i - j)
    }
  }
  sigma
}

make_cor_exchangeable <- function(m, rho = 0.3) {
  sigma <- matrix(rho, m, m)
  diag(sigma) <- 1
  sigma
}

ensure_pd <- function(sigma, tol = 1e-6) {
  eig <- eigen(sigma, symmetric = TRUE, only.values = TRUE)$values
  if (min(eig) < tol) {
    sigma <- nearPD(sigma, corr = TRUE)$mat
    sigma <- as.matrix(sigma)
  }
  sigma
}

# ==============================================================================
# HELPER FUNCTIONS: SIMULATION OF GENE Z-SCORES
# ==============================================================================

simulate_gene_z <- function(m, sigma, n_causal = 0, effect_size = 0,
                            causal_direction = "positive") {
  sigma_pd <- ensure_pd(sigma)
  z <- as.numeric(mvrnorm(1, mu = rep(0, m), Sigma = sigma_pd))

  if (n_causal > 0 && effect_size > 0) {
    causal_idx <- sample(seq_len(m), min(n_causal, m))

    if (causal_direction == "positive") {
      z[causal_idx] <- z[causal_idx] + effect_size
    } else if (causal_direction == "negative") {
      z[causal_idx] <- z[causal_idx] - effect_size
    } else if (causal_direction == "mixed") {
      n_pos <- ceiling(length(causal_idx) / 2)
      z[causal_idx[1:n_pos]] <- z[causal_idx[1:n_pos]] + effect_size
      if (n_pos < length(causal_idx)) {
        z[causal_idx[(n_pos + 1):length(causal_idx)]] <-
          z[causal_idx[(n_pos + 1):length(causal_idx)]] - effect_size
      }
    }
  }

  names(z) <- paste0("GENE", seq_len(m))
  z
}

z_to_p_onesided <- function(z) {
  p <- pnorm(z, lower.tail = FALSE)
  pmax(pmin(p, 1 - min_p), min_p)
}

z_to_p_twosided <- function(z) {
  p <- 2 * pnorm(-abs(z))
  pmax(pmin(p, 1 - min_p), min_p)
}

create_gene_results <- function(z_scores, two_sided = FALSE) {
  if (two_sided) {
    p_vals <- z_to_p_twosided(z_scores)
  } else {
    p_vals <- z_to_p_onesided(z_scores)
  }

  data.frame(
    GENE = names(z_scores),
    P = p_vals,
    ZSTAT = z_scores,
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# COMPONENT TEST FUNCTIONS (USING PROPER IMPLEMENTATIONS)
# ==============================================================================

compute_acat <- function(p_vals) {
  p_use <- p_vals
  if (exists("fix_p_for_acat", mode = "function")) {
    p_use <- fix_p_for_acat(p_use, min_p = min_p)
  } else {
    p_use <- pmax(pmin(p_use, 1 - min_p), min_p)
  }

  if (requireNamespace("ACAT", quietly = TRUE)) {
    p_acat <- ACAT::ACAT(Pvals = p_use)
  } else {
    tan_vals <- tan((0.5 - p_use) * pi)
    tstat <- mean(tan_vals)
    p_acat <- 0.5 - atan(tstat) / pi
  }

  max(min(p_acat, 1 - min_p), min_p)
}

compute_fisher <- function(p_vals) {
  p_use <- p_vals
  if (exists("fix_p_for_acat", mode = "function")) {
    p_use <- fix_p_for_acat(p_use, min_p = min_p)
  } else {
    p_use <- pmax(pmin(p_use, 1 - min_p), min_p)
  }

  m <- length(p_use)
  t_fisher <- -2 * sum(log(p_use))
  p_fisher <- pchisq(t_fisher, df = 2 * m, lower.tail = FALSE)
  max(min(p_fisher, 1 - min_p), min_p)
}

compute_stouffer <- function(z_vals, alternative = "greater") {
  m <- length(z_vals)
  z_stouffer <- sum(z_vals) / sqrt(m)

  if (alternative == "greater") {
    p <- pnorm(z_stouffer, lower.tail = FALSE)
  } else if (alternative == "less") {
    p <- pnorm(z_stouffer, lower.tail = TRUE)
  } else {
    p <- 2 * pnorm(-abs(z_stouffer))
  }
  max(min(p, 1 - min_p), min_p)
}

compute_minp <- function(p_vals) {
  p_use <- p_vals
  if (exists("fix_p_for_acat", mode = "function")) {
    p_use <- fix_p_for_acat(p_use, min_p = min_p)
  } else {
    p_use <- pmax(pmin(p_use, 1 - min_p), min_p)
  }

  m <- length(p_use)
  p_min <- min(p_use)

  if (requireNamespace("metap", quietly = TRUE)) {
    p_minp <- tryCatch(metap::minimump(p_use)$p, error = function(e) NA_real_)
    if (!is.finite(p_minp) || is.na(p_minp)) {
      p_minp <- 1 - (1 - p_min)^m
    }
  } else {
    p_minp <- 1 - (1 - p_min)^m
  }

  max(min(p_minp, 1 - min_p), min_p)
}

# TFisher using actual TFisher package (CORRECT implementation)
compute_tfisher <- function(p_vals, tg = tau_grid) {
  m <- length(p_vals)

  if (use_tfisher_pkg) {
    tg_use <- as.numeric(tg)
    tg_use <- tg_use[is.finite(tg_use) & !is.na(tg_use) & tg_use > 0]
    if (!length(tg_use)) return(NA_real_)
    tg_use <- sort(unique(tg_use), decreasing = FALSE)

    p_use <- p_vals
    if (exists("fix_p_for_acat", mode = "function")) {
      p_use <- fix_p_for_acat(p_use, min_p = min_p)
    } else {
      p_use <- pmax(pmin(p_use, 1 - min_p), min_p)
    }

    tryCatch({
      omni_out <- TFisher::stat.soft.omni(p = p_use, TAU1 = tg_use, M = NULL)
      Wo <- as.numeric(omni_out$omni)
      p_tfisher <- as.numeric(TFisher::p.soft.omni(q = Wo, n = m, TAU1 = tg_use, M = NULL))
      return(max(p_tfisher, min_p))
    }, error = function(e) {
      tryCatch({
        # Fallback to single-tau TFisher using TFisher's analytic null
        tau <- tg_use[1]
        stat <- TFisher::stat.soft(p = p_use, tau1 = tau)
        cdf <- TFisher::p.soft(q = stat, n = m, tau1 = tau, M = NULL)
        p_rt <- 1 - as.numeric(cdf)
        return(max(p_rt, min_p))
      }, error = function(e2) {
        return(0.5)
      })
    })
  } else {
    # Fallback approximation (less accurate)
    tfisher_ps <- sapply(tg, function(tau) {
      below_tau <- p_vals[p_vals < tau]
      if (length(below_tau) == 0) return(1)
      w_tau <- -2 * sum(log(below_tau))
      df_approx <- 2 * length(below_tau)
      pchisq(w_tau, df = df_approx, lower.tail = FALSE)
    })
    p_tfisher <- min(tfisher_ps) * length(tg)  # Bonferroni
    return(max(min(p_tfisher, 1 - min_p), min_p))
  }
}

compute_omnibus <- function(component_ps) {
  compute_acat(component_ps)
}

# Compute all component p-values analytically
compute_components_analytic <- function(gene_results, tg = tau_grid,
                                        stouffer_alt = "greater") {
  p_vals <- gene_results$P
  z_vals <- gene_results$ZSTAT

  p_acat <- compute_acat(p_vals)
  p_fisher <- compute_fisher(p_vals)
  p_stouffer <- compute_stouffer(z_vals, stouffer_alt)
  p_minp <- compute_minp(p_vals)
  p_tfisher <- compute_tfisher(p_vals, tg)

  component_ps <- c(p_acat, p_fisher, p_tfisher, p_minp, p_stouffer)
  p_omnibus <- compute_omnibus(component_ps)

  list(
    acat = p_acat,
    fisher = p_fisher,
    tfisher = p_tfisher,
    minp = p_minp,
    stouffer = p_stouffer,
    omnibus = p_omnibus
  )
}

#' Compute all component p-values with MVN calibration
compute_components_mvn <- function(gene_results, Sigma, B = b_perm,
                                   tg = tau_grid, two_sided = FALSE,
                                   stouffer_alt = "greater") {
  m <- nrow(gene_results)

  # Observed (analytic) p-values
  obs <- compute_components_analytic(gene_results, tg = tg, stouffer_alt = stouffer_alt)

  # Null simulations with correlation
  Sigma_pd <- ensure_pd(Sigma)
  null_stats <- matrix(NA, nrow = B, ncol = 6)
  colnames(null_stats) <- c("acat", "fisher", "tfisher", "minp", "stouffer", "omnibus")

  for (b in seq_len(B)) {
    Z_null <- as.numeric(mvrnorm(1, mu = rep(0, m), Sigma = Sigma_pd))
    P_null <- if (two_sided) z_to_p_twosided(Z_null) else z_to_p_onesided(Z_null)

    gene_null <- data.frame(
      GENE = paste0("GENE", seq_len(m)),
      P = P_null,
      ZSTAT = Z_null,
      stringsAsFactors = FALSE
    )

    null_res <- compute_components_analytic(gene_null, tg = tg, stouffer_alt = stouffer_alt)
    null_stats[b, ] <- unlist(null_res)
  }

  emp_p <- function(obs_val, null_vals) {
    (1 + sum(null_vals <= obs_val, na.rm = TRUE)) / (sum(!is.na(null_vals)) + 1)
  }

  list(
    acat     = emp_p(obs$acat, null_stats[, "acat"]),
    fisher   = emp_p(obs$fisher, null_stats[, "fisher"]),
    tfisher  = emp_p(obs$tfisher, null_stats[, "tfisher"]),
    minp     = emp_p(obs$minp, null_stats[, "minp"]),
    stouffer = emp_p(obs$stouffer, null_stats[, "stouffer"]),
    omnibus  = emp_p(obs$omnibus, null_stats[, "omnibus"])
  )
}

# ==============================================================================
# PRECOMPUTED NULL DISTRIBUTIONS (KEY SPEEDUP)
# ==============================================================================

# Precompute null distribution for a given (m, Sigma) condition
# This is computed ONCE and reused across all replicates
precompute_null_distribution <- function(m, sigma, b = 500L, tg = tau_grid,
                                         two_sided = FALSE, seed = NULL,
                                         stouffer_alt = "greater") {
  if (!is.null(seed)) set.seed(seed)

  sigma_pd <- ensure_pd(sigma)

  null_stats <- matrix(NA, nrow = b, ncol = 6)
  colnames(null_stats) <- c("acat", "fisher", "tfisher", "minp",
                            "stouffer", "omnibus")

  for (i in seq_len(b)) {
    z_null <- as.numeric(mvrnorm(1, mu = rep(0, m), Sigma = sigma_pd))

    if (two_sided) {
      p_null <- z_to_p_twosided(z_null)
      s_alt <- "two.sided"
    } else {
      p_null <- z_to_p_onesided(z_null)
      s_alt <- stouffer_alt
    }

    gene_null <- data.frame(
      GENE = paste0("GENE", seq_len(m)),
      P = p_null,
      ZSTAT = z_null,
      stringsAsFactors = FALSE
    )

    res_null <- compute_components_analytic(gene_null, tg, s_alt)
    null_stats[i, ] <- unlist(res_null)
  }

  null_stats
}

# Calibrate p-values using precomputed null distribution
calibrate_with_null <- function(obs_pvals, null_pvals) {
  calibrated <- sapply(names(obs_pvals), function(method) {
    obs <- obs_pvals[[method]]
    null <- null_pvals[, method]
    (1 + sum(null <= obs, na.rm = TRUE)) / (sum(!is.na(null)) + 1)
  })
  as.list(calibrated)
}

empirical_p_from_null <- function(obs_values, null_values) {
  null_clean <- sort(null_values[is.finite(null_values) & !is.na(null_values)])
  if (!length(null_clean)) {
    return(rep(NA_real_, length(obs_values)))
  }
  (1 + findInterval(obs_values, null_clean, rightmost.closed = TRUE)) /
    (length(null_clean) + 1)
}

calibrate_component_matrix <- function(obs_mat, null_mat) {
  out <- matrix(
    NA_real_,
    nrow = nrow(obs_mat),
    ncol = ncol(obs_mat),
    dimnames = dimnames(obs_mat)
  )

  common_cols <- intersect(colnames(obs_mat), colnames(null_mat))
  for (comp in common_cols) {
    out[, comp] <- empirical_p_from_null(obs_mat[, comp], null_mat[, comp])
  }

  out
}

# ==============================================================================
# HELPER FUNCTIONS: METRICS
# ==============================================================================

compute_lambda <- function(p_values) {
  p <- p_values[!is.na(p_values) & p_values > 0 & p_values < 1]
  if (length(p) < 10) return(NA_real_)
  median(qchisq(1 - p, df = 1)) / qchisq(0.5, df = 1)
}

compute_lambda_se <- function(p_values, n_boot = lambda_boot, seed = NULL) {
  p <- p_values[!is.na(p_values) & p_values > 0 & p_values < 1]
  if (length(p) < 10 || n_boot <= 0) return(NA_real_)
  if (!is.null(seed)) set.seed(seed)
  boot_vals <- replicate(n_boot, {
    p_boot <- sample(p, length(p), replace = TRUE)
    compute_lambda(p_boot)
  })
  sd(boot_vals, na.rm = TRUE)
}

compute_type1 <- function(p_values, alpha = 0.05) {
  p <- p_values[!is.na(p_values)]
  if (length(p) == 0) return(NA_real_)
  mean(p < alpha)
}

compute_prop_se <- function(p_values, alpha = 0.05) {
  p <- p_values[!is.na(p_values)]
  n <- length(p)
  if (n == 0) return(NA_real_)
  phat <- mean(p < alpha)
  sqrt(phat * (1 - phat) / n)
}

compute_power <- function(p_values, alpha = 0.05) {
  compute_type1(p_values, alpha)
}

# ==============================================================================
# BLOCK A: POWER ACROSS SIGNAL SHAPES
# ==============================================================================

run_block_a <- function(n_sims = 200L, b = 500L, seed = 12345L,
                        pathway_sizes_arg = pathway_sizes,
                        effect_sizes_arg = c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0),
                        verbose = TRUE) {

  set.seed(seed)

  signal_types <- list(
    "Sparse_Strong" = list(prop_causal = 0.05, direction = "positive"),
    "Sparse_Moderate" = list(prop_causal = 0.10, direction = "positive"),
    "Dense_Weak" = list(prop_causal = 0.30, direction = "positive"),
    "Dense_Strong" = list(prop_causal = 0.50, direction = "positive"),
    "Mixed_Direction" = list(prop_causal = 0.30, direction = "mixed")
  )

  cor_types <- list(
    "LD_moderate" = function(m) make_cor_block(m, block_size = 20L, rho = 0.3),
    "LD_strong" = function(m) make_cor_block(m, block_size = 20L, rho = 0.7),
    "LD_independent" = function(m) make_cor_independent(m)
  )

  results_list <- list()

  total_conditions <- length(pathway_sizes_arg) * length(cor_types) *
                      length(signal_types) * length(effect_sizes_arg)
  cond_idx <- 0

  if (verbose) cat("\n=== BLOCK A: Power Across Signal Shapes ===\n")
  if (verbose) cat("Precomputing null distributions...\n")

  null_cache <- list()
  for (m in pathway_sizes_arg) {
    for (cor_name in names(cor_types)) {
      sigma <- cor_types[[cor_name]](m)
      cache_key <- paste(m, cor_name, sep = "_")

      if (verbose) cat(sprintf("  Precomputing: m=%d, cor=%s\n", m, cor_name))

      null_cache[[paste0(cache_key, "_1sided")]] <-
        precompute_null_distribution(
          m, sigma, b = b, two_sided = FALSE,
          seed = seed + m * 1000
        )
      null_cache[[paste0(cache_key, "_2sided")]] <-
        precompute_null_distribution(
          m, sigma, b = b, two_sided = TRUE,
          seed = seed + m * 1000 + 1
        )
    }
  }

  if (verbose) cat("\nRunning power simulations...\n")
  if (verbose) cat("Total iterations:", total_conditions * n_sims, "\n\n")

  for (m in pathway_sizes_arg) {
    for (cor_name in names(cor_types)) {
      sigma <- cor_types[[cor_name]](m)
      cache_key <- paste(m, cor_name, sep = "_")

      for (sig_name in names(signal_types)) {
        sig_params <- signal_types[[sig_name]]
        n_causal <- max(1, round(m * sig_params$prop_causal))
        two_sided <- (sig_params$direction == "mixed")
        stouffer_alt <- if (two_sided) "two.sided" else "greater"

        null_key <- paste0(cache_key, ifelse(two_sided, "_2sided", "_1sided"))
        null_dist <- null_cache[[null_key]]

        for (delta in effect_sizes_arg) {
          cond_idx <- cond_idx + 1
          if (verbose && cond_idx %% 10 == 0) {
            cat(sprintf("  Condition %d / %d\n", cond_idx, total_conditions))
          }

          p_analytic <- matrix(NA, nrow = n_sims, ncol = 6)
          p_mvn <- matrix(NA, nrow = n_sims, ncol = 6)
          colnames(p_analytic) <- colnames(p_mvn) <-
            c("acat", "fisher", "tfisher", "minp", "stouffer", "omnibus")

          for (s in seq_len(n_sims)) {
            z <- simulate_gene_z(m, sigma, n_causal, delta, sig_params$direction)
            gene_res <- create_gene_results(z, two_sided = two_sided)

            res_ana <- compute_components_analytic(gene_res, stouffer_alt = stouffer_alt)
            p_analytic[s, ] <- unlist(res_ana)

            res_cal <- calibrate_with_null(res_ana, null_dist)
            p_mvn[s, ] <- unlist(res_cal)
          }

          for (method in colnames(p_analytic)) {
            results_list[[length(results_list) + 1]] <- data.frame(
              block = "A",
              pathway_size = m,
              cor_structure = cor_name,
              signal_type = sig_name,
              effect_size = delta,
              n_causal = n_causal,
              method = method,
              calibration = "analytic",
              power_05 = compute_power(p_analytic[, method], 0.05),
              power_01 = compute_power(p_analytic[, method], 0.01),
              power_05_se = compute_prop_se(p_analytic[, method], 0.05),
              power_01_se = compute_prop_se(p_analytic[, method], 0.01),
              lambda = compute_lambda(p_analytic[, method]),
              lambda_se = compute_lambda_se(p_analytic[, method]),
              stringsAsFactors = FALSE
            )

            results_list[[length(results_list) + 1]] <- data.frame(
              block = "A",
              pathway_size = m,
              cor_structure = cor_name,
              signal_type = sig_name,
              effect_size = delta,
              n_causal = n_causal,
              method = method,
              calibration = "mvn",
              power_05 = compute_power(p_mvn[, method], 0.05),
              power_01 = compute_power(p_mvn[, method], 0.01),
              power_05_se = compute_prop_se(p_mvn[, method], 0.05),
              power_01_se = compute_prop_se(p_mvn[, method], 0.01),
              lambda = compute_lambda(p_mvn[, method]),
              lambda_se = compute_lambda_se(p_mvn[, method]),
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }

  if (verbose) cat("Block A complete.\n\n")
  bind_rows(results_list)
}

# Block A null simulations
run_block_a_null <- function(n_sims = 500L, b = 500L, seed = 13345L,
                             pathway_sizes_arg = pathway_sizes,
                             verbose = TRUE) {

  set.seed(seed)

  cor_types <- list(
    "LD_moderate" = function(m) make_cor_block(m, block_size = 20L, rho = 0.3),
    "LD_strong" = function(m) make_cor_block(m, block_size = 20L, rho = 0.7),
    "LD_independent" = function(m) make_cor_independent(m)
  )

  results_list <- list()

  total_conditions <- length(pathway_sizes_arg) * length(cor_types)
  cond_idx <- 0

  if (verbose) cat("\n=== BLOCK A (Null): Calibration Under Null ===\n")
  if (verbose) cat("Precomputing null distributions...\n")

  null_cache <- list()
  for (m in pathway_sizes_arg) {
    for (cor_name in names(cor_types)) {
      sigma <- cor_types[[cor_name]](m)
      cache_key <- paste(m, cor_name, sep = "_")
      if (verbose) cat(sprintf("  Precomputing: m=%d, cor=%s\n", m, cor_name))
      null_cache[[cache_key]] <-
        precompute_null_distribution(m, sigma, b = b, two_sided = FALSE,
                                     seed = seed + m * 1000)
    }
  }

  if (verbose) cat("\nRunning null simulations...\n")
  if (verbose) cat("Total iterations:", total_conditions * n_sims, "\n\n")

  for (m in pathway_sizes_arg) {
    for (cor_name in names(cor_types)) {
      cond_idx <- cond_idx + 1
      if (verbose) {
        cat(sprintf("  Condition %d / %d: m=%d, cor=%s\n",
                    cond_idx, total_conditions, m, cor_name))
      }

      sigma <- cor_types[[cor_name]](m)
      cache_key <- paste(m, cor_name, sep = "_")
      null_dist <- null_cache[[cache_key]]

      component_cols <- c("acat", "fisher", "tfisher", "minp", "stouffer")
      combined_cols <- c(component_cols, "omnibus_combined")
      p_analytic <- matrix(NA, nrow = n_sims, ncol = length(combined_cols))
      p_mvn <- matrix(NA, nrow = n_sims, ncol = length(combined_cols))
      p_mvn_alone <- numeric(n_sims)
      colnames(p_analytic) <- colnames(p_mvn) <- combined_cols

      for (s in seq_len(n_sims)) {
        z <- simulate_gene_z(m, sigma, n_causal = 0, effect_size = 0)
        gene_res <- create_gene_results(z)

        res_ana <- compute_components_analytic(gene_res)
        p_analytic[s, component_cols] <- unlist(res_ana[component_cols])
        p_analytic[s, "omnibus_combined"] <- res_ana$omnibus

        res_cal <- calibrate_with_null(res_ana, null_dist)
        p_mvn[s, component_cols] <- unlist(res_cal[component_cols])
        p_mvn[s, "omnibus_combined"] <- compute_omnibus(unlist(res_cal[component_cols]))
        p_mvn_alone[s] <- res_cal$omnibus
      }

      for (method in colnames(p_analytic)) {
        results_list[[length(results_list) + 1]] <- data.frame(
          block = "A_null",
          pathway_size = m,
          cor_structure = cor_name,
          method = method,
          calibration = "analytic",
          lambda = compute_lambda(p_analytic[, method]),
          lambda_se = compute_lambda_se(p_analytic[, method]),
          type1_05 = compute_type1(p_analytic[, method], 0.05),
          type1_05_se = compute_prop_se(p_analytic[, method], 0.05),
          type1_01 = compute_type1(p_analytic[, method], 0.01),
          type1_01_se = compute_prop_se(p_analytic[, method], 0.01),
          mean_p = mean(p_analytic[, method], na.rm = TRUE),
          median_p = median(p_analytic[, method], na.rm = TRUE),
          stringsAsFactors = FALSE
        )

        results_list[[length(results_list) + 1]] <- data.frame(
          block = "A_null",
          pathway_size = m,
          cor_structure = cor_name,
          method = method,
          calibration = "mvn",
          lambda = compute_lambda(p_mvn[, method]),
          lambda_se = compute_lambda_se(p_mvn[, method]),
          type1_05 = compute_type1(p_mvn[, method], 0.05),
          type1_05_se = compute_prop_se(p_mvn[, method], 0.05),
          type1_01 = compute_type1(p_mvn[, method], 0.01),
          type1_01_se = compute_prop_se(p_mvn[, method], 0.01),
          mean_p = mean(p_mvn[, method], na.rm = TRUE),
          median_p = median(p_mvn[, method], na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }

      results_list[[length(results_list) + 1]] <- data.frame(
        block = "A_null",
        pathway_size = m,
        cor_structure = cor_name,
        method = "omnibus_alone",
        calibration = "mvn_alone",
        lambda = compute_lambda(p_mvn_alone),
        lambda_se = compute_lambda_se(p_mvn_alone),
        type1_05 = compute_type1(p_mvn_alone, 0.05),
        type1_05_se = compute_prop_se(p_mvn_alone, 0.05),
        type1_01 = compute_type1(p_mvn_alone, 0.01),
        type1_01_se = compute_prop_se(p_mvn_alone, 0.01),
        mean_p = mean(p_mvn_alone, na.rm = TRUE),
        median_p = median(p_mvn_alone, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }

  if (verbose) cat("Block A (null) complete.\n\n")
  bind_rows(results_list)
}

# ==============================================================================
# BLOCK B: BROKEN COMPONENT STRESS TEST (ALL COMPONENTS AFFECTED)
# ==============================================================================

run_block_b <- function(n_sims = 500L, b = 500L, seed = 22345L,
                        pathway_sizes_arg = pathway_sizes,
                        rho_values = c(0, 0.3, 0.7),
                        broken_components = c("stouffer", "acat", "minp"),
                        broken_power = 0.7,
                        verbose = TRUE) {

  set.seed(seed)

  results_list <- list()

  if (verbose) cat("\n=== BLOCK B: Analytic vs MVN Under Correlation ===\n")
  if (verbose) cat("Testing: Correlation ignored + extra broken components\n\n")

  if (verbose) cat("Precomputing null distributions...\n")
  null_cache <- list()
  for (m in pathway_sizes_arg) {
    for (rho in rho_values) {
      sigma_true <- make_cor_block(m, block_size = 20L, rho = rho)
      cache_key <- paste(m, rho, sep = "_")
      if (verbose) cat(sprintf("  Precomputing: m=%d, rho=%.1f\n", m, rho))
      null_cache[[cache_key]] <-
        precompute_null_distribution(m, sigma_true, b = b, two_sided = FALSE,
                                     seed = seed + m * 1000 + round(rho * 100))
    }
  }

  for (m in pathway_sizes_arg) {
    for (rho in rho_values) {
      if (verbose) cat(sprintf("  m=%d, rho=%.1f\n", m, rho))

      sigma_true <- make_cor_block(m, block_size = 20L, rho = rho)
      cache_key <- paste(m, rho, sep = "_")
      null_dist <- null_cache[[cache_key]]

      p_analytic <- matrix(NA, nrow = n_sims, ncol = 6)
      p_mvn <- matrix(NA, nrow = n_sims, ncol = 6)
      colnames(p_analytic) <- colnames(p_mvn) <-
        c("acat", "fisher", "tfisher", "minp", "stouffer", "omnibus")

      for (s in seq_len(n_sims)) {
        z <- simulate_gene_z(m, sigma_true, n_causal = 0, effect_size = 0)
        gene_res <- create_gene_results(z)

        res_ana <- compute_components_analytic(gene_res)
        if (length(broken_components) > 0) {
          res_ana_broken <- apply_component_breakage(
            res_ana = res_ana,
            components_to_break = broken_components,
            break_type = "anti",
            severity = broken_power
          )
          p_analytic[s, ] <- unlist(res_ana_broken)
        } else {
          p_analytic[s, ] <- unlist(res_ana)
        }

        res_cal <- calibrate_with_null(res_ana, null_dist)
        p_mvn[s, ] <- unlist(res_cal)
      }

      for (method in colnames(p_analytic)) {
        results_list[[length(results_list) + 1]] <- data.frame(
          block = "B",
          pathway_size = m,
          rho = rho,
          method = method,
          scenario = "analytic_broken",
          lambda = compute_lambda(p_analytic[, method]),
          lambda_se = compute_lambda_se(p_analytic[, method]),
          type1_05 = compute_type1(p_analytic[, method], 0.05),
          type1_05_se = compute_prop_se(p_analytic[, method], 0.05),
          type1_01 = compute_type1(p_analytic[, method], 0.01),
          type1_01_se = compute_prop_se(p_analytic[, method], 0.01),
          mean_p = mean(p_analytic[, method], na.rm = TRUE),
          stringsAsFactors = FALSE
        )

        results_list[[length(results_list) + 1]] <- data.frame(
          block = "B",
          pathway_size = m,
          rho = rho,
          method = method,
          scenario = "mvn",
          lambda = compute_lambda(p_mvn[, method]),
          lambda_se = compute_lambda_se(p_mvn[, method]),
          type1_05 = compute_type1(p_mvn[, method], 0.05),
          type1_05_se = compute_prop_se(p_mvn[, method], 0.05),
          type1_01 = compute_type1(p_mvn[, method], 0.01),
          type1_01_se = compute_prop_se(p_mvn[, method], 0.01),
          mean_p = mean(p_mvn[, method], na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (verbose) cat("Block B complete.\n\n")

  bind_rows(results_list)
}

# ==============================================================================
# BLOCK B EXTENDED: COMPREHENSIVE BROKEN COMPONENT STRESS TEST
# ==============================================================================
# Tests omnibus robustness when individual components are miscalibrated:
#   - Anti-conservative (inflated): p-values made artificially small
#   - Conservative (deflated): p-values made artificially large
#   - Varying severity levels
#   - Breaking 1, 2, 3, or 4 components simultaneously
# ==============================================================================

#' Break a p-value in anti-conservative direction (inflate type I error)
#' @param p original p-value
#' @param severity how much to inflate (0.3 = severe, 0.5 = moderate, 0.7 = mild)
break_anticonservative <- function(p, severity = 0.5) {
  # For 0 < p < 1, using an exponent > 1 shrinks p-values and inflates type I error.
  # Map smaller severities to stronger breakage while preserving the existing API.
  exponent <- 1 / severity
  pmax(pmin(p^exponent, 1 - min_p), min_p)
}

#' Break a p-value in conservative direction (deflate power)
#' @param p original p-value
#' @param severity how much to deflate (2 = mild, 5 = moderate, 10 = severe)
break_conservative <- function(p, severity = 2) {
  # Multiply p by factor > 1, cap at 1
  pmax(pmin(p * severity, 1 - min_p), min_p)
}

#' Apply breakage to specific components
#' @param res_ana list of analytic p-values from compute_components_analytic
#' @param components_to_break character vector of component names to break
#' @param break_type "anti" for anti-conservative, "cons" for conservative
#' @param severity severity parameter for the breakage
apply_component_breakage <- function(res_ana, components_to_break,
                                     break_type = "anti", severity = 0.5) {
  res_broken <- res_ana

  for (comp in components_to_break) {
    if (comp %in% names(res_broken) && comp != "omnibus") {
      p_raw <- res_broken[[comp]]
      if (break_type == "anti") {
        res_broken[[comp]] <- break_anticonservative(p_raw, severity)
      } else {
        res_broken[[comp]] <- break_conservative(p_raw, severity)
      }
    }
  }

  # Recompute omnibus from the (possibly broken) components
  component_ps <- c(res_broken$acat, res_broken$fisher, res_broken$tfisher,
                    res_broken$minp, res_broken$stouffer)
  res_broken$omnibus <- compute_omnibus(component_ps)

  res_broken
}

run_block_b_extended <- function(n_sims = 500L, b = 500L, seed = 32345L,
                                  pathway_sizes_arg = c(25L, 50L),
                                  rho = 0.3,
                                  verbose = TRUE) {

  set.seed(seed)
  results_list <- list()

  if (verbose) cat("\n=== BLOCK B EXTENDED: Comprehensive Broken Component Stress Test ===\n")

  # Define breakage scenarios
  all_components <- c("acat", "fisher", "tfisher", "minp", "stouffer")

  # Anti-conservative severities (power < 1 makes p smaller)
  anti_severities <- list(
    "anti_mild" = 0.7,
    "anti_moderate" = 0.5,
    "anti_severe" = 0.3
  )

  # Conservative severities (multiply by factor > 1)
  cons_severities <- list(
    "cons_mild" = 2,
    "cons_moderate" = 5,
    "cons_severe" = 10
  )

  # Number of components to break simultaneously
  n_broken_options <- c(1, 2, 3, 4)

  # Precompute null distributions
  if (verbose) cat("Precomputing null distributions...\n")
  null_cache <- list()
  for (m in pathway_sizes_arg) {
    sigma_true <- make_cor_block(m, block_size = 20L, rho = rho)
    cache_key <- as.character(m)
    if (verbose) cat(sprintf("  Precomputing: m=%d, rho=%.1f\n", m, rho))
    null_cache[[cache_key]] <-
      precompute_null_distribution(m, sigma_true, b = b, two_sided = FALSE,
                                   seed = seed + m * 1000)
  }

  # Run simulations for each pathway size
  for (m in pathway_sizes_arg) {
    sigma_true <- make_cor_block(m, block_size = 20L, rho = rho)
    cache_key <- as.character(m)
    null_dist <- null_cache[[cache_key]]

    if (verbose) cat(sprintf("\nPathway size m=%d:\n", m))

    # 1. Baseline: no breakage (control)
    if (verbose) cat("  Running baseline (no breakage)...\n")
    p_baseline <- matrix(NA, nrow = n_sims, ncol = 6)
    p_baseline_mvn <- matrix(NA, nrow = n_sims, ncol = 6)
    colnames(p_baseline) <- colnames(p_baseline_mvn) <-
      c("acat", "fisher", "tfisher", "minp", "stouffer", "omnibus")

    for (s in seq_len(n_sims)) {
      z <- simulate_gene_z(m, sigma_true, n_causal = 0, effect_size = 0)
      gene_res <- create_gene_results(z)
      res_ana <- compute_components_analytic(gene_res)
      p_baseline[s, ] <- unlist(res_ana)
      res_cal <- calibrate_with_null(res_ana, null_dist)
      p_baseline_mvn[s, ] <- unlist(res_cal)
    }

    # Store baseline results
    for (method in colnames(p_baseline)) {
      results_list[[length(results_list) + 1]] <- data.frame(
        block = "B_ext",
        pathway_size = m,
        break_type = "none",
        severity = "baseline",
        n_broken = 0,
        broken_components = "none",
        method = method,
        scenario = "analytic",
        lambda = compute_lambda(p_baseline[, method]),
        lambda_se = compute_lambda_se(p_baseline[, method]),
        type1_05 = compute_type1(p_baseline[, method], 0.05),
        type1_05_se = compute_prop_se(p_baseline[, method], 0.05),
        type1_01 = compute_type1(p_baseline[, method], 0.01),
        type1_01_se = compute_prop_se(p_baseline[, method], 0.01),
        stringsAsFactors = FALSE
      )
      results_list[[length(results_list) + 1]] <- data.frame(
        block = "B_ext",
        pathway_size = m,
        break_type = "none",
        severity = "baseline",
        n_broken = 0,
        broken_components = "none",
        method = method,
        scenario = "mvn",
        lambda = compute_lambda(p_baseline_mvn[, method]),
        lambda_se = compute_lambda_se(p_baseline_mvn[, method]),
        type1_05 = compute_type1(p_baseline_mvn[, method], 0.05),
        type1_05_se = compute_prop_se(p_baseline_mvn[, method], 0.05),
        type1_01 = compute_type1(p_baseline_mvn[, method], 0.01),
        type1_01_se = compute_prop_se(p_baseline_mvn[, method], 0.01),
        stringsAsFactors = FALSE
      )
    }

    # 2. Anti-conservative breakage scenarios
    for (sev_name in names(anti_severities)) {
      sev_val <- anti_severities[[sev_name]]

      for (n_broken in n_broken_options) {
        # Sample which components to break
        set.seed(seed + m + n_broken * 100)
        components_to_break <- sample(all_components, n_broken)
        comp_str <- paste(sort(components_to_break), collapse = "+")

        if (verbose) cat(sprintf("  %s, n_broken=%d (%s)...\n",
                                 sev_name, n_broken, comp_str))

        p_broken <- matrix(NA, nrow = n_sims, ncol = 6)
        p_broken_mvn <- matrix(NA, nrow = n_sims, ncol = 6)
        colnames(p_broken) <- colnames(p_broken_mvn) <-
          c("acat", "fisher", "tfisher", "minp", "stouffer", "omnibus")

        for (s in seq_len(n_sims)) {
          z <- simulate_gene_z(m, sigma_true, n_causal = 0, effect_size = 0)
          gene_res <- create_gene_results(z)
          res_ana <- compute_components_analytic(gene_res)

          # Apply breakage
          res_broken <- apply_component_breakage(res_ana, components_to_break,
                                                  break_type = "anti",
                                                  severity = sev_val)
          p_broken[s, ] <- unlist(res_broken)

          # MVN calibration on original (unbroken) analytic
          res_cal <- calibrate_with_null(res_ana, null_dist)
          p_broken_mvn[s, ] <- unlist(res_cal)
        }

        for (method in colnames(p_broken)) {
          results_list[[length(results_list) + 1]] <- data.frame(
            block = "B_ext",
            pathway_size = m,
            break_type = "anti",
            severity = sev_name,
            n_broken = n_broken,
            broken_components = comp_str,
            method = method,
            scenario = "analytic_broken",
            lambda = compute_lambda(p_broken[, method]),
            lambda_se = compute_lambda_se(p_broken[, method]),
            type1_05 = compute_type1(p_broken[, method], 0.05),
            type1_05_se = compute_prop_se(p_broken[, method], 0.05),
            type1_01 = compute_type1(p_broken[, method], 0.01),
            type1_01_se = compute_prop_se(p_broken[, method], 0.01),
            stringsAsFactors = FALSE
          )
          results_list[[length(results_list) + 1]] <- data.frame(
            block = "B_ext",
            pathway_size = m,
            break_type = "anti",
            severity = sev_name,
            n_broken = n_broken,
            broken_components = comp_str,
            method = method,
            scenario = "mvn",
            lambda = compute_lambda(p_broken_mvn[, method]),
            lambda_se = compute_lambda_se(p_broken_mvn[, method]),
            type1_05 = compute_type1(p_broken_mvn[, method], 0.05),
            type1_05_se = compute_prop_se(p_broken_mvn[, method], 0.05),
            type1_01 = compute_type1(p_broken_mvn[, method], 0.01),
            type1_01_se = compute_prop_se(p_broken_mvn[, method], 0.01),
            stringsAsFactors = FALSE
          )
        }
      }
    }

    # 3. Conservative breakage scenarios
    for (sev_name in names(cons_severities)) {
      sev_val <- cons_severities[[sev_name]]

      for (n_broken in n_broken_options) {
        set.seed(seed + m + n_broken * 100 + 500)
        components_to_break <- sample(all_components, n_broken)
        comp_str <- paste(sort(components_to_break), collapse = "+")

        if (verbose) cat(sprintf("  %s, n_broken=%d (%s)...\n",
                                 sev_name, n_broken, comp_str))

        p_broken <- matrix(NA, nrow = n_sims, ncol = 6)
        p_broken_mvn <- matrix(NA, nrow = n_sims, ncol = 6)
        colnames(p_broken) <- colnames(p_broken_mvn) <-
          c("acat", "fisher", "tfisher", "minp", "stouffer", "omnibus")

        for (s in seq_len(n_sims)) {
          z <- simulate_gene_z(m, sigma_true, n_causal = 0, effect_size = 0)
          gene_res <- create_gene_results(z)
          res_ana <- compute_components_analytic(gene_res)

          res_broken <- apply_component_breakage(res_ana, components_to_break,
                                                  break_type = "cons",
                                                  severity = sev_val)
          p_broken[s, ] <- unlist(res_broken)

          res_cal <- calibrate_with_null(res_ana, null_dist)
          p_broken_mvn[s, ] <- unlist(res_cal)
        }

        for (method in colnames(p_broken)) {
          results_list[[length(results_list) + 1]] <- data.frame(
            block = "B_ext",
            pathway_size = m,
            break_type = "cons",
            severity = sev_name,
            n_broken = n_broken,
            broken_components = comp_str,
            method = method,
            scenario = "analytic_broken",
            lambda = compute_lambda(p_broken[, method]),
            lambda_se = compute_lambda_se(p_broken[, method]),
            type1_05 = compute_type1(p_broken[, method], 0.05),
            type1_05_se = compute_prop_se(p_broken[, method], 0.05),
            type1_01 = compute_type1(p_broken[, method], 0.01),
            type1_01_se = compute_prop_se(p_broken[, method], 0.01),
            stringsAsFactors = FALSE
          )
          results_list[[length(results_list) + 1]] <- data.frame(
            block = "B_ext",
            pathway_size = m,
            break_type = "cons",
            severity = sev_name,
            n_broken = n_broken,
            broken_components = comp_str,
            method = method,
            scenario = "mvn",
            lambda = compute_lambda(p_broken_mvn[, method]),
            lambda_se = compute_lambda_se(p_broken_mvn[, method]),
            type1_05 = compute_type1(p_broken_mvn[, method], 0.05),
            type1_05_se = compute_prop_se(p_broken_mvn[, method], 0.05),
            type1_01 = compute_type1(p_broken_mvn[, method], 0.01),
            type1_01_se = compute_prop_se(p_broken_mvn[, method], 0.01),
            stringsAsFactors = FALSE
          )
        }
      }
    }

    # 4. Mixed breakage: some anti-conservative, some conservative
    if (verbose) cat("  Running mixed breakage scenarios...\n")

    mixed_scenarios <- list(
      list(anti = c("acat", "minp"), cons = c("fisher", "stouffer")),
      list(anti = c("tfisher"), cons = c("acat", "fisher", "minp")),
      list(anti = c("acat", "fisher"), cons = c("stouffer"))
    )

    for (i in seq_along(mixed_scenarios)) {
      scenario <- mixed_scenarios[[i]]
      anti_comps <- scenario$anti
      cons_comps <- scenario$cons
      comp_str <- paste0("anti:", paste(anti_comps, collapse = "+"),
                         "_cons:", paste(cons_comps, collapse = "+"))

      if (verbose) cat(sprintf("  mixed scenario %d (%s)...\n", i, comp_str))

      p_broken <- matrix(NA, nrow = n_sims, ncol = 6)
      p_broken_mvn <- matrix(NA, nrow = n_sims, ncol = 6)
      colnames(p_broken) <- colnames(p_broken_mvn) <-
        c("acat", "fisher", "tfisher", "minp", "stouffer", "omnibus")

      for (s in seq_len(n_sims)) {
        z <- simulate_gene_z(m, sigma_true, n_causal = 0, effect_size = 0)
        gene_res <- create_gene_results(z)
        res_ana <- compute_components_analytic(gene_res)

        # Apply anti-conservative breakage
        res_broken <- apply_component_breakage(res_ana, anti_comps,
                                                break_type = "anti",
                                                severity = 0.5)
        # Then apply conservative breakage
        res_broken <- apply_component_breakage(res_broken, cons_comps,
                                                break_type = "cons",
                                                severity = 5)
        p_broken[s, ] <- unlist(res_broken)

        res_cal <- calibrate_with_null(res_ana, null_dist)
        p_broken_mvn[s, ] <- unlist(res_cal)
      }

      for (method in colnames(p_broken)) {
        results_list[[length(results_list) + 1]] <- data.frame(
          block = "B_ext",
          pathway_size = m,
          break_type = "mixed",
          severity = "moderate",
          n_broken = length(anti_comps) + length(cons_comps),
          broken_components = comp_str,
          method = method,
          scenario = "analytic_broken",
          lambda = compute_lambda(p_broken[, method]),
          lambda_se = compute_lambda_se(p_broken[, method]),
          type1_05 = compute_type1(p_broken[, method], 0.05),
          type1_05_se = compute_prop_se(p_broken[, method], 0.05),
          type1_01 = compute_type1(p_broken[, method], 0.01),
          type1_01_se = compute_prop_se(p_broken[, method], 0.01),
          stringsAsFactors = FALSE
        )
        results_list[[length(results_list) + 1]] <- data.frame(
          block = "B_ext",
          pathway_size = m,
          break_type = "mixed",
          severity = "moderate",
          n_broken = length(anti_comps) + length(cons_comps),
          broken_components = comp_str,
          method = method,
          scenario = "mvn",
          lambda = compute_lambda(p_broken_mvn[, method]),
          lambda_se = compute_lambda_se(p_broken_mvn[, method]),
          type1_05 = compute_type1(p_broken_mvn[, method], 0.05),
          type1_05_se = compute_prop_se(p_broken_mvn[, method], 0.05),
          type1_01 = compute_type1(p_broken_mvn[, method], 0.01),
          type1_01_se = compute_prop_se(p_broken_mvn[, method], 0.01),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (verbose) cat("\nBlock B Extended complete.\n\n")
  bind_rows(results_list)
}

# Plotting function for Block B Extended
plot_block_b_extended <- function(results, output_dir = "simulation_results") {

  df <- results %>%
    filter(block == "B_ext", scenario == "analytic_broken" | scenario == "analytic") %>%
    mutate(
      lambda_plot = pmin(pmax(lambda, 0), lambda_cap),  # Clamp lambda to [0, cap]
      method_label = factor(toupper(method),
                            levels = c("ACAT", "FISHER", "TFISHER", "MINP", "STOUFFER", "OMNIBUS")),
      severity_label = factor(severity,
                              levels = c("baseline", "anti_mild", "anti_moderate", "anti_severe",
                                         "cons_mild", "cons_moderate", "cons_severe", "moderate")),
      break_label = case_when(
        break_type == "none" ~ "Baseline",
        break_type == "anti" ~ "Anti-conservative (Inflated)",
        break_type == "cons" ~ "Conservative (Deflated)",
        break_type == "mixed" ~ "Mixed"
      )
    )

  # Custom color palette with OMNIBUS in dark purple (not yellow)
  method_colors <- c("ACAT" = "#E41A1C", "FISHER" = "#377EB8", "TFISHER" = "#4DAF4A",
                     "MINP" = "#984EA3", "STOUFFER" = "#FF7F00", "OMNIBUS" = "#000000")

  # ---- Summary data for line plots ----
  df_summary <- df %>%
    filter(break_type %in% c("anti", "cons")) %>%
    group_by(pathway_size, break_type, n_broken, method_label) %>%
    summarise(
      mean_type1 = mean(type1_05, na.rm = TRUE),
      se_type1 = sd(type1_05, na.rm = TRUE) / sqrt(n()),
      mean_lambda = mean(lambda_plot, na.rm = TRUE),
      se_lambda = sd(lambda_plot, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )

  # Separate data for inflated and deflated
  df_inflated <- df_summary %>% filter(break_type == "anti")
  df_deflated <- df_summary %>% filter(break_type == "cons")

  # Common theme for line plots
  line_plot_theme <- sim_theme +
    theme(plot.title = element_text(size = 22, face = "bold"),
          legend.position = "top",
          legend.box = "horizontal")

  # ---- Panel A: Lambda - Inflated only ----
  pA <- ggplot(df_inflated, aes(x = n_broken, y = mean_lambda, color = method_label)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
    facet_wrap(~ pathway_size, labeller = labeller(pathway_size = function(x) paste0("m=", x))) +
    scale_color_manual(values = method_colors, name = "Method") +
    labs(title = "A. Genomic Control (λ) - Inflated Components",
         x = "Number of Broken Components",
         y = expression(lambda)) +
    line_plot_theme

  # ---- Panel B: Type I Error - Inflated only ----
  pB <- ggplot(df_inflated, aes(x = n_broken, y = mean_type1, color = method_label)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.8) +
    facet_wrap(~ pathway_size, labeller = labeller(pathway_size = function(x) paste0("m=", x))) +
    scale_color_manual(values = method_colors, name = "Method") +
    labs(title = "B. Type I Error - Inflated Components",
         x = "Number of Broken Components",
         y = "Type I Error Rate") +
    line_plot_theme +
    theme(legend.position = "none")

  # ---- Panel C: Lambda - Deflated only ----
  pC <- ggplot(df_deflated, aes(x = n_broken, y = mean_lambda, color = method_label)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
    facet_wrap(~ pathway_size, labeller = labeller(pathway_size = function(x) paste0("m=", x))) +
    scale_color_manual(values = method_colors, name = "Method") +
    labs(title = "C. Genomic Control (λ) - Deflated Components",
         x = "Number of Broken Components",
         y = expression(lambda)) +
    line_plot_theme +
    theme(legend.position = "none")

  # ---- Panel D: Type I Error - Deflated only ----
  pD <- ggplot(df_deflated, aes(x = n_broken, y = mean_type1, color = method_label)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.8) +
    facet_wrap(~ pathway_size, labeller = labeller(pathway_size = function(x) paste0("m=", x))) +
    scale_color_manual(values = method_colors, name = "Method") +
    labs(title = "D. Type I Error - Deflated Components",
         x = "Number of Broken Components",
         y = "Type I Error Rate") +
    line_plot_theme +
    theme(legend.position = "none")

  # ---- Prepare comparison data with method names for Best ----
  df_individuals <- df %>%
    filter(break_type %in% c("anti", "cons"), n_broken >= 2, method != "omnibus") %>%
    mutate(deviation = abs(type1_05 - 0.05),
           deviation_lambda = abs(lambda_plot - 1))

  # Best individual per scenario - KEEP THE METHOD NAME
  df_best_type1 <- df_individuals %>%
    group_by(pathway_size, break_type, severity, n_broken, broken_components) %>%
    slice_min(deviation, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(comparison_group = paste0("Best(", method_label, ")"))

  df_best_lambda <- df_individuals %>%
    group_by(pathway_size, break_type, severity, n_broken, broken_components) %>%
    slice_min(deviation_lambda, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(comparison_group = paste0("Best(", method_label, ")"))

  # Average of all individuals
  df_avg <- df_individuals %>%
    group_by(pathway_size, break_type, severity, n_broken, broken_components) %>%
    summarise(
      type1_05 = mean(type1_05, na.rm = TRUE),
      lambda_plot = mean(lambda_plot, na.rm = TRUE),
      deviation = mean(deviation, na.rm = TRUE),
      deviation_lambda = mean(deviation_lambda, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(comparison_group = "Avg Individual")

  # Omnibus
  df_omnibus <- df %>%
    filter(break_type %in% c("anti", "cons"), n_broken >= 2, method == "omnibus") %>%
    mutate(
      deviation = abs(type1_05 - 0.05),
      deviation_lambda = abs(lambda_plot - 1),
      comparison_group = "OMNIBUS"
    )

  # For Type I comparison
  df_compare_type1 <- bind_rows(
    df_best_type1 %>% dplyr::select(pathway_size, break_type, severity, n_broken,
                              type1_05, deviation, comparison_group),
    df_avg %>% dplyr::select(pathway_size, break_type, severity, n_broken,
                       type1_05, deviation, comparison_group),
    df_omnibus %>% dplyr::select(pathway_size, break_type, severity, n_broken,
                           type1_05, deviation, comparison_group)
  )

  # For Lambda comparison
  df_compare_lambda <- bind_rows(
    df_best_lambda %>% dplyr::select(pathway_size, break_type, severity, n_broken,
                               lambda_plot, deviation_lambda, comparison_group),
    df_avg %>% dplyr::select(pathway_size, break_type, severity, n_broken,
                       lambda_plot, deviation_lambda, comparison_group),
    df_omnibus %>% dplyr::select(pathway_size, break_type, severity, n_broken,
                           lambda_plot, deviation_lambda, comparison_group)
  )

  # Create color mapping for comparison groups
  comparison_colors <- c("OMNIBUS" = "#2166AC", "Avg Individual" = "#E41A1C")
  # Add colors for Best(METHOD)
  for (m in c("ACAT", "FISHER", "TFISHER", "MINP", "STOUFFER")) {
    comparison_colors[paste0("Best(", m, ")")] <- "#4DAF4A"  # All Best in green
  }

  # ---- Panel E: Lambda comparison - Omnibus vs Individual ----
  pE <- ggplot(df_compare_lambda, aes(x = comparison_group, y = lambda_plot, fill = comparison_group)) +
    geom_boxplot(alpha = 0.8, outlier.alpha = 0.5) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
    facet_grid(break_type ~ severity,
               labeller = labeller(
                 break_type = c("anti" = "Inflated", "cons" = "Deflated")
               )) +
    scale_fill_manual(values = comparison_colors, name = "Group") +
    labs(title = "E. Genomic Control (λ): Omnibus vs Individual Methods",
         subtitle = "Best(METHOD) = whichever method was closest to λ=1 in each scenario",
         x = "",
         y = expression(lambda)) +
    sim_theme +
    theme(legend.position = "top",
          strip.text = element_text(size = 14),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1))

  # ---- Panel F: Type I comparison - Omnibus vs Individual ----
  pF <- ggplot(df_compare_type1, aes(x = comparison_group, y = type1_05, fill = comparison_group)) +
    geom_boxplot(alpha = 0.8, outlier.alpha = 0.5) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.8) +
    facet_grid(break_type ~ severity,
               labeller = labeller(
                 break_type = c("anti" = "Inflated", "cons" = "Deflated")
               )) +
    scale_fill_manual(values = comparison_colors, name = "Group") +
    labs(title = "F. Type I Error: Omnibus vs Individual Methods",
         subtitle = "Best(METHOD) = whichever method was closest to α=0.05 in each scenario",
         x = "",
         y = "Type I Error Rate") +
    sim_theme +
    theme(legend.position = "none",
          strip.text = element_text(size = 14),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1))

  # Save individual plots
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  ggsave(file.path(output_dir, "block_b_ext_A_lambda_inflated.png"),
         pA, width = 14, height = 8, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, "block_b_ext_B_type1_inflated.png"),
         pB, width = 14, height = 8, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, "block_b_ext_C_lambda_deflated.png"),
         pC, width = 14, height = 8, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, "block_b_ext_D_type1_deflated.png"),
         pD, width = 14, height = 8, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, "block_b_ext_E_lambda_comparison.png"),
         pE, width = 16, height = 10, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, "block_b_ext_F_type1_comparison.png"),
         pF, width = 16, height = 10, dpi = 300, bg = "white")

  # ---- Combined main figure ----
  p_combined_titled <- (pA | pB) / (pC | pD) / (pE | pF) +
    plot_annotation(
      title = "Stress Test: Omnibus Robustness to Component Miscalibration",
      theme = theme(
        plot.title = element_text(size = 26, face = "bold", hjust = 0.5)
      )
    ) +
    plot_layout(heights = c(1, 1, 1.2))

  ggsave(file.path(output_dir, "block_b_ext_combined_main.png"),
         p_combined_titled, width = 20, height = 24, dpi = 300, bg = "white")

  # ---- Create summary CSV tables ----
  # For CSVs, we need a simpler comparison_group that aggregates Best(METHOD) into "Best (post-hoc)"
  df_compare_csv <- bind_rows(
    df_best_type1 %>%
      dplyr::select(pathway_size, break_type, severity, n_broken, type1_05, deviation) %>%
      mutate(comparison_group = "Best (post-hoc)"),
    df_avg %>%
      dplyr::select(pathway_size, break_type, severity, n_broken, type1_05, deviation) %>%
      mutate(comparison_group = "Avg Individual"),
    df_omnibus %>%
      dplyr::select(pathway_size, break_type, severity, n_broken, type1_05, deviation) %>%
      mutate(comparison_group = "OMNIBUS")
  )

  # Summary: Omnibus vs Best Individual comparison
  summary_table <- df_compare_csv %>%
    group_by(break_type, severity, comparison_group) %>%
    summarise(
      n_scenarios = n(),
      mean_type1 = round(mean(type1_05, na.rm = TRUE), 4),
      sd_type1 = round(sd(type1_05, na.rm = TRUE), 4),
      mean_deviation = round(mean(deviation, na.rm = TRUE), 4),
      sd_deviation = round(sd(deviation, na.rm = TRUE), 4),
      min_type1 = round(min(type1_05, na.rm = TRUE), 4),
      max_type1 = round(max(type1_05, na.rm = TRUE), 4),
      .groups = "drop"
    ) %>%
    arrange(break_type, severity, comparison_group)

  write.csv(summary_table,
            file.path(output_dir, "block_b_ext_summary_omnibus_vs_best.csv"),
            row.names = FALSE)

  # Detailed summary by pathway size and n_broken
  detailed_table <- df_compare_csv %>%
    group_by(pathway_size, break_type, severity, n_broken, comparison_group) %>%
    summarise(
      mean_type1 = round(mean(type1_05, na.rm = TRUE), 4),
      mean_deviation = round(mean(deviation, na.rm = TRUE), 4),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = comparison_group,
      values_from = c(mean_type1, mean_deviation),
      names_sep = "_"
    ) %>%
    mutate(
      omnibus_vs_best = `mean_deviation_OMNIBUS` - `mean_deviation_Best (post-hoc)`,
      omnibus_vs_avg = `mean_deviation_OMNIBUS` - `mean_deviation_Avg Individual`,
      omnibus_better_than_avg = omnibus_vs_avg < 0
    ) %>%
    arrange(break_type, severity, pathway_size, n_broken)

  write.csv(detailed_table,
            file.path(output_dir, "block_b_ext_detailed_comparison.csv"),
            row.names = FALSE)

  # Overall summary
  overall_summary <- df_compare_csv %>%
    group_by(comparison_group) %>%
    summarise(
      n_total = n(),
      mean_type1 = round(mean(type1_05, na.rm = TRUE), 4),
      sd_type1 = round(sd(type1_05, na.rm = TRUE), 4),
      mean_deviation = round(mean(deviation, na.rm = TRUE), 4),
      sd_deviation = round(sd(deviation, na.rm = TRUE), 4),
      pct_within_2x_nominal = round(mean(type1_05 <= 0.10) * 100, 1),
      .groups = "drop"
    )

  write.csv(overall_summary,
            file.path(output_dir, "block_b_ext_overall_summary.csv"),
            row.names = FALSE)

  cat("\nCSV tables saved:\n")
  cat("  - block_b_ext_summary_omnibus_vs_best.csv\n")
  cat("  - block_b_ext_detailed_comparison.csv\n")
  cat("  - block_b_ext_overall_summary.csv\n")

  list(
    panel_A_lambda_inflated = pA,
    panel_B_type1_inflated = pB,
    panel_C_lambda_deflated = pC,
    panel_D_type1_deflated = pD,
    panel_E_lambda_comparison = pE,
    panel_F_type1_comparison = pF,
    combined = p_combined_titled,
    summary_table = summary_table,
    detailed_table = detailed_table,
    overall_summary = overall_summary
  )
}

# ==============================================================================
# BLOCK B SUMMARY: AVERAGED BROKEN-COMPONENT STRESS TEST
# ==============================================================================

mode_string <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  tab <- sort(table(x), decreasing = TRUE)
  names(tab)[1]
}

run_block_b_summary <- function(n_sims = 500L, b = 500L, seed = 42345L,
                                pathway_sizes_arg = c(5L, 25L, 50L),
                                rho_values = c(0, 0.3, 0.7),
                                n_broken_options = 1:4,
                                broken_power = 0.7,
                                verbose = TRUE) {

  set.seed(seed)
  all_components <- c("acat", "fisher", "tfisher", "minp", "stouffer")
  results_list <- list()

  if (verbose) cat("\n=== BLOCK B SUMMARY: Averaged Broken-Component Stress Test ===\n")

  null_cache <- list()
  observed_cache <- list()
  mvn_cache <- list()

  if (verbose) cat("Precomputing null distributions and baseline null replicates...\n")

  for (m in pathway_sizes_arg) {
    for (rho in rho_values) {
      sigma_true <- make_cor_block(m, block_size = 20L, rho = rho)
      cache_key <- paste(m, rho, sep = "_")

      if (verbose) cat(sprintf("  Precomputing: m=%d, rho=%.1f\n", m, rho))

      null_dist <- precompute_null_distribution(
        m, sigma_true, b = b, two_sided = FALSE,
        seed = seed + m * 1000 + round(rho * 100)
      )
      null_cache[[cache_key]] <- null_dist

      obs_analytic <- matrix(NA_real_, nrow = n_sims, ncol = 6)
      obs_mvn <- matrix(NA_real_, nrow = n_sims, ncol = 6)
      colnames(obs_analytic) <- colnames(obs_mvn) <-
        c("acat", "fisher", "tfisher", "minp", "stouffer", "omnibus")

      for (s in seq_len(n_sims)) {
        z <- simulate_gene_z(m, sigma_true, n_causal = 0, effect_size = 0)
        gene_res <- create_gene_results(z)
        res_ana <- compute_components_analytic(gene_res)
        obs_analytic[s, ] <- unlist(res_ana)
        obs_mvn[s, ] <- unlist(calibrate_with_null(res_ana, null_dist))
      }

      observed_cache[[cache_key]] <- obs_analytic
      mvn_cache[[cache_key]] <- obs_mvn
    }
  }

  for (m in pathway_sizes_arg) {
    for (rho in rho_values) {
      cache_key <- paste(m, rho, sep = "_")
      obs_analytic <- observed_cache[[cache_key]]
      obs_mvn <- mvn_cache[[cache_key]]

      if (verbose) cat(sprintf("  Summarising: m=%d, rho=%.1f\n", m, rho))

      omnibus_mvn_lambda <- compute_lambda(obs_mvn[, "omnibus"])
      omnibus_mvn_type1 <- compute_type1(obs_mvn[, "omnibus"], 0.05)

      for (n_broken in n_broken_options) {
        combo_matrix <- utils::combn(all_components, n_broken)
        n_combos <- ncol(combo_matrix)

        lambda_best <- numeric(n_combos)
        type1_best <- numeric(n_combos)
        lambda_best_method <- character(n_combos)
        type1_best_method <- character(n_combos)

        for (combo_idx in seq_len(n_combos)) {
          broken_components <- combo_matrix[, combo_idx]
          broken_mat <- obs_analytic

          for (comp in broken_components) {
            broken_mat[, comp] <- break_anticonservative(broken_mat[, comp], broken_power)
          }

          broken_mat[, "omnibus"] <- apply(
            broken_mat[, c("acat", "fisher", "tfisher", "minp", "stouffer"), drop = FALSE],
            1,
            compute_omnibus
          )

          comp_lambda <- vapply(
            all_components,
            function(comp) compute_lambda(broken_mat[, comp]),
            numeric(1)
          )
          comp_type1 <- vapply(
            all_components,
            function(comp) compute_type1(broken_mat[, comp], 0.05),
            numeric(1)
          )

          best_lambda_idx <- which.min(abs(comp_lambda - 1))
          best_type1_idx <- which.min(abs(comp_type1 - 0.05))

          lambda_best[combo_idx] <- comp_lambda[best_lambda_idx]
          type1_best[combo_idx] <- comp_type1[best_type1_idx]
          lambda_best_method[combo_idx] <- all_components[best_lambda_idx]
          type1_best_method[combo_idx] <- all_components[best_type1_idx]
        }

        results_list[[length(results_list) + 1]] <- data.frame(
          block = "B_summary",
          pathway_size = m,
          rho = rho,
          n_broken = n_broken,
          strategy = "best_component",
          lambda = mean(lambda_best, na.rm = TRUE),
          lambda_se = stats::sd(lambda_best, na.rm = TRUE) / sqrt(n_combos),
          type1_05 = mean(type1_best, na.rm = TRUE),
          type1_05_se = stats::sd(type1_best, na.rm = TRUE) / sqrt(n_combos),
          winner_label_lambda = mode_string(lambda_best_method),
          winner_label_type1 = mode_string(type1_best_method),
          n_combos = n_combos,
          stringsAsFactors = FALSE
        )

        results_list[[length(results_list) + 1]] <- data.frame(
          block = "B_summary",
          pathway_size = m,
          rho = rho,
          n_broken = n_broken,
          strategy = "omnibus_mvn",
          lambda = omnibus_mvn_lambda,
          lambda_se = 0,
          type1_05 = omnibus_mvn_type1,
          type1_05_se = 0,
          winner_label_lambda = "omnibus",
          winner_label_type1 = "omnibus",
          n_combos = n_combos,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  bind_rows(results_list)
}

plot_block_b_summary <- function(results) {
  df <- results %>%
    filter(block == "B_summary") %>%
    mutate(
      cor_structure = dplyr::case_when(
        rho == 0 ~ "LD_independent",
        rho == 0.3 ~ "LD_moderate",
        rho == 0.7 ~ "LD_strong",
        TRUE ~ paste0("rho=", rho)
      ),
      strategy_label = dplyr::case_when(
        strategy == "best_component" ~ "Best Broken Component",
        strategy == "omnibus_mvn" ~ "Omnibus MVN",
        TRUE ~ strategy
      ),
      strategy_label = factor(strategy_label,
                              levels = c("Best Broken Component", "Omnibus MVN"))
    )

  strategy_colors <- c(
    "Best Broken Component" = "#E41A1C",
    "Omnibus MVN" = "#377EB8"
  )

  p_lambda <- ggplot(df, aes(x = n_broken, y = lambda, color = strategy_label)) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2.4) +
    geom_errorbar(aes(ymin = lambda - lambda_se, ymax = lambda + lambda_se),
                  width = 0.12, linewidth = 0.5) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.7) +
    facet_grid(pathway_size ~ cor_structure,
               labeller = labeller(pathway_size = function(x) paste0("m=", x))) +
    scale_color_manual(values = strategy_colors, labels = c(
      "Best Broken Component",
      "Omnibus MVN"
    )) +
    scale_x_continuous(breaks = 1:4) +
    labs(title = "Block B Summary: Broken Components vs Omnibus MVN",
         subtitle = "Best component averaged over all combinations of broken tests of size k",
         x = "Number of Broken Components",
         y = expression(lambda),
         color = "Strategy") +
    sim_theme

  p_type1 <- ggplot(df, aes(x = n_broken, y = type1_05, color = strategy_label)) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2.4) +
    geom_errorbar(aes(ymin = pmax(type1_05 - type1_05_se, 0),
                      ymax = type1_05 + type1_05_se),
                  width = 0.12, linewidth = 0.5) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.7) +
    facet_grid(pathway_size ~ cor_structure,
               labeller = labeller(pathway_size = function(x) paste0("m=", x))) +
    scale_color_manual(values = strategy_colors, labels = c(
      "Best Broken Component",
      "Omnibus MVN"
    )) +
    scale_x_continuous(breaks = 1:4) +
    labs(title = "Block B Summary: Type I Error Under Component Breakage",
         subtitle = "Best component averaged over all combinations of broken tests of size k",
         x = "Number of Broken Components",
         y = "Type I Error Rate",
         color = "Strategy") +
    sim_theme

  list(lambda = p_lambda, type1 = p_type1)
}

# ==============================================================================
# BLOCK B PROFILE: COMPONENT-WISE STRESS TEST WITH OMNIBUS RE-CALIBRATION
# ==============================================================================

apply_component_breakage_matrix <- function(p_mat, components_to_break,
                                            broken_power = 0.7) {
  broken_mat <- p_mat
  for (comp in intersect(components_to_break, colnames(broken_mat))) {
    broken_mat[, comp] <- break_anticonservative(broken_mat[, comp], broken_power)
  }
  broken_mat
}

run_block_b_profile <- function(n_sims = 500L, b = 500L, seed = 52345L,
                                pathway_sizes_arg = c(5L, 25L, 50L),
                                rho_values = c(0, 0.3, 0.7),
                                n_broken_options = 1:4,
                                broken_power = 0.7,
                                verbose = TRUE) {

  set.seed(seed)
  all_components <- c("acat", "fisher", "tfisher", "minp", "stouffer")
  method_order <- c(
    all_components,
    "omnibus_analytic",
    "omnibus_mvn_combined",
    "omnibus_mvn_alone"
  )
  results_list <- list()

  if (verbose) cat("\n=== BLOCK B PROFILE: Component-Wise Broken-Component Stress Test ===\n")
  if (verbose) cat("Precomputing baseline null replicates and null distributions...\n")

  observed_cache <- list()
  null_cache <- list()

  for (m in pathway_sizes_arg) {
    for (rho in rho_values) {
      sigma_true <- make_cor_block(m, block_size = 20L, rho = rho)
      cache_key <- paste(m, rho, sep = "_")

      if (verbose) cat(sprintf("  Precomputing: m=%d, rho=%.1f\n", m, rho))

      null_cache[[cache_key]] <- precompute_null_distribution(
        m, sigma_true, b = b, two_sided = FALSE,
        seed = seed + m * 1000 + round(rho * 100)
      )

      obs_analytic <- matrix(NA_real_, nrow = n_sims, ncol = 6)
      colnames(obs_analytic) <- c("acat", "fisher", "tfisher", "minp", "stouffer", "omnibus")

      for (s in seq_len(n_sims)) {
        z <- simulate_gene_z(m, sigma_true, n_causal = 0, effect_size = 0)
        gene_res <- create_gene_results(z)
        obs_analytic[s, ] <- unlist(compute_components_analytic(gene_res))
      }

      observed_cache[[cache_key]] <- obs_analytic
    }
  }

  for (m in pathway_sizes_arg) {
    for (rho in rho_values) {
      cache_key <- paste(m, rho, sep = "_")
      obs_analytic <- observed_cache[[cache_key]]
      null_dist <- null_cache[[cache_key]]

      if (verbose) cat(sprintf("  Profiling: m=%d, rho=%.1f\n", m, rho))

      for (n_broken in n_broken_options) {
        combo_matrix <- utils::combn(all_components, n_broken)
        n_combos <- ncol(combo_matrix)

        metric_lambda <- matrix(NA_real_, nrow = n_combos, ncol = length(method_order))
        metric_type1 <- matrix(NA_real_, nrow = n_combos, ncol = length(method_order))
        colnames(metric_lambda) <- colnames(metric_type1) <- method_order

        for (combo_idx in seq_len(n_combos)) {
          broken_components <- combo_matrix[, combo_idx]

          broken_obs <- apply_component_breakage_matrix(
            obs_analytic[, all_components, drop = FALSE],
            broken_components,
            broken_power = broken_power
          )
          broken_obs_omni <- apply(broken_obs, 1, compute_omnibus)

          broken_null <- apply_component_breakage_matrix(
            null_dist[, all_components, drop = FALSE],
            broken_components,
            broken_power = broken_power
          )
          broken_null_omni <- apply(broken_null, 1, compute_omnibus)
          broken_obs_cal <- calibrate_component_matrix(broken_obs, broken_null)
          omnibus_mvn_combined <- apply(broken_obs_cal, 1, compute_omnibus)
          omnibus_mvn_alone <- empirical_p_from_null(broken_obs_omni, broken_null_omni)

          for (comp in all_components) {
            metric_lambda[combo_idx, comp] <- compute_lambda(broken_obs[, comp])
            metric_type1[combo_idx, comp] <- compute_type1(broken_obs[, comp], 0.05)
          }

          metric_lambda[combo_idx, "omnibus_analytic"] <- compute_lambda(broken_obs_omni)
          metric_type1[combo_idx, "omnibus_analytic"] <- compute_type1(broken_obs_omni, 0.05)

          metric_lambda[combo_idx, "omnibus_mvn_combined"] <- compute_lambda(omnibus_mvn_combined)
          metric_type1[combo_idx, "omnibus_mvn_combined"] <- compute_type1(omnibus_mvn_combined, 0.05)

          metric_lambda[combo_idx, "omnibus_mvn_alone"] <- compute_lambda(omnibus_mvn_alone)
          metric_type1[combo_idx, "omnibus_mvn_alone"] <- compute_type1(omnibus_mvn_alone, 0.05)
        }

        for (method in method_order) {
          results_list[[length(results_list) + 1]] <- data.frame(
            block = "B_profile",
            pathway_size = m,
            rho = rho,
            n_broken = n_broken,
            method = method,
            lambda = mean(metric_lambda[, method], na.rm = TRUE),
            lambda_se = stats::sd(metric_lambda[, method], na.rm = TRUE) / sqrt(n_combos),
            type1_05 = mean(metric_type1[, method], na.rm = TRUE),
            type1_05_se = stats::sd(metric_type1[, method], na.rm = TRUE) / sqrt(n_combos),
            n_combos = n_combos,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  bind_rows(results_list)
}

plot_block_b_profile <- function(results, include_analytic_omnibus = FALSE) {
  method_levels <- c("acat", "fisher", "tfisher", "minp", "stouffer",
                     "omnibus_mvn_combined", "omnibus_mvn_alone")
  if (include_analytic_omnibus) {
    method_levels <- c("acat", "fisher", "tfisher", "minp", "stouffer",
                       "omnibus_analytic", "omnibus_mvn_combined", "omnibus_mvn_alone")
  }

  method_labels <- c(
    "acat" = "ACAT",
    "fisher" = "Fisher",
    "tfisher" = "TFisher",
    "minp" = "minP",
    "stouffer" = "Stouffer",
    "omnibus_analytic" = "Omnibus (analytic)",
    "omnibus_mvn_combined" = "Omnibus Combined",
    "omnibus_mvn_alone" = "Omnibus Alone"
  )

  method_colors <- c(
    "ACAT" = "#E41A1C",
    "Fisher" = "#377EB8",
    "TFisher" = "#4DAF4A",
    "minP" = "#984EA3",
    "Stouffer" = "#FF7F00",
    "Omnibus (analytic)" = "#666666",
    "Omnibus Combined" = "#000000",
    "Omnibus Alone" = "#8C564B"
  )

  df <- results %>%
    filter(block == "B_profile",
           method %in% method_levels) %>%
    mutate(
      cor_structure = dplyr::case_when(
        rho == 0 ~ "LD_independent",
        rho == 0.3 ~ "LD_moderate",
        rho == 0.7 ~ "LD_strong",
        TRUE ~ paste0("rho=", rho)
      ),
      method = factor(method, levels = method_levels),
      method_label = factor(method_labels[as.character(method)],
                            levels = unname(method_labels[method_levels]))
    )

  p_lambda <- ggplot(df, aes(x = n_broken, y = lambda,
                             color = method_label, group = method_label)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.8) +
    geom_errorbar(aes(ymin = pmax(lambda - lambda_se, 0),
                      ymax = lambda + lambda_se),
                  width = 0.08, linewidth = 0.45, alpha = 0.9) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.7) +
    facet_grid(pathway_size ~ cor_structure,
               labeller = labeller(pathway_size = function(x) paste0("m=", x))) +
    scale_color_manual(values = method_colors) +
    scale_x_continuous(breaks = 1:4) +
    labs(title = "Block B: Component Breakage Stress Test",
         subtitle = "Average over all combinations of k broken components; compare component-wise MVN recombination vs direct omnibus calibration",
         x = "Number of Broken Components",
         y = expression(lambda),
         color = "Method") +
    sim_theme

  p_type1 <- ggplot(df, aes(x = n_broken, y = type1_05,
                            color = method_label, group = method_label)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.8) +
    geom_errorbar(aes(ymin = pmax(type1_05 - type1_05_se, 0),
                      ymax = type1_05 + type1_05_se),
                  width = 0.08, linewidth = 0.45, alpha = 0.9) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.7) +
    facet_grid(pathway_size ~ cor_structure,
               labeller = labeller(pathway_size = function(x) paste0("m=", x))) +
    scale_color_manual(values = method_colors) +
    scale_x_continuous(breaks = 1:4) +
    labs(title = "Block B: Type I Error Under Component Breakage",
         subtitle = "Average over all combinations of k broken components; omnibus recalibrated under matching broken MVN null",
         x = "Number of Broken Components",
         y = "Type I Error Rate",
         color = "Method") +
    sim_theme

  list(lambda = p_lambda, type1 = p_type1)
}

# ==============================================================================
# BLOCK B': INPUT MISCALIBRATION STRESS TEST
# ==============================================================================
# This is a more realistic stress test: instead of breaking individual component
# tests, we inflate/deflate the GENE p-values (the input to all tests) and see
# how each component test responds. Different tests have different sensitivities.

run_block_b_input <- function(n_sims = 300L, b = 500L, seed = 54321L,
                               pathway_sizes_arg = c(5L, 25L, 50L),
                               rho = 0.3,
                               verbose = TRUE) {

  set.seed(seed)
  results_list <- list()

  if (verbose) {
    cat("\n")
    cat(strrep("=", 70), "\n")
    cat("BLOCK B': INPUT MISCALIBRATION STRESS TEST\n")
    cat(strrep("=", 70), "\n")
    cat("Testing: How do component tests respond when GENE p-values are miscalibrated?\n")
    cat("  - Anti-conservative: gene p-values are too small (inflated significance)\
")
    cat("  - Conservative: gene p-values are too large (deflated significance)\n")
    cat("  - All tests receive the SAME distorted input\n\n")
  }

  # Severity levels for gene p-value distortion
  # p^severity transformation:
  #   - severity > 1: p^2 makes p SMALLER (e.g., 0.1^2 = 0.01) → Anti-conservative (inflated)
  #   - severity < 1: p^0.5 makes p LARGER (e.g., 0.01^0.5 = 0.1) → Conservative (deflated)
  severities <- list(
    baseline = 1.0,      # No distortion
    cons_mild = 0.8,     # Mild: p larger → conservative
    cons_moderate = 0.5, # Moderate: p larger → conservative
    cons_severe = 0.3,   # Severe: p larger → conservative
    anti_mild = 1.25,    # Mild: p smaller → anti-conservative
    anti_moderate = 2.0, # Moderate: p smaller → anti-conservative
    anti_severe = 3.0    # Severe: p smaller → anti-conservative
  )

  # Precompute null distributions for MVN calibration
  if (verbose) cat("Precomputing null distributions for calibration...\n")
  null_cache <- list()
  sigma_cache <- list()

  for (m in pathway_sizes_arg) {
    sigma <- make_cor_block(m, block_size = max(5L, m %/% 5), rho = rho)
    sigma_cache[[as.character(m)]] <- sigma
    null_cache[[as.character(m)]] <- precompute_null_distribution(
      m, sigma, b = b, two_sided = FALSE, seed = seed + m * 100
    )
  }

  # Run simulations
  for (m in pathway_sizes_arg) {
    if (verbose) cat(sprintf("\nPathway size m=%d\n", m))

    sigma <- sigma_cache[[as.character(m)]]
    null_dist <- null_cache[[as.character(m)]]

    for (sev_name in names(severities)) {
      sev_val <- severities[[sev_name]]

      # Determine break type
      if (sev_name == "baseline") {
        break_type <- "none"
      } else if (startsWith(sev_name, "anti")) {
        break_type <- "anti"
      } else {
        break_type <- "cons"
      }

      if (verbose) cat(sprintf("  Severity: %s (exponent=%.2f)...\n", sev_name, sev_val))

      # Storage for p-values
      p_analytic <- matrix(NA, nrow = n_sims, ncol = 6)
      p_mvn <- matrix(NA, nrow = n_sims, ncol = 6)
      colnames(p_analytic) <- colnames(p_mvn) <-
        c("acat", "fisher", "tfisher", "minp", "stouffer", "omnibus")

      for (s in seq_len(n_sims)) {
        # Generate null z-scores from MVN
        z <- simulate_gene_z(m, sigma, n_causal = 0, effect_size = 0)

        # Convert to p-values
        gene_p <- 2 * pnorm(-abs(z))

        # Apply distortion to gene p-values: p^severity
        # For anti-conservative (severity < 1): p^0.5 makes p smaller
        # For conservative (severity > 1): p^2 makes p larger
        gene_p_distorted <- pmin(pmax(gene_p^sev_val, 1e-300), 1 - 1e-10)

        # Create gene results with distorted p-values
        gene_res <- data.frame(
          GENE = paste0("G", seq_len(m)),
          ZSTAT = qnorm(gene_p_distorted / 2, lower.tail = FALSE),  # Back to z for some methods
          P = gene_p_distorted,
          stringsAsFactors = FALSE
        )

        # Compute all component tests on distorted input
        res_ana <- compute_components_analytic(gene_res)
        p_analytic[s, ] <- unlist(res_ana)

        # Also compute MVN-calibrated version (using original null distribution)
        # Note: This tests if MVN calibration can rescue distorted input
        res_mvn <- calibrate_with_null(res_ana, null_dist)
        p_mvn[s, ] <- unlist(res_mvn)
      }

      # Store results for each method
      for (method in colnames(p_analytic)) {
        results_list[[length(results_list) + 1]] <- data.frame(
          block = "B_input",
          pathway_size = m,
          break_type = break_type,
          severity = sev_name,
          severity_value = sev_val,
          method = method,
          scenario = "analytic",
          lambda = compute_lambda(p_analytic[, method]),
          lambda_se = compute_lambda_se(p_analytic[, method]),
          type1_05 = compute_type1(p_analytic[, method], 0.05),
          type1_05_se = compute_prop_se(p_analytic[, method], 0.05),
          type1_01 = compute_type1(p_analytic[, method], 0.01),
          type1_01_se = compute_prop_se(p_analytic[, method], 0.01),
          stringsAsFactors = FALSE
        )

        results_list[[length(results_list) + 1]] <- data.frame(
          block = "B_input",
          pathway_size = m,
          break_type = break_type,
          severity = sev_name,
          severity_value = sev_val,
          method = method,
          scenario = "mvn_calibrated",
          lambda = compute_lambda(p_mvn[, method]),
          lambda_se = compute_lambda_se(p_mvn[, method]),
          type1_05 = compute_type1(p_mvn[, method], 0.05),
          type1_05_se = compute_prop_se(p_mvn[, method], 0.05),
          type1_01 = compute_type1(p_mvn[, method], 0.01),
          type1_01_se = compute_prop_se(p_mvn[, method], 0.01),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (verbose) {
    cat("\n")
    cat(strrep("=", 70), "\n")
    cat("Block B' (Input Miscalibration) complete.\n")
    cat(strrep("=", 70), "\n\n")
  }

  bind_rows(results_list)
}

# Plotting function for Block B' Input Miscalibration
plot_block_b_input <- function(results, output_dir = "simulation_results") {

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Filter to analytic scenario (the main results)
  df <- results %>%
    filter(block == "B_input", scenario == "analytic") %>%
    mutate(
      method_label = factor(toupper(method),
                            levels = c("ACAT", "FISHER", "TFISHER", "MINP", "STOUFFER", "OMNIBUS")),
      severity_label = factor(severity,
                              levels = c("baseline", "anti_mild", "anti_moderate", "anti_severe",
                                         "cons_mild", "cons_moderate", "cons_severe")),
      break_label = case_when(
        break_type == "none" ~ "Baseline",
        break_type == "anti" ~ "Anti-conservative\n(Inflated)",
        break_type == "cons" ~ "Conservative\n(Deflated)"
      )
    )

  # Color palette - OMNIBUS in dark blue for visibility
  method_colors <- c("ACAT" = "#E41A1C", "FISHER" = "#377EB8", "TFISHER" = "#4DAF4A",
                     "MINP" = "#984EA3", "STOUFFER" = "#FF7F00", "OMNIBUS" = "#000000")

  # Common theme
  plot_theme <- sim_theme +
    theme(
      legend.position = "top",
      legend.box = "horizontal",
      strip.text = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )

  # ---- Panel A: Type I Error by Severity (faceted by pathway size) ----
  pA <- ggplot(df %>% filter(break_type != "none"),
               aes(x = severity_label, y = type1_05, color = method_label, group = method_label)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.8) +
    facet_grid(break_type ~ pathway_size,
               scales = "free_y",
               labeller = labeller(
                 pathway_size = function(x) paste0("m = ", x),
                 break_type = c("anti" = "Inflated Input", "cons" = "Deflated Input")
               )) +
    scale_color_manual(values = method_colors, name = "Method") +
    scale_y_continuous(limits = c(0, NA)) +
    labs(title = "A. Type I Error Rate by Input Distortion Severity",
         subtitle = "Gene p-values distorted as p^severity before running all tests",
         x = "Severity Level",
         y = "Type I Error Rate (α = 0.05)") +
    plot_theme

  # ---- Panel B: Lambda by Severity ----
  pB <- ggplot(df %>% filter(break_type != "none"),
               aes(x = severity_label, y = lambda, color = method_label, group = method_label)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
    facet_grid(break_type ~ pathway_size,
               scales = "free_y",
               labeller = labeller(
                 pathway_size = function(x) paste0("m = ", x),
                 break_type = c("anti" = "Inflated Input", "cons" = "Deflated Input")
               )) +
    scale_color_manual(values = method_colors, name = "Method") +
    labs(title = "B. Genomic Control (λ) by Input Distortion Severity",
         subtitle = "λ = 1 indicates proper calibration; λ > 1 indicates inflation",
         x = "Severity Level",
         y = expression(lambda)) +
    plot_theme +
    theme(legend.position = "none")

  # ---- Panel C: Boxplot comparison across all scenarios ----
  df_compare <- df %>%
    filter(break_type != "none") %>%
    mutate(is_omnibus = ifelse(method == "omnibus", "OMNIBUS", "Individual"))

  pC <- ggplot(df_compare,
               aes(x = method_label, y = type1_05, fill = method_label)) +
    geom_boxplot(alpha = 0.8, outlier.alpha = 0.5) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.8) +
    facet_wrap(~ break_type,
               labeller = labeller(break_type = c("anti" = "Inflated Input",
                                                   "cons" = "Deflated Input"))) +
    scale_fill_manual(values = method_colors) +
    labs(title = "C. Type I Error Distribution by Method",
         subtitle = "Across all pathway sizes and severity levels",
         x = "Method",
         y = "Type I Error Rate") +
    plot_theme +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 0, hjust = 0.5))

  # ---- Panel D: Lambda boxplot ----
  pD <- ggplot(df_compare,
               aes(x = method_label, y = lambda, fill = method_label)) +
    geom_boxplot(alpha = 0.8, outlier.alpha = 0.5) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
    facet_wrap(~ break_type,
               labeller = labeller(break_type = c("anti" = "Inflated Input",
                                                   "cons" = "Deflated Input"))) +
    scale_fill_manual(values = method_colors) +
    labs(title = "D. Genomic Control (λ) Distribution by Method",
         subtitle = "Across all pathway sizes and severity levels",
         x = "Method",
         y = expression(lambda)) +
    plot_theme +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 0, hjust = 0.5))

  # Save individual plots
  ggsave(file.path(output_dir, "block_b_input_A_type1_by_severity.png"),
         pA, width = 14, height = 10, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, "block_b_input_B_lambda_by_severity.png"),
         pB, width = 14, height = 10, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, "block_b_input_C_type1_boxplot.png"),
         pC, width = 12, height = 8, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, "block_b_input_D_lambda_boxplot.png"),
         pD, width = 12, height = 8, dpi = 300, bg = "white")

  # ---- Combined figure ----
  p_combined <- (pA / pB) | (pC / pD)
  p_combined <- p_combined +
    plot_annotation(
      title = "Stress Test: Component Response to Gene P-value Miscalibration",
      subtitle = "All tests receive the same distorted gene p-values (p^severity)",
      theme = theme(
        plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5)
      )
    )

  ggsave(file.path(output_dir, "block_b_input_combined.png"),
         p_combined, width = 22, height = 16, dpi = 300, bg = "white")

  # ---- Summary tables ----
  # Summary by method and break type
  summary_by_method <- df %>%
    filter(break_type != "none") %>%
    group_by(method_label, break_type) %>%
    summarise(
      n_scenarios = n(),
      mean_type1 = round(mean(type1_05, na.rm = TRUE), 4),
      sd_type1 = round(sd(type1_05, na.rm = TRUE), 4),
      mean_lambda = round(mean(lambda, na.rm = TRUE), 4),
      sd_lambda = round(sd(lambda, na.rm = TRUE), 4),
      max_type1 = round(max(type1_05, na.rm = TRUE), 4),
      .groups = "drop"
    ) %>%
    arrange(break_type, desc(mean_type1))

  write.csv(summary_by_method,
            file.path(output_dir, "block_b_input_summary_by_method.csv"),
            row.names = FALSE)

  # Detailed results
  detailed_table <- df %>%
    filter(break_type != "none") %>%
    dplyr::select(pathway_size, break_type, severity, severity_value,
                  method_label, type1_05, lambda) %>%
    arrange(break_type, severity, pathway_size, method_label)

  write.csv(detailed_table,
            file.path(output_dir, "block_b_input_detailed.csv"),
            row.names = FALSE)

  cat("\nPlots saved:\n")
  cat("  - block_b_input_A_type1_by_severity.png\n")
  cat("  - block_b_input_B_lambda_by_severity.png\n")
  cat("  - block_b_input_C_type1_boxplot.png\n")
  cat("  - block_b_input_D_lambda_boxplot.png\n")
  cat("  - block_b_input_combined.png\n")
  cat("\nCSV tables saved:\n")
  cat("  - block_b_input_summary_by_method.csv\n")
  cat("  - block_b_input_detailed.csv\n")

  list(
    panel_A = pA,
    panel_B = pB,
    panel_C = pC,
    panel_D = pD,
    combined = p_combined,
    summary_by_method = summary_by_method,
    detailed_table = detailed_table
  )
}

# ==============================================================================
# BLOCK C: MISSING CORRELATION SCENARIOS
# ==============================================================================

damage_correlation <- function(sigma, missing_frac = 0.3) {
  m <- nrow(sigma)
  sigma_damaged <- sigma

  off_diag <- which(upper.tri(sigma), arr.ind = TRUE)
  n_damage <- round(nrow(off_diag) * missing_frac)

  if (n_damage > 0) {
    damage_idx <- sample(seq_len(nrow(off_diag)), n_damage)

    for (k in damage_idx) {
      i <- off_diag[k, 1]
      j <- off_diag[k, 2]
      sigma_damaged[i, j] <- 0
      sigma_damaged[j, i] <- 0
    }
  }

  ensure_pd(sigma_damaged)
}

run_block_c <- function(n_sims = 500L, b = 500L, seed = 32345L,
                        pathway_sizes_arg = pathway_sizes,
                        missing_fracs = c(0, 0.1, 0.3, 0.5, 0.7),
                        verbose = TRUE) {

  set.seed(seed)

  results_list <- list()

  if (verbose) cat("\n=== BLOCK C: Missing Correlation Scenarios ===\n")
  if (verbose) cat("Testing: Impact of incomplete correlation matrices\n\n")

  if (verbose) cat("Precomputing null distributions...\n")
  null_cache_true <- list()
  null_cache_damaged <- list()
  sigma_damaged_cache <- list()

  for (m in pathway_sizes_arg) {
    sigma_true <- make_cor_block(m, block_size = 20L, rho = 0.5)
    null_cache_true[[as.character(m)]] <-
      precompute_null_distribution(m, sigma_true, b = b, two_sided = FALSE,
                                   seed = seed + m * 1000)

    for (miss_frac in missing_fracs) {
      sigma_damaged <- damage_correlation(sigma_true, miss_frac)
      cache_key <- paste(m, miss_frac, sep = "_")
      sigma_damaged_cache[[cache_key]] <- sigma_damaged
      null_cache_damaged[[cache_key]] <-
        precompute_null_distribution(m, sigma_damaged, b = b, two_sided = FALSE,
                                     seed = seed + m * 1000 + round(miss_frac * 100))
    }
  }

  for (m in pathway_sizes_arg) {
    sigma_true <- make_cor_block(m, block_size = 20L, rho = 0.5)
    null_true <- null_cache_true[[as.character(m)]]

    for (miss_frac in missing_fracs) {
      if (verbose) cat(sprintf("  m=%d, missing=%.0f%%\n", m, miss_frac * 100))

      cache_key <- paste(m, miss_frac, sep = "_")
      sigma_damaged <- sigma_damaged_cache[[cache_key]]
      null_damaged <- null_cache_damaged[[cache_key]]

      component_cols <- c("acat", "fisher", "tfisher", "minp", "stouffer")
      combined_cols <- c(component_cols, "omnibus_combined")
      p_analytic <- matrix(NA, nrow = n_sims, ncol = length(combined_cols))
      p_true_mvn <- matrix(NA, nrow = n_sims, ncol = length(combined_cols))
      p_imputed_mvn <- matrix(NA, nrow = n_sims, ncol = length(combined_cols))
      p_true_mvn_alone <- numeric(n_sims)
      p_imputed_mvn_alone <- numeric(n_sims)
      colnames(p_analytic) <- colnames(p_true_mvn) <- colnames(p_imputed_mvn) <- combined_cols

      for (s in seq_len(n_sims)) {
        z <- simulate_gene_z(m, sigma_true, n_causal = 0, effect_size = 0)
        gene_res <- create_gene_results(z)

        res_ana <- compute_components_analytic(gene_res)
        p_analytic[s, component_cols] <- unlist(res_ana[component_cols])
        p_analytic[s, "omnibus_combined"] <- res_ana$omnibus

        res_true <- calibrate_with_null(res_ana, null_true)
        p_true_mvn[s, component_cols] <- unlist(res_true[component_cols])
        p_true_mvn[s, "omnibus_combined"] <- compute_omnibus(unlist(res_true[component_cols]))
        p_true_mvn_alone[s] <- res_true$omnibus

        res_imputed <- calibrate_with_null(res_ana, null_damaged)
        p_imputed_mvn[s, component_cols] <- unlist(res_imputed[component_cols])
        p_imputed_mvn[s, "omnibus_combined"] <- compute_omnibus(unlist(res_imputed[component_cols]))
        p_imputed_mvn_alone[s] <- res_imputed$omnibus
      }

      for (method in colnames(p_analytic)) {
        results_list[[length(results_list) + 1]] <- data.frame(
          block = "C",
          pathway_size = m,
          missing_frac = miss_frac,
          method = method,
          strategy = "analytic_fallback",
          lambda = compute_lambda(p_analytic[, method]),
          lambda_se = compute_lambda_se(p_analytic[, method]),
          type1_05 = compute_type1(p_analytic[, method], 0.05),
          type1_05_se = compute_prop_se(p_analytic[, method], 0.05),
          type1_01 = compute_type1(p_analytic[, method], 0.01),
          type1_01_se = compute_prop_se(p_analytic[, method], 0.01),
          stringsAsFactors = FALSE
        )

        results_list[[length(results_list) + 1]] <- data.frame(
          block = "C",
          pathway_size = m,
          missing_frac = miss_frac,
          method = method,
          strategy = "mvn_true_cor",
          lambda = compute_lambda(p_true_mvn[, method]),
          lambda_se = compute_lambda_se(p_true_mvn[, method]),
          type1_05 = compute_type1(p_true_mvn[, method], 0.05),
          type1_05_se = compute_prop_se(p_true_mvn[, method], 0.05),
          type1_01 = compute_type1(p_true_mvn[, method], 0.01),
          type1_01_se = compute_prop_se(p_true_mvn[, method], 0.01),
          stringsAsFactors = FALSE
        )

        results_list[[length(results_list) + 1]] <- data.frame(
          block = "C",
          pathway_size = m,
          missing_frac = miss_frac,
          method = method,
          strategy = "mvn_imputed_cor",
          lambda = compute_lambda(p_imputed_mvn[, method]),
          lambda_se = compute_lambda_se(p_imputed_mvn[, method]),
          type1_05 = compute_type1(p_imputed_mvn[, method], 0.05),
          type1_05_se = compute_prop_se(p_imputed_mvn[, method], 0.05),
          type1_01 = compute_type1(p_imputed_mvn[, method], 0.01),
          type1_01_se = compute_prop_se(p_imputed_mvn[, method], 0.01),
          stringsAsFactors = FALSE
        )
      }

      results_list[[length(results_list) + 1]] <- data.frame(
        block = "C",
        pathway_size = m,
        missing_frac = miss_frac,
        method = "omnibus_alone",
        strategy = "mvn_true_cor",
        lambda = compute_lambda(p_true_mvn_alone),
        lambda_se = compute_lambda_se(p_true_mvn_alone),
        type1_05 = compute_type1(p_true_mvn_alone, 0.05),
        type1_05_se = compute_prop_se(p_true_mvn_alone, 0.05),
        type1_01 = compute_type1(p_true_mvn_alone, 0.01),
        type1_01_se = compute_prop_se(p_true_mvn_alone, 0.01),
        stringsAsFactors = FALSE
      )

      results_list[[length(results_list) + 1]] <- data.frame(
        block = "C",
        pathway_size = m,
        missing_frac = miss_frac,
        method = "omnibus_alone",
        strategy = "mvn_imputed_cor",
        lambda = compute_lambda(p_imputed_mvn_alone),
        lambda_se = compute_lambda_se(p_imputed_mvn_alone),
        type1_05 = compute_type1(p_imputed_mvn_alone, 0.05),
        type1_05_se = compute_prop_se(p_imputed_mvn_alone, 0.05),
        type1_01 = compute_type1(p_imputed_mvn_alone, 0.01),
        type1_01_se = compute_prop_se(p_imputed_mvn_alone, 0.01),
        stringsAsFactors = FALSE
      )
    }
  }

  if (verbose) cat("Block C complete.\n\n")
  bind_rows(results_list)
}

# ==============================================================================
# BLOCK D: ADAPTIVE OMNIBUS WITH TRAIN/TEST SPLIT
# ==============================================================================

# Learn which components to drop based on training null simulations
learn_component_mask <- function(training_null_pvals, lambda_threshold = 1.2,
                                 type1_threshold = 0.10) {
  methods <- c("acat", "fisher", "tfisher", "minp", "stouffer")
  keep_mask <- setNames(rep(TRUE, length(methods)), methods)

  for (method in methods) {
    lambda <- compute_lambda(training_null_pvals[, method])
    type1 <- compute_type1(training_null_pvals[, method], 0.05)

    if (!is.na(lambda) && lambda > lambda_threshold) {
      keep_mask[method] <- FALSE
    }
    if (!is.na(type1) && type1 > type1_threshold) {
      keep_mask[method] <- FALSE
    }
  }

  keep_mask
}

# Compute adaptive omnibus using only "good" components
compute_adaptive_omnibus <- function(gene_results, keep_mask, tg = tau_grid) {
  p_vals <- gene_results$P
  z_vals <- gene_results$ZSTAT

  component_ps <- c()

  if (keep_mask["acat"]) {
    component_ps <- c(component_ps, compute_acat(p_vals))
  }
  if (keep_mask["fisher"]) {
    component_ps <- c(component_ps, compute_fisher(p_vals))
  }
  if (keep_mask["tfisher"]) {
    component_ps <- c(component_ps, compute_tfisher(p_vals, tg))
  }
  if (keep_mask["minp"]) {
    component_ps <- c(component_ps, compute_minp(p_vals))
  }
  if (keep_mask["stouffer"]) {
    component_ps <- c(component_ps, compute_stouffer(z_vals))
  }

  if (length(component_ps) == 0) {
    return(0.5)
  }

  compute_omnibus(component_ps)
}

run_block_d <- function(n_train = 200L, n_test = 300L, b = 500L, seed = 42345L,
                        pathway_sizes_arg = pathway_sizes,
                        rho_values = c(0, 0.3, 0.7),
                        verbose = TRUE) {

  set.seed(seed)

  results_list <- list()

  if (verbose) cat("\n=== BLOCK D: Adaptive Omnibus with Train/Test Split ===\n")
  if (verbose) cat("Learning which components to drop from training null data\n\n")

  for (m in pathway_sizes_arg) {
    for (rho in rho_values) {
      if (verbose) cat(sprintf("  m=%d, rho=%.1f\n", m, rho))

      sigma_true <- make_cor_block(m, block_size = 20L, rho = rho)

      # ===== TRAINING PHASE =====
      # Generate training null data with BROKEN approach (analytic only)
      train_p_broken <- matrix(NA, nrow = n_train, ncol = 5)
      colnames(train_p_broken) <- c("acat", "fisher", "tfisher", "minp", "stouffer")

      for (s in seq_len(n_train)) {
        z <- simulate_gene_z(m, sigma_true, n_causal = 0, effect_size = 0)
        gene_res <- create_gene_results(z)
        res <- compute_components_analytic(gene_res)
        train_p_broken[s, ] <- c(res$acat, res$fisher, res$tfisher,
                                 res$minp, res$stouffer)
      }

      # Learn which components are miscalibrated
      keep_mask <- learn_component_mask(train_p_broken,
                                        lambda_threshold = 1.15,
                                        type1_threshold = 0.08)

      dropped <- names(keep_mask)[!keep_mask]
      drop_label <- if (length(dropped) > 0) {
        paste0("Drop: ", paste(dropped, collapse = ", "))
      } else {
        "Drop: none"
      }

      if (verbose) {
        if (length(dropped) > 0) {
          cat(sprintf("    Dropping: %s\n", paste(dropped, collapse = ", ")))
        } else {
          cat("    Keeping all components\n")
        }
      }

      # ===== TEST PHASE =====
      null_dist_correct <- precompute_null_distribution(
        m, sigma_true, b = b, two_sided = FALSE, seed = seed + m * 1000
      )

      p_naive_omni <- numeric(n_test)
      p_calibrated_omni_combined <- numeric(n_test)
      p_calibrated_omni_alone <- numeric(n_test)
      p_adaptive_omni <- numeric(n_test)
      p_adaptive_mvn_combined <- numeric(n_test)
      p_adaptive_mvn_alone <- numeric(n_test)

      # Precompute MVN null for adaptive statistic (fixed mask)
      null_adaptive <- numeric(b)
      for (i in seq_len(b)) {
        z_null <- simulate_gene_z(m, sigma_true, n_causal = 0, effect_size = 0)
        gene_null <- create_gene_results(z_null)
        null_adaptive[i] <- compute_adaptive_omnibus(gene_null, keep_mask)
      }

      for (s in seq_len(n_test)) {
        z <- simulate_gene_z(m, sigma_true, n_causal = 0, effect_size = 0)
        gene_res <- create_gene_results(z)

        # Naive omnibus (all analytic)
        res_ana <- compute_components_analytic(gene_res)
        p_naive_omni[s] <- res_ana$omnibus

        # MVN component calibration then recombination
        res_cal <- calibrate_with_null(res_ana, null_dist_correct)
        cal_components <- unlist(res_cal[c("acat", "fisher", "tfisher", "minp", "stouffer")])
        p_calibrated_omni_combined[s] <- compute_omnibus(cal_components)

        # Direct MVN calibration of the omnibus statistic itself
        p_calibrated_omni_alone[s] <- res_cal$omnibus

        # Adaptive omnibus (drop bad components; analytic p-value)
        p_adaptive_omni[s] <- compute_adaptive_omnibus(gene_res, keep_mask)

        keep_components <- names(keep_mask)[keep_mask]
        if (length(keep_components) > 0) {
          p_adaptive_mvn_combined[s] <- compute_omnibus(cal_components[keep_components])
        } else {
          p_adaptive_mvn_combined[s] <- 0.5
        }

        # Direct MVN calibration of the adaptive omnibus statistic
        p_adaptive_mvn_alone[s] <- empirical_p_from_null(p_adaptive_omni[s], null_adaptive)
      }

      results_list[[length(results_list) + 1]] <- data.frame(
        block = "D",
        pathway_size = m,
        rho = rho,
        method = "omnibus_analytical",
        n_components = 5,
        drop_label = drop_label,
        lambda = compute_lambda(p_naive_omni),
        lambda_se = compute_lambda_se(p_naive_omni),
        type1_05 = compute_type1(p_naive_omni, 0.05),
        type1_05_se = compute_prop_se(p_naive_omni, 0.05),
        type1_01 = compute_type1(p_naive_omni, 0.01),
        type1_01_se = compute_prop_se(p_naive_omni, 0.01),
        stringsAsFactors = FALSE
      )

      results_list[[length(results_list) + 1]] <- data.frame(
        block = "D",
        pathway_size = m,
        rho = rho,
        method = "omnibus_mvn_combined",
        n_components = 5,
        drop_label = drop_label,
        lambda = compute_lambda(p_calibrated_omni_combined),
        lambda_se = compute_lambda_se(p_calibrated_omni_combined),
        type1_05 = compute_type1(p_calibrated_omni_combined, 0.05),
        type1_05_se = compute_prop_se(p_calibrated_omni_combined, 0.05),
        type1_01 = compute_type1(p_calibrated_omni_combined, 0.01),
        type1_01_se = compute_prop_se(p_calibrated_omni_combined, 0.01),
        stringsAsFactors = FALSE
      )

      results_list[[length(results_list) + 1]] <- data.frame(
        block = "D",
        pathway_size = m,
        rho = rho,
        method = "omnibus_mvn_alone",
        n_components = 5,
        drop_label = drop_label,
        lambda = compute_lambda(p_calibrated_omni_alone),
        lambda_se = compute_lambda_se(p_calibrated_omni_alone),
        type1_05 = compute_type1(p_calibrated_omni_alone, 0.05),
        type1_05_se = compute_prop_se(p_calibrated_omni_alone, 0.05),
        type1_01 = compute_type1(p_calibrated_omni_alone, 0.01),
        type1_01_se = compute_prop_se(p_calibrated_omni_alone, 0.01),
        stringsAsFactors = FALSE
      )

      results_list[[length(results_list) + 1]] <- data.frame(
        block = "D",
        pathway_size = m,
        rho = rho,
        method = "omnibus_adaptive",
        n_components = sum(keep_mask),
        drop_label = drop_label,
        lambda = compute_lambda(p_adaptive_omni),
        lambda_se = compute_lambda_se(p_adaptive_omni),
        type1_05 = compute_type1(p_adaptive_omni, 0.05),
        type1_05_se = compute_prop_se(p_adaptive_omni, 0.05),
        type1_01 = compute_type1(p_adaptive_omni, 0.01),
        type1_01_se = compute_prop_se(p_adaptive_omni, 0.01),
        stringsAsFactors = FALSE
      )

      results_list[[length(results_list) + 1]] <- data.frame(
        block = "D",
        pathway_size = m,
        rho = rho,
        method = "omnibus_adaptive_mvn_combined",
        n_components = sum(keep_mask),
        drop_label = drop_label,
        lambda = compute_lambda(p_adaptive_mvn_combined),
        lambda_se = compute_lambda_se(p_adaptive_mvn_combined),
        type1_05 = compute_type1(p_adaptive_mvn_combined, 0.05),
        type1_05_se = compute_prop_se(p_adaptive_mvn_combined, 0.05),
        type1_01 = compute_type1(p_adaptive_mvn_combined, 0.01),
        type1_01_se = compute_prop_se(p_adaptive_mvn_combined, 0.01),
        stringsAsFactors = FALSE
      )

      results_list[[length(results_list) + 1]] <- data.frame(
        block = "D",
        pathway_size = m,
        rho = rho,
        method = "omnibus_adaptive_mvn_alone",
        n_components = sum(keep_mask),
        drop_label = drop_label,
        lambda = compute_lambda(p_adaptive_mvn_alone),
        lambda_se = compute_lambda_se(p_adaptive_mvn_alone),
        type1_05 = compute_type1(p_adaptive_mvn_alone, 0.05),
        type1_05_se = compute_prop_se(p_adaptive_mvn_alone, 0.05),
        type1_01 = compute_type1(p_adaptive_mvn_alone, 0.01),
        type1_01_se = compute_prop_se(p_adaptive_mvn_alone, 0.01),
        stringsAsFactors = FALSE
      )
    }
  }

  if (verbose) cat("Block D complete.\n\n")
  bind_rows(results_list)
}

# ==============================================================================
# BLOCK E: LEAVE-ONE-OUT OMNIBUS (MVN CALIBRATION)
# ==============================================================================

run_block_e <- function(n_sims = n_sims_null, b = b_perm, seed = 52345L,
                        pathway_sizes_arg = pathway_sizes,
                        verbose = TRUE) {

  set.seed(seed)

  comp_cols <- c("acat", "fisher", "tfisher", "minp", "stouffer")
  loo_variants <- c(
    "omnibus_minus_acat",
    "omnibus_minus_fisher",
    "omnibus_minus_tfisher",
    "omnibus_minus_minp",
    "omnibus_minus_stouffer",
    "omnibus_all"
  )
  mvn_modes <- c("combined", "alone")

  cor_types <- list(
    "LD_moderate" = function(m) make_cor_block(m, block_size = 20L, rho = 0.3),
    "LD_strong" = function(m) make_cor_block(m, block_size = 20L, rho = 0.7),
    "LD_independent" = function(m) make_cor_independent(m)
  )

  results_list <- list()

  if (verbose) cat("\n=== BLOCK E: Leave-one-out Omnibus (MVN) ===\n")
  if (verbose) cat("Precomputing null distributions...\n")

  null_cache <- list()
  for (m in pathway_sizes_arg) {
    for (cor_name in names(cor_types)) {
      sigma <- cor_types[[cor_name]](m)
      cache_key <- paste(m, cor_name, sep = "_")
      if (verbose) cat(sprintf("  Precomputing: m=%d, cor=%s\n", m, cor_name))

      null_stats <- precompute_null_distribution(
        m, sigma, b = b, two_sided = FALSE, seed = seed + m * 1000
      )
      null_comps <- null_stats[, comp_cols, drop = FALSE]

      null_loo <- list()
      null_loo[["omnibus_all"]] <- null_stats[, "omnibus"]
      for (comp in comp_cols) {
        keep_comps <- setdiff(comp_cols, comp)
        null_loo[[paste0("omnibus_minus_", comp)]] <-
          apply(null_comps[, keep_comps, drop = FALSE], 1, compute_omnibus)
      }
      null_cache[[cache_key]] <- list(
        null_stats = null_stats,
        null_loo = null_loo
      )
    }
  }

  for (m in pathway_sizes_arg) {
    for (cor_name in names(cor_types)) {
      if (verbose) cat(sprintf("  m=%d, cor=%s\n", m, cor_name))

      sigma <- cor_types[[cor_name]](m)
      cache_key <- paste(m, cor_name, sep = "_")
      null_bundle <- null_cache[[cache_key]]
      null_stats <- null_bundle$null_stats
      null_loo <- null_bundle$null_loo

      p_combined <- matrix(NA, nrow = n_sims, ncol = length(loo_variants))
      p_alone <- matrix(NA, nrow = n_sims, ncol = length(loo_variants))
      colnames(p_combined) <- colnames(p_alone) <- loo_variants

      for (s in seq_len(n_sims)) {
        z <- simulate_gene_z(m, sigma, n_causal = 0, effect_size = 0)
        gene_res <- create_gene_results(z)

        res_ana <- compute_components_analytic(gene_res)
        comp_vals <- unlist(res_ana[comp_cols])
        res_cal <- calibrate_with_null(res_ana, null_stats)
        comp_cal <- unlist(res_cal[comp_cols])

        obs_loo <- list(
          omnibus_all = res_ana$omnibus,
          omnibus_minus_acat = compute_omnibus(comp_vals[comp_cols != "acat"]),
          omnibus_minus_fisher = compute_omnibus(comp_vals[comp_cols != "fisher"]),
          omnibus_minus_tfisher = compute_omnibus(comp_vals[comp_cols != "tfisher"]),
          omnibus_minus_minp = compute_omnibus(comp_vals[comp_cols != "minp"]),
          omnibus_minus_stouffer = compute_omnibus(comp_vals[comp_cols != "stouffer"])
        )

        obs_loo_combined <- list(
          omnibus_all = compute_omnibus(comp_cal),
          omnibus_minus_acat = compute_omnibus(comp_cal[comp_cols != "acat"]),
          omnibus_minus_fisher = compute_omnibus(comp_cal[comp_cols != "fisher"]),
          omnibus_minus_tfisher = compute_omnibus(comp_cal[comp_cols != "tfisher"]),
          omnibus_minus_minp = compute_omnibus(comp_cal[comp_cols != "minp"]),
          omnibus_minus_stouffer = compute_omnibus(comp_cal[comp_cols != "stouffer"])
        )

        for (method in loo_variants) {
          obs_val <- obs_loo[[method]]
          null_vals <- null_loo[[method]]
          p_alone[s, method] <- empirical_p_from_null(obs_val, null_vals)
          p_combined[s, method] <- obs_loo_combined[[method]]
        }
      }

      for (method in loo_variants) {
        results_list[[length(results_list) + 1]] <- data.frame(
          block = "E",
          pathway_size = m,
          cor_structure = cor_name,
          method = method,
          calibration = "combined",
          lambda = compute_lambda(p_combined[, method]),
          lambda_se = compute_lambda_se(p_combined[, method]),
          type1_05 = compute_type1(p_combined[, method], 0.05),
          type1_05_se = compute_prop_se(p_combined[, method], 0.05),
          type1_01 = compute_type1(p_combined[, method], 0.01),
          type1_01_se = compute_prop_se(p_combined[, method], 0.01),
          stringsAsFactors = FALSE
        )

        results_list[[length(results_list) + 1]] <- data.frame(
          block = "E",
          pathway_size = m,
          cor_structure = cor_name,
          method = method,
          calibration = "alone",
          lambda = compute_lambda(p_alone[, method]),
          lambda_se = compute_lambda_se(p_alone[, method]),
          type1_05 = compute_type1(p_alone[, method], 0.05),
          type1_05_se = compute_prop_se(p_alone[, method], 0.05),
          type1_01 = compute_type1(p_alone[, method], 0.01),
          type1_01_se = compute_prop_se(p_alone[, method], 0.01),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (verbose) cat("Block E complete.\n\n")
  bind_rows(results_list)
}

# ==============================================================================
# PLOTTING FUNCTIONS
# ==============================================================================

plot_block_a_null <- function(results) {
  method_order <- c(
    "ACAT", "FISHER", "MINP", "STOUFFER", "TFISHER",
    "OMNIBUS\nCombined", "OMNIBUS\nAlone"
  )
  cor_order <- c("LD_moderate", "LD_strong", "LD_independent")
  df <- results %>%
    filter(block == "A_null") %>%
    mutate(
      lambda_plot = pmin(lambda, lambda_cap),
      cor_structure = factor(cor_structure, levels = cor_order),
      method_label = factor(dplyr::case_when(
        method == "acat" ~ "ACAT",
        method == "fisher" ~ "FISHER",
        method == "minp" ~ "MINP",
        method == "stouffer" ~ "STOUFFER",
        method == "tfisher" ~ "TFISHER",
        method == "omnibus_combined" ~ "OMNIBUS\nCombined",
        method == "omnibus_alone" ~ "OMNIBUS\nAlone",
        TRUE ~ toupper(method)
      ), levels = method_order),
      fill_calibration = dplyr::if_else(calibration == "analytic", "analytic", "mvn"),
      label_lambda = ifelse(lambda > lambda_cap,
                            paste0(">", lambda_cap),
                            sprintf("%.2f", lambda_plot)),
      label_type1 = sprintf("%.3f", type1_05),
      lambda_lower = pmax(lambda_plot - lambda_se, 0),
      lambda_upper = pmin(lambda_plot + lambda_se, lambda_cap),
      type1_lower = pmax(type1_05 - type1_05_se, 0),
      type1_upper = pmin(type1_05 + type1_05_se, 1)
    )

  p_lambda <- ggplot(df, aes(x = method_label, y = lambda_plot, fill = fill_calibration)) +
    geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
    geom_errorbar(aes(ymin = lambda_lower, ymax = lambda_upper),
                  position = position_dodge(0.8), width = 0.2) +
    geom_text(aes(label = label_lambda),
              position = position_dodge(0.8),
              angle = 90, vjust = -0.2, size = 2.3) +
    geom_vline(xintercept = 6.5, linetype = "dotted", color = "grey40", linewidth = 0.5) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.7) +
    facet_grid(pathway_size ~ cor_structure,
               labeller = labeller(pathway_size = function(x) paste0("m=", x))) +
    scale_fill_manual(values = c("analytic" = "#E41A1C", "mvn" = "#377EB8"),
                      labels = c("Analytic", "MVN Calibrated")) +
    labs(title = "Block A: Lambda Under Null",
         x = "Method", y = expression(lambda), fill = "Calibration") +
    sim_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

  p_type1 <- ggplot(df, aes(x = method_label, y = type1_05, fill = fill_calibration)) +
    geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
    geom_errorbar(aes(ymin = type1_lower, ymax = type1_upper),
                  position = position_dodge(0.8), width = 0.2) +
    geom_text(aes(label = label_type1),
              position = position_dodge(0.8),
              angle = 90, vjust = -0.2, size = 2.3) +
    geom_vline(xintercept = 6.5, linetype = "dotted", color = "grey40", linewidth = 0.5) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.7) +
    facet_grid(pathway_size ~ cor_structure,
               labeller = labeller(pathway_size = function(x) paste0("m=", x))) +
    scale_fill_manual(values = c("analytic" = "#E41A1C", "mvn" = "#377EB8"),
                      labels = c("Analytic", "MVN Calibrated")) +
    labs(title = "Block A: Type I Error Under Null",
         x = "Method", y = "Type I Error Rate", fill = "Calibration") +
    sim_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

  list(lambda = p_lambda, type1 = p_type1)
}

# ==============================================================================
# BLOCK F: RANDOM TWO-COMPONENT DROPOUT (MVN CALIBRATION)
# ==============================================================================

run_block_f <- function(n_sims = n_sims_null, b = b_perm, seed = 62345L,
                        pathway_sizes_arg = pathway_sizes,
                        verbose = TRUE) {

  set.seed(seed)

  comp_cols <- c("acat", "fisher", "tfisher", "minp", "stouffer")

  cor_types <- list(
    "LD_moderate" = function(m) make_cor_block(m, block_size = 20L, rho = 0.3),
    "LD_strong" = function(m) make_cor_block(m, block_size = 20L, rho = 0.7),
    "LD_independent" = function(m) make_cor_independent(m)
  )

  results_list <- list()

  if (verbose) cat("\n=== BLOCK F: Random Two-Component Dropout (MVN) ===\n")
  if (verbose) cat("Precomputing null distributions...\n")

  null_cache <- list()
  for (m in pathway_sizes_arg) {
    for (cor_name in names(cor_types)) {
      sigma <- cor_types[[cor_name]](m)
      cache_key <- paste(m, cor_name, sep = "_")
      if (verbose) cat(sprintf("  Precomputing: m=%d, cor=%s\n", m, cor_name))

      null_stats <- precompute_null_distribution(
        m, sigma, b = b, two_sided = FALSE, seed = seed + m * 1000
      )
      null_comps <- null_stats[, comp_cols, drop = FALSE]

      set.seed(seed + m * 2000 + nchar(cor_name))
      null_random2drop <- vapply(seq_len(nrow(null_comps)), function(i) {
        drop_pair <- sample(comp_cols, 2, replace = FALSE)
        keep <- setdiff(comp_cols, drop_pair)
        compute_omnibus(null_comps[i, keep])
      }, numeric(1))

      null_cache[[cache_key]] <- list(
        omnibus_all = null_stats[, "omnibus"],
        omnibus_random2drop = null_random2drop
      )
    }
  }

  for (m in pathway_sizes_arg) {
    for (cor_name in names(cor_types)) {
      if (verbose) cat(sprintf("  m=%d, cor=%s\n", m, cor_name))

      sigma <- cor_types[[cor_name]](m)
      cache_key <- paste(m, cor_name, sep = "_")
      null_loo <- null_cache[[cache_key]]

      p_mvn <- matrix(NA, nrow = n_sims, ncol = 2)
      colnames(p_mvn) <- c("omnibus_all", "omnibus_random2drop")

      for (s in seq_len(n_sims)) {
        z <- simulate_gene_z(m, sigma, n_causal = 0, effect_size = 0)
        gene_res <- create_gene_results(z)

        res_ana <- compute_components_analytic(gene_res)
        comp_vals <- unlist(res_ana[comp_cols])

        drop_pair <- sample(comp_cols, 2, replace = FALSE)
        keep <- setdiff(comp_cols, drop_pair)
        obs_random2drop <- compute_omnibus(comp_vals[keep])

        p_mvn[s, "omnibus_all"] <-
          (1 + sum(null_loo$omnibus_all <= res_ana$omnibus, na.rm = TRUE)) /
          (sum(!is.na(null_loo$omnibus_all)) + 1)

        p_mvn[s, "omnibus_random2drop"] <-
          (1 + sum(null_loo$omnibus_random2drop <= obs_random2drop, na.rm = TRUE)) /
          (sum(!is.na(null_loo$omnibus_random2drop)) + 1)
      }

      for (method in colnames(p_mvn)) {
        results_list[[length(results_list) + 1]] <- data.frame(
          block = "F",
          pathway_size = m,
          cor_structure = cor_name,
          method = method,
          calibration = "mvn",
          lambda = compute_lambda(p_mvn[, method]),
          lambda_se = compute_lambda_se(p_mvn[, method]),
          type1_05 = compute_type1(p_mvn[, method], 0.05),
          type1_05_se = compute_prop_se(p_mvn[, method], 0.05),
          type1_01 = compute_type1(p_mvn[, method], 0.01),
          type1_01_se = compute_prop_se(p_mvn[, method], 0.01),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (verbose) cat("Block F complete.\n\n")
  bind_rows(results_list)
}

plot_block_f <- function(results) {
  cor_order <- c("LD_moderate", "LD_strong", "LD_independent")
  method_order <- c("omnibus_all", "omnibus_random2drop")
  method_labels <- c("Omnibus (all)", "Omnibus (drop 2)")

  df <- results %>%
    filter(block == "F") %>%
    mutate(
      cor_structure = factor(cor_structure, levels = cor_order),
      method = factor(method, levels = method_order),
      method_label = factor(method_labels[match(as.character(method), method_order)],
                            levels = method_labels),
      lambda_plot = pmin(lambda, lambda_cap),
      label_lambda = ifelse(lambda > lambda_cap,
                            paste0(">", lambda_cap),
                            sprintf("%.2f", lambda_plot)),
      label_type1 = sprintf("%.3f", type1_05),
      lambda_lower = pmax(lambda_plot - lambda_se, 0),
      lambda_upper = pmin(lambda_plot + lambda_se, lambda_cap),
      type1_lower = pmax(type1_05 - type1_05_se, 0),
      type1_upper = pmin(type1_05 + type1_05_se, 1)
    )

  p_lambda <- ggplot(df, aes(x = method_label, y = lambda_plot, fill = method_label)) +
    geom_col(width = 0.7, alpha = 0.85) +
    geom_errorbar(aes(ymin = lambda_lower, ymax = lambda_upper), width = 0.2) +
    geom_text(aes(label = label_lambda),
              angle = 90, vjust = -0.2, size = 2.3) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.7) +
    facet_grid(pathway_size ~ cor_structure,
               labeller = labeller(pathway_size = function(x) paste0("m=", x))) +
    labs(title = "Block F: Random Two-Component Dropout (MVN)",
         x = "Omnibus Variant", y = expression(lambda)) +
    sim_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          legend.position = "none")

  p_type1 <- ggplot(df, aes(x = method_label, y = type1_05, fill = method_label)) +
    geom_col(width = 0.7, alpha = 0.85) +
    geom_errorbar(aes(ymin = type1_lower, ymax = type1_upper), width = 0.2) +
    geom_text(aes(label = label_type1),
              angle = 90, vjust = -0.2, size = 2.3) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.7) +
    facet_grid(pathway_size ~ cor_structure,
               labeller = labeller(pathway_size = function(x) paste0("m=", x))) +
    labs(title = "Block F: Type I Error (Random Two-Component Dropout)",
         x = "Omnibus Variant", y = "Type I Error Rate") +
    sim_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          legend.position = "none")

  list(lambda = p_lambda, type1 = p_type1)
}

plot_block_a_power <- function(results) {
  method_order <- c("ACAT", "FISHER", "MINP", "STOUFFER", "TFISHER", "OMNIBUS")
  df <- results %>%
    filter(block == "A", calibration == "mvn") %>%
    mutate(
      method_label = factor(toupper(method), levels = method_order)
    )

  ggplot(df, aes(x = effect_size, y = power_05,
                 color = method_label, linetype = method_label)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    facet_grid(pathway_size ~ signal_type,
               labeller = labeller(pathway_size = function(x) paste0("m=", x))) +
    scale_color_manual(values = c(
      "ACAT" = "#E41A1C", "FISHER" = "#377EB8", "TFISHER" = "#4DAF4A",
      "MINP" = "#984EA3", "STOUFFER" = "#FF7F00", "OMNIBUS" = "#000000"
    )) +
    labs(title = "Block A: Power (MVN Calibration)",
         x = "Effect Size", y = "Power",
         color = "Method", linetype = "Method") +
    sim_theme +
    theme(legend.position = "right")
}

plot_block_b <- function(results) {
  method_order <- c("ACAT", "FISHER", "MINP", "STOUFFER", "TFISHER", "OMNIBUS")
  df <- results %>%
    filter(block == "B") %>%
    mutate(
      lambda_plot = pmin(lambda, lambda_cap),
      method_label = factor(toupper(method), levels = method_order),
      label_lambda = ifelse(lambda > lambda_cap,
                            paste0(">", lambda_cap),
                            sprintf("%.2f", lambda_plot)),
      label_type1 = sprintf("%.3f", type1_05),
      lambda_lower = pmax(lambda_plot - lambda_se, 0),
      lambda_upper = pmin(lambda_plot + lambda_se, lambda_cap),
      type1_lower = pmax(type1_05 - type1_05_se, 0),
      type1_upper = pmin(type1_05 + type1_05_se, 1),
      rho_label = dplyr::case_when(
        rho == 0   ~ "LD_independent",
        rho == 0.3 ~ "LD_moderate",
        rho == 0.7 ~ "LD_strong",
        TRUE       ~ paste0("rho = ", rho)
      ),
      rho_label = factor(rho_label, levels = c("LD_moderate", "LD_strong", "LD_independent"))
    )

  p_lambda <- ggplot(df, aes(x = method_label, y = lambda_plot, fill = scenario)) +
    geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
    geom_errorbar(aes(ymin = lambda_lower, ymax = lambda_upper),
                  position = position_dodge(0.8), width = 0.2) +
    geom_text(aes(label = label_lambda),
              position = position_dodge(0.8),
              angle = 90, vjust = -0.2, size = 2.3) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.7) +
    facet_grid(pathway_size ~ rho_label,
               labeller = labeller(pathway_size = function(x) paste0("m=", x))) +
    scale_fill_manual(
      values = c("analytic_broken" = "#E41A1C", "mvn" = "#377EB8"),
      labels = c("Analytic (ACAT+minP+Stouffer broken)", "MVN Calibrated")
    ) +
    labs(title = "Block B: Lambda - Analytic vs MVN",
         x = "Method", y = expression(lambda), fill = "Scenario") +
    sim_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

  p_type1 <- ggplot(df, aes(x = method_label, y = type1_05, fill = scenario)) +
    geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
    geom_errorbar(aes(ymin = type1_lower, ymax = type1_upper),
                  position = position_dodge(0.8), width = 0.2) +
    geom_text(aes(label = label_type1),
              position = position_dodge(0.8),
              angle = 90, vjust = -0.2, size = 2.3) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.7) +
    facet_grid(pathway_size ~ rho_label,
               labeller = labeller(pathway_size = function(x) paste0("m=", x))) +
    scale_fill_manual(
      values = c("analytic_broken" = "#E41A1C", "mvn" = "#377EB8"),
      labels = c("Analytic (ACAT+TFisher broken)", "MVN Calibrated")
    ) +
    labs(title = "Block B: Type I Error - Analytic vs MVN",
         x = "Method", y = "Type I Error Rate", fill = "Scenario") +
    sim_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

  list(lambda = p_lambda, type1 = p_type1)
}

plot_block_c <- function(results) {
  df <- results %>%
    filter(block == "C") %>%
    mutate(
      lambda_plot = pmin(lambda, lambda_cap),
      omni_variant = factor(dplyr::case_when(
        method == "omnibus_combined" ~ "Combined",
        method == "omnibus_alone" ~ "Alone",
        TRUE ~ NA_character_
      ), levels = c("Combined", "Alone")),
      pathway_size = as.integer(pathway_size),
      facet_label = factor(paste0("m=", pathway_size),
                           levels = paste0("m=", sort(unique(pathway_size))))
    )

  df_omni <- df %>%
    filter(method %in% c("omnibus_combined", "omnibus_alone")) %>%
    mutate(
      lambda_lower = pmax(lambda_plot - lambda_se, 0),
      lambda_upper = pmin(lambda_plot + lambda_se, lambda_cap),
      type1_lower = pmax(type1_05 - type1_05_se, 0),
      type1_upper = pmin(type1_05 + type1_05_se, 1)
    )

  p_lambda <- ggplot(df_omni, aes(x = missing_frac, y = lambda_plot, color = strategy,
                                  linetype = omni_variant, shape = omni_variant)) +
    geom_ribbon(aes(ymin = lambda_lower, ymax = lambda_upper, fill = strategy),
                alpha = 0.15, show.legend = FALSE) +
    geom_point(size = 2.8) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
    geom_ribbon(aes(ymin = 0.95, ymax = 1.05), alpha = 0.1, fill = "gray", color = NA) +
    facet_wrap(~ facet_label, ncol = 3) +
    scale_color_manual(
      values = c("analytic_fallback" = "#E41A1C",
                 "mvn_true_cor" = "#377EB8",
                 "mvn_imputed_cor" = "#4DAF4A"),
      labels = c("Analytic Fallback", "MVN True Cor", "MVN Imputed (0)")
    ) +
    scale_linetype_manual(values = c("Combined" = "solid", "Alone" = "dashed")) +
    scale_shape_manual(values = c("Combined" = 16, "Alone" = 17)) +
    labs(title = "Block C: Omnibus Lambda vs Missing Correlations",
         x = "Fraction Missing", y = expression(lambda),
         color = "Strategy", linetype = "Omnibus", shape = "Omnibus") +
    sim_theme

  p_type1 <- ggplot(df_omni, aes(x = missing_frac, y = type1_05, color = strategy,
                                 linetype = omni_variant, shape = omni_variant)) +
    geom_ribbon(aes(ymin = type1_lower, ymax = type1_upper, fill = strategy),
                alpha = 0.15, show.legend = FALSE) +
    geom_point(size = 2.8) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "black", linewidth = 0.5) +
    facet_wrap(~ facet_label, ncol = 3) +
    scale_color_manual(
      values = c("analytic_fallback" = "#E41A1C",
                 "mvn_true_cor" = "#377EB8",
                 "mvn_imputed_cor" = "#4DAF4A"),
      labels = c("Analytic Fallback", "MVN True Cor", "MVN Imputed (0)")
    ) +
    scale_linetype_manual(values = c("Combined" = "solid", "Alone" = "dashed")) +
    scale_shape_manual(values = c("Combined" = 16, "Alone" = 17)) +
    labs(title = "Block C: Omnibus Type I Error vs Missing Correlations",
         x = "Fraction Missing", y = "Type I Error Rate",
         color = "Strategy", linetype = "Omnibus", shape = "Omnibus") +
    sim_theme

  list(lambda = p_lambda, type1 = p_type1)
}

plot_block_d <- function(results) {
  df <- results %>%
    filter(block == "D") %>%
    mutate(
      rho_label = dplyr::case_when(
        rho == 0   ~ "LD_independent",
        rho == 0.3 ~ "LD_moderate",
        rho == 0.7 ~ "LD_strong",
        TRUE       ~ paste0("rho = ", rho)
      ),
      rho_label = factor(rho_label, levels = c("LD_moderate", "LD_strong", "LD_independent")),
      method_label = dplyr::case_when(
        method == "omnibus_analytical"    ~ "Analytical",
        method == "omnibus_mvn_combined"  ~ "MVN Combined",
        method == "omnibus_mvn_alone"     ~ "MVN Alone",
        method == "omnibus_adaptive"      ~ "Adaptive+Analytical",
        method == "omnibus_adaptive_mvn_combined" ~ "Adaptive+Combined",
        method == "omnibus_adaptive_mvn_alone"    ~ "Adaptive+Alone",
        TRUE                              ~ method
      ),
      method_label = factor(method_label,
                            levels = c("Analytical", "MVN Combined", "MVN Alone",
                                       "Adaptive+Analytical", "Adaptive+Combined",
                                       "Adaptive+Alone")),
      lambda_plot = pmin(lambda, lambda_cap),
      label_lambda = ifelse(lambda > lambda_cap,
                            paste0(">", lambda_cap),
                            sprintf("%.2f", lambda_plot)),
      label_type1 = sprintf("%.3f", type1_05),
      lambda_lower = pmax(lambda_plot - lambda_se, 0),
      lambda_upper = pmin(lambda_plot + lambda_se, lambda_cap),
      type1_lower = pmax(type1_05 - type1_05_se, 0),
      type1_upper = pmin(type1_05 + type1_05_se, 1)
    )

  drop_labels <- df %>%
    filter(method == "omnibus_adaptive") %>%
    group_by(pathway_size, rho_label) %>%
    summarise(
      drop_label = dplyr::first(drop_label),
      .groups = "drop"
    )

  p_comparison <- ggplot(df, aes(x = method_label, y = lambda_plot, fill = method_label)) +
    geom_col(width = 0.7, alpha = 0.8) +
    geom_errorbar(aes(ymin = lambda_lower, ymax = lambda_upper), width = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.7) +
    geom_text(aes(label = label_lambda),
              angle = 90, vjust = -0.2, size = 2.3) +
    geom_text(data = drop_labels,
              aes(x = 2, y = Inf, label = drop_label),
              inherit.aes = FALSE, size = 2.3, vjust = 1.2) +
    facet_grid(pathway_size ~ rho_label,
               labeller = labeller(pathway_size = function(x) paste0("m=", x))) +
    scale_fill_manual(
      values = c("Analytical" = "#E41A1C",
                 "MVN Combined" = "#377EB8",
                 "MVN Alone" = "#1F78B4",
                 "Adaptive+Analytical" = "#4DAF4A",
                 "Adaptive+Combined" = "#984EA3",
                 "Adaptive+Alone" = "#A65628")
    ) +
    labs(title = "Block D: Adaptive Omnibus Comparison",
         x = "Method", y = expression(lambda)) +
    sim_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          legend.position = "none") +
    coord_cartesian(clip = "off")

  p_type1 <- ggplot(df, aes(x = method_label, y = type1_05, fill = method_label)) +
    geom_col(width = 0.7, alpha = 0.8) +
    geom_errorbar(aes(ymin = type1_lower, ymax = type1_upper), width = 0.2) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.7) +
    geom_text(aes(label = label_type1),
              angle = 90, vjust = -0.2, size = 2.3) +
    geom_text(data = drop_labels,
              aes(x = 2, y = 0.31, label = drop_label),
              inherit.aes = FALSE, size = 2.3) +
    facet_grid(pathway_size ~ rho_label,
               labeller = labeller(pathway_size = function(x) paste0("m=", x))) +
    scale_fill_manual(
      values = c("Analytical" = "#E41A1C",
                 "MVN Combined" = "#377EB8",
                 "MVN Alone" = "#1F78B4",
                 "Adaptive+Analytical" = "#4DAF4A",
                 "Adaptive+Combined" = "#984EA3",
                 "Adaptive+Alone" = "#A65628")
    ) +
    labs(title = "Block D: Type I Error Comparison",
         x = "Method", y = "Type I Error Rate") +
    sim_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          legend.position = "none")

  list(comparison = p_comparison, type1 = p_type1)
}

plot_block_e <- function(results) {
  cor_order <- c("LD_moderate", "LD_strong", "LD_independent")
  method_order <- c(
    "omnibus_minus_acat",
    "omnibus_minus_fisher",
    "omnibus_minus_tfisher",
    "omnibus_minus_minp",
    "omnibus_minus_stouffer",
    "omnibus_all"
  )
  method_labels <- c(
    "Omnibus - ACAT",
    "Omnibus - Fisher",
    "Omnibus - TFisher",
    "Omnibus - minP",
    "Omnibus - Stouffer",
    "Omnibus (all)"
  )

  df <- results %>%
    filter(block == "E") %>%
    mutate(
      cor_structure = factor(cor_structure, levels = cor_order),
      method = factor(method, levels = method_order),
      method_label = paste0(
        method_labels[match(as.character(method), method_order)],
        ifelse(calibration == "combined", " (Combined)", " (Alone)")
      ),
      method_label = factor(
        method_label,
        levels = c(
          paste0(method_labels, " (Combined)"),
          paste0(method_labels, " (Alone)")
        )
      ),
      lambda_plot = pmin(lambda, lambda_cap),
      label_lambda = ifelse(lambda > lambda_cap,
                            paste0(">", lambda_cap),
                            sprintf("%.2f", lambda_plot)),
      label_type1 = sprintf("%.3f", type1_05),
      lambda_lower = pmax(lambda_plot - lambda_se, 0),
      lambda_upper = pmin(lambda_plot + lambda_se, lambda_cap),
      type1_lower = pmax(type1_05 - type1_05_se, 0),
      type1_upper = pmin(type1_05 + type1_05_se, 1)
    )

  p_lambda <- ggplot(df, aes(x = method_label, y = lambda_plot, fill = calibration)) +
    geom_col(width = 0.7, alpha = 0.85) +
    geom_errorbar(aes(ymin = lambda_lower, ymax = lambda_upper), width = 0.2) +
    geom_text(aes(label = label_lambda),
              angle = 90, vjust = -0.2, size = 2.3) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.7) +
    facet_grid(pathway_size ~ cor_structure,
               labeller = labeller(pathway_size = function(x) paste0("m=", x))) +
    scale_fill_manual(values = c("combined" = "#377EB8", "alone" = "#8C564B")) +
    labs(title = "Block E: Leave-one-out Omnibus Lambda (MVN)",
         x = "Omnibus Variant", y = expression(lambda), fill = "MVN mode") +
    sim_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

  p_type1 <- ggplot(df, aes(x = method_label, y = type1_05, fill = calibration)) +
    geom_col(width = 0.7, alpha = 0.85) +
    geom_errorbar(aes(ymin = type1_lower, ymax = type1_upper), width = 0.2) +
    geom_text(aes(label = label_type1),
              angle = 90, vjust = -0.2, size = 2.3) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.7) +
    facet_grid(pathway_size ~ cor_structure,
               labeller = labeller(pathway_size = function(x) paste0("m=", x))) +
    scale_fill_manual(values = c("combined" = "#4DAF4A", "alone" = "#A65628")) +
    labs(title = "Block E: Leave-one-out Omnibus Type I Error (MVN)",
         x = "Omnibus Variant", y = "Type I Error Rate", fill = "MVN mode") +
    sim_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

  list(lambda = p_lambda, type1 = p_type1)
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

run_all_simulations <- function(output_dir = "simulation_results",
                                run_block_a_flag = TRUE,
                                run_block_b_flag = TRUE,
                                run_block_c_flag = TRUE,
                                run_block_d_flag = TRUE,
                                run_block_e_flag = TRUE,
                                run_block_f_flag = TRUE,
                                reduced = FALSE) {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Set parameters based on reduced flag
  if (reduced) {
    n_null <- 50L
    n_alt <- 30L
    b_sim <- 100L
    pathway_sizes_run <- c(20L, 50L)
    effect_sizes_run <- c(1.0, 2.0)
  } else {
    n_null <- n_sims_null
    n_alt <- n_sims_alt
    b_sim <- b_perm
    pathway_sizes_run <- pathway_sizes
    effect_sizes_run <- effect_sizes
  }

  results <- list()

  # Block A
  if (run_block_a_flag) {
    cat("\n", strrep("=", 60), "\n")
    cat("RUNNING BLOCK A: POWER ACROSS SIGNAL SHAPES\n")
    cat(strrep("=", 60), "\n")

    results$a_null <- run_block_a_null(
      n_sims = n_null, b = b_sim, seed = seed_base + 1000,
      pathway_sizes_arg = pathway_sizes_run
    )

    results$a_power <- run_block_a(
      n_sims = n_alt, b = b_sim, seed = seed_base,
      pathway_sizes_arg = pathway_sizes_run,
      effect_sizes_arg = effect_sizes_run
    )

    saveRDS(results$a_null, file.path(output_dir, "block_a_null.rds"))
    saveRDS(results$a_power, file.path(output_dir, "block_a_power.rds"))

    plots_a_null <- plot_block_a_null(results$a_null)
    ggsave(file.path(output_dir, "block_a_null_lambda.png"),
           plots_a_null$lambda, width = 14, height = 10, dpi = 300, bg = "white")
    ggsave(file.path(output_dir, "block_a_null_type1.png"),
           plots_a_null$type1, width = 14, height = 10, dpi = 300, bg = "white")

    plot_a_power <- plot_block_a_power(results$a_power)
    ggsave(file.path(output_dir, "block_a_power.png"),
           plot_a_power, width = 16, height = 12, dpi = 300, bg = "white")

    cat("Block A complete.\n")
  }

  # Block B
  if (run_block_b_flag) {
    cat("\n", strrep("=", 60), "\n")
    cat("RUNNING BLOCK B: ANALYTIC VS MVN STRESS TEST\n")
    cat(strrep("=", 60), "\n")

    results$b <- run_block_b(
      n_sims = n_null, b = b_sim, seed = seed_base + 2000,
      pathway_sizes_arg = if (reduced) c(20L, 50L) else pathway_sizes_run,
      rho_values = c(0, 0.3, 0.7)
    )

    saveRDS(results$b, file.path(output_dir, "block_b.rds"))

    plots_b <- plot_block_b(results$b)
    ggsave(file.path(output_dir, "block_b_lambda.png"),
           plots_b$lambda, width = 14, height = 10, dpi = 300, bg = "white")
    ggsave(file.path(output_dir, "block_b_type1.png"),
           plots_b$type1, width = 14, height = 10, dpi = 300, bg = "white")

    cat("Block B complete.\n")
  }

  # Block C
  if (run_block_c_flag) {
    cat("\n", strrep("=", 60), "\n")
    cat("RUNNING BLOCK C: MISSING CORRELATION SCENARIOS\n")
    cat(strrep("=", 60), "\n")

    results$c <- run_block_c(
      n_sims = n_null, b = b_sim, seed = seed_base + 3000,
      pathway_sizes_arg = if (reduced) c(10L, 30L, 50L) else pathway_sizes_run,
      missing_fracs = c(0, 0.1, 0.3, 0.5, 0.7)
    )

    saveRDS(results$c, file.path(output_dir, "block_c.rds"))

    plots_c <- plot_block_c(results$c)
    ggsave(file.path(output_dir, "block_c_lambda.png"),
           plots_c$lambda, width = 10, height = 8, dpi = 300, bg = "white")
    ggsave(file.path(output_dir, "block_c_type1.png"),
           plots_c$type1, width = 10, height = 8, dpi = 300, bg = "white")

    cat("Block C complete.\n")
  }

  # Block D
  if (run_block_d_flag) {
    cat("\n", strrep("=", 60), "\n")
    cat("RUNNING BLOCK D: ADAPTIVE OMNIBUS\n")
    cat(strrep("=", 60), "\n")

    results$d <- run_block_d(
      n_train = if (reduced) 50L else 200L,
      n_test = if (reduced) 100L else 300L,
      b = b_sim, seed = seed_base + 4000,
      pathway_sizes_arg = if (reduced) c(10L, 30L, 50L) else pathway_sizes_run,
      rho_values = c(0, 0.3, 0.7)
    )

    saveRDS(results$d, file.path(output_dir, "block_d.rds"))

    plots_d <- plot_block_d(results$d)
    ggsave(file.path(output_dir, "block_d_comparison.png"),
           plots_d$comparison, width = 12, height = 10, dpi = 300, bg = "white")
    ggsave(file.path(output_dir, "block_d_type1.png"),
           plots_d$type1, width = 12, height = 10, dpi = 300, bg = "white")

    cat("Block D complete.\n")
  }

  # Block E
  if (run_block_e_flag) {
    cat("\n", strrep("=", 60), "\n")
    cat("RUNNING BLOCK E: LEAVE-ONE-OUT OMNIBUS (MVN)\n")
    cat(strrep("=", 60), "\n")

    results$e <- run_block_e(
      n_sims = n_null, b = b_sim, seed = seed_base + 5000,
      pathway_sizes_arg = if (reduced) c(20L, 50L) else pathway_sizes_run
    )

    saveRDS(results$e, file.path(output_dir, "block_e.rds"))

    plots_e <- plot_block_e(results$e)
    ggsave(file.path(output_dir, "block_e_lambda.png"),
           plots_e$lambda, width = 14, height = 10, dpi = 300, bg = "white")
    ggsave(file.path(output_dir, "block_e_type1.png"),
           plots_e$type1, width = 14, height = 10, dpi = 300, bg = "white")

    cat("Block E complete.\n")
  }

  # Block F
  if (run_block_f_flag) {
    cat("\n", strrep("=", 60), "\n")
    cat("RUNNING BLOCK F: RANDOM TWO-COMPONENT DROPOUT (MVN)\n")
    cat(strrep("=", 60), "\n")

    results$f <- run_block_f(
      n_sims = n_null, b = b_sim, seed = seed_base + 6000,
      pathway_sizes_arg = if (reduced) c(20L, 50L) else pathway_sizes_run
    )

    saveRDS(results$f, file.path(output_dir, "block_f.rds"))

    plots_f <- plot_block_f(results$f)
    ggsave(file.path(output_dir, "block_f_lambda.png"),
           plots_f$lambda, width = 14, height = 10, dpi = 300, bg = "white")
    ggsave(file.path(output_dir, "block_f_type1.png"),
           plots_f$type1, width = 14, height = 10, dpi = 300, bg = "white")

    cat("Block F complete.\n")
  }

  # Summary
  cat("\n", strrep("=", 60), "\n")
  cat("SIMULATION SUMMARY\n")
  cat(strrep("=", 60), "\n")

  if (!is.null(results$a_null)) {
    cat("\n--- Block A (Null) Summary ---\n")
    summ_a <- results$a_null %>%
      group_by(cor_structure, calibration) %>%
      summarise(
        mean_lambda = mean(lambda, na.rm = TRUE),
        mean_type1 = mean(type1_05, na.rm = TRUE),
        .groups = "drop"
      )
    print(summ_a)
  }

  if (!is.null(results$b)) {
    cat("\n--- Block B Summary (Analytic vs MVN) ---\n")
    summ_b <- results$b %>%
      group_by(rho, scenario, method) %>%
      summarise(
        mean_lambda = mean(lambda, na.rm = TRUE),
        mean_type1 = mean(type1_05, na.rm = TRUE),
        .groups = "drop"
      )
    print(summ_b)
  }

  if (!is.null(results$c)) {
    cat("\n--- Block C Summary (Omnibus) ---\n")
    summ_c <- results$c %>%
      filter(method == "omnibus") %>%
      group_by(missing_frac, strategy) %>%
      summarise(
        mean_lambda = mean(lambda, na.rm = TRUE),
        mean_type1 = mean(type1_05, na.rm = TRUE),
        .groups = "drop"
      )
    print(summ_c)
  }

  if (!is.null(results$d)) {
    cat("\n--- Block D Summary (Adaptive Omnibus) ---\n")
    summ_d <- results$d %>%
      group_by(rho, method) %>%
      summarise(
        mean_lambda = mean(lambda, na.rm = TRUE),
        mean_type1 = mean(type1_05, na.rm = TRUE),
        .groups = "drop"
      )
    print(summ_d)
  }

  if (!is.null(results$e)) {
    cat("\n--- Block E Summary (Leave-one-out Omnibus) ---\n")
    summ_e <- results$e %>%
      group_by(cor_structure, method) %>%
      summarise(
        mean_lambda = mean(lambda, na.rm = TRUE),
        mean_type1 = mean(type1_05, na.rm = TRUE),
        .groups = "drop"
      )
    print(summ_e)
  }

  if (!is.null(results$f)) {
    cat("\n--- Block F Summary (Random Two-Component Dropout) ---\n")
    summ_f <- results$f %>%
      group_by(cor_structure, method) %>%
      summarise(
        mean_lambda = mean(lambda, na.rm = TRUE),
        mean_type1 = mean(type1_05, na.rm = TRUE),
        .groups = "drop"
      )
    print(summ_f)
  }

  cat("\n", strrep("=", 60), "\n")
  cat("ALL SIMULATIONS COMPLETE\n")
  cat("Results saved to:", output_dir, "\n")
  cat(strrep("=", 60), "\n")

  invisible(results)
}

# ==============================================================================
# USAGE
# ==============================================================================

# Uncomment to run all simulations:
# results <- run_all_simulations(reduced = FALSE)

# Or run specific blocks:
# results <- run_all_simulations(reduced = FALSE,
#         run_block_a_flag = FALSE,
#         run_block_b_flag = FALSE,
#         run_block_c_flag = FALSE,
#         run_block_d_flag = TRUE,
#         run_block_e_flag = FALSE,
#         run_block_f_flag = FALSE
# )


cat("\n")
cat("=======================================================\n")
cat("SIMULATION FRAMEWORK v2 LOADED\n")
cat("=======================================================\n")
cat("\n")
cat("Key improvements from v1:\n")
cat("  - Uses TFisher package for correct null distribution\n")
cat("  - Precomputes null distributions (massive speedup)\n")
cat("  - Block B includes stress tests for broken component combinations\n")
cat("  - Block D implements adaptive omnibus with train/test split\n")
cat("  - Uses two-sided p-values for mixed direction scenarios\n")
cat("\n")
cat("To run simulations:\n")
cat("  results <- run_all_simulations(reduced = TRUE)   # Quick test\n")
cat("  results <- run_all_simulations(reduced = FALSE)  # Full run\n")
cat("\n")
cat("Results will be saved to 'simulation_results/' directory\n")
cat("=======================================================\n")
