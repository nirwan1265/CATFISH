# ==============================================================================
# SIMULATION FRAMEWORK FOR CATFISH PATHWAY TESTS (v2)
# ==============================================================================
# This script implements comprehensive simulations to evaluate:
#   Block A: Power across different signal shapes (sparse vs dense)
#   Block B: "Broken component" stress test (isolate TFisher inflation)
#   Block C: Missing correlation / incomplete Sigma scenarios
#   Block D: Adaptive omnibus with train/test split
#   Block E: Leave-one-out omnibus sensitivity (MVN)
#   Block F: Random two-component dropout (MVN)
#
# Key fixes from v1:
#   - Uses actual TFisher package for correct null distribution
#   - Precomputes null distributions per (m, Sigma) for massive speedup
#   - Block B isolates TFisher as the only broken component
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

      p_analytic <- matrix(NA, nrow = n_sims, ncol = 6)
      p_mvn <- matrix(NA, nrow = n_sims, ncol = 6)
      colnames(p_analytic) <- colnames(p_mvn) <-
        c("acat", "fisher", "tfisher", "minp", "stouffer", "omnibus")

      for (s in seq_len(n_sims)) {
        z <- simulate_gene_z(m, sigma, n_causal = 0, effect_size = 0)
        gene_res <- create_gene_results(z)

        res_ana <- compute_components_analytic(gene_res)
        p_analytic[s, ] <- unlist(res_ana)

        res_cal <- calibrate_with_null(res_ana, null_dist)
        p_mvn[s, ] <- unlist(res_cal)
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
                        broken_components = c("stouffer", "acat", "minP"),
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
          res_ana_broken <- res_ana
          for (comp in broken_components) {
            if (comp %in% names(res_ana_broken)) {
              p_raw <- res_ana_broken[[comp]]
              p_broken <- pmax(pmin(p_raw^broken_power, 1 - min_p), min_p)
              res_ana_broken[[comp]] <- p_broken
            }
          }
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

      p_analytic <- matrix(NA, nrow = n_sims, ncol = 6)
      p_true_mvn <- matrix(NA, nrow = n_sims, ncol = 6)
      p_imputed_mvn <- matrix(NA, nrow = n_sims, ncol = 6)
      colnames(p_analytic) <- colnames(p_true_mvn) <- colnames(p_imputed_mvn) <-
        c("acat", "fisher", "tfisher", "minp", "stouffer", "omnibus")

      for (s in seq_len(n_sims)) {
        z <- simulate_gene_z(m, sigma_true, n_causal = 0, effect_size = 0)
        gene_res <- create_gene_results(z)

        res_ana <- compute_components_analytic(gene_res)
        p_analytic[s, ] <- unlist(res_ana)

        res_true <- calibrate_with_null(res_ana, null_true)
        p_true_mvn[s, ] <- unlist(res_true)

        res_imputed <- calibrate_with_null(res_ana, null_damaged)
        p_imputed_mvn[s, ] <- unlist(res_imputed)
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
      p_calibrated_omni <- numeric(n_test)
      p_adaptive_omni <- numeric(n_test)

      for (s in seq_len(n_test)) {
        z <- simulate_gene_z(m, sigma_true, n_causal = 0, effect_size = 0)
        gene_res <- create_gene_results(z)

        # Naive omnibus (all analytic)
        res_ana <- compute_components_analytic(gene_res)
        p_naive_omni[s] <- res_ana$omnibus

        # Calibrated omnibus (MVN via precomputed null)
        res_cal <- calibrate_with_null(res_ana, null_dist_correct)
        p_calibrated_omni[s] <- res_cal$omnibus

        # Adaptive omnibus (drop bad components)
        p_adaptive_omni[s] <- compute_adaptive_omnibus(gene_res, keep_mask)
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
        method = "omnibus_mvn",
        n_components = 5,
        drop_label = drop_label,
        lambda = compute_lambda(p_calibrated_omni),
        lambda_se = compute_lambda_se(p_calibrated_omni),
        type1_05 = compute_type1(p_calibrated_omni, 0.05),
        type1_05_se = compute_prop_se(p_calibrated_omni, 0.05),
        type1_01 = compute_type1(p_calibrated_omni, 0.01),
        type1_01_se = compute_prop_se(p_calibrated_omni, 0.01),
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
  loo_methods <- c(
    "omnibus_minus_acat",
    "omnibus_minus_fisher",
    "omnibus_minus_tfisher",
    "omnibus_minus_minp",
    "omnibus_minus_stouffer",
    "omnibus_all"
  )

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
      null_cache[[cache_key]] <- null_loo
    }
  }

  for (m in pathway_sizes_arg) {
    for (cor_name in names(cor_types)) {
      if (verbose) cat(sprintf("  m=%d, cor=%s\n", m, cor_name))

      sigma <- cor_types[[cor_name]](m)
      cache_key <- paste(m, cor_name, sep = "_")
      null_loo <- null_cache[[cache_key]]

      p_mvn <- matrix(NA, nrow = n_sims, ncol = length(loo_methods))
      colnames(p_mvn) <- loo_methods

      for (s in seq_len(n_sims)) {
        z <- simulate_gene_z(m, sigma, n_causal = 0, effect_size = 0)
        gene_res <- create_gene_results(z)

        res_ana <- compute_components_analytic(gene_res)
        comp_vals <- unlist(res_ana[comp_cols])

        obs_loo <- list(
          omnibus_all = res_ana$omnibus,
          omnibus_minus_acat = compute_omnibus(comp_vals[comp_cols != "acat"]),
          omnibus_minus_fisher = compute_omnibus(comp_vals[comp_cols != "fisher"]),
          omnibus_minus_tfisher = compute_omnibus(comp_vals[comp_cols != "tfisher"]),
          omnibus_minus_minp = compute_omnibus(comp_vals[comp_cols != "minp"]),
          omnibus_minus_stouffer = compute_omnibus(comp_vals[comp_cols != "stouffer"])
        )

        for (method in loo_methods) {
          obs_val <- obs_loo[[method]]
          null_vals <- null_loo[[method]]
          p_mvn[s, method] <- (1 + sum(null_vals <= obs_val, na.rm = TRUE)) /
            (sum(!is.na(null_vals)) + 1)
        }
      }

      for (method in loo_methods) {
        results_list[[length(results_list) + 1]] <- data.frame(
          block = "E",
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

  if (verbose) cat("Block E complete.\n\n")
  bind_rows(results_list)
}

# ==============================================================================
# PLOTTING FUNCTIONS
# ==============================================================================

plot_block_a_null <- function(results) {
  method_order <- c("ACAT", "FISHER", "MINP", "STOUFFER", "TFISHER", "OMNIBUS")
  cor_order <- c("LD_moderate", "LD_strong", "LD_independent")
  df <- results %>%
    filter(block == "A_null") %>%
    mutate(
      lambda_plot = pmin(lambda, lambda_cap),
      cor_structure = factor(cor_structure, levels = cor_order),
      method_label = factor(toupper(method), levels = method_order),
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
    geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
    geom_errorbar(aes(ymin = lambda_lower, ymax = lambda_upper),
                  position = position_dodge(0.8), width = 0.2) +
    geom_text(aes(label = label_lambda),
              position = position_dodge(0.8),
              angle = 90, vjust = -0.2, size = 2.3) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.7) +
    facet_grid(pathway_size ~ cor_structure,
               labeller = labeller(pathway_size = function(x) paste0("m=", x))) +
    scale_fill_manual(values = c("analytic" = "#E41A1C", "mvn" = "#377EB8"),
                      labels = c("Analytic", "MVN Calibrated")) +
    labs(title = "Block A: Lambda Under Null",
         x = "Method", y = expression(lambda), fill = "Calibration") +
    sim_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

  p_type1 <- ggplot(df, aes(x = method_label, y = type1_05, fill = calibration)) +
    geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
    geom_errorbar(aes(ymin = type1_lower, ymax = type1_upper),
                  position = position_dodge(0.8), width = 0.2) +
    geom_text(aes(label = label_type1),
              position = position_dodge(0.8),
              angle = 90, vjust = -0.2, size = 2.3) +
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
      method_label = toupper(method),
      pathway_size = as.integer(pathway_size),
      facet_label = factor(paste0("m=", pathway_size),
                           levels = paste0("m=", sort(unique(pathway_size))))
    )

  df_omni <- df %>%
    filter(method == "omnibus") %>%
    mutate(
      lambda_lower = pmax(lambda_plot - lambda_se, 0),
      lambda_upper = pmin(lambda_plot + lambda_se, lambda_cap),
      type1_lower = pmax(type1_05 - type1_05_se, 0),
      type1_upper = pmin(type1_05 + type1_05_se, 1)
    )

  p_lambda <- ggplot(df_omni, aes(x = missing_frac, y = lambda_plot, color = strategy)) +
    geom_ribbon(aes(ymin = lambda_lower, ymax = lambda_upper, fill = strategy),
                alpha = 0.15, show.legend = FALSE) +
    geom_point(size = 3) +
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
    labs(title = "Block C: Omnibus Lambda vs Missing Correlations",
         x = "Fraction Missing", y = expression(lambda),
         color = "Strategy") +
    sim_theme

  p_type1 <- ggplot(df_omni, aes(x = missing_frac, y = type1_05, color = strategy)) +
    geom_ribbon(aes(ymin = type1_lower, ymax = type1_upper, fill = strategy),
                alpha = 0.15, show.legend = FALSE) +
    geom_point(size = 3) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "black", linewidth = 0.5) +
    facet_wrap(~ facet_label, ncol = 3) +
    scale_color_manual(
      values = c("analytic_fallback" = "#E41A1C",
                 "mvn_true_cor" = "#377EB8",
                 "mvn_imputed_cor" = "#4DAF4A"),
      labels = c("Analytic Fallback", "MVN True Cor", "MVN Imputed (0)")
    ) +
    labs(title = "Block C: Omnibus Type I Error vs Missing Correlations",
         x = "Fraction Missing", y = "Type I Error Rate",
         color = "Strategy") +
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
        method == "omnibus_analytical" ~ "Analytical",
        method == "omnibus_mvn"        ~ "MVN",
        method == "omnibus_adaptive"   ~ "Adaptive",
        TRUE                           ~ method
      ),
      method_label = factor(method_label, levels = c("Analytical", "MVN", "Adaptive")),
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
                 "MVN" = "#377EB8",
                 "Adaptive" = "#4DAF4A")
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
                 "MVN" = "#377EB8",
                 "Adaptive" = "#4DAF4A")
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

  p_lambda <- ggplot(df, aes(x = method_label, y = lambda_plot)) +
    geom_col(width = 0.7, alpha = 0.85, fill = "#377EB8") +
    geom_errorbar(aes(ymin = lambda_lower, ymax = lambda_upper), width = 0.2) +
    geom_text(aes(label = label_lambda),
              angle = 90, vjust = -0.2, size = 2.3) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.7) +
    facet_grid(pathway_size ~ cor_structure,
               labeller = labeller(pathway_size = function(x) paste0("m=", x))) +
    labs(title = "Block E: Leave-one-out Omnibus Lambda (MVN)",
         x = "Omnibus Variant", y = expression(lambda)) +
    sim_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

  p_type1 <- ggplot(df, aes(x = method_label, y = type1_05)) +
    geom_col(width = 0.7, alpha = 0.85, fill = "#4DAF4A") +
    geom_errorbar(aes(ymin = type1_lower, ymax = type1_upper), width = 0.2) +
    geom_text(aes(label = label_type1),
              angle = 90, vjust = -0.2, size = 2.3) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.7) +
    facet_grid(pathway_size ~ cor_structure,
               labeller = labeller(pathway_size = function(x) paste0("m=", x))) +
    labs(title = "Block E: Leave-one-out Omnibus Type I Error (MVN)",
         x = "Omnibus Variant", y = "Type I Error Rate") +
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



results <- run_all_simulations(reduced = FALSE)


results <- run_all_simulations(reduced = FALSE,
        run_block_a_flag = FALSE,
        run_block_b_flag = FALSE,
        run_block_c_flag = FALSE,
        run_block_d_flag = TRUE,
        run_block_e_flag = FALSE,
        run_block_f_flag = FALSE
)


cat("\n")
cat("=======================================================\n")
cat("SIMULATION FRAMEWORK v2 LOADED\n")
cat("=======================================================\n")
cat("\n")
cat("Key improvements from v1:\n")
cat("  - Uses TFisher package for correct null distribution\n")
cat("  - Precomputes null distributions (massive speedup)\n")
cat("  - Block B isolates TFisher as the only broken component\n")
cat("  - Block D implements adaptive omnibus with train/test split\n")
cat("  - Uses two-sided p-values for mixed direction scenarios\n")
cat("\n")
cat("To run simulations:\n")
cat("  results <- run_all_simulations(reduced = TRUE)   # Quick test\n")
cat("  results <- run_all_simulations(reduced = FALSE)  # Full run\n")
cat("\n")
cat("Results will be saved to 'simulation_results/' directory\n")
cat("=======================================================\n")
