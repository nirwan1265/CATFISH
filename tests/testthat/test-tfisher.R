# Tests for TFisher-related functions

# Helper: create mock gene results
mock_gene_results <- function(n = 10) {
  data.frame(
    GENE = paste0("gene", seq_len(n)),
    P = runif(n, 0.001, 0.999),
    ZSTAT = rnorm(n),
    NSNPS = sample(10:100, n, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

# ============================================================
# Tests for catfish_soft_tfisher_adaptive_pathways
# ============================================================

test_that("catfish_soft_tfisher_adaptive_pathways requires TFisher package", {
  skip_if_not_installed("TFisher")

  gene_results <- mock_gene_results(10)
  pathways <- list(pwy1 = c("gene1", "gene2", "gene3"))

  result <- catfish_soft_tfisher_adaptive_pathways(
    gene_results = gene_results,
    pathways = pathways
  )

  expect_s3_class(result, "data.frame")
})

test_that("catfish_soft_tfisher_adaptive_pathways returns correct columns", {
  skip_if_not_installed("TFisher")

  gene_results <- mock_gene_results(10)
  pathways <- list(pwy1 = c("gene1", "gene2", "gene3"))

  result <- catfish_soft_tfisher_adaptive_pathways(
    gene_results = gene_results,
    pathways = pathways
  )

  expected_cols <- c("pathway_id", "pathway_name", "n_genes", "gene_names",
                     "tau_hat", "tfisher_stat_hat", "tfisher_p_omni")
  expect_true(all(expected_cols %in% names(result)))
})

test_that("catfish_soft_tfisher_adaptive_pathways uses tau_grid", {
  skip_if_not_installed("TFisher")

  gene_results <- mock_gene_results(10)
  pathways <- list(pwy1 = c("gene1", "gene2", "gene3"))

  custom_tau <- c(0.1, 0.05, 0.01)

  result <- catfish_soft_tfisher_adaptive_pathways(
    gene_results = gene_results,
    pathways = pathways,
    tau_grid = custom_tau
  )

  # tau_hat should be one of the grid values
  expect_true(result$tau_hat[1] %in% custom_tau)
})

test_that("adaptive TFisher excludes tau = 1 consistently", {
  skip_if_not_installed("TFisher")

  gene_results <- data.frame(
    GENE = paste0("gene", 1:5),
    P = c(1e-4, 2e-3, 0.02, 0.15, 0.4),
    stringsAsFactors = FALSE
  )
  pathways <- list(pwy1 = gene_results$GENE)
  tau_with_fisher <- c(1, 0.1, 0.05, 0.01)

  result <- catfish_soft_tfisher_adaptive_pathways(
    gene_results = gene_results,
    pathways = pathways,
    tau_grid = tau_with_fisher
  )

  p_vec <- gene_results$P
  helper_p <- CATFISH:::.catfish_tfisher_adaptive_p(
    p = p_vec,
    tau_grid = tau_with_fisher
  )

  expect_false(isTRUE(result$tau_hat[1] == 1))
  expect_equal(result$tfisher_p_omni[1], helper_p, tolerance = 1e-12)
})

test_that("catfish_soft_tfisher_adaptive_pathways supports tau_grid = 'auto'", {
  skip_if_not_installed("TFisher")

  set.seed(1)
  gene_results <- mock_gene_results(200)
  pathways <- list(pwy1 = paste0("gene", 1:20), pwy2 = paste0("gene", 21:40))

  result <- catfish_soft_tfisher_adaptive_pathways(
    gene_results = gene_results,
    pathways = pathways,
    tau_grid = "auto"
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_true(all(result$tfisher_p_omni >= 0 & result$tfisher_p_omni <= 1, na.rm = TRUE))
  expect_true(all(result$tau_hat > 0, na.rm = TRUE))
})

test_that(".catfish_auto_tau_grid always returns adaptive taus below 1", {
  set.seed(11)
  g_auto <- CATFISH:::.catfish_auto_tau_grid(runif(500))
  expect_true(length(g_auto) >= 2)
  expect_true(all(g_auto > 0 & g_auto < 1))
  expect_equal(g_auto, sort(g_auto, decreasing = FALSE))
})

# ============================================================
# Tests for the data-driven (auto) tau grid selector
# ============================================================

test_that(".catfish_auto_tau_grid adapts to signal strength", {
  # under a uniform (null) p distribution, the quantile-based grid tracks the
  # classic thresholds (e.g. the 5% quantile of Uniform(0,1) is ~0.05)
  set.seed(2)
  p_null <- runif(5000)
  g_null <- CATFISH:::.catfish_auto_tau_grid(p_null)
  expect_true(length(g_null) >= 2)
  expect_true(all(g_null > 0 & g_null < 1))
  expect_true(any(abs(g_null - 0.05) < 0.03))

  # under strong signal (tiny p-values), the grid shrinks toward the extreme tail
  p_strong <- 10^(-runif(5000, 2, 12))   # p between 1e-12 and 1e-2
  g_strong <- CATFISH:::.catfish_auto_tau_grid(p_strong)
  expect_true(max(g_strong) < max(g_null))

  # too few genes -> classic fixed fallback grid
  g_small <- CATFISH:::.catfish_auto_tau_grid(runif(10))
  expect_equal(sort(g_small), sort(c(0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.001)))
})
