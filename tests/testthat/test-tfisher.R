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
# Tests for magcat_tfisher_pathways (truncated Fisher / TPM)
# ============================================================

test_that("magcat_tfisher_pathways returns correct columns", {
  gene_results <- mock_gene_results(10)
  pathways <- list(
    pwy1 = c("gene1", "gene2", "gene3"),
    pwy2 = c("gene4", "gene5")
  )

  result <- magcat_tfisher_pathways(
    gene_results = gene_results,
    pathways = pathways,
    B_perm = 0  # Skip permutations for speed
  )

  expected_cols <- c("pathway_id", "pathway_name", "n_genes", "gene_names",
                     "tpm_stat", "nincl", "tpm_p_analytic")
  expect_true(all(expected_cols %in% names(result)))
})

test_that("magcat_tfisher_pathways returns correct number of rows", {
  gene_results <- mock_gene_results(10)
  pathways <- list(
    pwy1 = c("gene1", "gene2"),
    pwy2 = c("gene3", "gene4"),
    pwy3 = c("gene5", "gene6", "gene7")
  )

  result <- magcat_tfisher_pathways(
    gene_results = gene_results,
    pathways = pathways,
    B_perm = 0
  )

  expect_equal(nrow(result), 3)
})

test_that("magcat_tfisher_pathways errors without pathways or species", {
  gene_results <- mock_gene_results(10)

  expect_error(
    magcat_tfisher_pathways(gene_results = gene_results),
    "must provide either"
  )
})

test_that("magcat_tfisher_pathways errors with both pathways and species", {
  gene_results <- mock_gene_results(10)
  pathways <- list(pwy1 = c("gene1", "gene2"))

  expect_error(
    magcat_tfisher_pathways(
      gene_results = gene_results,
      pathways = pathways,
      species = "maize"
    ),
    "not both"
  )
})

test_that("magcat_tfisher_pathways respects ptrunc parameter", {
  gene_results <- data.frame(
    GENE = paste0("gene", 1:5),
    P = c(0.01, 0.03, 0.06, 0.1, 0.5),  # 2 below 0.05, 3 below 0.1
    stringsAsFactors = FALSE
  )

  pathways <- list(pwy1 = paste0("gene", 1:5))

  result_05 <- magcat_tfisher_pathways(
    gene_results = gene_results,
    pathways = pathways,
    ptrunc = 0.05,
    B_perm = 0
  )

  result_10 <- magcat_tfisher_pathways(
    gene_results = gene_results,
    pathways = pathways,
    ptrunc = 0.1,
    B_perm = 0
  )

  # nincl should differ based on ptrunc
  expect_equal(result_05$nincl[1], 2)
  expect_equal(result_10$nincl[1], 3)
})

test_that("magcat_tfisher_pathways permutations produce valid p-values", {
  set.seed(42)
  gene_results <- mock_gene_results(20)
  pathways <- list(pwy1 = paste0("gene", 1:5))

  result <- magcat_tfisher_pathways(
    gene_results = gene_results,
    pathways = pathways,
    B_perm = 100,  # Small number for speed
    seed = 123
  )

  expect_true(!is.na(result$tpm_p_perm[1]))
  expect_true(result$tpm_p_perm[1] >= 0 & result$tpm_p_perm[1] <= 1)
})

# ============================================================
# Tests for magcat_soft_tfisher_pathways
# ============================================================

test_that("magcat_soft_tfisher_pathways requires TFisher package", {
  skip_if_not_installed("TFisher")

  gene_results <- mock_gene_results(10)
  pathways <- list(pwy1 = c("gene1", "gene2", "gene3"))

  result <- magcat_soft_tfisher_pathways(
    gene_results = gene_results,
    pathways = pathways,
    B_perm = 0
  )

  expect_s3_class(result, "data.frame")
})

test_that("magcat_soft_tfisher_pathways returns correct columns", {
  skip_if_not_installed("TFisher")

  gene_results <- mock_gene_results(10)
  pathways <- list(pwy1 = c("gene1", "gene2", "gene3"))

  result <- magcat_soft_tfisher_pathways(
    gene_results = gene_results,
    pathways = pathways,
    B_perm = 0
  )

  expected_cols <- c("pathway_id", "pathway_name", "n_genes", "gene_names",
                     "tfisher_stat", "tfisher_p_analytic")
  expect_true(all(expected_cols %in% names(result)))
})

# ============================================================
# Tests for magcat_soft_tfisher_adaptive_pathways
# ============================================================

test_that("magcat_soft_tfisher_adaptive_pathways requires TFisher package", {
  skip_if_not_installed("TFisher")

  gene_results <- mock_gene_results(10)
  pathways <- list(pwy1 = c("gene1", "gene2", "gene3"))

  result <- magcat_soft_tfisher_adaptive_pathways(
    gene_results = gene_results,
    pathways = pathways
  )

  expect_s3_class(result, "data.frame")
})

test_that("magcat_soft_tfisher_adaptive_pathways returns correct columns", {
  skip_if_not_installed("TFisher")

  gene_results <- mock_gene_results(10)
  pathways <- list(pwy1 = c("gene1", "gene2", "gene3"))

  result <- magcat_soft_tfisher_adaptive_pathways(
    gene_results = gene_results,
    pathways = pathways
  )

  expected_cols <- c("pathway_id", "pathway_name", "n_genes", "gene_names",
                     "tau_hat", "tfisher_stat_hat", "tfisher_p_omni")
  expect_true(all(expected_cols %in% names(result)))
})

test_that("magcat_soft_tfisher_adaptive_pathways uses tau_grid", {
  skip_if_not_installed("TFisher")

  gene_results <- mock_gene_results(10)
  pathways <- list(pwy1 = c("gene1", "gene2", "gene3"))

  custom_tau <- c(0.1, 0.05, 0.01)

  result <- magcat_soft_tfisher_adaptive_pathways(
    gene_results = gene_results,
    pathways = pathways,
    tau_grid = custom_tau
  )

  # tau_hat should be one of the grid values
  expect_true(result$tau_hat[1] %in% custom_tau)
})
