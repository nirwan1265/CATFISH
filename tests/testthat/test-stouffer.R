# Tests for Stouffer-related functions

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
# Tests for magcat_stoufferZ_pathways
# ============================================================

test_that("magcat_stoufferZ_pathways returns correct columns", {
  gene_results <- mock_gene_results(10)
  pathways <- list(
    pwy1 = c("gene1", "gene2", "gene3"),
    pwy2 = c("gene4", "gene5")
  )

  result <- magcat_stoufferZ_pathways(
    gene_results = gene_results,
    pathways = pathways
  )

  expected_cols <- c("pathway_id", "pathway_name", "n_genes", "gene_names",
                     "gene_zvals", "stouffer_Z", "stouffer_p")
  expect_true(all(expected_cols %in% names(result)))
})

test_that("magcat_stoufferZ_pathways returns correct number of rows", {
  gene_results <- mock_gene_results(10)
  pathways <- list(
    pwy1 = c("gene1", "gene2"),
    pwy2 = c("gene3", "gene4"),
    pwy3 = c("gene5", "gene6", "gene7")
  )

  result <- magcat_stoufferZ_pathways(
    gene_results = gene_results,
    pathways = pathways
  )

  expect_equal(nrow(result), 3)
})

test_that("magcat_stoufferZ_pathways handles data.frame pathways", {
  gene_results <- mock_gene_results(10)
  pathways_df <- data.frame(
    pathway_id = c("pwy1", "pwy1", "pwy2", "pwy2"),
    gene_id = c("gene1", "gene2", "gene3", "gene4"),
    pathway_name = c("Pathway 1", "Pathway 1", "Pathway 2", "Pathway 2"),
    stringsAsFactors = FALSE
  )

  result <- magcat_stoufferZ_pathways(
    gene_results = gene_results,
    pathways = pathways_df
  )

  expect_equal(nrow(result), 2)
})

test_that("magcat_stoufferZ_pathways errors without pathways or species", {
  gene_results <- mock_gene_results(10)

  expect_error(
    magcat_stoufferZ_pathways(gene_results = gene_results),
    "must provide either"
  )
})

test_that("magcat_stoufferZ_pathways errors with both pathways and species", {
  gene_results <- mock_gene_results(10)
  pathways <- list(pwy1 = c("gene1", "gene2"))

  expect_error(
    magcat_stoufferZ_pathways(
      gene_results = gene_results,
      pathways = pathways,
      species = "maize"
    ),
    "not both"
  )
})

test_that("magcat_stoufferZ_pathways validates required columns", {
  gene_results <- data.frame(
    GENE = c("a", "b"),
    WRONG_Z = c(1.5, -0.5)  # Missing ZSTAT
  )
  pathways <- list(pwy1 = c("a", "b"))

  expect_error(
    magcat_stoufferZ_pathways(gene_results = gene_results, pathways = pathways),
    "must contain columns"
  )
})

test_that("magcat_stoufferZ_pathways p-values are in valid range", {
  set.seed(42)
  gene_results <- mock_gene_results(20)
  pathways <- list(
    pwy1 = paste0("gene", 1:5),
    pwy2 = paste0("gene", 6:10)
  )

  result <- magcat_stoufferZ_pathways(
    gene_results = gene_results,
    pathways = pathways
  )

  p_vals <- result$stouffer_p[!is.na(result$stouffer_p)]
  expect_true(all(p_vals >= 0 & p_vals <= 1))
})

test_that("magcat_stoufferZ_pathways is case-insensitive", {
  gene_results <- data.frame(
    GENE = c("GENE1", "GENE2", "GENE3"),
    ZSTAT = c(2.5, 1.5, -0.5),
    stringsAsFactors = FALSE
  )

  pathways <- list(pwy1 = c("gene1", "gene2"))

  result <- magcat_stoufferZ_pathways(
    gene_results = gene_results,
    pathways = pathways
  )

  expect_equal(result$n_genes[1], 2)
})

test_that("magcat_stoufferZ_pathways results are sorted by p-value", {
  set.seed(123)
  gene_results <- mock_gene_results(20)
  pathways <- list(
    pwy1 = paste0("gene", 1:5),
    pwy2 = paste0("gene", 6:10),
    pwy3 = paste0("gene", 11:15)
  )

  result <- magcat_stoufferZ_pathways(
    gene_results = gene_results,
    pathways = pathways
  )

  p_vals <- result$stouffer_p[!is.na(result$stouffer_p)]
  expect_true(all(diff(p_vals) >= 0))
})

test_that("magcat_stoufferZ_pathways supports different alternatives", {
  set.seed(42)
  gene_results <- mock_gene_results(10)
  pathways <- list(pwy1 = paste0("gene", 1:5))

  result_greater <- magcat_stoufferZ_pathways(
    gene_results = gene_results,
    pathways = pathways,
    alternative = "greater"
  )

  result_two <- magcat_stoufferZ_pathways(
    gene_results = gene_results,
    pathways = pathways,
    alternative = "two.sided"
  )

  result_less <- magcat_stoufferZ_pathways(
    gene_results = gene_results,
    pathways = pathways,
    alternative = "less"
  )

  # All should return valid p-values
  expect_true(result_greater$stouffer_p[1] >= 0 & result_greater$stouffer_p[1] <= 1)
  expect_true(result_two$stouffer_p[1] >= 0 & result_two$stouffer_p[1] <= 1)
  expect_true(result_less$stouffer_p[1] >= 0 & result_less$stouffer_p[1] <= 1)

  # Z scores should be identical
  expect_equal(result_greater$stouffer_Z[1], result_two$stouffer_Z[1])
  expect_equal(result_greater$stouffer_Z[1], result_less$stouffer_Z[1])
})

test_that("magcat_stoufferZ_pathways supports weighted analysis", {
  set.seed(42)
  gene_results <- data.frame(
    GENE = paste0("gene", 1:10),
    ZSTAT = rnorm(10),
    NSNPS = sample(10:100, 10),
    stringsAsFactors = FALSE
  )

  pathways <- list(pwy1 = paste0("gene", 1:5))

  result_unweighted <- magcat_stoufferZ_pathways(
    gene_results = gene_results,
    pathways = pathways,
    weight_col = NULL
  )

  result_weighted <- magcat_stoufferZ_pathways(
    gene_results = gene_results,
    pathways = pathways,
    weight_col = "NSNPS"
  )

  # Both should produce valid results
  expect_true(!is.na(result_unweighted$stouffer_Z[1]))
  expect_true(!is.na(result_weighted$stouffer_Z[1]))

  # Z scores may differ due to weighting
  # (unless by chance they're the same)
})
