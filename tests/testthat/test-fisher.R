# Tests for Fisher-related functions

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
# Tests for magcat_fisher_pathways
# ============================================================

test_that("magcat_fisher_pathways returns correct columns", {
  gene_results <- mock_gene_results(10)
  pathways <- list(
    pwy1 = c("gene1", "gene2", "gene3"),
    pwy2 = c("gene4", "gene5")
  )

  result <- magcat_fisher_pathways(
    gene_results = gene_results,
    pathways = pathways
  )

  expected_cols <- c("pathway_id", "pathway_name", "n_genes", "gene_names", "fisher_p")
  expect_true(all(expected_cols %in% names(result)))
})

test_that("magcat_fisher_pathways returns correct number of rows", {
  gene_results <- mock_gene_results(10)
  pathways <- list(
    pwy1 = c("gene1", "gene2"),
    pwy2 = c("gene3", "gene4"),
    pwy3 = c("gene5", "gene6", "gene7")
  )

  result <- magcat_fisher_pathways(
    gene_results = gene_results,
    pathways = pathways
  )

  expect_equal(nrow(result), 3)
})

test_that("magcat_fisher_pathways handles data.frame pathways", {
  gene_results <- mock_gene_results(10)
  pathways_df <- data.frame(
    pathway_id = c("pwy1", "pwy1", "pwy2", "pwy2"),
    gene_id = c("gene1", "gene2", "gene3", "gene4"),
    pathway_name = c("Pathway 1", "Pathway 1", "Pathway 2", "Pathway 2"),
    stringsAsFactors = FALSE
  )

  result <- magcat_fisher_pathways(
    gene_results = gene_results,
    pathways = pathways_df
  )

  expect_equal(nrow(result), 2)
})

test_that("magcat_fisher_pathways errors without pathways or species", {
  gene_results <- mock_gene_results(10)

  expect_error(
    magcat_fisher_pathways(gene_results = gene_results),
    "must provide either"
  )
})

test_that("magcat_fisher_pathways errors with both pathways and species", {
  gene_results <- mock_gene_results(10)
  pathways <- list(pwy1 = c("gene1", "gene2"))

  expect_error(
    magcat_fisher_pathways(
      gene_results = gene_results,
      pathways = pathways,
      species = "maize"
    ),
    "not both"
  )
})

test_that("magcat_fisher_pathways validates required columns", {
  gene_results <- data.frame(WRONG_COL = c("a", "b"), P = c(0.1, 0.2))
  pathways <- list(pwy1 = c("a", "b"))

  expect_error(
    magcat_fisher_pathways(gene_results = gene_results, pathways = pathways),
    "must contain columns"
  )
})

test_that("magcat_fisher_pathways p-values are in valid range", {
  set.seed(42)
  gene_results <- mock_gene_results(20)
  pathways <- list(
    pwy1 = paste0("gene", 1:5),
    pwy2 = paste0("gene", 6:10)
  )

  result <- magcat_fisher_pathways(
    gene_results = gene_results,
    pathways = pathways
  )

  p_vals <- result$fisher_p[!is.na(result$fisher_p)]
  expect_true(all(p_vals >= 0 & p_vals <= 1))
})

test_that("magcat_fisher_pathways is case-insensitive", {
  gene_results <- data.frame(
    GENE = c("GENE1", "GENE2", "GENE3"),
    P = c(0.01, 0.05, 0.1),
    stringsAsFactors = FALSE
  )

  pathways <- list(pwy1 = c("gene1", "gene2"))

  result <- magcat_fisher_pathways(
    gene_results = gene_results,
    pathways = pathways
  )

  expect_equal(result$n_genes[1], 2)
})

test_that("magcat_fisher_pathways handles empty pathway matches", {
  gene_results <- mock_gene_results(5)
  pathways <- list(
    pwy1 = c("gene1", "gene2"),
    pwy_empty = c("nonexistent1", "nonexistent2")
  )

  result <- magcat_fisher_pathways(
    gene_results = gene_results,
    pathways = pathways
  )

  expect_equal(nrow(result), 2)
})

test_that("magcat_fisher_pathways results are sorted by p-value", {
  set.seed(123)
  gene_results <- mock_gene_results(20)
  pathways <- list(
    pwy1 = paste0("gene", 1:5),
    pwy2 = paste0("gene", 6:10),
    pwy3 = paste0("gene", 11:15)
  )

  result <- magcat_fisher_pathways(
    gene_results = gene_results,
    pathways = pathways
  )

  p_vals <- result$fisher_p[!is.na(result$fisher_p)]
  expect_true(all(diff(p_vals) >= 0))
})
