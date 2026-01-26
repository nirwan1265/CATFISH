# Tests for ACAT-related functions

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

# Helper: create mock pathways
mock_pathways <- function(gene_ids, n_pathways = 3, genes_per_pathway = 3) {
  pathways <- list()
  for (i in seq_len(n_pathways)) {
    pathway_genes <- sample(gene_ids, min(genes_per_pathway, length(gene_ids)))
    pathways[[paste0("pathway", i)]] <- pathway_genes
  }
  pathways
}

# ============================================================
# Tests for fix_p_for_acat
# ============================================================

test_that("fix_p_for_acat handles normal p-values", {
  p <- c(0.01, 0.05, 0.5, 0.95)
  result <- fix_p_for_acat(p)

expect_equal(length(result), 4)
  expect_true(all(result > 0))
  expect_true(all(result < 1))
})

test_that("fix_p_for_acat handles extreme values", {
  p <- c(0, 1e-20, 0.5, 1, 1.5)
  result <- fix_p_for_acat(p)

  # Should cap extremes
  expect_true(all(result > 0))
  expect_true(all(result < 1))
})

test_that("fix_p_for_acat removes NAs", {
  p <- c(0.1, NA, 0.5, NA, 0.9)
  result <- fix_p_for_acat(p)

  expect_equal(length(result), 3)
  expect_false(anyNA(result))
})

test_that("fix_p_for_acat handles empty input", {
  result <- fix_p_for_acat(numeric(0))
  expect_equal(length(result), 0)
})

test_that("fix_p_for_acat handles all-NA input", {
  result <- fix_p_for_acat(c(NA, NA, NA))
  expect_equal(length(result), 0)
})

# ============================================================
# Tests for magcat_load_pathways
# ============================================================

test_that("magcat_load_pathways validates species argument", {
  expect_error(
    magcat_load_pathways(species = "invalid_species"),
    "should be one of"
  )
})

test_that("magcat_load_pathways returns correct structure", {
  skip_if_not_installed("MAGCAT")

  # This will only work if package is installed with data
  result <- tryCatch(
    magcat_load_pathways(species = "maize"),
    error = function(e) NULL
  )

  if (!is.null(result)) {
    expect_s3_class(result, "data.frame")
    expect_true("pathway_id" %in% names(result))
    expect_true("pathway_name" %in% names(result))
    expect_true("gene_id" %in% names(result))
  }
})

# ============================================================
# Tests for magcat_acat_pathways
# ============================================================

test_that("magcat_acat_pathways requires ACAT package", {
  skip_if_not_installed("ACAT")

  gene_results <- mock_gene_results(10)
  pathways <- mock_pathways(gene_results$GENE)

  result <- magcat_acat_pathways(
    gene_results = gene_results,
    pathways = pathways
  )

  expect_s3_class(result, "data.frame")
})

test_that("magcat_acat_pathways returns correct columns", {
  skip_if_not_installed("ACAT")

  gene_results <- mock_gene_results(10)
  pathways <- mock_pathways(gene_results$GENE)

  result <- magcat_acat_pathways(
    gene_results = gene_results,
    pathways = pathways
  )

  expected_cols <- c("pathway_id", "pathway_name", "n_genes", "gene_names", "acat_p")
  expect_true(all(expected_cols %in% names(result)))
})

test_that("magcat_acat_pathways handles list pathways", {
  skip_if_not_installed("ACAT")

  gene_results <- mock_gene_results(10)
  pathways <- list(
    pwy1 = c("gene1", "gene2", "gene3"),
    pwy2 = c("gene4", "gene5"),
    pwy3 = c("gene6", "gene7", "gene8", "gene9")
  )

  result <- magcat_acat_pathways(
    gene_results = gene_results,
    pathways = pathways
  )

  expect_equal(nrow(result), 3)
})

test_that("magcat_acat_pathways handles data.frame pathways", {
  skip_if_not_installed("ACAT")

  gene_results <- mock_gene_results(10)
  pathways_df <- data.frame(
    pathway_id = c("pwy1", "pwy1", "pwy2", "pwy2", "pwy2"),
    gene_id = c("gene1", "gene2", "gene3", "gene4", "gene5"),
    pathway_name = c("Pathway 1", "Pathway 1", "Pathway 2", "Pathway 2", "Pathway 2"),
    stringsAsFactors = FALSE
  )

  result <- magcat_acat_pathways(
    gene_results = gene_results,
    pathways = pathways_df
  )

  expect_equal(nrow(result), 2)
})

test_that("magcat_acat_pathways errors without pathways or species", {
  gene_results <- mock_gene_results(10)

  expect_error(
    magcat_acat_pathways(gene_results = gene_results),
    "must provide either"
  )
})

test_that("magcat_acat_pathways errors with both pathways and species", {
  gene_results <- mock_gene_results(10)
  pathways <- list(pwy1 = c("gene1", "gene2"))

  expect_error(
    magcat_acat_pathways(
      gene_results = gene_results,
      pathways = pathways,
      species = "maize"
    ),
    "not both"
  )
})

test_that("magcat_acat_pathways validates required columns", {
  skip_if_not_installed("ACAT")

  gene_results <- data.frame(WRONG_COL = c("a", "b"), P = c(0.1, 0.2))
  pathways <- list(pwy1 = c("a", "b"))

  expect_error(
    magcat_acat_pathways(gene_results = gene_results, pathways = pathways),
    "missing required column"
  )
})

test_that("magcat_acat_pathways handles pathways with no matching genes", {
  skip_if_not_installed("ACAT")

  gene_results <- mock_gene_results(5)
  pathways <- list(
    pwy1 = c("gene1", "gene2"),
    pwy_no_match = c("nonexistent1", "nonexistent2")
  )

  result <- magcat_acat_pathways(
    gene_results = gene_results,
    pathways = pathways
  )

  # Should still return 2 rows, but one with NA p-value
  expect_equal(nrow(result), 2)
  expect_true(any(result$n_genes == 0) || any(is.na(result$acat_p)))
})

test_that("magcat_acat_pathways is case-insensitive for gene matching", {
  skip_if_not_installed("ACAT")

  gene_results <- data.frame(
    GENE = c("GENE1", "GENE2", "GENE3"),
    P = c(0.01, 0.05, 0.1),
    stringsAsFactors = FALSE
  )

  pathways <- list(pwy1 = c("gene1", "gene2"))  # lowercase

  result <- magcat_acat_pathways(
    gene_results = gene_results,
    pathways = pathways
  )

  expect_equal(result$n_genes[1], 2)
})

test_that("magcat_acat_pathways results are sorted by p-value", {
  skip_if_not_installed("ACAT")

  set.seed(123)
  gene_results <- data.frame(
    GENE = paste0("gene", 1:20),
    P = runif(20, 0.001, 0.999),
    stringsAsFactors = FALSE
  )

  pathways <- list(
    pwy1 = paste0("gene", 1:5),
    pwy2 = paste0("gene", 6:10),
    pwy3 = paste0("gene", 11:15),
    pwy4 = paste0("gene", 16:20)
  )

  result <- magcat_acat_pathways(
    gene_results = gene_results,
    pathways = pathways
  )

  # Remove NAs for comparison
  p_vals <- result$acat_p[!is.na(result$acat_p)]
  expect_true(all(diff(p_vals) >= 0))  # Should be sorted ascending
})
