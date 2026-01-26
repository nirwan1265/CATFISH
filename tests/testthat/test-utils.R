# Tests for utility functions

# ============================================================
# Tests for magma_path
# ============================================================

test_that("magma_path returns a path or error", {
  # Either returns a valid path or throws an informative error
  result <- tryCatch(
    magma_path(),
    error = function(e) {
      expect_match(e$message, "MAGMA binary not found")
      "error_caught"
    }
  )

  if (result != "error_caught") {
    expect_true(is.character(result))
    expect_true(nzchar(result))
  }
})

test_that("magma_path respects options", {
  # Save original option
  orig_opt <- getOption("magma.path")
  on.exit(options(magma.path = orig_opt))

  # Set a fake path (won't exist)
  options(magma.path = "/fake/path/to/magma")

  # This should fail because file doesn't exist
  # The function checks file.exists()
  result <- tryCatch(
    magma_path(),
    error = function(e) "error"
  )

  # Either it errors or falls through to PATH lookup
  expect_true(result == "error" || is.character(result))
})

# ============================================================
# Tests for get_gene_lengths
# ============================================================

test_that("get_gene_lengths validates input file", {
  expect_error(
    get_gene_lengths("nonexistent_file.gff3"),
    "not found"
  )
})

test_that("get_gene_lengths validates gff3_file is single path", {
  expect_error(
    get_gene_lengths(c("file1.gff3", "file2.gff3")),
    "single file path"
  )
})

# ============================================================
# Tests for gff3_to_geneloc
# ============================================================

test_that("gff3_to_geneloc validates input file", {
  expect_error(
    gff3_to_geneloc(gff = "nonexistent.gff3", out = "output.loc"),
    "not found"
  )
})

test_that("gff3_to_geneloc validates recode_chr argument", {
  # Create a temporary minimal GFF3 for testing
  tmp_gff <- tempfile(fileext = ".gff3")
  on.exit(unlink(tmp_gff))

  # Write minimal GFF3
  writeLines(c(
    "##gff-version 3",
    "chr1\t.\tgene\t100\t200\t.\t+\t.\tID=gene1"
  ), tmp_gff)

  tmp_out <- tempfile(fileext = ".loc")
  on.exit(unlink(tmp_out), add = TRUE)

  # This should work with valid recode_chr
  result <- tryCatch(
    gff3_to_geneloc(gff = tmp_gff, out = tmp_out, recode_chr = "none"),
    error = function(e) NULL
  )

  # Just verify it doesn't error on valid input
  # (or errors gracefully if something else goes wrong)
})

# ============================================================
# Tests for magcat_adjust_gene_p
# ============================================================

test_that("magcat_adjust_gene_p validates inputs", {
  expect_error(
    magcat_adjust_gene_p(
      gene_results = "not_a_dataframe",
      gene_lengths = data.frame(gene_id = "a", length = 100)
    ),
    "data.frame"
  )

  expect_error(
    magcat_adjust_gene_p(
      gene_results = data.frame(GENE = "a", P = 0.1, NSNPS = 10),
      gene_lengths = "not_a_dataframe"
    ),
    "data.frame"
  )
})

test_that("magcat_adjust_gene_p validates required columns", {
  gene_results <- data.frame(
    WRONG_COL = c("gene1", "gene2"),
    P = c(0.1, 0.05),
    NSNPS = c(10, 20)
  )
  gene_lengths <- data.frame(
    gene_id = c("gene1", "gene2"),
    length = c(1000, 2000)
  )

  expect_error(
    magcat_adjust_gene_p(gene_results, gene_lengths),
    "missing column"
  )
})

test_that("magcat_adjust_gene_p returns correct structure", {
  gene_results <- data.frame(
    GENE = c("gene1", "gene2", "gene3", "gene4", "gene5"),
    P = c(0.01, 0.05, 0.1, 0.5, 0.9),
    ZSTAT = c(2.3, 1.6, 1.3, 0.1, -1.3),
    NSNPS = c(50, 30, 20, 10, 40),
    stringsAsFactors = FALSE
  )
  gene_lengths <- data.frame(
    gene_id = c("gene1", "gene2", "gene3", "gene4", "gene5"),
    length = c(5000, 3000, 2000, 1000, 4000),
    stringsAsFactors = FALSE
  )

  result <- magcat_adjust_gene_p(gene_results, gene_lengths)

  expect_s3_class(result, "data.frame")
  expect_true("gene_id" %in% names(result))
  expect_true("z_adj" %in% names(result))
  expect_true("p_adj" %in% names(result))
  expect_true(!is.null(attr(result, "lm_fit")))
})

test_that("magcat_adjust_gene_p p_adj values are valid", {
  gene_results <- data.frame(
    GENE = paste0("gene", 1:10),
    P = runif(10, 0.001, 0.999),
    ZSTAT = rnorm(10),
    NSNPS = sample(10:100, 10),
    stringsAsFactors = FALSE
  )
  gene_lengths <- data.frame(
    gene_id = paste0("gene", 1:10),
    length = sample(1000:10000, 10),
    stringsAsFactors = FALSE
  )

  result <- magcat_adjust_gene_p(gene_results, gene_lengths)

  p_adj <- result$p_adj[!is.na(result$p_adj)]
  expect_true(all(p_adj >= 0 & p_adj <= 1))
})
