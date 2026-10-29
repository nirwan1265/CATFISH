# Tests for the Westfall-Young engine (R/westfall_young.R) and its use as the
# omnibus calibration in catfish_omni2_pathways(perm_mode = "wy").

# Small synthetic GWAS: 3 chromosomes, banded correlation between neighbouring
# genes, a few associated genes (one far below 1e-300), overlapping pathways,
# a duplicated gene row and a missing p-value.
wy_mock <- function(seed = 99) {
  set.seed(seed)
  NG <- 300
  chr  <- rep(c("1", "2", "3"), each = NG / 3)
  gene <- as.character(100000 + seq_len(NG))
  pairs <- list()
  for (c0 in unique(chr)) {
    ii <- which(chr == c0)
    for (k in 1:3) {
      a <- ii[1:(length(ii) - k)]
      pairs[[length(pairs) + 1]] <- data.frame(gene1 = gene[a], gene2 = gene[a + k],
                                               r = 0.45 * 0.6^(k - 1))
    }
  }
  cp <- do.call(rbind, pairs)
  R <- diag(NG)
  R[cbind(match(cp$gene1, gene), match(cp$gene2, gene))] <- cp$r
  R[cbind(match(cp$gene2, gene), match(cp$gene1, gene))] <- cp$r
  z <- as.vector(t(chol(R)) %*% rnorm(NG))
  z[c(5, 6, 110)] <- z[c(5, 6, 110)] + c(40, 8, 6)
  g <- data.frame(GENE = gene, CHR = chr, P = pnorm(z, lower.tail = FALSE),
                  ZSTAT = z, stringsAsFactors = FALSE)
  g <- rbind(g, g[17, ])
  g$P[33] <- NA
  pw <- lapply(1:30, function(i) sample(c(gene, "999999999"), sample(10:30, 1)))
  names(pw) <- sprintf("PW_%02d", 1:30)
  list(g = g, cp = cp, pw = pw)
}

test_that("log-space ACAT matches ACAT::ACAT and handles extreme / negative cases", {
  skip_if_not_installed("ACAT")
  p <- c(0.2, 0.03, 0.6, 0.9)
  expect_equal(exp(.wy_acat_log(log(p))), ACAT::ACAT(p), tolerance = 1e-10)
  # p far below double precision stays finite and ordered
  a <- .wy_acat_log(c(-900, log(0.5)))
  b <- .wy_acat_log(c(-800, log(0.5)))
  expect_true(is.finite(a) && is.finite(b) && a < b)
  # all genes near 1: T is negative; result is a valid p close to 1, not 0.5
  q <- exp(.wy_acat_log(log(c(0.99, 0.995, 0.98))))
  expect_true(q > 0.5 && q <= 1)
})

test_that("log-space minP keeps resolution when the Sidak p underflows", {
  m <- 50
  lp <- c(log(1e-12), rep(log(0.5), m - 1))
  expect_equal(.wy_minp_log(lp), log(m) + log(1e-12), tolerance = 1e-6)
  expect_true(is.finite(.wy_minp_log(c(-2000, rep(log(0.5), 9)))))
})

test_that("TFisher component follows the drop rule and the 1e-16 cap", {
  skip_if_not_installed("TFisher")
  tg <- c(0.01, 0.05, 0.5)
  # no gene below any threshold -> undefined (NA), so it leaves the omnibus
  expect_true(is.na(.wy_tfisher_log(log(c(0.8, 0.9, 0.95)), tg)))
  # very strong signal is capped at 1e-16
  expect_equal(.wy_tfisher_log(c(-500, log(0.4), log(0.7)), tg), log(1e-16))
})

test_that("catfish_wy_pathways returns valid, internally consistent p-values", {
  d <- wy_mock()
  res <- suppressMessages(catfish_wy_pathways(
    d$g, pathways = d$pw, cor_pairs = d$cp, B = 60L, seed = 1L,
    min_genes = 10L, max_genes = 200L, verbose = FALSE))
  expect_true(nrow(res) > 0)
  B <- 60L
  for (cl in c("p_calibrated", "p_fwer_wy", "p_fwer_std")) {
    v <- res[[cl]][is.finite(res[[cl]])]
    expect_true(all(v >= 1 / (B + 1) - 1e-12 & v <= 1))
  }
  # family-wise p is never smaller than the pathway's own calibrated p
  expect_true(all(res$p_fwer_wy >= res$p_calibrated - 1e-12))
  expect_equal(attr(res, "B"), B)
  expect_true(is.finite(attr(res, "wy_threshold_05")))
  # same seed -> identical result
  res2 <- suppressMessages(catfish_wy_pathways(
    d$g, pathways = d$pw, cor_pairs = d$cp, B = 60L, seed = 1L,
    min_genes = 10L, max_genes = 200L, verbose = FALSE))
  expect_identical(res, res2)
})

test_that("catfish_wy_pathways rejects a data-adaptive tau grid", {
  d <- wy_mock()
  expect_error(catfish_wy_pathways(d$g, pathways = d$pw, cor_pairs = d$cp,
                                   B = 10L, tau_grid = "auto", verbose = FALSE),
               "fixed numeric grid")
})

test_that("omni2 perm_mode = 'wy' reproduces catfish_wy_pathways exactly", {
  d <- wy_mock()
  sa <- suppressMessages(catfish_wy_pathways(
    d$g, pathways = d$pw, cor_pairs = d$cp, B = 60L, seed = 7L,
    min_genes = 10L, max_genes = 200L, verbose = FALSE))
  om <- suppressWarnings(suppressMessages(catfish_omni2_pathways(
    d$g, pathways = d$pw, gene_col = "GENE", p_raw_col = "P", chr_col = "CHR",
    min_genes = 10L, max_genes = 200L, perm_mode = "wy", B_perm = 60L,
    magma_cor_pairs = d$cp, seed = 7L)))
  m <- match(sa$pathway_id, om$pathway_id)
  expect_false(anyNA(m))
  expect_equal(om$omni_p_cal[m],         sa$p_calibrated)
  expect_equal(om$omni_p_fwer_wy[m],     sa$p_fwer_wy)
  expect_equal(om$omni_p_fwer_wy_std[m], sa$p_fwer_std)
  expect_equal(om$omni_z_std[m],         sa$z_std)
  expect_equal(om$omni_p_wy_stat[m],     sa$p_omnibus)
  expect_equal(om$omni_log10_p_wy_stat[m], sa$log10_p_omnibus)
  # the reported p is the WY-calibrated one
  expect_equal(om$omni_p_final, om$omni_p_cal)
  expect_true(all(om$calib_mode == "wy"))
  expect_equal(attr(om, "wy_threshold_05"), attr(sa, "wy_threshold_05"))
})

test_that("omni2 Westfall-Young needs chromosomes and no longer offers MVN", {
  d <- wy_mock()
  expect_error(catfish_omni2_pathways(d$g, pathways = d$pw, perm_mode = "wy",
                                      B_perm = 10L, magma_cor_pairs = d$cp),
               "chr_col")
  expect_error(catfish_omni2_pathways(d$g, pathways = d$pw, perm_mode = "mvn",
                                      B_perm = 10L, magma_cor_pairs = d$cp))
})

test_that("omni2 analytic path is unchanged by the Westfall-Young columns", {
  d <- wy_mock()
  a <- suppressWarnings(catfish_omni2_pathways(d$g, pathways = d$pw, min_genes = 10L))
  expect_true(all(is.na(a$omni_p_cal)))
  expect_equal(a$omni_p_final, a$omni_p_analytic)
  expect_true(all(a$calib_mode == "analytic"))
})
