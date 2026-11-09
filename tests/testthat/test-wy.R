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

test_that("genomic control: shift removed, invariant to a genome-wide shift, idempotent", {
  skip_if_not_installed("TFisher")
  m <- wy_mock()
  ## keep every p representable in double precision so that shifted copies of
  ## the input are exact (the mock's z = 40 gene underflows to 0 otherwise)
  m$g$P <- pmax(m$g$P, 1e-200)
  run <- function(g) suppressMessages(catfish_wy_pathways(
    g, pathways = m$pw, cor_pairs = m$cp, B = 200, min_genes = 10, seed = 3,
    verbose = FALSE))
  r0 <- run(m$g)
  ## the removed shift is the median gene z over the genes analysed
  ok <- is.finite(m$g$P) & !duplicated(m$g$GENE)
  z  <- qnorm(log(m$g$P[ok]), lower.tail = FALSE, log.p = TRUE)
  expect_equal(attr(r0, "gc_shift"), median(z))
  expect_true(all(r0$gc_shift == attr(r0, "gc_shift")))
  ## shifting every gene z by the same constant changes nothing but gc_shift
  g1 <- m$g; g1$P <- pnorm(qnorm(m$g$P, lower.tail = FALSE) + 0.5, lower.tail = FALSE)
  r1 <- run(g1)
  expect_equal(r1$pathway_id, r0$pathway_id)
  expect_equal(r1$z_std, r0$z_std, tolerance = 1e-6)
  expect_equal(r1$p_fwer_std, r0$p_fwer_std)
  expect_equal(attr(r1, "gc_shift"), attr(r0, "gc_shift") + 0.5, tolerance = 1e-6)
  ## centred input: nothing left to remove
  g2 <- m$g; g2$P[ok] <- pnorm(z - median(z), lower.tail = FALSE)
  r2 <- run(g2)
  expect_lt(abs(attr(r2, "gc_shift")), 1e-8)
  expect_equal(r2$z_std, r0$z_std, tolerance = 1e-6)
})

test_that("gc = FALSE turns genomic control off and reports the median it left in", {
  skip_if_not_installed("TFisher")
  m <- wy_mock(); m$g$P <- pmax(m$g$P, 1e-200)
  run <- function(g, gc) suppressMessages(catfish_wy_pathways(
    g, pathways = m$pw, cor_pairs = m$cp, B = 200, min_genes = 10, seed = 3,
    gc = gc, verbose = FALSE))
  r_on <- run(m$g, TRUE)
  ok <- is.finite(m$g$P) & !duplicated(m$g$GENE)
  z  <- qnorm(log(m$g$P[ok]), lower.tail = FALSE, log.p = TRUE)
  ## off on centred input == on on raw input
  g2 <- m$g; g2$P[ok] <- pnorm(z - median(z), lower.tail = FALSE)
  r_off <- run(g2, FALSE)
  expect_equal(r_off$pathway_id, r_on$pathway_id)
  expect_equal(r_off$z_std, r_on$z_std, tolerance = 1e-6)
  expect_equal(attr(r_off, "gc_shift"), 0)
  expect_true(all(r_off$gc_shift == 0))
  ## off on raw input keeps the shift, and says so
  r_raw <- run(m$g, FALSE)
  expect_equal(attr(r_raw, "median_gene_z"), median(z))
  expect_equal(attr(r_raw, "gc_shift"), 0)
  expect_message(catfish_wy_pathways(m$g, pathways = m$pw, cor_pairs = m$cp, B = 50,
                                     min_genes = 10, seed = 3, gc = FALSE),
                 "genomic control off")
})

test_that("a large median gene z prints the competitive-test note", {
  skip_if_not_installed("TFisher")
  m <- wy_mock(); m$g$P <- pmax(m$g$P, 1e-200)
  g1 <- m$g; g1$P <- pnorm(qnorm(m$g$P, lower.tail = FALSE) + 1, lower.tail = FALSE)
  expect_message(catfish_wy_pathways(g1, pathways = m$pw, cor_pairs = m$cp, B = 50,
                                     min_genes = 10, seed = 3),
                 "closer to a competitive test")
  expect_false(grepl("note", .catfish_gc_message(0.1, TRUE, "WY")))
  expect_true(grepl("deflated", .catfish_gc_message(-0.8, TRUE, "WY")))
})

test_that("omni2 applies the same genomic control", {
  m <- wy_mock()
  om <- suppressWarnings(catfish_omni2_pathways(
    m$g, pathways = m$pw, gene_col = "GENE", p_raw_col = "P", min_genes = 10,
    perm_mode = "none", tau_grid = c(0.01, 0.05, 0.5)))
  ok <- is.finite(m$g$P) & !duplicated(m$g$GENE)
  expect_equal(attr(om, "gc_shift"),
               median(qnorm(log(m$g$P[ok]), lower.tail = FALSE, log.p = TRUE)))
  om0 <- suppressMessages(suppressWarnings(catfish_omni2_pathways(
    m$g, pathways = m$pw, gene_col = "GENE", p_raw_col = "P", min_genes = 10,
    perm_mode = "none", tau_grid = c(0.01, 0.05, 0.5), gc = FALSE)))
  expect_equal(attr(om0, "gc_shift"), 0)
  expect_equal(attr(om0, "median_gene_z"), attr(om, "gc_shift"))
})

test_that("competitive = FALSE leaves the output unchanged", {
  skip_if_not_installed("TFisher")
  m <- wy_mock()
  run <- function(...) suppressMessages(catfish_wy_pathways(
    m$g, pathways = m$pw, cor_pairs = m$cp, B = 100, min_genes = 10, seed = 3,
    verbose = FALSE, ...))
  r0 <- run(); r1 <- run(competitive = TRUE, B_comp = 50)
  expect_false(any(c("z_std_comp", "p_fwer_comp", "p_calibrated_comp") %in% names(r0)))
  expect_null(attr(r0, "B_comp"))
  ## the self-contained columns do not depend on the competitive run
  cols <- setdiff(names(r0), "pathway_name")
  expect_equal(r1[, cols], r0[, cols])
})

test_that("competitive test: valid p-values, replicates capped at genes - 1", {
  skip_if_not_installed("TFisher")
  m <- wy_mock()
  m$g$START <- ave(seq_len(nrow(m$g)), m$g$CHR, FUN = seq_along) * 1000
  r <- suppressMessages(catfish_wy_pathways(
    m$g, pathways = m$pw, cor_pairs = m$cp, B = 100, min_genes = 10, seed = 3,
    competitive = TRUE, pos_col = "START", B_comp = 5000, verbose = FALSE))
  ng <- sum(is.finite(m$g$P) & !duplicated(m$g$GENE))
  expect_equal(attr(r, "B_comp"), ng - 1L)
  for (cl in c("p_fwer_comp", "p_calibrated_comp")) {
    expect_true(all(r[[cl]] > 0 & r[[cl]] <= 1, na.rm = TRUE))
    expect_true(all(r[[cl]] >= 1 / ng))
  }
  ## family-wise p is never below the pathway's own calibrated p
  expect_true(all(r$p_fwer_comp >= r$p_calibrated_comp - 1e-12))
  expect_error(catfish_wy_pathways(m$g, pathways = m$pw, cor_pairs = m$cp, B = 10,
                                   competitive = TRUE, pos_col = "NOPE"),
               "missing column")
})

test_that("competitive test ignores a genome-wide shift that the self-contained test reads as signal", {
  skip_if_not_installed("TFisher")
  set.seed(11)
  NG <- 600
  g <- data.frame(GENE = paste0("g", 1:NG), CHR = rep(1:3, each = NG / 3),
                  START = rep(seq_len(NG / 3), 3) * 1000)
  z <- rnorm(NG) + 1                       # every gene shifted: polygenic background
  z[1:12] <- z[1:12] + 4                   # one real pathway on top of it
  g$P <- pnorm(z, lower.tail = FALSE)
  pw <- c(list(REAL = g$GENE[seq(1, 12)]),
          setNames(lapply(1:40, function(i) sample(g$GENE[-(1:12)], 12)), paste0("RND", 1:40)))
  r <- suppressWarnings(suppressMessages(catfish_wy_pathways(
    g, pathways = pw, B = 300, min_genes = 10, seed = 5, gc = FALSE,
    competitive = TRUE, pos_col = "START", B_comp = 300, verbose = FALSE)))
  rnd <- r[r$pathway_id != "REAL", ]
  ## self-contained, uncentred: the shift makes random pathways significant
  expect_gt(mean(rnd$p_fwer_std <= 0.05), 0.5)
  ## competitive: random pathways are null, the real one is found
  expect_lte(mean(rnd$p_calibrated_comp <= 0.05), 0.15)
  expect_lte(r$p_fwer_comp[r$pathway_id == "REAL"], 0.05)
})

test_that("circular shift keeps a pathway's genes together", {
  ## a pathway of adjacent genes gets, in every replicate, a block of adjacent
  ## scores from somewhere else in the genome
  lp <- log(seq(0.01, 0.99, length.out = 50))
  ord <- c(26:50, 1:25)                     # genome order differs from input order
  paths <- list(A = 30:34)
  pos <- integer(50); pos[ord] <- seq_len(50)
  lp_ord <- lp[ord]
  for (s in c(1L, 7L, 49L)) {
    lpn <- lp_ord[((pos - 1L + s) %% 50) + 1L]
    got <- match(lpn[paths$A], lp_ord)
    expect_true(all(diff(got) == 1L | diff(got) == -49L))
  }
})
