## ============================================================
##  MAGCAT omnibus (fresh, fast, 6-method, with global + MVN resampling)
##
##  Methods (component p-values per pathway):
##    1) ACAT (gene-level)
##    2) Fisher (gene-level)
##    3) Adaptive soft TFisher (analytic p.soft across tau_grid)
##    4) minP (Sidak across genes)
##    5) Stouffer Z (from gene-level Z)
##    6) MAGMA competitive (optional; you pass the data.frame output "out")
##
##  Omnibus across METHODS (per pathway): ACAT or minP (Sidak across methods)
##
##  Permutation / resampling (omnibus only):
##    - global_resampling: resample genes from a global pool (fast heuristic)
##    - mvn_resampling: simulate correlated Z ~ MVN(0, R_S) using gene-gene correlations
##      (requires magma_cor_file or magma_cor_pairs)
##
#' MAGCAT omnibus pathway p-values (6-method framework; optional global + MVN calibration)
#'
#' Computes pathway-level p-values from gene-level statistics using up to six component
#' methods and combines them into a single omnibus p-value per pathway. Optionally,
#' calibrates the omnibus p-values using either (i) global gene resampling or (ii)
#' LD-aware MVN resampling based on MAGMA gene–gene correlations.
#'
#' @details
#' \strong{Component methods (per pathway):}
#' \enumerate{
#'   \item \strong{ACAT (gene-level):} ACAT combine of gene p-values in the pathway.
#'   \item \strong{Fisher (gene-level):} Fisher's method on gene p-values in the pathway.
#'   \item \strong{Adaptive soft TFisher (gene-level):} For each \code{tau} in \code{tau_grid},
#'         computes \code{TFisher::stat.soft()} and its analytic p-value via \code{TFisher::p.soft()},
#'         then takes the minimum p-value over \code{tau_grid}.
#'   \item \strong{minP (gene-level):} Sidak-adjusted minimum gene p-value in the pathway.
#'   \item \strong{Stouffer Z (gene-level):} Stouffer combination of gene Z-scores (optionally weighted).
#'   \item \strong{MAGMA competitive (optional):} Uses your precomputed MAGMA competitive
#'         gene-set p-value (provided via \code{magma_out}).
#' }
#'
#' \strong{Omnibus across methods (per pathway):}
#' Combines the component method p-values using either:
#' \itemize{
#'   \item \code{omnibus="ACAT"}: ACAT across method p-values, or
#'   \item \code{omnibus="minP"}: Sidak-adjusted minP across method p-values.
#' }
#'
#' \strong{Resampling / calibration (omnibus only):}
#' \itemize{
#'   \item \code{perm_mode="global"}: resamples gene p-values (and optionally Z's) from a global pool
#'         to build a null distribution of the omnibus statistic.
#'   \item \code{perm_mode="mvn"}: simulates correlated gene Z's
#'         \eqn{Z \sim \mathcal{N}(0, R_S)} using a pathway-specific correlation matrix \eqn{R_S}
#'         built from MAGMA gene–gene correlations. P-values for p-based methods are derived from the
#'         simulated Z's.
#' }
#'
#' \strong{MVN p-value marginals:}
#' When \code{perm_mode} includes \code{"mvn"}, p-values for p-based component methods are obtained from
#' the simulated MVN Z's in one of two ways:
#' \itemize{
#'   \item \code{mvn_marginal="uniform"} (recommended): Gaussian-copula mapping
#'         \eqn{U=\Phi(Z)}, \eqn{p=2\min(U,1-U)} so each marginal is Uniform(0,1).
#'   \item \code{mvn_marginal="empirical"} (advanced): Gaussian copula for dependence, but map
#'         \eqn{U=\Phi(Z)} to an empirical null p-value distribution \eqn{F_0^{-1}(U)} derived from
#'         \code{pool_p} (or an automatically constructed pool).
#' }
#'
#' @param gene_results A data.frame of gene-level results (typically one row per gene).
#'   Must contain the gene identifier column specified by \code{gene_col} and at least one
#'   p-value column (\code{p_raw_col} and/or \code{p_adj_col}). If you want Stouffer, it must
#'   also contain \code{z_col}. If you want weighted Stouffer, it must contain \code{weight_col}.
#'
#' @param pathways Pathway definitions. Either:
#'   \itemize{
#'     \item a data.frame with columns \code{pathway_id} and \code{gene_id} (and optional \code{pathway_name}), or
#'     \item a named list where each element is a character vector of gene IDs (names are pathway IDs).
#'   }
#'   Provide \code{pathways} \strong{or} \code{species} (not both).
#'
#' @param species Character. If \code{pathways} is \code{NULL}, a species keyword
#'   (one of "maize", "sorghum", "arabidopsis", "plant", "fly") used to auto-load built-in
#'   pathway sets. Provide \code{species} \strong{or} \code{pathways} (not both).
#'
#' @param pmn_gene_col Character. When \code{species} is used, the name of the gene identifier column
#'   in the built-in PMN pathway table (e.g., \code{"Gene-name"}). Ignored if \code{pathways} is provided.
#'
#' @param gene_col Character. Column name in \code{gene_results} containing gene IDs.
#'   Default \code{"GENE"}.
#'
#' @param p_raw_col Character. Column name in \code{gene_results} containing raw gene-level p-values.
#'   Default \code{"P"}.
#'
#' @param p_adj_col Character or \code{NULL}. Optional column name for adjusted/preferred gene p-values.
#'   If provided and present, this column is used for the p-based component methods (ACAT/Fisher/TFisher/minP).
#'   If \code{NULL}, \code{p_raw_col} is used.
#'
#' @param z_col Character or \code{NULL}. Optional column name in \code{gene_results} containing gene-level Z
#'   statistics used for Stouffer. If \code{NULL}, Stouffer is not computed (set to \code{NA}).
#'
#' @param weight_col Character or \code{NULL}. Optional column name in \code{gene_results} containing weights
#'   for Stouffer (e.g., gene size or other reliability weights). If \code{NULL}, equal weights are used.
#'
#' @param tau_grid Numeric vector of \eqn{\tau} values for adaptive soft TFisher. The method computes analytic
#'   p-values for each \eqn{\tau} and selects the minimum. Default is a typical decreasing grid.
#'
#' @param tau_cap Numeric. Maximum allowed \eqn{\tau} value (sanity bound). Default \code{1}.
#'
#' @param min_p Numeric. Lower bound for p-value clipping to avoid numerical issues (e.g., \code{log(0)} in Fisher,
#'   \code{tan()} blowups in ACAT). P-values are clipped to \code{[min_p, 1-min_p]}. Default \code{1e-15}.
#'
#' @param min_genes Integer. Minimum number of genes required for a pathway to be evaluated. Pathways with fewer
#'   genes are dropped. Default \code{2L}.
#'
#' @param magma_out Optional data.frame of MAGMA competitive gene-set results. Must include \code{pathway_id} and
#'   \code{magma_pvalue}. If \code{NULL}, MAGMA competitive is not included.
#'
#' @param include_magma_in_omni Logical. If \code{TRUE} and \code{magma_out} is provided, include MAGMA competitive
#'   p-values as an additional component in the omnibus across methods. Default \code{TRUE}.
#'
#' @param include_magma_in_perm Logical. If \code{TRUE}, also include MAGMA competitive within the resampling null.
#'   This is generally \strong{not recommended} because there is no cheap principled null generator for competitive
#'   MAGMA; default \code{FALSE}.
#'
#' @param omnibus Character. How to combine the component method p-values into one omnibus p-value per pathway:
#'   \code{"ACAT"} (ACAT across methods) or \code{"minP"} (Sidak-adjusted minimum across methods).
#'
#' @param B_perm Integer or \code{NULL}. Convenience parameter for the number of resampling iterations. When provided,
#'   it is mapped to \code{B_global} and/or \code{B_mvn} according to \code{perm_mode}.
#'
#' @param perm_mode Character. Which resampling/calibration mode to run for the omnibus:
#'   \code{"none"} (analytic only), \code{"global"} (global gene resampling), \code{"mvn"} (LD-aware MVN resampling),
#'   or \code{"both"} (run both; final p-value prefers MVN if available).
#'
#' @param mvn_marginal Character. Only used when \code{perm_mode} includes \code{"mvn"}.
#'   Controls how simulated MVN Z's are converted to p-values for p-based component methods:
#'   \code{"uniform"} (Gaussian copula; Uniform(0,1) marginals) or \code{"empirical"} (map copula uniforms to an
#'   empirical null p distribution).
#'
#' @param mvn_pool Character. Only used when \code{mvn_marginal="empirical"} and \code{pool_p} is \code{NULL}.
#'   Determines which gene-level p-values are used to construct the empirical null pool:
#'   \code{"use"} (the p column used for observed p-based methods; typically \code{p_adj_col} if present) or
#'   \code{"raw"} (\code{p_raw_col}).
#'
#' @param pool_p Numeric vector or \code{NULL}. Optional explicit empirical null pool of p-values used when
#'   \code{mvn_marginal="empirical"}. If \code{NULL}, an empirical pool is constructed automatically (excluding
#'   pathway genes to reduce leakage).
#'
#' @param B_global Integer. Number of global resampling iterations (back-compat / low-level control). Ignored if
#'   \code{B_perm} is provided (then \code{B_perm} + \code{perm_mode} determine \code{B_global}).
#'
#' @param B_mvn Integer. Number of MVN resampling iterations (back-compat / low-level control). Ignored if
#'   \code{B_perm} is provided.
#'
#' @param perm_pool Character. Only affects global resampling. Which p-value pool to resample from:
#'   \code{"obs"} uses the same p-values used for observed p-based methods (typically \code{p_adj_col} if present),
#'   while \code{"raw"} uses \code{p_raw_col}.
#'
#' @param magma_cor_file Character or \code{NULL}. Path to a 3-column gene correlation pairs file:
#'   \code{gene1 gene2 r}. Required for \code{perm_mode="mvn"} unless \code{magma_cor_pairs} is provided.
#'
#' @param magma_cor_pairs data.frame or \code{NULL}. In-memory correlation pairs table with columns
#'   \code{gene1}, \code{gene2}, \code{r}. If provided, it overrides \code{magma_cor_file}.
#'
#' @param make_PD Logical. If \code{TRUE}, force the pathway correlation matrix \eqn{R_S} to be positive definite
#'   (using \code{Matrix::nearPD()}) before Cholesky. Strongly recommended. Default \code{TRUE}.
#'
#' @param seed Integer or \code{NULL}. Random seed for reproducible resampling.
#'
#' @param output Logical. If \code{TRUE}, write results to CSV in \code{out_dir} and attach the file path as an attribute.
#'
#' @param out_dir Character. Output directory used when \code{output=TRUE}. Created if it does not exist.
#'
#' @return
#' A data.frame with one row per pathway, including:
#' \itemize{
#'   \item \code{acat_p}, \code{fisher_p}, \code{tfisher_p_analytic}, \code{minp_p_analytic}, \code{stouffer_p_analytic}
#'   \item \code{magma_pvalue} (if provided)
#'   \item \code{omni_p_analytic} and BH-adjusted version
#'   \item \code{omni_p_global} (if global resampling run) and BH-adjusted version
#'   \item \code{omni_p_mvn} (if MVN resampling run) and BH-adjusted version
#'   \item \code{omni_p_final} (priority: MVN > global > analytic) and BH-adjusted version
#' }
#'
# ============================================================
# MAGCAT omnibus (fast, 6-method; global + MVN calibration for OMNIBUS only)
#   - Component methods are computed by YOUR MAGCAT wrappers (analytic only)
#   - Permutations/resampling are ONLY for the OMNIBUS
# ============================================================

# ----------------------------
# Helper: safely grab MAGCAT functions (exported wrappers)
# Note: These functions are all in the same package, so we just use get()
# ----------------------------
.magcat_get_fun <- function(fname) {
  if (exists(fname, mode = "function", envir = parent.frame())) {
    return(get(fname, mode = "function", envir = parent.frame()))
  }
  if (exists(fname, mode = "function")) {
    return(get(fname, mode = "function"))
  }
  NULL
}

# ----------------------------
# MVN: Z -> p helpers (uniform vs empirical marginals)
# ----------------------------
.mvn_z_to_p <- function(Z, min_p = 1e-15) {
  U <- stats::pnorm(Z)
  P <- 2 * pmin(U, 1 - U)
  P <- pmax(pmin(P, 1 - min_p), min_p)
  P
}

.empirical_qfun <- function(pool_p) {
  pool_p <- suppressWarnings(as.numeric(pool_p))
  pool_p <- pool_p[is.finite(pool_p) & !is.na(pool_p) & pool_p > 0 & pool_p < 1]
  pool_p <- sort(pool_p)
  n <- length(pool_p)
  if (n < 1000) warning("pool_p is small (n=", n, "); empirical mapping may be noisy.", call. = FALSE)

  stats::approxfun(
    x = (1:n) / (n + 1),
    y = pool_p,
    method = "linear",
    rule = 2,
    ties = "ordered"
  )
}

.mvn_z_to_empirical_p <- function(Z, pool_p, min_p = 1e-15) {
  q0 <- .empirical_qfun(pool_p)
  U  <- stats::pnorm(Z)
  P  <- matrix(q0(as.vector(U)), nrow = nrow(Z), ncol = ncol(Z))
  P  <- pmax(pmin(P, 1 - min_p), min_p)
  P
}

# ----------------------------
# p-fixing (uses your fix_p_for_acat if available)
# ----------------------------
.magcat_fix_p <- function(p, min_p = 1e-15, do_fix = TRUE) {
  p <- suppressWarnings(as.numeric(p))
  p <- p[is.finite(p) & !is.na(p)]
  if (!length(p)) return(numeric(0))

  if (isTRUE(do_fix)) {
    # fix_p_for_acat is defined in this package (ACAT_wrappers.R)
    if (exists("fix_p_for_acat", mode = "function")) {
      fp <- get("fix_p_for_acat", mode = "function")
      p <- fp(p, min_p = min_p)
      p <- p[is.finite(p) & !is.na(p)]
      if (!length(p)) return(numeric(0))
      return(p)
    }
  }

  # fallback clip
  p[p <= 0] <- min_p
  p[p >= 1] <- 1 - min_p
  p
}

.magcat_bh <- function(p) {
  p <- as.numeric(p)
  out <- rep(NA_real_, length(p))
  idx <- which(is.finite(p) & !is.na(p))
  if (length(idx)) out[idx] <- stats::p.adjust(p[idx], method = "BH")
  out
}

## ----------------------------
## NEW tests
## ----------------------------

.magcat_qvalue <- function(p) {
  p <- as.numeric(p)
  out <- rep(NA_real_, length(p))
  idx <- which(is.finite(p) & !is.na(p))

  if (!length(idx)) return(out)

  if (requireNamespace("qvalue", quietly = TRUE)) {
    qobj <- qvalue::qvalue(p[idx])
    out[idx] <- as.numeric(qobj$qvalues)
  } else {
    warning("Package 'qvalue' not installed; falling back to BH.", call. = FALSE)
    out[idx] <- stats::p.adjust(p[idx], method = "BH")
  }
  out
}

# BKY / two-stage BH (Benjamini-Krieger-Yekutieli style)
.magcat_bky <- function(p, q = 0.10) {
  p <- as.numeric(p)
  out <- rep(NA_real_, length(p))
  idx <- which(is.finite(p) & !is.na(p))
  if (!length(idx)) return(out)

  p0 <- p[idx]
  m  <- length(p0)

  # stage 1
  q1 <- q / (1 + q)
  adj1 <- stats::p.adjust(p0, method = "BH")
  R1   <- sum(adj1 <= q1, na.rm = TRUE)
  m0hat <- max(1, m - R1)

  # stage 2 (adaptive BH level inflation)
  q2 <- q * m / m0hat
  out[idx] <- stats::p.adjust(p0, method = "BH") / q * q2  # rescale BH adj p's to target q2
  out[idx] <- pmin(out[idx], 1)
  out
}

.magcat_ihw <- function(p, covar) {
  p <- as.numeric(p)
  covar <- as.numeric(covar)
  out <- rep(NA_real_, length(p))
  idx <- which(is.finite(p) & !is.na(p) & is.finite(covar) & !is.na(covar))

  if (!length(idx)) return(out)

  if (requireNamespace("IHW", quietly = TRUE)) {
    ihw_res <- IHW::ihw(p[idx], covariates = covar[idx], alpha = 0.1)
    out[idx] <- IHW::adj_pvalues(ihw_res)
  } else {
    warning("Package 'IHW' not installed; returning NA for IHW q-values.", call. = FALSE)
  }
  out
}

## ----------------------------
## NEW tests END
## ----------------------------


# ----------------------------
# Omnibus across METHODS (ACAT or Sidak-minP)
# ----------------------------
.magcat_acat_combine <- function(p, min_p = 1e-15, do_fix = TRUE) {
  p <- .magcat_fix_p(p, min_p = min_p, do_fix = do_fix)
  p <- p[p > 0 & p < 1]
  if (!length(p)) return(NA_real_)
  tstat <- mean(tan((0.5 - p) * pi))
  out <- 0.5 - atan(tstat) / pi
  out <- max(min(out, 1 - min_p), min_p)
  as.numeric(out)
}

.magcat_minp_sidak <- function(p, min_p = 1e-15, do_fix = TRUE) {
  p <- .magcat_fix_p(p, min_p = min_p, do_fix = do_fix)
  p <- p[p > 0 & p < 1]
  k <- length(p)
  if (k < 1L) return(NA_real_)
  pmin <- min(p)
  1 - (1 - pmin)^k
}

.magcat_omni_methods <- function(p_methods, omnibus = c("ACAT", "minP"), min_p = 1e-15, do_fix = TRUE) {
  omnibus <- match.arg(omnibus)
  pv <- .magcat_fix_p(p_methods, min_p = min_p, do_fix = do_fix)
  pv <- pv[pv > 0 & pv < 1]
  if (!length(pv)) return(NA_real_)

  if (omnibus == "ACAT") {
    .magcat_acat_combine(pv, min_p = min_p, do_fix = do_fix)
  } else {
    k <- length(pv)
    pmin <- min(pv)
    1 - (1 - pmin)^k
  }
}


####################################################################################
####### NEW FUNCTIONS FOR INDIVIDUAL SIMULATION
####################################################################################

# empirical calibration for "p-like" scores where smaller = more extreme
.emp_p_obs_lower <- function(p_null, p_obs) {
  p_null <- as.numeric(p_null)
  p_obs  <- as.numeric(p_obs)
  p_null <- p_null[is.finite(p_null) & !is.na(p_null)]
  B <- length(p_null)
  if (!is.finite(p_obs) || is.na(p_obs) || B < 1L) return(NA_real_)
  (1 + sum(p_null <= p_obs)) / (B + 1)
}

# OLD VERSION
# # for each null draw b, return its own calibrated p_b = (1 + #{b': p_{b'} <= p_b})/(B+1)
# .emp_p_null_lower <- function(p_null) {
#   p_null <- as.numeric(p_null)
#   B <- length(p_null)
#   if (B < 1L) return(numeric(0))
#   # NA/inf become 1 (uninformative)
#   p2 <- p_null
#   bad <- !is.finite(p2) | is.na(p2)
#   if (any(bad)) p2[bad] <- 1
#   (1 + rank(p2, ties.method = "max")) / (B + 1)
# }


# NEW VERSION
.emp_p_null_lower <- function(p_null) {
  p2 <- as.numeric(p_null)
  B <- length(p2)
  if (B < 1L) return(numeric(0))
  bad <- !is.finite(p2) | is.na(p2)
  if (any(bad)) p2[bad] <- 1

  # leave-one-out empirical p for each null draw:
  # allows minimum 1/(B+1), matching observed Monte Carlo support
  rank_max <- rank(p2, ties.method = "max")
  as.numeric(rank_max) / (B + 1)
}



####################################################################################
####### END
####################################################################################


# ----------------------------
# MAGMA correlation helpers (MVN)
# ----------------------------
.magma_read_gene_cor_pairs <- function(magma_cor_file = NULL, magma_cor_pairs = NULL) {
  if (!is.null(magma_cor_pairs)) {
    tab <- magma_cor_pairs
    if (!is.data.frame(tab) || ncol(tab) < 3) stop("magma_cor_pairs must be a data.frame with >=3 cols.", call. = FALSE)
    tab <- tab[, 1:3]
    names(tab) <- c("gene1", "gene2", "r")
    tab$gene1 <- as.character(tab$gene1)
    tab$gene2 <- as.character(tab$gene2)
    tab$r <- suppressWarnings(as.numeric(tab$r))
    tab <- tab[is.finite(tab$r) & !is.na(tab$r), , drop = FALSE]
    return(tab)
  }

  if (is.null(magma_cor_file)) return(NULL)
  if (!file.exists(magma_cor_file)) stop("magma_cor_file does not exist: ", magma_cor_file, call. = FALSE)

  if (requireNamespace("data.table", quietly = TRUE)) {
    tab <- data.table::fread(magma_cor_file, header = FALSE, data.table = FALSE, fill = TRUE, showProgress = FALSE)
  } else {
    tab <- utils::read.table(magma_cor_file, header = FALSE, stringsAsFactors = FALSE, fill = TRUE)
  }
  if (ncol(tab) < 3) stop("magma_cor_file must have 3 columns: gene1 gene2 r", call. = FALSE)

  tab <- tab[, 1:3]
  names(tab) <- c("gene1", "gene2", "r")
  tab$gene1 <- as.character(tab$gene1)
  tab$gene2 <- as.character(tab$gene2)
  tab$r <- suppressWarnings(as.numeric(tab$r))
  tab <- tab[is.finite(tab$r) & !is.na(tab$r), , drop = FALSE]
  tab
}

# Old without shrinkage
# .magma_build_R_from_pairs <- function(genes, cor_pairs, make_PD = TRUE) {
#   genes <- as.character(genes)
#   genes <- genes[!is.na(genes) & genes != ""]
#   d <- length(genes)
#   if (d < 2L) return(NULL)

#   R <- diag(1, d)
#   rownames(R) <- genes
#   colnames(R) <- genes

#   if (!is.null(cor_pairs) && nrow(cor_pairs)) {
#     pos <- stats::setNames(seq_len(d), genes)

#     g1 <- as.character(cor_pairs$gene1)
#     g2 <- as.character(cor_pairs$gene2)
#     r  <- suppressWarnings(as.numeric(cor_pairs$r))

#     keep <- is.finite(r) & !is.na(r) & (g1 %in% genes) & (g2 %in% genes)
#     if (any(keep)) {
#       g1 <- g1[keep]; g2 <- g2[keep]; r <- r[keep]
#       i <- unname(pos[g1]); j <- unname(pos[g2])
#       ok <- is.finite(i) & is.finite(j) & !is.na(i) & !is.na(j) & (i != j)
#       i <- i[ok]; j <- j[ok]; r <- r[ok]

#       if (length(i)) {
#         r <- pmax(pmin(r, 0.999), -0.999)
#         R[cbind(i, j)] <- r
#         R[cbind(j, i)] <- r
#       }
#     }
#   }

#   if (isTRUE(make_PD)) {
#     if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix package required for make_PD=TRUE.", call. = FALSE)
#     npd <- Matrix::nearPD(R, corr = TRUE)
#     R <- as.matrix(npd$mat)
#   }

#   R
# }

# New with shinkage
.magma_build_R_from_pairs <- function(genes, cor_pairs, make_PD = TRUE, eps = 1e-6) {
  genes <- as.character(genes)
  genes <- genes[!is.na(genes) & genes != ""]
  d <- length(genes)
  if (d < 2L) return(NULL)

  R <- diag(1, d)
  rownames(R) <- genes
  colnames(R) <- genes

  if (!is.null(cor_pairs) && nrow(cor_pairs)) {
    pos <- stats::setNames(seq_len(d), genes)

    g1 <- as.character(cor_pairs$gene1)
    g2 <- as.character(cor_pairs$gene2)
    r  <- suppressWarnings(as.numeric(cor_pairs$r))

    keep <- is.finite(r) & !is.na(r) & (g1 %in% genes) & (g2 %in% genes)
    if (any(keep)) {
      g1 <- g1[keep]; g2 <- g2[keep]; r <- r[keep]
      i <- unname(pos[g1]); j <- unname(pos[g2])
      ok <- is.finite(i) & is.finite(j) & !is.na(i) & !is.na(j) & (i != j)
      i <- i[ok]; j <- j[ok]; r <- r[ok]

      if (length(i)) {
        r <- pmax(pmin(r, 0.999), -0.999)
        R[cbind(i, j)] <- r
        R[cbind(j, i)] <- r
      }
    }
  }

  if (isTRUE(make_PD)) {
    if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix package required for make_PD=TRUE.", call. = FALSE)
    npd <- Matrix::nearPD(R, corr = TRUE)
    R <- as.matrix(npd$mat)
  }

  # NEW: tiny diagonal shrinkage/jitter so chol() is stable
  eps <- as.numeric(eps)
  if (is.finite(eps) && eps > 0) {
    R <- (1 - eps) * R + eps * diag(1, ncol(R))
  }

  R
}


.magma_simulate_Z_from_R <- function(R, B) {
  B <- as.integer(B)
  d <- ncol(R)
  if (B < 1L || d < 2L) return(NULL)

  U <- tryCatch(chol(R), error = function(e) NULL)
  if (is.null(U)) stop("chol(R) failed even after PD-fix. Check correlation input.", call. = FALSE)

  Z0 <- matrix(stats::rnorm(B * d), nrow = B, ncol = d)
  Z  <- Z0 %*% U
  colnames(Z) <- colnames(R)
  Z
}

# ----------------------------
# Global sampler helper (WITHOUT replacement per permutation)
# ----------------------------
.magcat_sample_idx_mat <- function(n_pool, d, B) {
  n_pool <- as.integer(n_pool)
  d      <- as.integer(d)
  B      <- as.integer(B)
  if (n_pool < 1L || d < 1L || B < 1L) return(NULL)

  if (n_pool >= d) {
    idx <- replicate(B, sample.int(n_pool, size = d, replace = FALSE))
    t(idx) # B x d
  } else {
    matrix(sample.int(n_pool, size = B * d, replace = TRUE), nrow = B, ncol = d)
  }
}

# ----------------------------
# prepare (fast reuse)
# ----------------------------
magcat_omni2_prepare <- function(gene_results,
                                pathways = NULL,
                                species = NULL,
                                pmn_gene_col = NULL,
                                gene_col = "GENE",
                                p_raw_col = "P",
                                p_adj_col = NULL,
                                z_col = NULL,
                                weight_col = NULL,
                                tau_grid = c(0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.001),
                                tau_cap = 1,
                                min_p = 1e-15,
                                min_genes = 2L,
                                seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  if (!is.data.frame(gene_results)) stop("gene_results must be a data.frame.", call. = FALSE)
  if (!(gene_col %in% names(gene_results))) stop("gene_results missing gene_col: ", gene_col, call. = FALSE)

  if (!is.null(pathways) && !is.null(species)) stop("Provide either pathways OR species, not both.", call. = FALSE)
  if (is.null(pathways) && is.null(species)) stop("Need pathways (list/df) or species.", call. = FALSE)

  # load pathways if species (use wrapper helper from this package)
  if (is.null(pathways)) {
    # magcat_load_pathways is defined in this package (ACAT_wrappers.R)
    if (exists("magcat_load_pathways", mode = "function")) {
      pathways <- magcat_load_pathways(species = species, gene_col = pmn_gene_col)
    } else {
      stop("Could not find magcat_load_pathways(). This should not happen.", call. = FALSE)
    }
  }

  # choose p column for p-based methods: prefer adjusted if provided and present
  p_use_col <- p_raw_col
  if (!is.null(p_adj_col) && (p_adj_col %in% names(gene_results))) p_use_col <- p_adj_col
  if (!(p_use_col %in% names(gene_results))) stop("Chosen p column not found in gene_results: ", p_use_col, call. = FALSE)
  if (!(p_raw_col %in% names(gene_results))) p_raw_col <- NULL

  # gene ids
  g_raw  <- as.character(gene_results[[gene_col]])
  g_norm <- tolower(trimws(g_raw))
  okg <- !is.na(g_norm) & g_norm != "" & g_norm != "unknown"
  g_raw  <- g_raw[okg]
  g_norm <- g_norm[okg]

  # p (use)
  p_use <- suppressWarnings(as.numeric(gene_results[[p_use_col]]))[okg]
  okp <- is.finite(p_use) & !is.na(p_use)
  g_raw  <- g_raw[okp]
  g_norm <- g_norm[okp]
  p_use  <- p_use[okp]

  # raw p aligned
  p_raw_vec <- NULL
  if (!is.null(p_raw_col) && (p_raw_col %in% names(gene_results))) {
    p_raw_vec <- suppressWarnings(as.numeric(gene_results[[p_raw_col]]))[okg][okp]
  }

  # z aligned
  z_use <- NULL
  if (!is.null(z_col)) {
    if (!(z_col %in% names(gene_results))) stop("z_col not found: ", z_col, call. = FALSE)
    z_use <- suppressWarnings(as.numeric(gene_results[[z_col]]))[okg][okp]
  }

  # weights aligned
  w_use <- NULL
  if (!is.null(weight_col)) {
    if (!(weight_col %in% names(gene_results))) stop("weight_col not found: ", weight_col, call. = FALSE)
    w_use <- suppressWarnings(as.numeric(gene_results[[weight_col]]))[okg][okp]
  }

  # dedup by normalized gene id (keep first)
  keep1 <- !duplicated(g_norm)
  g_raw  <- g_raw[keep1]
  g_norm <- g_norm[keep1]
  p_use  <- p_use[keep1]
  if (!is.null(p_raw_vec)) p_raw_vec <- p_raw_vec[keep1]
  if (!is.null(z_use)) z_use <- z_use[keep1]
  if (!is.null(w_use)) w_use <- w_use[keep1]

  gene_map <- stats::setNames(g_raw, g_norm)
  gene_p <- stats::setNames(p_use, g_norm)
  gene_z <- if (!is.null(z_use)) stats::setNames(z_use, g_norm) else NULL
  gene_w <- if (!is.null(w_use)) stats::setNames(w_use, g_norm) else NULL

  # coerce pathways to long
  if (is.data.frame(pathways)) {
    if (!all(c("pathway_id","gene_id") %in% names(pathways))) stop("pathways df must have pathway_id and gene_id.", call. = FALSE)
    if (!("pathway_name" %in% names(pathways))) pathways$pathway_name <- pathways$pathway_id
    pw_long <- pathways[, c("pathway_id","pathway_name","gene_id"), drop = FALSE]
    pw_long$gene_id <- as.character(pw_long$gene_id)
  } else if (is.list(pathways)) {
    pid <- names(pathways); if (is.null(pid)) pid <- paste0("PWY_", seq_along(pathways))
    pw_long <- data.frame(
      pathway_id   = rep(pid, lengths(pathways)),
      pathway_name = rep(pid, lengths(pathways)),
      gene_id      = unlist(pathways, use.names = FALSE),
      stringsAsFactors = FALSE
    )
  } else {
    stop("pathways must be a data.frame or list.", call. = FALSE)
  }

  pw_long$gene_id   <- gsub("^gene:", "", pw_long$gene_id, ignore.case = TRUE)
  pw_long$gene_norm <- tolower(trimws(pw_long$gene_id))
  pw_long <- pw_long[!is.na(pw_long$gene_norm) & pw_long$gene_norm != "" & pw_long$gene_norm != "unknown", , drop = FALSE]

  g_list_norm <- split(pw_long$gene_norm, pw_long$pathway_id)
  pw_name     <- tapply(pw_long$pathway_name, pw_long$pathway_id, function(x) x[1])

  # intersect with available genes
  g_list_norm <- lapply(g_list_norm, function(g) unique(g[g %in% names(gene_p)]))

  # drop small pathways early
  min_genes <- as.integer(min_genes)
  if (!is.na(min_genes) && min_genes >= 1L) {
    keep_pw <- lengths(g_list_norm) >= min_genes
    g_list_norm <- g_list_norm[keep_pw]
  }

  pid <- names(g_list_norm)

  # per pathway vectors
  p_list <- lapply(g_list_norm, function(g) unname(gene_p[g]))
  z_list <- if (!is.null(gene_z)) lapply(g_list_norm, function(g) unname(gene_z[g])) else NULL
  w_list <- if (!is.null(gene_w)) lapply(g_list_norm, function(g) unname(gene_w[g])) else NULL
  g_list_raw <- lapply(g_list_norm, function(g) unname(gene_map[g]))
  names(g_list_raw) <- names(g_list_norm)  

  # tau sanitize
  tau_grid <- suppressWarnings(as.numeric(tau_grid))
  tau_grid <- tau_grid[is.finite(tau_grid) & !is.na(tau_grid) & tau_grid > 0 & tau_grid <= tau_cap]
  tau_grid <- sort(unique(tau_grid), decreasing = TRUE)
  if (!length(tau_grid)) tau_grid <- 0.05

  # minimal gene_results to feed your wrapper functions (already deduped/cleaned)
  n <- length(g_raw)
  if (n == 0L) stop("No valid genes left after filtering gene_results.", call. = FALSE)

  # IMPORTANT: create the data.frame with n rows immediately
  gene_results_core <- data.frame(
    setNames(list(g_raw), gene_col),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # then add remaining columns
  gene_results_core[[p_use_col]] <- p_use
  if (!is.null(z_col))      gene_results_core[[z_col]] <- z_use
  if (!is.null(weight_col)) gene_results_core[[weight_col]] <- w_use

  # pathway list to feed your wrapper functions
  pathways_list_raw <- g_list_raw
  names(pathways_list_raw) <- pid

  structure(
    list(
      pathway_id   = pid,
      pathway_name = as.character(pw_name[pid]),
      g_list_norm  = g_list_norm,
      g_list_raw   = g_list_raw,
      p_list       = p_list,
      z_list       = z_list,
      w_list       = w_list,

      gene_p_all_use = p_use,
      gene_p_all_raw = p_raw_vec,
      gene_z_all     = z_use,
      gene_w_all     = w_use,

      gene_col   = gene_col,
      p_use_col  = p_use_col,
      p_raw_col  = p_raw_col,
      p_adj_col  = p_adj_col,
      z_col      = z_col,
      weight_col = weight_col,
      tau_grid   = tau_grid,
      tau_cap    = tau_cap,
      min_p      = min_p,
      min_genes  = min_genes,
      seed       = seed,

      gene_results_core  = gene_results_core,
      pathways_list_raw  = pathways_list_raw
    ),
    class = "magcat_omni2_prep"
  )
}

# ----------------------------
# Run (components from YOUR wrappers; perms only for omnibus)
# ----------------------------
magcat_omni2_run <- function(prep,
                            omnibus = c("ACAT","minP"),
                            # component wrapper params
                            do_fix = TRUE,
                            stouffer_min_abs_w = 1e-8,
                            stouffer_alternative = c("greater","two.sided","less"),
                            # magma competitive
                            magma_out = NULL,
                            include_magma_in_omni = TRUE,
                            include_magma_in_perm = FALSE,
                            # omnibus calibration
                            B_global = 0L,
                            B_mvn = 0L,
                            perm_pool = c("raw","obs"),
                            mvn_marginal = c("uniform","empirical"),
                            mvn_pool     = c("use","raw"),
                            pool_p       = NULL,
                            magma_cor_file = NULL,
                            magma_cor_pairs = NULL,
                             mvn_calibrate_components = FALSE,   # <<<<<< ADD THIS LINE
                            make_PD = TRUE,
                            output = FALSE,
                            out_dir = "magcat_omni2") {

  if (!inherits(prep, "magcat_omni2_prep")) stop("prep must come from magcat_omni2_prepare().", call. = FALSE)
  omnibus   <- match.arg(omnibus)
  perm_pool <- match.arg(perm_pool)
  mvn_marginal <- match.arg(mvn_marginal)
  mvn_pool     <- match.arg(mvn_pool)
  stouffer_alternative <- match.arg(stouffer_alternative)

  B_global <- as.integer(B_global)
  B_mvn    <- as.integer(B_mvn)
  do_perm  <- (B_global > 0L || B_mvn > 0L)

  pid <- prep$pathway_id
  n_genes <- lengths(prep$g_list_norm)

  res <- data.frame(
    pathway_id   = pid,
    pathway_name = prep$pathway_name,
    n_genes      = as.integer(n_genes),
    stringsAsFactors = FALSE
  )

  # ============================================================
  # OBSERVED component p-values: CALL YOUR MAGCAT WRAPPERS (analytic only)
  # ============================================================
  acat_fun    <- .magcat_get_fun("magcat_acat_pathways")
  fisher_fun  <- .magcat_get_fun("magcat_fisher_pathways")
  minp_fun    <- .magcat_get_fun("magcat_minp_pathways")
  tfisher_fun <- .magcat_get_fun("magcat_soft_tfisher_adaptive_pathways")
  stouf_fun   <- .magcat_get_fun("magcat_stoufferZ_pathways")

  if (is.null(acat_fun) || is.null(fisher_fun) || is.null(minp_fun) || is.null(tfisher_fun)) {
    stop("Could not find required MAGCAT wrappers. Load MAGCAT or source your wrapper files first.", call. = FALSE)
  }

  gene_df <- prep$gene_results_core
  pw_list <- prep$pathways_list_raw

  
  acat_tab <- acat_fun(
  gene_results = gene_df,
  pathways     = pw_list,
  gene_col     = prep$gene_col,
  p_col        = prep$p_use_col,
  min_p        = prep$min_p,
  do_fix       = isTRUE(do_fix),
  output       = FALSE,
  out_dir      = tempdir()
)
res$acat_p <- acat_tab$acat_p[match(res$pathway_id, acat_tab$pathway_id)]

  fish_tab <- fisher_fun(
    gene_results = gene_df,
    pathways     = pw_list,
    gene_col     = prep$gene_col,
    p_col        = prep$p_use_col,
    min_p        = prep$min_p,
    do_fix       = isTRUE(do_fix),
    output       = FALSE
  )
  res$fisher_p <- fish_tab$fisher_p[match(res$pathway_id, fish_tab$pathway_id)]

  tf_tab <- tfisher_fun(
  gene_results = gene_df,
  pathways     = pw_list,
  gene_col     = prep$gene_col,
  p_col        = prep$p_use_col,
  tau_grid     = prep$tau_grid,
  min_p        = prep$min_p,
  do_fix       = isTRUE(do_fix),
  output       = FALSE
)

res$tfisher_p_analytic <- tf_tab$tfisher_p_omni[match(res$pathway_id, tf_tab$pathway_id)]
res$tau_hat            <- tf_tab$tau_hat[match(res$pathway_id, tf_tab$pathway_id)]
res$tfisher_stat_hat   <- tf_tab$tfisher_stat_hat[match(res$pathway_id, tf_tab$pathway_id)]


minp_tab <- minp_fun(
  gene_results     = gene_df,
  pathways         = pw_list,
  gene_col         = prep$gene_col,
  p_col            = prep$p_use_col,
  min_p            = prep$min_p,
  do_fix           = isTRUE(do_fix),
  output           = FALSE
)
res$minp_p_analytic <- minp_tab$minp_p[match(res$pathway_id, minp_tab$pathway_id)]


res$stouffer_p_analytic <- NA_real_
if (!is.null(prep$z_col)) {
  if (is.null(stouf_fun)) stop("z_col provided but magcat_stoufferZ_pathways() not found.", call. = FALSE)

  st_tab <- stouf_fun(
    gene_results   = gene_df,
    pathways       = pw_list,
    gene_col       = prep$gene_col,
    z_col          = prep$z_col,
    weight_col     = prep$weight_col,
    min_abs_w      = stouffer_min_abs_w,
    alternative    = stouffer_alternative,
    output         = FALSE
  )

  # your simplified version returns stouffer_p (not stouffer_p_analytic)
  if ("stouffer_p" %in% names(st_tab)) {
    res$stouffer_p_analytic <- st_tab$stouffer_p[match(res$pathway_id, st_tab$pathway_id)]
  } else {
    res$stouffer_p_analytic <- st_tab$stouffer_p_analytic[match(res$pathway_id, st_tab$pathway_id)]
  }
}




  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ##### OLD
  # # ACAT
  # acat_tab <- acat_fun(
  #   gene_results = gene_df,
  #   pathways     = pw_list,
  #   gene_col     = prep$gene_col,
  #   p_col        = prep$p_use_col,
  #   min_p        = prep$min_p,
  #   do_fix       = isTRUE(do_fix),
  #   B            = 0L,
  #   seed         = prep$seed,
  #   output       = FALSE
  # )
  # res$acat_p <- acat_tab$acat_p[match(res$pathway_id, acat_tab$pathway_id)]

  # # Fisher
  # fish_tab <- fisher_fun(
  #   gene_results = gene_df,
  #   pathways     = pw_list,
  #   gene_col     = prep$gene_col,
  #   p_col        = prep$p_use_col,
  #   min_p        = prep$min_p,
  #   do_fix       = isTRUE(do_fix),
  #   output       = FALSE
  # )
  # res$fisher_p <- fish_tab$fisher_p[match(res$pathway_id, fish_tab$pathway_id)]

  # # Adaptive soft TFisher (analytic only; NO per-method permutations)
  # tf_tab <- tfisher_fun(
  #   gene_results     = gene_df,
  #   pathways         = pw_list,
  #   gene_col         = prep$gene_col,
  #   p_col            = prep$p_use_col,
  #   tau_grid         = prep$tau_grid,
  #   min_p            = prep$min_p,
  #   do_fix           = isTRUE(do_fix),
  #   B_perm           = 0L,
  #   perm_mode        = "resample_global",  # irrelevant when B_perm=0
  #   magma_genes_out  = NULL,
  #   magma_cor_file   = NULL,
  #   make_PD          = make_PD,
  #   seed             = prep$seed,
  #   analytic_logical = TRUE,
  #   output           = FALSE
  # )
  # res$tfisher_p_analytic <- tf_tab$tfisher_p_analytic[match(res$pathway_id, tf_tab$pathway_id)]
  # res$tau_hat            <- tf_tab$tau_hat[match(res$pathway_id, tf_tab$pathway_id)]
  # res$tfisher_stat_hat   <- tf_tab$tfisher_stat_hat[match(res$pathway_id, tf_tab$pathway_id)]

  # # minP (analytic only; NO per-method permutations)
  # minp_tab <- minp_fun(
  #   gene_results     = gene_df,
  #   pathways         = pw_list,
  #   gene_col         = prep$gene_col,
  #   p_col            = prep$p_use_col,
  #   min_p            = prep$min_p,
  #   do_fix           = isTRUE(do_fix),
  #   B_perm           = 0L,
  #   seed             = prep$seed,
  #   analytic_logical = TRUE,
  #   output           = FALSE
  # )
  # res$minp_p_analytic <- minp_tab$minp_p_analytic[match(res$pathway_id, minp_tab$pathway_id)]

  # # StoufferZ (analytic only) if z_col provided
  # res$stouffer_p_analytic <- NA_real_
  # if (!is.null(prep$z_col)) {
  #   if (is.null(stouf_fun)) stop("z_col provided but magcat_stoufferZ_pathways() not found.", call. = FALSE)

  #   st_tab <- stouf_fun(
  #     gene_results    = gene_df,
  #     pathways        = pw_list,
  #     gene_col        = prep$gene_col,
  #     z_col           = prep$z_col,
  #     weight_col      = prep$weight_col,
  #     min_abs_w       = stouffer_min_abs_w,
  #     B_perm          = 0L,
  #     perm_mode       = "resample_global", # irrelevant when B_perm=0
  #     magma_genes_out = NULL,
  #     magma_cor_file  = NULL,
  #     make_PD         = make_PD,
  #     seed            = prep$seed,
  #     alternative     = stouffer_alternative,
  #     output          = FALSE
  #   )
  #   res$stouffer_p_analytic <- st_tab$stouffer_p_analytic[match(res$pathway_id, st_tab$pathway_id)]
  # }

  # MAGMA competitive (optional)
  res$magma_pvalue <- NA_real_
  if (!is.null(magma_out)) {
    if (!is.data.frame(magma_out)) stop("magma_out must be a data.frame (MAGMA competitive output).", call. = FALSE)
    if (!all(c("pathway_id","magma_pvalue") %in% names(magma_out))) {
      stop("magma_out must have columns pathway_id and magma_pvalue.", call. = FALSE)
    }
    res$magma_pvalue <- magma_out$magma_pvalue[match(res$pathway_id, magma_out$pathway_id)]
  }

  # If permutations requested but MAGMA not included in perm, omit it from omni_p_* (keep as separate column)
  include_magma_eff <- isTRUE(include_magma_in_omni) && any(is.finite(res$magma_pvalue))
  if (do_perm && include_magma_eff && !isTRUE(include_magma_in_perm)) {
    warning("Permutations requested but include_magma_in_perm=FALSE; MAGMA competitive is kept in 'magma_pvalue' but omitted from omni_p_* so the permuted omnibus is well-defined.", call. = FALSE)
    include_magma_eff <- FALSE
  }

  # omnibus analytic (across methods)
  res$omni_p_analytic <- vapply(seq_len(nrow(res)), function(i) {
    comps <- c(res$acat_p[i],
               res$fisher_p[i],
               res$tfisher_p_analytic[i],
               res$minp_p_analytic[i],
               res$stouffer_p_analytic[i])
    if (isTRUE(include_magma_eff) && is.finite(res$magma_pvalue[i])) comps <- c(comps, res$magma_pvalue[i])
    .magcat_omni_methods(comps, omnibus = omnibus, min_p = prep$min_p, do_fix = isTRUE(do_fix))
  }, numeric(1))
  res$omni_p_analytic_BH <- .magcat_bh(res$omni_p_analytic)

  # ============================================================
  # OMNIBUS calibration: global resampling (ONLY omnibus)
  # ============================================================
  res$omni_p_global <- NA_real_
  res$omni_p_global_BH <- NA_real_

  if (B_global > 0L) {

    if (!is.null(prep$seed)) set.seed(prep$seed)

    pool_p <- NULL
    if (perm_pool == "raw" && !is.null(prep$p_raw_col) && !is.null(prep$gene_p_all_raw)) {
      pool_p <- suppressWarnings(as.numeric(prep$gene_p_all_raw))
    } else {
      pool_p <- suppressWarnings(as.numeric(prep$gene_p_all_use))
    }

    pool_z <- NULL
    pool_w <- NULL
    if (!is.null(prep$z_col) && !is.null(prep$gene_z_all)) {
      pool_z <- suppressWarnings(as.numeric(prep$gene_z_all))
      if (!is.null(prep$weight_col) && !is.null(prep$gene_w_all)) {
        pool_w <- suppressWarnings(as.numeric(prep$gene_w_all))
      }
    }

    ok <- is.finite(pool_p) & !is.na(pool_p) & pool_p > 0 & pool_p < 1
    if (!is.null(pool_z)) ok <- ok & is.finite(pool_z) & !is.na(pool_z)
    if (!is.null(pool_w)) ok <- ok & is.finite(pool_w) & !is.na(pool_w)

    pool_p <- pool_p[ok]
    if (!is.null(pool_z)) pool_z <- pool_z[ok]
    if (!is.null(pool_w)) pool_w <- pool_w[ok]

    pool_p <- .magcat_fix_p(pool_p, min_p = prep$min_p, do_fix = isTRUE(do_fix))
    if (length(pool_p) < 2L) stop("Global resampling pool has <2 valid p-values.", call. = FALSE)

    n_pool <- length(pool_p)

    # helper for stouffer p from Z matrix
    .stouffer_p_from_Zmat <- function(Zmat, w = NULL) {
      if (is.null(Zmat) || !ncol(Zmat)) return(rep(NA_real_, nrow(Zmat)))
      if (is.null(w)) {
        Zs <- rowSums(Zmat) / sqrt(ncol(Zmat))
      } else {
        w <- suppressWarnings(as.numeric(w))
        bad <- !is.finite(w) | is.na(w) | abs(w) <= stouffer_min_abs_w
        if (any(bad)) {
          w_ok <- abs(w[!bad])
          repl <- if (length(w_ok)) stats::median(w_ok) else 1
          w[bad] <- repl
        }
        denom <- sqrt(sum(w^2))
        Zs <- as.numeric(Zmat %*% w) / denom
      }

      if (stouffer_alternative == "greater") {
        return(stats::pnorm(Zs, lower.tail = FALSE))
      } else if (stouffer_alternative == "less") {
        return(stats::pnorm(Zs, lower.tail = TRUE))
      } else {
        return(2 * stats::pnorm(-abs(Zs)))
      }
    }

    for (i in seq_len(nrow(res))) {

      d <- res$n_genes[i]
      if (!is.finite(d) || d < prep$min_genes) next

      omni_obs <- res$omni_p_analytic[i]
      if (!is.finite(omni_obs) || is.na(omni_obs)) next

      idx_mat <- .magcat_sample_idx_mat(n_pool = n_pool, d = d, B = B_global)
      if (is.null(idx_mat)) next

      P <- matrix(pool_p[idx_mat], nrow = B_global, ncol = d)

      # ACAT across genes (fast)
      pA <- vapply(seq_len(B_global), function(b) .magcat_acat_combine(P[b, ], min_p = prep$min_p, do_fix = isTRUE(do_fix)), numeric(1))
      # Fisher across genes (fast)
      pF <- stats::pchisq(-2 * rowSums(log(P)), df = 2 * d, lower.tail = FALSE)

      # TFisher adaptive (slower; but only in omnibus null)
      pT <- rep(NA_real_, B_global)
      if (!requireNamespace("TFisher", quietly = TRUE)) stop("TFisher required for omnibus TFisher null.", call. = FALSE)
      for (b in seq_len(B_global)) {
        pb <- .magcat_fix_p(P[b, ], min_p = prep$min_p, do_fix = isTRUE(do_fix))
        if (length(pb) < 2L) next

        # compute adaptive p (min over tau) using TFisher directly
        best <- Inf
        for (tau in prep$tau_grid) {
          st <- TFisher::stat.soft(p = pb, tau1 = tau)
          pv <- 1 - as.numeric(TFisher::p.soft(q = st, n = length(pb), tau1 = tau, M = NULL))
          if (is.finite(pv) && pv < best) best <- pv
        }
        pT[b] <- best
      }

      # minP Sidak across genes
      pM <- {
        pmin <- apply(P, 1, min)
        1 - (1 - pmin)^d
      }

      # Stouffer null if we have Z pool
      pS <- rep(NA_real_, B_global)
      if (!is.null(pool_z)) {
        Zmat <- matrix(pool_z[idx_mat], nrow = B_global, ncol = d)
        wvec <- NULL
        if (!is.null(pool_w)) wvec <- as.numeric(pool_w[idx_mat[1, ]]) # (approx) keep typical scale; optional
        pS <- .stouffer_p_from_Zmat(Zmat, w = NULL)  # weights in global null typically omitted
      }

      omni_null <- vapply(seq_len(B_global), function(b) {
        comps <- c(pA[b], pF[b], pT[b], pM[b], pS[b])
        if (isTRUE(include_magma_eff) && is.finite(res$magma_pvalue[i])) comps <- c(comps, res$magma_pvalue[i])
        .magcat_omni_methods(comps, omnibus = omnibus, min_p = prep$min_p, do_fix = isTRUE(do_fix))
      }, numeric(1))

      res$omni_p_global[i] <- (1 + sum(omni_null <= omni_obs, na.rm = TRUE)) / (B_global + 1)
    }

    res$omni_p_global_BH <- .magcat_bh(res$omni_p_global)
  }

  # ============================================================
  # OLD VERSION
  # ============================================================
  # ============================================================
  # OMNIBUS calibration: MVN resampling (LD-aware; ONLY omnibus)
  # ============================================================
  # res$omni_p_mvn <- NA_real_
  # res$omni_p_mvn_BH <- NA_real_

  # if (B_mvn > 0L) {

  #   cor_pairs <- .magma_read_gene_cor_pairs(magma_cor_file = magma_cor_file, magma_cor_pairs = magma_cor_pairs)
  #   if (is.null(cor_pairs) || !nrow(cor_pairs)) stop("MVN resampling requires magma_cor_file or magma_cor_pairs.", call. = FALSE)

  #   if (!is.null(prep$seed)) set.seed(prep$seed)

  #   for (i in seq_len(nrow(res))) {

  #     d <- res$n_genes[i]
  #     if (!is.finite(d) || d < prep$min_genes) next

  #     omni_obs <- res$omni_p_analytic[i]
  #     if (!is.finite(omni_obs) || is.na(omni_obs)) next

  #     genes_S <- prep$g_list_raw[[res$pathway_id[i]]]
  #     genes_S <- as.character(genes_S)
  #     genes_S <- genes_S[!is.na(genes_S) & genes_S != ""]
  #     if (length(genes_S) < prep$min_genes) next

  #     R <- .magma_build_R_from_pairs(genes = genes_S, cor_pairs = cor_pairs, make_PD = make_PD)
  #     if (is.null(R)) next

  #     Z <- .magma_simulate_Z_from_R(R, B = B_mvn)
  #     if (is.null(Z)) next

  #     if (mvn_marginal == "uniform") {
  #       P <- .mvn_z_to_p(Z, min_p = prep$min_p)
  #     } else {
  #       if (!is.null(pool_p)) {
  #         pool_this <- pool_p
  #       } else {
  #         # empirical pool from observed p-values (excluding pathway genes is ideal; minimal version here)
  #         pool_this <- prep$gene_p_all_use
  #       }
  #       P <- .mvn_z_to_empirical_p(Z, pool_p = pool_this, min_p = prep$min_p)
  #     }

  #     # ACAT/Fisher/minP from P
  #     pA <- vapply(seq_len(B_mvn), function(b) .magcat_acat_combine(P[b, ], min_p = prep$min_p, do_fix = isTRUE(do_fix)), numeric(1))
  #     pF <- stats::pchisq(-2 * rowSums(log(P)), df = 2 * ncol(P), lower.tail = FALSE)

  #     # TFisher adaptive from P
  #     pT <- rep(NA_real_, B_mvn)
  #     if (!requireNamespace("TFisher", quietly = TRUE)) stop("TFisher required for omnibus TFisher null.", call. = FALSE)
  #     for (b in seq_len(B_mvn)) {
  #       pb <- .magcat_fix_p(P[b, ], min_p = prep$min_p, do_fix = isTRUE(do_fix))
  #       if (length(pb) < 2L) next

  #       best <- Inf
  #       for (tau in prep$tau_grid) {
  #         st <- TFisher::stat.soft(p = pb, tau1 = tau)
  #         pv <- 1 - as.numeric(TFisher::p.soft(q = st, n = length(pb), tau1 = tau, M = NULL))
  #         if (is.finite(pv) && pv < best) best <- pv
  #       }
  #       pT[b] <- best
  #     }

  #     pM <- {
  #       k <- ncol(P)
  #       pmin <- apply(P, 1, min)
  #       1 - (1 - pmin)^k
  #     }

  #     # Stouffer from MVN Z (use observed pathway weights if provided)
  #     pS <- rep(NA_real_, B_mvn)
  #     if (!is.null(prep$z_col) && !is.null(prep$z_list)) {
  #       w_i <- NULL
  #       if (!is.null(prep$w_list)) w_i <- prep$w_list[[res$pathway_id[i]]]

  #       # clean weights like wrapper
  #       if (!is.null(w_i)) {
  #         w_i <- suppressWarnings(as.numeric(w_i))
  #         bad <- !is.finite(w_i) | is.na(w_i) | abs(w_i) <= stouffer_min_abs_w
  #         if (any(bad)) {
  #           w_ok <- abs(w_i[!bad])
  #           repl <- if (length(w_ok)) stats::median(w_ok) else 1
  #           w_i[bad] <- repl
  #         }
  #       }

  #       Zs <- if (is.null(w_i)) {
  #         rowSums(Z) / sqrt(ncol(Z))
  #       } else {
  #         as.numeric(Z %*% w_i) / sqrt(sum(w_i^2))
  #       }

  #       if (stouffer_alternative == "greater") {
  #         pS <- stats::pnorm(Zs, lower.tail = FALSE)
  #       } else if (stouffer_alternative == "less") {
  #         pS <- stats::pnorm(Zs, lower.tail = TRUE)
  #       } else {
  #         pS <- 2 * stats::pnorm(-abs(Zs))
  #       }
  #     }

  #     omni_null <- vapply(seq_len(B_mvn), function(b) {
  #       comps <- c(pA[b], pF[b], pT[b], pM[b], pS[b])
  #       if (isTRUE(include_magma_eff) && is.finite(res$magma_pvalue[i])) comps <- c(comps, res$magma_pvalue[i])
  #       .magcat_omni_methods(comps, omnibus = omnibus, min_p = prep$min_p, do_fix = isTRUE(do_fix))
  #     }, numeric(1))

  #     res$omni_p_mvn[i] <- (1 + sum(omni_null <= omni_obs, na.rm = TRUE)) / (B_mvn + 1)
  #   }

  #   res$omni_p_mvn_BH <- .magcat_bh(res$omni_p_mvn)
  # }

  # # final p priority: MVN > global > analytic
  # res$omni_p_final <- res$omni_p_analytic
  # if (B_global > 0L) res$omni_p_final <- res$omni_p_global
  # if (B_mvn > 0L)    res$omni_p_final <- res$omni_p_mvn
  # res$omni_p_final_BH <- .magcat_bh(res$omni_p_final)

  # res <- res[order(res$omni_p_final, decreasing = FALSE, na.last = TRUE), , drop = FALSE]

  # if (isTRUE(output)) {
  #   if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  #   tag <- paste0("omni_", tolower(omnibus),
  #                 if (B_mvn > 0L) "_mvn" else if (B_global > 0L) "_global" else "_analytic")
  #   out_path <- file.path(out_dir, paste0(tag, ".csv"))
  #   utils::write.csv(res, out_path, row.names = FALSE)
  #   attr(res, "file") <- out_path
  # }

  # res

  # ============================================================
  # NEW VERSION- INDIVIDUAL SIMULATION
  # ============================================================

  # ============================================================
  # MVN resampling (LD-aware): option to calibrate COMPONENTS first
  # ============================================================
  res$omni_p_mvn <- NA_real_
  res$omni_p_mvn_BH <- NA_real_

  # NEW: store MVN-calibrated component p-values if requested
  res$acat_p_mvn_cal     <- NA_real_
  res$fisher_p_mvn_cal   <- NA_real_
  res$tfisher_p_mvn_cal  <- NA_real_
  res$minp_p_mvn_cal     <- NA_real_
  res$stouffer_p_mvn_cal <- NA_real_

  # (optional) store the omnibus built from component-calibrated p's (not yet calibrated)
  res$omni_p_mvn_compcal <- NA_real_

  if (B_mvn > 0L) {

    cor_pairs <- .magma_read_gene_cor_pairs(magma_cor_file = magma_cor_file, magma_cor_pairs = magma_cor_pairs)
    if (is.null(cor_pairs) || !nrow(cor_pairs)) stop("MVN resampling requires magma_cor_file or magma_cor_pairs.", call. = FALSE)

    if (!is.null(prep$seed)) set.seed(prep$seed)

    for (i in seq_len(nrow(res))) {

      d <- res$n_genes[i]
      if (!is.finite(d) || d < prep$min_genes) next

      # observed component p's (analytic) already computed above
      pA_obs <- res$acat_p[i]
      pF_obs <- res$fisher_p[i]
      pT_obs <- res$tfisher_p_analytic[i]
      pM_obs <- res$minp_p_analytic[i]
      pS_obs <- res$stouffer_p_analytic[i]

      # if we have literally nothing finite, skip
      if (!any(is.finite(c(pA_obs, pF_obs, pT_obs, pM_obs, pS_obs)))) next

      genes_S <- prep$g_list_raw[[res$pathway_id[i]]]
      genes_S <- as.character(genes_S)
      genes_S <- genes_S[!is.na(genes_S) & genes_S != ""]
      if (length(genes_S) < prep$min_genes) next

      R <- .magma_build_R_from_pairs(genes = genes_S, cor_pairs = cor_pairs, make_PD = make_PD)
      if (is.null(R)) next

      Z <- .magma_simulate_Z_from_R(R, B = B_mvn)
      if (is.null(Z)) next

      # OLD one
      # map MVN Z -> p matrix for p-based methods
      # if (mvn_marginal == "uniform") {
      #   P <- .mvn_z_to_p(Z, min_p = prep$min_p)
      # } else {
      #   pool_this <- if (!is.null(pool_p)) pool_p else prep$gene_p_all_use
      #   P <- .mvn_z_to_empirical_p(Z, pool_p = pool_this, min_p = prep$min_p)
      # }


      # map MVN Z -> p matrix for p-based methods
      if (mvn_marginal == "uniform") {
        P <- .mvn_z_to_p(Z, min_p = prep$min_p)
      } else {

        if (!is.null(pool_p)) {
          # user provided a pool explicitly (no gene IDs to exclude from)
          pool_this <- pool_p
        } else {
          # build a pool from ALL genes EXCEPT the pathway genes (prevents leakage)
          gene_all <- tolower(trimws(as.character(prep$gene_results_core[[prep$gene_col]])))
          genes_S2 <- tolower(trimws(as.character(genes_S)))

          keep <- !(gene_all %in% genes_S2)
          pool_this <- prep$gene_p_all_use[keep]

          # safety: if exclusion leaves you with too few genes, fall back
          if (length(pool_this) < 1000L) {
            pool_this <- prep$gene_p_all_use
          }
        }

        P <- .mvn_z_to_empirical_p(Z, pool_p = pool_this, min_p = prep$min_p)
      }


      # --- NULL component p-values (same “analytic formulas” you already use) ---
      pA_null <- vapply(seq_len(B_mvn), function(b) {
        .magcat_acat_combine(P[b, ], min_p = prep$min_p, do_fix = isTRUE(do_fix))
      }, numeric(1))

      pF_null <- stats::pchisq(-2 * rowSums(log(P)), df = 2 * ncol(P), lower.tail = FALSE)

      pT_null <- rep(NA_real_, B_mvn)
      if (!requireNamespace("TFisher", quietly = TRUE)) stop("TFisher required for omnibus TFisher null.", call. = FALSE)
      for (b in seq_len(B_mvn)) {
        pb <- .magcat_fix_p(P[b, ], min_p = prep$min_p, do_fix = isTRUE(do_fix))
        if (length(pb) < 2L) next
        best <- Inf
        for (tau in prep$tau_grid) {
          st <- TFisher::stat.soft(p = pb, tau1 = tau)
          pv <- 1 - as.numeric(TFisher::p.soft(q = st, n = length(pb), tau1 = tau, M = NULL))
          if (is.finite(pv) && pv < best) best <- pv
        }
        pT_null[b] <- best
      }

      pM_null <- {
        k <- ncol(P)
        pmin <- apply(P, 1, min)
        1 - (1 - pmin)^k
      }

      # stouffer null p-values from MVN Z (if stouffer is enabled)
      pS_null <- rep(NA_real_, B_mvn)
      if (!is.null(prep$z_col) && !is.null(prep$z_list)) {
        w_i <- NULL
        if (!is.null(prep$w_list)) w_i <- prep$w_list[[res$pathway_id[i]]]

        if (!is.null(w_i)) {
          w_i <- suppressWarnings(as.numeric(w_i))
          bad <- !is.finite(w_i) | is.na(w_i) | abs(w_i) <= stouffer_min_abs_w
          if (any(bad)) {
            w_ok <- abs(w_i[!bad])
            repl <- if (length(w_ok)) stats::median(w_ok) else 1
            w_i[bad] <- repl
          }
        }

        Zs <- if (is.null(w_i)) {
          rowSums(Z) / sqrt(ncol(Z))
        } else {
          as.numeric(Z %*% w_i) / sqrt(sum(w_i^2))
        }

        if (stouffer_alternative == "greater") {
          pS_null <- stats::pnorm(Zs, lower.tail = FALSE)
        } else if (stouffer_alternative == "less") {
          pS_null <- stats::pnorm(Zs, lower.tail = TRUE)
        } else {
          pS_null <- 2 * stats::pnorm(-abs(Zs))
        }
      }

      # ============================================================
      # NEW OPTION: MVN-calibrate components first, then build omnibus
      # ============================================================
      if (isTRUE(mvn_calibrate_components)) {

        # observed component -> MVN-calibrated component p (empirical vs null)
        pA_cal_obs <- if (is.finite(pA_obs)) .emp_p_obs_lower(pA_null, pA_obs) else NA_real_
        pF_cal_obs <- if (is.finite(pF_obs)) .emp_p_obs_lower(pF_null, pF_obs) else NA_real_
        pT_cal_obs <- if (is.finite(pT_obs)) .emp_p_obs_lower(pT_null, pT_obs) else NA_real_
        pM_cal_obs <- if (is.finite(pM_obs)) .emp_p_obs_lower(pM_null, pM_obs) else NA_real_
        pS_cal_obs <- if (is.finite(pS_obs)) .emp_p_obs_lower(pS_null, pS_obs) else NA_real_

        res$acat_p_mvn_cal[i]     <- pA_cal_obs
        res$fisher_p_mvn_cal[i]   <- pF_cal_obs
        res$tfisher_p_mvn_cal[i]  <- pT_cal_obs
        res$minp_p_mvn_cal[i]     <- pM_cal_obs
        res$stouffer_p_mvn_cal[i] <- pS_cal_obs

        comps_obs <- c(pA_cal_obs, pF_cal_obs, pT_cal_obs, pM_cal_obs, pS_cal_obs)
        comps_obs <- comps_obs[is.finite(comps_obs) & !is.na(comps_obs)]
        if (!length(comps_obs)) next

        omni_obs <- .magcat_omni_methods(comps_obs, omnibus = omnibus, min_p = prep$min_p, do_fix = isTRUE(do_fix))
        res$omni_p_mvn_compcal[i] <- omni_obs

        # per-draw component-calibrated p vectors (keeps cross-method dependence automatically)
        pA_cal_null <- .emp_p_null_lower(pA_null)
        pF_cal_null <- .emp_p_null_lower(pF_null)
        pT_cal_null <- .emp_p_null_lower(pT_null)
        pM_cal_null <- .emp_p_null_lower(pM_null)
        pS_cal_null <- .emp_p_null_lower(pS_null)

        # build omnibus null from the per-draw calibrated component p’s
        omni_null <- rep(NA_real_, B_mvn)
        for (b in seq_len(B_mvn)) {
          comps_b <- c(pA_cal_null[b], pF_cal_null[b], pT_cal_null[b], pM_cal_null[b], pS_cal_null[b])
          comps_b <- comps_b[is.finite(comps_b) & !is.na(comps_b)]
          if (!length(comps_b)) next
          omni_null[b] <- .magcat_omni_methods(comps_b, omnibus = omnibus, min_p = prep$min_p, do_fix = isTRUE(do_fix))
        }

        res$omni_p_mvn[i] <- (1 + sum(omni_null <= omni_obs, na.rm = TRUE)) / (B_mvn + 1)

      } else {

        # ORIGINAL behavior: calibrate omnibus directly from null analytic component p’s
        omni_obs <- res$omni_p_analytic[i]
        if (!is.finite(omni_obs) || is.na(omni_obs)) next

        omni_null <- vapply(seq_len(B_mvn), function(b) {
          comps <- c(pA_null[b], pF_null[b], pT_null[b], pM_null[b], pS_null[b])
          comps <- comps[is.finite(comps) & !is.na(comps)]
          if (!length(comps)) return(NA_real_)
          .magcat_omni_methods(comps, omnibus = omnibus, min_p = prep$min_p, do_fix = isTRUE(do_fix))
        }, numeric(1))

        res$omni_p_mvn[i] <- (1 + sum(omni_null <= omni_obs, na.rm = TRUE)) / (B_mvn + 1)
      }
    }

    res$omni_p_mvn_BH <- .magcat_bh(res$omni_p_mvn)
  }

  # ============================================================
  # FINALIZE + RETURN (YOU NEED THIS)
  # ============================================================
  res$omni_p_final <- res$omni_p_analytic
  if (B_global > 0L) res$omni_p_final <- res$omni_p_global
  if (B_mvn    > 0L) res$omni_p_final <- res$omni_p_mvn
  res$omni_p_final_BH <- .magcat_bh(res$omni_p_final)

  # more lenient options
  res$omni_p_final_q_storey <- .magcat_qvalue(res$omni_p_final)
  res$omni_p_final_q_bky10  <- .magcat_bky(res$omni_p_final, q = 0.10)

  # IHW using pathway size as a covariate (log scale is typical)
  res$omni_p_final_q_ihw10  <- .magcat_ihw(res$omni_p_final, covar = log1p(res$n_genes))

  

  # ============================================================
  # ADD: genes + gene pvals used (semicolon-separated)  [FIXED]
  # ============================================================

  res$genes_used <- vapply(res$pathway_id, function(pid) {
    g <- prep$g_list_raw[[pid]]
    if (is.null(g) || !length(g)) return(NA_character_)
    paste(g, collapse = ";")
  }, character(1))

  res$gene_pvals_used <- vapply(res$pathway_id, function(pid) {
    p <- prep$p_list[[pid]]   # <-- THIS is the correct per-pathway p vector
    if (is.null(p) || !length(p)) return(NA_character_)
    paste(formatC(p, format = "e", digits = 3), collapse = ";")
  }, character(1))

  
  
  res <- res[order(res$omni_p_final, decreasing = FALSE, na.last = TRUE), , drop = FALSE]

  if (isTRUE(output)) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    tag <- paste0("omni_", tolower(omnibus),
                  if (B_mvn > 0L) "_mvn" else if (B_global > 0L) "_global" else "_analytic")
    out_path <- file.path(out_dir, paste0(tag, ".csv"))
    utils::write.csv(res, out_path, row.names = FALSE)
    attr(res, "file") <- out_path
  }

  # ============================================================
  # DIAGNOSTICS (dominant component, agreement, calibration impact)
  # Add AFTER sorting and BEFORE write.csv()
  # ============================================================

  # Which component p-value columns exist?
  comp_cols <- c("acat_p", "fisher_p", "tfisher_p_analytic", "minp_p_analytic", "stouffer_p_analytic")
  comp_cols <- intersect(comp_cols, names(res))

  # Optional: include MAGMA as a "component" in these diagnostics (ONLY if you want it)
  # (Recommended only if MAGMA is actually being used as a method in your omnibus)
  # if ("magma_pvalue" %in% names(res)) comp_cols <- c(comp_cols, "magma_pvalue")

  # --- dominant component: which method has the smallest p-value (among available methods)
  res$dominant_component <- NA_character_
  if (length(comp_cols) > 0L) {
    pmat <- as.matrix(res[, comp_cols, drop = FALSE])
    res$dominant_component <- apply(pmat, 1, function(x) {
      if (all(is.na(x))) return(NA_character_)
      comp_cols[which.min(replace(x, is.na(x), Inf))]
    })
  }

  # --- agreement score: fraction of available component methods with p < alpha
  alpha <- 0.05
  res$agreement_score <- NA_real_
  res$agreement_k <- NA_integer_   # how many methods significant
  res$agreement_m <- NA_integer_   # how many methods available (non-NA)

  if (length(comp_cols) > 0L) {
    pmat <- as.matrix(res[, comp_cols, drop = FALSE])
    sig  <- (pmat < alpha)

    m_avail <- rowSums(!is.na(pmat))
    k_sig   <- rowSums(sig, na.rm = TRUE)

    res$agreement_m <- as.integer(m_avail)
    res$agreement_k <- as.integer(k_sig)
    res$agreement_score <- ifelse(m_avail > 0, k_sig / m_avail, NA_real_)
  }

  # --- calibration impact: compare analytic omnibus vs the calibrated omnibus you actually report
  # Use omni_p_final (because it's MVN if available, else global, else analytic)
  res$calibration_impact <- NA_character_
  if ("omni_p_analytic" %in% names(res) && "omni_p_final" %in% names(res)) {
    a <- res$omni_p_analytic
    c <- res$omni_p_final

    res$calibration_impact <- ifelse(is.na(a) | is.na(c), NA_character_,
      ifelse(a < alpha & c >= alpha, "Significant_to_nonsig",
        ifelse(a >= alpha & c < alpha, "Nonsig_to_significant",
          ifelse(a < alpha & c < alpha, "Remains_significant", "Remains_nonsignificant")
        )
      )
    )
  }

  return(res)


  # ============================================================
  # NEW VERSION- END
  # ============================================================

}


#' MAGCAT omnibus pathway p-values (6-method framework; optional global + MVN calibration)
#'
#' Compute pathway-level p-values from gene-level results using up to six component
#' methods (ACAT, Fisher, adaptive soft TFisher, minP, Stouffer Z, optional MAGMA
#' competitive) and combine them into a single omnibus p-value per pathway.
#'
#' Resampling / calibration is applied **only to the omnibus** (not to each component),
#' using either:
#' \itemize{
#'   \item \strong{Global resampling} (fast heuristic): resample genes from a global pool.
#'   \item \strong{MVN resampling} (LD-aware): simulate \eqn{Z \sim N(0, R_S)} using pathway-specific
#'         gene–gene correlation matrices built from MAGMA correlation pairs, then convert simulated
#'         Z's into p-values for p-based methods.
#' }
#'
#' @details
#' \strong{Component methods per pathway}
#' \enumerate{
#'   \item \strong{ACAT:} ACAT combine of gene p-values.
#'   \item \strong{Fisher:} Fisher's method on gene p-values.
#'   \item \strong{Adaptive soft TFisher:} computes analytic p-values on a \code{tau_grid} and takes the minimum.
#'   \item \strong{minP:} Sidak-adjusted minimum gene p-value.
#'   \item \strong{Stouffer Z:} combines gene Z's (optionally weighted; supports \code{alternative}).
#'   \item \strong{MAGMA competitive (optional):} uses precomputed competitive gene-set p-values you provide.
#' }
#'
#' \strong{Omnibus across methods per pathway}
#' \itemize{
#'   \item \code{omnibus="ACAT"}: ACAT across method p-values, or
#'   \item \code{omnibus="minP"}: Sidak-adjusted minP across method p-values.
#' }
#'
#' \strong{Which gene-level p-values are used?}
#' P-based methods use \code{p_adj_col} if you provide it \emph{and} it exists in \code{gene_results};
#' otherwise they use \code{p_raw_col}.
#'
#' \strong{MVN marginal mapping (only when \code{perm_mode} includes \code{"mvn"})}
#' \itemize{
#'   \item \code{mvn_marginal="uniform"} (recommended): Gaussian-copula mapping
#'         \eqn{U=\Phi(Z)}, \eqn{p=2\min(U,1-U)} so each marginal is Uniform(0,1).
#'   \item \code{mvn_marginal="empirical"}: maps copula uniforms to an empirical null p distribution
#'         derived internally from gene-level p-values (excluding the pathway's own genes to reduce leakage),
#'         unless you supply \code{pool_p}.
#' }
#'
#' @param gene_results data.frame of gene-level results (one row per gene recommended).
#'   Must contain \code{gene_col} and at least one p-value column (\code{p_raw_col} and/or \code{p_adj_col}).
#'   If you want Stouffer, it must contain \code{z_col}. If you want weighted Stouffer, it must contain
#'   \code{weight_col}.
#'
#' @param pathways Pathway definitions (provide \code{pathways} OR \code{species}).
#'   Either a data.frame with columns \code{pathway_id}, \code{gene_id} (optional \code{pathway_name}),
#'   or a named list mapping \code{pathway_id -> character vector of genes}.
#'
#' @param species Optional character (provide \code{species} OR \code{pathways}), one of
#'   \code{"maize"}, \code{"sorghum"}, \code{"arabidopsis"}, \code{"plant"}.
#'   When provided, pathways are loaded via \code{magcat_load_pathways()}.
#'
#' @param pmn_gene_col Optional character. When \code{species} is used, passed to
#'   \code{magcat_load_pathways(species=..., gene_col=pmn_gene_col)}.
#'
#' @param gene_col Character. Gene ID column in \code{gene_results}. Default \code{"GENE"}.
#' @param p_raw_col Character. Raw p-value column in \code{gene_results}. Default \code{"P"}.
#' @param p_adj_col Character or NULL. Preferred/adjusted p-value column. If supplied and present, used
#'   for ACAT/Fisher/TFisher/minP. Default NULL.
#' @param z_col Character or NULL. Gene Z column for Stouffer. If NULL, Stouffer is skipped. Default NULL.
#' @param weight_col Character or NULL. Optional weights for Stouffer. Default NULL (equal weights).
#'
#' @param tau_grid Numeric vector of \eqn{\tau} values for adaptive soft TFisher. Default is a typical grid.
#' @param tau_cap Numeric. Upper bound for \code{tau_grid} (sanity bound). Default 1.
#' @param min_p Numeric. P-values are clipped to \code{[min_p, 1-min_p]} to avoid numeric issues. Default 1e-15.
#' @param min_genes Integer. Minimum genes required per pathway. Default 2.
#'
#' @param do_fix Logical. If TRUE, uses your \code{fix_p_for_acat()} if available (else clips). Default TRUE.
#' @param stouffer_min_abs_w Numeric. Non-finite / <=0 Stouffer weights replaced by this. Default 1e-8.
#' @param stouffer_alternative One of \code{"greater"}, \code{"two.sided"}, \code{"less"} for Stouffer. Default "greater".
#'
#' @param magma_out Optional data.frame of MAGMA competitive results with columns \code{pathway_id}, \code{magma_pvalue}.
#' @param include_magma_in_omni Logical. Include MAGMA p-values as a component in the omnibus (if provided). Default TRUE.
#' @param include_magma_in_perm Logical. Include MAGMA p-values inside the resampling null. Generally not recommended.
#'   If FALSE and permutations are run, MAGMA is omitted from the permuted omnibus. Default FALSE.
#'
#' @param omnibus Character. \code{"ACAT"} or \code{"minP"} for combining method p-values. Default \code{"ACAT"}.
#'
#' @param B_perm Integer or NULL. Convenience parameter for number of resampling iterations. If provided,
#'   it overrides \code{B_global}/\code{B_mvn} according to \code{perm_mode}.
#' @param perm_mode Character. One of \code{"none"}, \code{"global"}, \code{"mvn"}, \code{"both"}.
#'
#' @param mvn_marginal Character. Only for MVN mode: \code{"uniform"} or \code{"empirical"}.
#' @param mvn_pool Character. Only for MVN + empirical mode when \code{pool_p} is NULL.
#'   \code{"use"} uses the same p-values used for observed p-based methods; \code{"raw"} prefers \code{p_raw_col}.
#'   (If raw p-values are unavailable, falls back to \code{"use"}.)
#' @param pool_p Optional numeric vector of p-values to define the empirical null pool (advanced).
#'   If supplied, it is used directly (no per-pathway exclusion).
#'
#' @param B_global Integer. Number of global-resampling iterations (used when \code{B_perm} is NULL).
#' @param B_mvn Integer. Number of MVN iterations (used when \code{B_perm} is NULL).
#' @param perm_pool Character. For global resampling only: \code{"raw"} (use raw p-values) or \code{"obs"} (use p used by methods).
#'
#' @param magma_cor_file Character or NULL. Path to 3-col correlation pairs file: \code{gene1 gene2 r}.
#'   Required for MVN unless \code{magma_cor_pairs} is provided.
#' @param magma_cor_pairs data.frame or NULL. In-memory correlation pairs with columns \code{gene1}, \code{gene2}, \code{r}.
#' @param make_PD Logical. If TRUE, forces correlation matrices to be positive definite before Cholesky. Default TRUE.
#' @param mvn_calibrate_components Logical. If TRUE, MVN-calibrates each component first, then builds the omnibus
#'   (and calibrates that omnibus). Default FALSE.
#'
#' @param seed Integer or NULL. Random seed.
#' @param output Logical. If TRUE, writes CSV to \code{out_dir} and attaches \code{"file"} attribute. Default FALSE.
#' @param out_dir Character. Output directory. Default \code{"magcat_omni2"}.
#'
#' @return data.frame with one row per pathway, including:
#' \itemize{
#'   \item component p-values (analytic): \code{acat_p}, \code{fisher_p}, \code{tfisher_p_analytic}, \code{minp_p_analytic}, \code{stouffer_p_analytic}
#'   \item optional \code{magma_pvalue}
#'   \item \code{omni_p_analytic}, \code{omni_p_global}, \code{omni_p_mvn}, and \code{omni_p_final} (+ BH/qvalue extras if you kept them in \code{magcat_omni2_run})
#' }
#'
#' @examples
#' \dontrun{
#' # Load gene results from MAGMA output
#' gene_results <- read.delim("magma_output.genes.out")
#'
#' # Run omnibus test with analytic p-values only (no permutation)
#' omni_res <- magcat_omni2_pathways(
#'   gene_results = gene_results,
#'   species = "maize",
#'   gene_col = "GENE",
#'   p_raw_col = "P",
#'   z_col = "ZSTAT",
#'   perm_mode = "none"
#' )
#' head(omni_res)
#'
#' # With MVN-based calibration (requires gene-gene correlations)
#' omni_mvn <- magcat_omni2_pathways(
#'   gene_results = gene_results,
#'   species = "maize",
#'   perm_mode = "mvn",
#'   B_perm = 10000,
#'   magma_cor_file = "magma_output.genes.raw"
#' )
#'
#' # With global resampling
#' omni_global <- magcat_omni2_pathways(
#'   gene_results = gene_results,
#'   species = "maize",
#'   perm_mode = "global",
#'   B_perm = 5000
#' )
#' }
#'
#' @seealso \code{\link{magcat_acat_pathways}}, \code{\link{magcat_fisher_pathways}},
#'   \code{\link{magcat_minp_pathways}}, \code{\link{magcat_stoufferZ_pathways}}
#' @export
magcat_omni2_pathways <- function(gene_results,
                                 pathways = NULL,
                                 species = NULL,
                                 pmn_gene_col = NULL,
                                 gene_col = "GENE",
                                 p_raw_col = "P",
                                 p_adj_col = NULL,
                                 z_col = NULL,
                                 weight_col = NULL,
                                 tau_grid = c(0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.001),
                                 tau_cap = 1,
                                 min_p = 1e-15,
                                 min_genes = 2L,
                                 # component knobs
                                 do_fix = TRUE,
                                 stouffer_min_abs_w = 1e-8,
                                 stouffer_alternative = c("greater","two.sided","less"),
                                 # MAGMA competitive (optional)
                                 magma_out = NULL,
                                 include_magma_in_omni = TRUE,
                                 include_magma_in_perm = FALSE,
                                 # omnibus
                                 omnibus = c("ACAT","minP"),
                                 # resampling
                                 B_perm = NULL,
                                 perm_mode = c("none","global","mvn","both"),
                                 mvn_marginal = c("uniform","empirical"),
                                 mvn_pool     = c("use","raw"),
                                 pool_p       = NULL,
                                 mvn_calibrate_components = FALSE,
                                 # back-compat / low-level
                                 B_global = 0L,
                                 B_mvn = 0L,
                                 perm_pool = c("raw","obs"),
                                 magma_cor_file = NULL,
                                 magma_cor_pairs = NULL,
                                 make_PD = TRUE,
                                 seed = NULL,
                                 output = FALSE,
                                 out_dir = "magcat_omni2") {

  # ---- basic arg matching / validation ----
  omnibus   <- match.arg(omnibus)
  perm_mode <- match.arg(perm_mode)
  perm_pool <- match.arg(perm_pool)
  mvn_marginal <- match.arg(mvn_marginal)
  mvn_pool     <- match.arg(mvn_pool)
  stouffer_alternative <- match.arg(stouffer_alternative)

  if (!is.null(species)) {
    species <- match.arg(species, c("maize","sorghum","arabidopsis","plant","fly"))
  }

  min_p <- as.numeric(min_p)
  if (!is.finite(min_p) || is.na(min_p) || min_p <= 0 || min_p >= 0.5) {
    stop("min_p must be in (0, 0.5).", call. = FALSE)
  }

  min_genes <- as.integer(min_genes)
  if (!is.finite(min_genes) || is.na(min_genes) || min_genes < 1L) {
    stop("min_genes must be a positive integer.", call. = FALSE)
  }

  # ---- map B_perm + perm_mode -> B_global/B_mvn (B_perm overrides) ----
  if (!is.null(B_perm)) {
    B_perm <- as.integer(B_perm)
    if (!is.finite(B_perm) || is.na(B_perm) || B_perm < 0L) stop("B_perm must be >= 0.", call. = FALSE)

    if (perm_mode == "none" || B_perm == 0L) {
      B_global <- 0L; B_mvn <- 0L
    } else if (perm_mode == "global") {
      B_global <- B_perm; B_mvn <- 0L
    } else if (perm_mode == "mvn") {
      B_global <- 0L; B_mvn <- B_perm
    } else { # both
      B_global <- B_perm; B_mvn <- B_perm
    }
  } else {
    # respect perm_mode even if user sets B_global/B_mvn
    B_global <- as.integer(B_global); if (is.na(B_global) || B_global < 0L) B_global <- 0L
    B_mvn    <- as.integer(B_mvn);    if (is.na(B_mvn)    || B_mvn    < 0L) B_mvn <- 0L

    if (perm_mode == "none") {
      B_global <- 0L; B_mvn <- 0L
    } else if (perm_mode == "global") {
      B_mvn <- 0L
    } else if (perm_mode == "mvn") {
      B_global <- 0L
    } else {
      # both: keep as-is
    }
  }

  # MVN requirements
  if (B_mvn > 0L) {
    if (is.null(magma_cor_pairs) && is.null(magma_cor_file)) {
      stop("perm_mode includes 'mvn' (B_mvn > 0) but no gene-correlation input was provided.\n",
           "Provide magma_cor_file (path to 3-col gene1 gene2 r) or magma_cor_pairs (data.frame).",
           call. = FALSE)
    }
  }

  # Stouffer requirements (skip silently if z_col is NULL)
  if (is.null(z_col) && !is.null(weight_col)) {
    warning("weight_col provided but z_col is NULL; Stouffer is skipped and weight_col is ignored.", call. = FALSE)
    weight_col <- NULL
  }

  # ---- prepare (standardizes gene_results + pathways) ----
  if (!exists("magcat_omni2_prepare", mode = "function")) {
    stop("magcat_omni2_prepare() not found. Source your MAGCAT omnibus code first.", call. = FALSE)
  }
  if (!exists("magcat_omni2_run", mode = "function")) {
    stop("magcat_omni2_run() not found. Source your MAGCAT omnibus code first.", call. = FALSE)
  }

  prep <- magcat_omni2_prepare(
    gene_results = gene_results,
    pathways = pathways,
    species = species,
    pmn_gene_col = pmn_gene_col,
    gene_col = gene_col,
    p_raw_col = p_raw_col,
    p_adj_col = p_adj_col,
    z_col = z_col,
    weight_col = weight_col,
    tau_grid = tau_grid,
    tau_cap = tau_cap,
    min_p = min_p,
    min_genes = min_genes,
    seed = seed
  )

  # ---- MVN empirical pool choice (only when pool_p is NULL) ----
  # If you want "raw" empirical pool and raw is available, we swap the internal pool used by MVN mapping.
  # This preserves the per-pathway exclusion logic inside magcat_omni2_run().
  if (B_mvn > 0L && mvn_marginal == "empirical" && is.null(pool_p) && mvn_pool == "raw") {
    if (!is.null(prep$gene_p_all_raw)) {
      prep$gene_p_all_use <- prep$gene_p_all_raw
    } else {
      warning("mvn_pool='raw' requested but raw p-values are unavailable; falling back to mvn_pool='use'.", call. = FALSE)
    }
  }

  # ---- run (components analytic; permutations ONLY for omnibus) ----
  magcat_omni2_run(
    prep = prep,
    omnibus = omnibus,
    do_fix = do_fix,
    stouffer_min_abs_w = stouffer_min_abs_w,
    stouffer_alternative = stouffer_alternative,
    magma_out = magma_out,
    include_magma_in_omni = include_magma_in_omni,
    include_magma_in_perm = include_magma_in_perm,
    B_global = B_global,
    B_mvn = B_mvn,
    perm_pool = perm_pool,
    mvn_marginal = mvn_marginal,
    mvn_pool = mvn_pool,
    pool_p = pool_p,
    mvn_calibrate_components = mvn_calibrate_components,
    magma_cor_file = magma_cor_file,
    magma_cor_pairs = magma_cor_pairs,
    make_PD = make_PD,
    output = output,
    out_dir = out_dir
  )
}