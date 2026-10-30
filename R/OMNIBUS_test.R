## ============================================================
##  CATFISH omnibus (6-method, with global resampling + Westfall-Young)
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
##  Calibration (omnibus only):
##    - global resampling: resample genes from a global pool (fast heuristic)
##    - Westfall-Young: one genome-wide null bank of correlated gene z's
##      (MVN by chromosome from MAGMA gene-gene correlations), read per pathway
##      (calibrated p) and across pathways (family-wise min-P adjusted p).
##      Engine in R/westfall_young.R.
##
#' CATFISH omnibus pathway p-values (internal run step)
#'
#' Computes pathway-level p-values from a \code{catfish_omni2_prepare()} object,
#' combines components into an omnibus, and optionally calibrates the omnibus by
#' global resampling and/or Westfall--Young. See \code{catfish_omni2_pathways()}
#' for the full description of arguments and output columns.
#'
#' @keywords internal
#' @noRd
# ============================================================
# CATFISH omnibus (fast, 6-method; global + Westfall-Young calibration for OMNIBUS only)
#   - Component methods are computed by YOUR CATFISH wrappers (analytic only)
#   - Permutations/resampling are ONLY for the OMNIBUS
# ============================================================

# ----------------------------
# Helper: safely grab CATFISH functions (exported wrappers)
# Note: These functions are all in the same package, so we just use get()
# ----------------------------
.catfish_get_fun <- function(fname) {
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
.mvn_z_to_p <- function(Z, min_p = 1e-15, sided = c("one", "two")) {
  sided <- match.arg(sided)
  if (sided == "two") {
    U <- stats::pnorm(Z)
    P <- 2 * pmin(U, 1 - U)
  } else {
    # one-sided upper tail: matches MAGMA gene-level p-values (large Z = strong signal)
    P <- stats::pnorm(Z, lower.tail = FALSE)
  }
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

.mvn_z_to_empirical_p <- function(Z, pool_p, min_p = 1e-15, sided = c("one", "two")) {
  sided <- match.arg(sided)
  q0 <- .empirical_qfun(pool_p)
  # one-sided upper tail so large Z maps to small p (consistent with MAGMA);
  # two-sided folds the copula uniform as before.
  U  <- if (sided == "two") stats::pnorm(Z) else stats::pnorm(Z, lower.tail = FALSE)
  P  <- matrix(q0(as.vector(U)), nrow = nrow(Z), ncol = ncol(Z))
  P  <- pmax(pmin(P, 1 - min_p), min_p)
  P
}

# ----------------------------
# p-fixing (uses your fix_p_for_acat if available)
# ----------------------------
.catfish_fix_p <- function(p, min_p = 1e-50, do_fix = TRUE) {
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

.catfish_bh <- function(p) {
  p <- as.numeric(p)
  out <- rep(NA_real_, length(p))
  idx <- which(is.finite(p) & !is.na(p))
  if (length(idx)) out[idx] <- stats::p.adjust(p[idx], method = "BH")
  out
}

## ----------------------------
## NEW tests
## ----------------------------

.catfish_qvalue <- function(p) {
  p <- as.numeric(p)
  out <- rep(NA_real_, length(p))
  idx <- which(is.finite(p) & !is.na(p))

  if (!length(idx)) return(out)

  if (requireNamespace("qvalue", quietly = TRUE)) {
    qobj <- tryCatch(
      qvalue::qvalue(p[idx]),
      error = function(e) NULL
    )
    if (!is.null(qobj)) {
      out[idx] <- as.numeric(qobj$qvalues)
    } else {
      # qvalue failed (e.g., pi0 <= 0); fall back to BH
      out[idx] <- stats::p.adjust(p[idx], method = "BH")
    }
  } else {
    warning("Package 'qvalue' not installed; falling back to BH.", call. = FALSE)
    out[idx] <- stats::p.adjust(p[idx], method = "BH")
  }
  out
}

# BKY / two-stage BH (Benjamini-Krieger-Yekutieli style)
.catfish_bky <- function(p, q = 0.10) {
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

.catfish_ihw <- function(p, covar) {
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


# ----------------------------
# Omnibus across METHODS (ACAT or Sidak-minP)
# ----------------------------
.catfish_acat_combine <- function(p, min_p = 1e-50, do_fix = TRUE) {
  p <- .catfish_fix_p(p, min_p = min_p, do_fix = do_fix)
  p <- p[p > 0 & p < 1]
  if (!length(p)) return(NA_real_)
  tstat <- mean(tan((0.5 - p) * pi))
  out <- 0.5 - atan(tstat) / pi
  out <- max(min(out, 1 - min_p), min_p)
  as.numeric(out)
}

.catfish_minp_sidak <- function(p, min_p = 1e-50, do_fix = TRUE) {
  p <- .catfish_fix_p(p, min_p = min_p, do_fix = do_fix)
  p <- p[p > 0 & p < 1]
  k <- length(p)
  if (k < 1L) return(NA_real_)
  pmin <- min(p)
  1 - (1 - pmin)^k
}

.catfish_sanitize_tau_grid <- function(tau_grid,
                                       p_pool = NULL,
                                       tau_cap = 1,
                                       allow_tau_one = FALSE,
                                       fallback = NULL,
                                       sort_decreasing = FALSE) {
  if (is.character(tau_grid) && length(tau_grid) == 1L && tolower(tau_grid) == "auto") {
    tau_grid <- .catfish_auto_tau_grid(
      p_pool = p_pool,
      tau_cap = tau_cap,
      allow_tau_one = allow_tau_one
    )
  }

  tg <- suppressWarnings(as.numeric(tau_grid))
  tg <- tg[is.finite(tg) & !is.na(tg)]
  if (isTRUE(allow_tau_one)) {
    tg <- tg[tg > 0 & tg <= tau_cap]
  } else {
    tg <- tg[tg > 0 & tg < 1 & tg <= tau_cap]
  }
  tg <- sort(unique(tg), decreasing = sort_decreasing)

  if (!length(tg) && !is.null(fallback)) {
    fb <- suppressWarnings(as.numeric(fallback))
    fb <- fb[is.finite(fb) & !is.na(fb)]
    if (isTRUE(allow_tau_one)) {
      fb <- fb[fb > 0 & fb <= tau_cap]
    } else {
      fb <- fb[fb > 0 & fb < 1 & fb <= tau_cap]
    }
    tg <- sort(unique(fb), decreasing = sort_decreasing)
  }

  tg
}

.catfish_tfisher_adaptive_p <- function(p, tau_grid, min_p = 1e-50, do_fix = TRUE) {
  if (!requireNamespace("TFisher", quietly = TRUE)) return(NA_real_)

  pb <- .catfish_fix_p(p, min_p = min_p, do_fix = do_fix)
  pb <- pb[is.finite(pb) & !is.na(pb)]
  if (length(pb) < 2L) return(NA_real_)

  tg <- .catfish_sanitize_tau_grid(
    tau_grid = tau_grid,
    p_pool = pb,
    tau_cap = 1,
    allow_tau_one = FALSE,
    sort_decreasing = FALSE
  )
  if (!length(tg)) return(NA_real_)

  out <- tryCatch({
    omni_out <- TFisher::stat.soft.omni(p = pb, TAU1 = tg, M = NULL)
    Wo <- as.numeric(omni_out$omni)
    as.numeric(TFisher::p.soft.omni(q = Wo, n = length(pb), TAU1 = tg, M = NULL))
  }, error = function(e) {
    NA_real_
  })

  if (!is.finite(out) || is.na(out)) return(NA_real_)
  max(min(out, 1 - min_p), min_p)
}

.catfish_omni_methods <- function(p_methods, omnibus = c("ACAT", "minP"), min_p = 1e-50, do_fix = TRUE) {
  omnibus <- match.arg(omnibus)
  pv <- .catfish_fix_p(p_methods, min_p = min_p, do_fix = do_fix)
  pv <- pv[pv > 0 & pv < 1]
  if (!length(pv)) return(NA_real_)

  if (omnibus == "ACAT") {
    .catfish_acat_combine(pv, min_p = min_p, do_fix = do_fix)
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


# NEW VERSION
.emp_p_null_lower <- function(p_null) {
  p2 <- as.numeric(p_null)
  B <- length(p2)
  if (B < 1L) return(numeric(0))
  bad <- !is.finite(p2) | is.na(p2)
  if (any(bad)) p2[bad] <- Inf

  # Leave-one-out plus-one empirical p for each null draw on the same
  # tail-probability scale used by .emp_p_obs_lower().
  rank_max <- rank(p2, ties.method = "max")
  as.numeric(rank_max) / B
}


# ----------------------------
# GPD tail extrapolation helpers (Knijnenburg et al. 2009)
# ----------------------------

#' Convert p-value to score for GPD tail fitting
#' @noRd
.tail_score_from_p <- function(p, floor_p = 1e-300) {
  -log10(pmax(as.numeric(p), floor_p))
}

#' Hybrid empirical + GPD tail p-value estimator
#'
#' Uses standard plus-one Monte Carlo when exceedances are sufficient,
#' otherwise fits a GPD to the extreme tail for extrapolation.
#'
#' @param p_null Numeric vector of null p-values
#' @param p_obs Observed p-value
#' @param method "empirical" (standard) or "hybrid_gpd" (GPD tail extrapolation)
#' @param switch_exceed If >= this many null draws beat observed, use empirical
#' @param gpd_k Target number of extreme null draws for GPD fit
#' @param min_B Minimum null draws required for GPD

#' @param min_tail Minimum tail points required for GPD fit
#' @param floor_p Floor for p-values to avoid log(0)
#' @return Calibrated p-value
#' @noRd
.emp_or_gpd_p_lower <- function(p_null, p_obs,
                                method = c("empirical", "hybrid_gpd"),
                                switch_exceed = 10L,
                                gpd_k = 250L,
                                min_B = 10000L,
                                min_tail = 50L,
                                floor_p = 1e-300) {
  method <- match.arg(method)

  p_null <- as.numeric(p_null)
  p_obs  <- as.numeric(p_obs)

  p_null <- p_null[is.finite(p_null) & !is.na(p_null)]
  B <- length(p_null)

  if (!is.finite(p_obs) || is.na(p_obs) || B < 1L) return(NA_real_)

  # always compute the standard plus-one empirical p-value
  n_exc <- sum(p_null <= p_obs, na.rm = TRUE)
  p_emp <- (1 + n_exc) / (B + 1)

  if (method == "empirical") return(p_emp)

  # only switch to GPD when the empirical tail is too sparse
  if (B < min_B || n_exc >= switch_exceed) return(p_emp)

  # requires a GPD fitter; use evir package
  if (!requireNamespace("evir", quietly = TRUE)) {
    warning("Package 'evir' not installed; falling back to empirical p-value.", call. = FALSE)
    return(p_emp)
  }

  s_null <- .tail_score_from_p(p_null, floor_p = floor_p)
  s_obs  <- .tail_score_from_p(p_obs,  floor_p = floor_p)

  s_null <- s_null[is.finite(s_null) & !is.na(s_null)]
  if (length(s_null) < min_B) return(p_emp)

  # choose threshold so that roughly the top-k null scores are exceedances
  s_ord <- sort(s_null, decreasing = TRUE)
  k <- min(as.integer(gpd_k), floor(0.10 * B))
  k <- max(k, min_tail)
  if ((k + 1L) >= B) return(p_emp)

  u <- s_ord[k + 1L]
  n_u <- sum(s_null > u)
  if (n_u < min_tail) return(p_emp)

  # if observed score is not beyond the threshold, empirical is enough
  if (s_obs <= u) return(p_emp)

  fit <- tryCatch(
    evir::gpd(s_null, threshold = u, method = "ml"),
    error = function(e) NULL
  )
  if (is.null(fit)) return(p_emp)

  xi   <- unname(fit$par.ests["xi"])
  beta <- unname(fit$par.ests["beta"])

  if (!is.finite(xi) || !is.finite(beta) || beta <= 0) return(p_emp)

  y <- s_obs - u

  tail_cond <- if (abs(xi) < 1e-8) {
    exp(-y / beta)
  } else {
    z <- 1 + xi * y / beta
    if (z <= 0) return(p_emp)
    z^(-1 / xi)
  }

  p_u   <- n_u / B
  p_gpd <- p_u * tail_cond

  p_gpd <- max(min(p_gpd, 1), floor_p)

  # use the extrapolated tail estimate
  p_gpd
}



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



# Function with shinkage
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

  # Safety: ensure R has no NA/Inf

  R[!is.finite(R)] <- 0
  diag(R) <- 1

  U <- tryCatch(chol(R), error = function(e) NULL)
  if (is.null(U)) stop("chol(R) failed even after PD-fix. Check correlation input.", call. = FALSE)

  Z0 <- matrix(stats::rnorm(B * d), nrow = B, ncol = d)
  Z  <- Z0 %*% U
  # Safety: clamp any Inf Z values
  Z[!is.finite(Z)] <- 0
  colnames(Z) <- colnames(R)
  Z
}

# ----------------------------
# Global sampler helper (WITHOUT replacement per permutation)
# ----------------------------
.catfish_sample_idx_mat <- function(n_pool, d, B) {
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
# Stouffer p from Z matrix helper (for global/mvn_global resampling)
# ----------------------------
.stouffer_p_from_Zmat_helper <- function(Zmat, w = NULL, alternative = "greater",
                                         min_abs_w = 1e-8) {
  if (is.null(Zmat) || !ncol(Zmat)) return(rep(NA_real_, nrow(Zmat)))
  if (is.null(w)) {
    Zs <- rowSums(Zmat) / sqrt(ncol(Zmat))
  } else if (is.matrix(w)) {
    if (!all(dim(w) == dim(Zmat))) {
      stop("When w is a matrix, it must have the same dimensions as Zmat.", call. = FALSE)
    }
    W <- matrix(suppressWarnings(as.numeric(w)), nrow = nrow(w), ncol = ncol(w))
    bad <- !is.finite(W) | is.na(W) | abs(W) <= min_abs_w
    if (any(bad)) {
      W_abs <- abs(W)
      W_abs[bad] <- NA_real_
      row_repl <- apply(W_abs, 1, function(x) {
        x <- x[is.finite(x) & !is.na(x)]
        if (length(x)) stats::median(x) else 1
      })
      W[bad] <- row_repl[row(W)][bad]
    }
    denom <- sqrt(rowSums(W^2))
    denom[!is.finite(denom) | is.na(denom) | denom <= 0] <- 1
    Zs <- rowSums(Zmat * W) / denom
  } else {
    w <- suppressWarnings(as.numeric(w))
    bad <- !is.finite(w) | is.na(w) | abs(w) <= min_abs_w
    if (any(bad)) {
      w_ok <- abs(w[!bad])
      repl <- if (length(w_ok)) stats::median(w_ok) else 1
      w[bad] <- repl
    }
    denom <- sqrt(sum(w^2))
    Zs <- as.numeric(Zmat %*% w) / denom
  }

  if (alternative == "greater") {
    p_out <- stats::pnorm(Zs, lower.tail = FALSE)
  } else if (alternative == "less") {
    p_out <- stats::pnorm(Zs, lower.tail = TRUE)
  } else {
    p_out <- 2 * stats::pnorm(-abs(Zs))
  }
  pmax(pmin(p_out, 1 - 1e-15), 1e-15)
}

# ----------------------------
# Single-pathway global resampling helper (for mvn_global fallback)
# ----------------------------
.do_global_resampling_pathway <- function(d, B, pool_p, pool_z, pool_w,
                                          omni_obs, tau_grid, min_p, do_fix,
                                          omnibus, stouffer_alternative,
                                          stouffer_min_abs_w, include_magma,
                                          magma_pval = NA_real_,
                                          tail_mode = "empirical",
                                          tail_switch_exceed = 10L,
                                          tail_gpd_k = 250L,
                                          tail_min_B = 10000L,
                                          tail_min_tail = 50L) {
  n_pool <- length(pool_p)
  idx_mat <- .catfish_sample_idx_mat(n_pool = n_pool, d = d, B = B)
  if (is.null(idx_mat)) return(NA_real_)

  P <- matrix(pool_p[idx_mat], nrow = B, ncol = d)

  # ACAT across genes
  pA <- vapply(seq_len(B), function(b) {
    .catfish_acat_combine(P[b, ], min_p = min_p, do_fix = do_fix)
  }, numeric(1))

  # Fisher across genes
  pF <- stats::pchisq(-2 * rowSums(log(P)), df = 2 * d, lower.tail = FALSE)

  # TFisher adaptive
  pT <- rep(NA_real_, B)
  if (requireNamespace("TFisher", quietly = TRUE)) {
    for (b in seq_len(B)) {
      pT[b] <- .catfish_tfisher_adaptive_p(
        P[b, ],
        tau_grid = tau_grid,
        min_p = min_p,
        do_fix = do_fix
      )
    }
  }

  # minP Sidak
  pM <- {
    pmin_vec <- apply(P, 1, min)
    1 - (1 - pmin_vec)^d
  }

  # Stouffer
  pS <- rep(NA_real_, B)
  if (!is.null(pool_z)) {
    Zmat <- matrix(pool_z[idx_mat], nrow = B, ncol = d)
    Wmat <- NULL
    if (!is.null(pool_w)) Wmat <- matrix(pool_w[idx_mat], nrow = B, ncol = d)
    pS <- .stouffer_p_from_Zmat_helper(Zmat, w = Wmat,
                                       alternative = stouffer_alternative,
                                       min_abs_w = stouffer_min_abs_w)
  }

  # Build omnibus null
  omni_null <- vapply(seq_len(B), function(b) {
    comps <- c(pA[b], pF[b], pT[b], pM[b], pS[b])
    if (isTRUE(include_magma) && is.finite(magma_pval)) comps <- c(comps, magma_pval)
    .catfish_omni_methods(comps, omnibus = omnibus, min_p = min_p, do_fix = do_fix)
  }, numeric(1))

  # Empirical or GPD-extrapolated p-value
  .emp_or_gpd_p_lower(
    p_null = omni_null,
    p_obs  = omni_obs,
    method = tail_mode,
    switch_exceed = tail_switch_exceed,
    gpd_k = tail_gpd_k,
    min_B = tail_min_B,
    min_tail = tail_min_tail
  )
}

# ----------------------------
# Data-driven (auto) tau grid for adaptive soft TFisher
# ----------------------------
# Picks soft-TFisher thresholds from the empirical quantiles of the pooled
# gene-level p-values. Under a uniform (null) p distribution this reproduces the
# classic c(0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.001) grid; under strong signal
# the quantiles shrink toward the extreme tail, giving a strict-like grid with no
# manual toggle. The grid is estimated ONCE from the global pool and then held
# fixed, so it acts as a fixed transform for the MVN null calibration.
.catfish_auto_tau_grid <- function(p_pool,
                                   probs = c(0.20, 0.10, 0.05, 0.02, 0.01, 0.005, 0.001),
                                   tau_floor = 1e-8,
                                   tau_cap = 0.5,
                                   allow_tau_one = FALSE,
                                   min_n = 50L) {
  p <- suppressWarnings(as.numeric(p_pool))
  p <- p[is.finite(p) & !is.na(p) & p > 0 & p < 1]
  tau_upper <- if (isTRUE(allow_tau_one)) tau_cap else min(tau_cap, 1 - sqrt(.Machine$double.eps))

  # too few genes to estimate quantiles: fall back to the classic fixed grid
  if (length(p) < min_n) {
    g <- c(0.20, 0.10, 0.05, 0.02, 0.01, 0.005, 0.001)
    g <- g[g > 0 & g <= tau_upper]
    return(sort(unique(g), decreasing = FALSE))
  }

  q <- as.numeric(stats::quantile(p, probs = probs, na.rm = TRUE, names = FALSE))
  q <- q[is.finite(q)]
  q <- pmin(pmax(q, tau_floor), tau_upper)   # keep numerically safe + valid
  q <- sort(unique(q), decreasing = FALSE)

  # guarantee at least two distinct thresholds
  if (length(q) < 2L) {
    hi <- min(tau_upper, max(q, 0.05))
    lo <- max(tau_floor, hi / 100)
    q  <- sort(unique(c(hi, lo)), decreasing = FALSE)
  }
  q
}

# ----------------------------
# prepare (fast reuse)
# ----------------------------
catfish_omni2_prepare <- function(gene_results,
                                pathways = NULL,
                                species = NULL,
                                pmn_gene_col = NULL,
                                gene_col = "GENE",
                                p_raw_col = "P",
                                p_adj_col = NULL,
                                z_col = NULL,
                                weight_col = NULL,
                                chr_col = NULL,
                                tau_grid = "auto",
                                tau_cap = 1,
                                min_p = 1e-50,
                                min_genes = 2L,
                                max_genes = Inf,
                                seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  if (!is.data.frame(gene_results)) stop("gene_results must be a data.frame.", call. = FALSE)
  if (!(gene_col %in% names(gene_results))) stop("gene_results missing gene_col: ", gene_col, call. = FALSE)

  if (!is.null(pathways) && !is.null(species)) stop("Provide either pathways OR species, not both.", call. = FALSE)
  if (is.null(pathways) && is.null(species)) stop("Need pathways (list/df) or species.", call. = FALSE)

  # load pathways if species (use wrapper helper from this package)
  if (is.null(pathways)) {
    # catfish_load_pathways is defined in this package (ACAT_wrappers.R)
    if (exists("catfish_load_pathways", mode = "function")) {
      pathways <- catfish_load_pathways(species = species, gene_col = pmn_gene_col)
    } else {
      stop("Could not find catfish_load_pathways(). This should not happen.", call. = FALSE)
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

  # chromosome aligned (needed by the Westfall-Young null: genes are drawn
  # jointly within a chromosome and independently across chromosomes)
  chr_use <- NULL
  if (!is.null(chr_col)) {
    if (!(chr_col %in% names(gene_results))) stop("chr_col not found: ", chr_col, call. = FALSE)
    chr_use <- as.character(gene_results[[chr_col]])[okg][okp]
  }

  # dedup by normalized gene id (keep first)
  keep1 <- !duplicated(g_norm)
  g_raw  <- g_raw[keep1]
  g_norm <- g_norm[keep1]
  p_use  <- p_use[keep1]
  if (!is.null(p_raw_vec)) p_raw_vec <- p_raw_vec[keep1]
  if (!is.null(z_use)) z_use <- z_use[keep1]
  if (!is.null(w_use)) w_use <- w_use[keep1]
  if (!is.null(chr_use)) chr_use <- chr_use[keep1]

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
  max_genes <- suppressWarnings(as.numeric(max_genes))
  if (length(max_genes) == 1L && is.finite(max_genes)) {
    g_list_norm <- g_list_norm[lengths(g_list_norm) <= max_genes]
  }
  if (!length(g_list_norm)) stop("No pathway passes the min_genes/max_genes filter.", call. = FALSE)

  pid <- names(g_list_norm)

  # per pathway vectors
  p_list <- lapply(g_list_norm, function(g) unname(gene_p[g]))
  z_list <- if (!is.null(gene_z)) lapply(g_list_norm, function(g) unname(gene_z[g])) else NULL
  w_list <- if (!is.null(gene_w)) lapply(g_list_norm, function(g) unname(gene_w[g])) else NULL
  g_list_raw <- lapply(g_list_norm, function(g) unname(gene_map[g]))
  names(g_list_raw) <- names(g_list_norm)  

  # Adaptive TFisher excludes tau = 1 because Fisher is already a separate
  # CATFISH component; keeping both duplicates the untruncated Fisher corner.
  tau_grid <- .catfish_sanitize_tau_grid(
    tau_grid = tau_grid,
    p_pool = p_use,
    tau_cap = tau_cap,
    allow_tau_one = FALSE,
    fallback = 0.05,
    sort_decreasing = FALSE
  )

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
      gene_id_all    = g_norm,
      gene_chr_all   = chr_use,

      gene_col   = gene_col,
      p_use_col  = p_use_col,
      p_raw_col  = p_raw_col,
      p_adj_col  = p_adj_col,
      z_col      = z_col,
      weight_col = weight_col,
      chr_col    = chr_col,
      tau_grid   = tau_grid,
      tau_cap    = tau_cap,
      min_p      = min_p,
      min_genes  = min_genes,
      max_genes  = max_genes,
      seed       = seed,

      gene_results_core  = gene_results_core,
      pathways_list_raw  = pathways_list_raw
    ),
    class = "catfish_omni2_prep"
  )
}

# ----------------------------
# Run (components from YOUR wrappers; perms only for omnibus)
# ----------------------------
catfish_omni2_run <- function(prep,
                            omnibus = c("ACAT","minP"),
                            # component wrapper params
                            use_tfisher = TRUE,
                            do_fix = TRUE,
                            stouffer_min_abs_w = 1e-8,
                            stouffer_alternative = c("greater","two.sided","less"),
                            # magma competitive
                            magma_out = NULL,
                            include_magma_in_omni = TRUE,
                            include_magma_in_perm = FALSE,
                            # omnibus calibration
                            B_global = 0L,
                            B_wy = 0L,
                            perm_mode = c("none","global","wy","both"),
                            perm_pool = c("raw","obs"),
                            wy_tau_grid = c(0.01, 0.05, 0.5),
                            magma_cor_file = NULL,
                            magma_cor_pairs = NULL,
                            # GPD tail extrapolation (Knijnenburg et al. 2009)
                            tail_mode = c("empirical", "hybrid_gpd"),
                            tail_switch_exceed = 10L,
                            tail_gpd_k = 250L,
                            tail_min_B = 10000L,
                            tail_min_tail = 50L,
                            output = FALSE,
                            out_dir = "catfish_omni2") {

  if (!inherits(prep, "catfish_omni2_prep")) stop("prep must come from catfish_omni2_prepare().", call. = FALSE)
  omnibus   <- match.arg(omnibus)
  perm_mode <- match.arg(perm_mode)
  perm_pool <- match.arg(perm_pool)
  stouffer_alternative <- match.arg(stouffer_alternative)
  tail_mode <- match.arg(tail_mode)

  B_global <- as.integer(B_global)
  B_wy     <- as.integer(B_wy)
  if (is.na(B_wy) || B_wy < 0L) B_wy <- 0L
  if (!(perm_mode %in% c("wy", "both"))) B_wy <- 0L
  do_perm  <- (B_global > 0L || B_wy > 0L)

  pid <- prep$pathway_id
  n_genes <- lengths(prep$g_list_norm)

  res <- data.frame(
    pathway_id   = pid,
    pathway_name = prep$pathway_name,
    n_genes      = as.integer(n_genes),
    stringsAsFactors = FALSE
  )

  # Track which calibration mode was used per pathway
  res$calib_mode <- NA_character_

  # ============================================================
  # OBSERVED component p-values: CALL YOUR CATFISH WRAPPERS (analytic only)
  # ============================================================
  acat_fun    <- .catfish_get_fun("catfish_acat_pathways")
  fisher_fun  <- .catfish_get_fun("catfish_fisher_pathways")
  minp_fun    <- .catfish_get_fun("catfish_minp_pathways")
  tfisher_fun <- .catfish_get_fun("catfish_soft_tfisher_adaptive_pathways")
  stouf_fun   <- .catfish_get_fun("catfish_stoufferZ_pathways")

  if (is.null(acat_fun) || is.null(fisher_fun) || is.null(minp_fun) || is.null(tfisher_fun)) {
    stop("Could not find required CATFISH wrappers. Load CATFISH or source your wrapper files first.", call. = FALSE)
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

  # TFisher (optional - can be disabled with use_tfisher = FALSE)
  res$tfisher_p_analytic <- NA_real_
  res$tfisher_p_min      <- NA_real_
  res$tau_hat            <- NA_real_
  res$tfisher_stat_hat   <- NA_real_

  if (isTRUE(use_tfisher)) {
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
    res$tfisher_p_min      <- tf_tab$tfisher_p_min[match(res$pathway_id, tf_tab$pathway_id)]
    res$tau_hat            <- tf_tab$tau_hat[match(res$pathway_id, tf_tab$pathway_id)]
    res$tfisher_stat_hat   <- tf_tab$tfisher_stat_hat[match(res$pathway_id, tf_tab$pathway_id)]
  }

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
  if (is.null(stouf_fun)) stop("z_col provided but catfish_stoufferZ_pathways() not found.", call. = FALSE)

  st_tab <- stouf_fun(
    gene_results   = gene_df,
    pathways       = pw_list,
    gene_col       = prep$gene_col,
    z_col          = prep$z_col,
    weight_col     = prep$weight_col,
    min_abs_w      = stouffer_min_abs_w,
    min_p          = prep$min_p,
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
               res$minp_p_analytic[i],
               res$stouffer_p_analytic[i])
    # Include TFisher only if enabled
    if (isTRUE(use_tfisher) && is.finite(res$tfisher_p_analytic[i])) {
      comps <- c(comps, res$tfisher_p_analytic[i])
    }
    if (isTRUE(include_magma_eff) && is.finite(res$magma_pvalue[i])) comps <- c(comps, res$magma_pvalue[i])
    .catfish_omni_methods(comps, omnibus = omnibus, min_p = prep$min_p, do_fix = isTRUE(do_fix))
  }, numeric(1))
  res$omni_p_analytic_BH <- .catfish_bh(res$omni_p_analytic)

  # ============================================================
  # OMNIBUS calibration: global resampling (ONLY omnibus)
  # ============================================================
  res$omni_p_global <- NA_real_
  res$omni_p_global_BH <- NA_real_

  if (B_global > 0L && perm_mode %in% c("global", "both")) {

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

    pool_p <- .catfish_fix_p(pool_p, min_p = prep$min_p, do_fix = isTRUE(do_fix))
    if (length(pool_p) < 2L) stop("Global resampling pool has <2 valid p-values.", call. = FALSE)

    n_pool <- length(pool_p)

    # helper for stouffer p from Z matrix
    .stouffer_p_from_Zmat <- function(Zmat, w = NULL) {
      .stouffer_p_from_Zmat_helper(
        Zmat = Zmat,
        w = w,
        alternative = stouffer_alternative,
        min_abs_w = stouffer_min_abs_w
      )
    }

    n_pathways <- nrow(res)
    for (i in seq_len(n_pathways)) {
      # Progress message
      if (i == 1L || i %% 10L == 0L || i == n_pathways) {
        message(sprintf("[Global] Pathway %d/%d: %s", i, n_pathways, res$pathway_id[i]))
      }

      d <- res$n_genes[i]
      if (!is.finite(d) || d < prep$min_genes) next

      omni_obs <- res$omni_p_analytic[i]
      if (!is.finite(omni_obs) || is.na(omni_obs)) next

      idx_mat <- .catfish_sample_idx_mat(n_pool = n_pool, d = d, B = B_global)
      if (is.null(idx_mat)) next

      P <- matrix(pool_p[idx_mat], nrow = B_global, ncol = d)

      # ACAT across genes (fast)
      pA <- vapply(seq_len(B_global), function(b) .catfish_acat_combine(P[b, ], min_p = prep$min_p, do_fix = isTRUE(do_fix)), numeric(1))
      # Fisher across genes (fast)
      pF <- stats::pchisq(-2 * rowSums(log(P)), df = 2 * d, lower.tail = FALSE)

      # TFisher adaptive (slower; but only in omnibus null)
      pT <- rep(NA_real_, B_global)
      if (!requireNamespace("TFisher", quietly = TRUE)) stop("TFisher required for omnibus TFisher null.", call. = FALSE)
      for (b in seq_len(B_global)) {
        pT[b] <- .catfish_tfisher_adaptive_p(
          P[b, ],
          tau_grid = prep$tau_grid,
          min_p = prep$min_p,
          do_fix = isTRUE(do_fix)
        )
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
        Wmat <- NULL
        if (!is.null(pool_w)) Wmat <- matrix(pool_w[idx_mat], nrow = B_global, ncol = d)
        pS <- .stouffer_p_from_Zmat(Zmat, w = Wmat)
      }

      omni_null <- vapply(seq_len(B_global), function(b) {
        comps <- c(pA[b], pF[b], pT[b], pM[b], pS[b])
        if (isTRUE(include_magma_eff) && is.finite(res$magma_pvalue[i])) comps <- c(comps, res$magma_pvalue[i])
        .catfish_omni_methods(comps, omnibus = omnibus, min_p = prep$min_p, do_fix = isTRUE(do_fix))
      }, numeric(1))

      res$omni_p_global[i] <- .emp_or_gpd_p_lower(
        p_null = omni_null,
        p_obs  = omni_obs,
        method = tail_mode,
        switch_exceed = tail_switch_exceed,
        gpd_k = tail_gpd_k,
        min_B = tail_min_B,
        min_tail = tail_min_tail
      )
      # Track the active calibration in modes that actually run the global block.
      if (perm_mode %in% c("global", "both")) {
        res$calib_mode[i] <- "global"
      }
    }

    res$omni_p_global_BH <- .catfish_bh(res$omni_p_global)
  }


  # ============================================================
  # Westfall-Young calibration (LD-aware)
  # One genome-wide null bank: each of B_wy replicates draws every gene's z
  # jointly from N(0, R) (R from the MAGMA gene-gene correlations, block-
  # diagonal by chromosome) and re-scores every pathway on that same draw.
  # Column-wise the bank gives each pathway a calibrated omnibus p; row-wise
  # its minimum gives the single-step min-P Westfall-Young family-wise p.
  # The omnibus here is ACAT-O over ACAT, Fisher, Stouffer, minP and
  # (optionally) soft TFisher, all computed in log space from gene p-values.
  # ============================================================
  res$omni_p_wy_stat     <- NA_real_   # observed log-space omnibus (as a p)
  res$omni_log10_p_wy_stat <- NA_real_ # same on log10 scale (no underflow)
  res$omni_z_std         <- NA_real_   # standardised against own null
  res$omni_p_cal         <- NA_real_   # per-pathway calibrated p
  res$omni_p_cal_BH      <- NA_real_
  res$omni_p_fwer_wy     <- NA_real_   # WY family-wise p, raw scale
  res$omni_p_fwer_wy_std <- NA_real_   # WY family-wise p, standardised scale
  wy_info <- NULL

  if (B_wy > 0L) {
    if (is.null(prep$gene_chr_all))
      stop("perm_mode = 'wy' needs chromosome labels: pass chr_col to catfish_omni2_prepare().",
           call. = FALSE)
    if (is.character(wy_tau_grid))
      stop("wy_tau_grid must be a fixed numeric grid; a data-adaptive grid breaks ",
           "the strong family-wise control of the Westfall-Young procedure.", call. = FALSE)
    wy_tau_grid <- sort(unique(as.numeric(wy_tau_grid)))
    if (!length(wy_tau_grid) || any(!is.finite(wy_tau_grid)) ||
        any(wy_tau_grid <= 0 | wy_tau_grid >= 1))
      stop("wy_tau_grid values must lie in (0, 1).", call. = FALSE)
    if (omnibus != "ACAT")
      warning("The Westfall-Young omnibus always combines components with ACAT; ",
              "omnibus = '", omnibus, "' applies only to omni_p_analytic/omni_p_global.",
              call. = FALSE)
    if (any(is.finite(res$magma_pvalue)))
      warning("MAGMA competitive p-values cannot be re-evaluated on null draws and ",
              "are not part of the Westfall-Young omnibus.", call. = FALSE)

    cp_wy <- .wy_read_cor_pairs(magma_cor_pairs, magma_cor_file)
    if (is.null(cp_wy) || !nrow(cp_wy))
      warning("No gene correlations supplied: every gene is treated as independent.",
              call. = FALSE)

    gid_wy <- prep$gene_id_all
    glp_wy <- log(pmin(pmax(suppressWarnings(as.numeric(prep$gene_p_all_use)), 1e-320), 1))
    idx_of <- stats::setNames(seq_along(gid_wy), gid_wy)
    paths_wy <- lapply(prep$g_list_norm[res$pathway_id], function(g) as.integer(idx_of[g]))
    names(paths_wy) <- res$pathway_id

    eng <- .wy_engine(gid_wy, glp_wy, prep$gene_chr_all, paths_wy, cp_wy,
                      B_wy, isTRUE(use_tfisher), wy_tau_grid, prep$seed, TRUE)

    res$omni_p_wy_stat     <- exp(eng$obs_log)
    res$omni_log10_p_wy_stat <- eng$obs_log / log(10)
    res$omni_z_std         <- eng$z_std
    res$omni_p_cal         <- eng$p_cal
    res$omni_p_cal_BH      <- .catfish_bh(res$omni_p_cal)
    res$omni_p_fwer_wy     <- eng$p_fwer_wy
    res$omni_p_fwer_wy_std <- eng$p_fwer_wy_std
    res$calib_mode[is.finite(res$omni_p_cal)] <- "wy"

    wy_info <- list(
      B = B_wy,
      tau_grid = wy_tau_grid,
      threshold_05 = exp(stats::quantile(eng$minp_null, 0.05, names = FALSE)),
      repaired_chr = eng$chol$repaired
    )
  }

  # ============================================================
  # FINALIZE + RETURN (YOU NEED THIS)
  # ============================================================
  res$omni_p_final <- res$omni_p_analytic
  if (B_global > 0L) {
    idx_g <- is.finite(res$omni_p_global) & !is.na(res$omni_p_global)
    if (any(idx_g)) res$omni_p_final[idx_g] <- res$omni_p_global[idx_g]
  }
  if (B_wy > 0L) {
    idx_w <- is.finite(res$omni_p_cal) & !is.na(res$omni_p_cal)
    if (any(idx_w)) res$omni_p_final[idx_w] <- res$omni_p_cal[idx_w]
  }
  res$omni_p_final_BH <- .catfish_bh(res$omni_p_final)

  # Fill in calib_mode for pathways that didn't get resampled
  res$calib_mode[is.na(res$calib_mode)] <- "analytic"

  # more lenient options
  res$omni_p_final_q_storey <- .catfish_qvalue(res$omni_p_final)
  res$omni_p_final_q_bky10  <- .catfish_bky(res$omni_p_final, q = 0.10)

  # IHW using pathway size as a covariate (log scale is typical)
  res$omni_p_final_q_ihw10  <- .catfish_ihw(res$omni_p_final, covar = log1p(res$n_genes))

  

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
                  if (B_wy > 0L) "_wy" else if (B_global > 0L) "_global" else "_analytic")
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
  # Use omni_p_final (WY-calibrated if available, else global, else analytic)
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

  if (!is.null(wy_info)) {
    attr(res, "wy_threshold_05") <- wy_info$threshold_05
    attr(res, "wy_B")            <- wy_info$B
    attr(res, "wy_tau_grid")     <- wy_info$tau_grid
  }
  return(res)

}


#' CATFISH omnibus pathway p-values (analytic, global resampling, or Westfall--Young)
#'
#' Compute pathway-level p-values from gene-level results using up to six component
#' methods (ACAT, Fisher, adaptive soft TFisher, minP, Stouffer Z, optional MAGMA
#' competitive) and combine them into a single omnibus p-value per pathway.
#'
#' Calibration is applied **only to the omnibus** (not to each component), using either:
#' \itemize{
#'   \item \strong{Global resampling} (fast heuristic): resample genes from a global pool.
#'   \item \strong{Westfall--Young} (LD-aware, family-wise): one genome-wide null bank. Each of
#'         \code{B_wy} replicates draws every gene's \eqn{z} jointly from \eqn{N(0, R)}, with \eqn{R}
#'         built from MAGMA gene--gene correlations (block-diagonal by chromosome), and re-scores
#'         every pathway on that same draw. Read column-wise, the bank gives each pathway a
#'         calibrated omnibus p-value; read row-wise, its minimum gives the single-step min-P
#'         Westfall--Young family-wise p-value. See \code{\link{catfish_wy_pathways}} for the
#'         stand-alone version.
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
#' otherwise they use \code{p_raw_col}. The Westfall--Young omnibus uses the same column.
#'
#' \strong{Westfall--Young omnibus (\code{perm_mode = "wy"} or \code{"both"})}
#' The statistic calibrated by the null bank is ACAT-O across ACAT, Fisher, Stouffer
#' (from gene p-values), minP and, if \code{use_tfisher = TRUE}, soft TFisher
#' \eqn{W_o} on the fixed grid \code{wy_tau_grid}, all evaluated in log space so
#' that gene p-values far below \eqn{10^{-300}} keep their order. TFisher is
#' omitted from a pathway's omnibus when no gene falls below any threshold and is
#' capped at \eqn{10^{-16}}. MAGMA competitive p-values are not part of this
#' statistic, because they cannot be re-evaluated on null draws. Family-wise
#' control holds in the strong sense under the multivariate normal model, because
#' each pathway's statistic depends only on its own genes; this is why
#' \code{wy_tau_grid} must be fixed rather than estimated from the data.
#'
#' @param gene_results data.frame of gene-level results (one row per gene recommended).
#'   Must contain \code{gene_col} and at least one p-value column (\code{p_raw_col} and/or \code{p_adj_col}).
#'   If you want Stouffer, it must contain \code{z_col}. If you want weighted Stouffer, it must contain
#'   \code{weight_col}. For Westfall--Young it must also contain \code{chr_col}.
#'
#' @param pathways Pathway definitions (provide \code{pathways} OR \code{species}).
#'   Either a data.frame with columns \code{pathway_id}, \code{gene_id} (optional \code{pathway_name}),
#'   or a named list mapping \code{pathway_id -> character vector of genes}.
#'
#' @param species Optional character (provide \code{species} OR \code{pathways}), one of
#'   \code{"maize"}, \code{"sorghum"}, \code{"arabidopsis"}, \code{"plant"}, \code{"fly"},
#'   \code{"human"} (HumanCyc) or \code{"human_reactome"} (Reactome, Entrez IDs).
#'   When provided, pathways are loaded via \code{catfish_load_pathways()}.
#'
#' @param pmn_gene_col Optional character. When \code{species} is used, passed to
#'   \code{catfish_load_pathways(species=..., gene_col=pmn_gene_col)}.
#'
#' @param gene_col Character. Gene ID column in \code{gene_results}. Default \code{"GENE"}.
#' @param p_raw_col Character. Raw p-value column in \code{gene_results}. Default \code{"P"}.
#' @param p_adj_col Character or NULL. Preferred/adjusted p-value column. If supplied and present, used
#'   for ACAT/Fisher/TFisher/minP. Default NULL.
#' @param z_col Character or NULL. Gene Z column for Stouffer. If NULL, Stouffer is skipped. Default NULL.
#' @param weight_col Character or NULL. Optional weights for Stouffer. Default NULL (equal weights).
#' @param chr_col Character or NULL. Chromosome column in \code{gene_results}. Required for
#'   \code{perm_mode = "wy"} / \code{"both"}; genes are drawn jointly within a chromosome and
#'   independently across chromosomes. Default NULL.
#'
#' @param tau_grid Adaptive soft-TFisher thresholds for the analytic component. Either a numeric vector of
#'   \eqn{\tau} values, or \code{"auto"} (default) to select the grid from the
#'   empirical quantiles of the pooled gene-level p-values.
#' @param tau_cap Numeric. Upper bound for \code{tau_grid} (sanity bound). Default 1.
#' @param min_p Numeric. P-values are clipped to \code{[min_p, 1-min_p]} to avoid numeric issues. Default 1e-50.
#' @param min_genes Integer. Minimum genes required per pathway. Default 2.
#' @param max_genes Numeric. Maximum genes per pathway (counted over genes present). Default \code{Inf}
#'   (no upper limit).
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
#'   The Westfall--Young omnibus always uses ACAT.
#'
#' @param B_perm Integer or NULL. Convenience parameter for number of resampling iterations. If provided,
#'   it overrides \code{B_global}/\code{B_wy} according to \code{perm_mode}.
#' @param perm_mode Character. One of \code{"none"} (analytic only), \code{"global"} (global gene
#'   resampling), \code{"wy"} (Westfall--Young) or \code{"both"} (global and Westfall--Young;
#'   \code{omni_p_final} prefers the Westfall--Young calibrated p-value).
#' @param wy_tau_grid Numeric. Fixed soft-TFisher grid used inside the Westfall--Young omnibus.
#'   Default \code{c(0.01, 0.05, 0.5)}: the grid of Zhang et al. (2020) without its \eqn{\tau = 1}
#'   corner, which is Fisher.
#'
#' @param B_global Integer. Number of global-resampling iterations (used when \code{B_perm} is NULL).
#' @param B_wy Integer. Number of genome-wide Westfall--Young null replicates (used when \code{B_perm} is NULL).
#' @param perm_pool Character. For global resampling only: \code{"raw"} (use raw p-values) or \code{"obs"} (use p used by methods).
#'
#' @param magma_cor_file Character or NULL. Path to 3-col correlation pairs file: \code{gene1 gene2 r}.
#'   Used by Westfall--Young; without it (or \code{magma_cor_pairs}) genes are treated as independent.
#' @param magma_cor_pairs data.frame or NULL. In-memory correlation pairs with columns \code{gene1}, \code{gene2}, \code{r}.
#'   If provided, it overrides \code{magma_cor_file}.
#'
#' @param seed Integer or NULL. Random seed.
#' @param output Logical. If TRUE, writes CSV to \code{out_dir} and attaches \code{"file"} attribute. Default FALSE.
#' @param out_dir Character. Output directory. Default \code{"catfish_omni2"}.
#' @param use_tfisher Logical. Include the adaptive soft-TFisher component. Default TRUE.
#' @param tail_mode Character. \code{"empirical"} or \code{"hybrid_gpd"} (GPD tail
#'   extrapolation for global-resampling p-values below the Monte Carlo floor). Default \code{"empirical"}.
#' @param tail_switch_exceed Integer. Switch to GPD extrapolation when the number of null
#'   exceedances falls below this. Default 10.
#' @param tail_gpd_k Integer. Number of extreme null draws used to fit the GPD. Default 250.
#' @param tail_min_B Integer. Minimum number of null draws required before GPD is attempted.
#'   Default 10000.
#' @param tail_min_tail Integer. Minimum tail observations required for a stable GPD fit.
#'   Default 50.
#'
#' @return data.frame with one row per pathway, including:
#' \itemize{
#'   \item component p-values (analytic): \code{acat_p}, \code{fisher_p}, \code{tfisher_p_analytic}, \code{minp_p_analytic}, \code{stouffer_p_analytic}
#'   \item optional \code{magma_pvalue}
#'   \item \code{omni_p_analytic}, \code{omni_p_global}
#'   \item Westfall--Young columns (when run): \code{omni_p_wy_stat} (observed log-space omnibus)
#'     and \code{omni_log10_p_wy_stat} (its log10, which does not underflow),
#'     \code{omni_z_std} (standardised against the pathway's own null), \code{omni_p_cal}
#'     (per-pathway calibrated p) and \code{omni_p_cal_BH}, \code{omni_p_fwer_wy_std}
#'     (Westfall--Young family-wise p on the standardised scale; the headline family-wise result)
#'     and \code{omni_p_fwer_wy} (raw scale)
#'   \item \code{omni_p_final} (priority: Westfall--Young calibrated > global > analytic) + BH/qvalue extras
#' }
#' When Westfall--Young is run, the attributes \code{wy_threshold_05}, \code{wy_B} and
#' \code{wy_tau_grid} are attached.
#'
#' @examples
#' ## Toy GWAS: 300 genes on 3 chromosomes, neighbouring genes correlated
#' ## (r = 0.3), 6 associated genes inside pathway "PWY_signal".
#' set.seed(1)
#' genes <- data.frame(GENE = paste0("g", 1:300), CHR = rep(1:3, each = 100))
#' R <- diag(100)
#' R[cbind(1:99, 2:100)] <- R[cbind(2:100, 1:99)] <- 0.3
#' z <- as.vector(replicate(3, t(chol(R)) %*% rnorm(100)))
#' z[1:6] <- z[1:6] + 4
#' genes$ZSTAT <- z
#' genes$P     <- pnorm(z, lower.tail = FALSE)
#' cor_pairs <- do.call(rbind, lapply(split(genes$GENE, genes$CHR), function(g)
#'   data.frame(gene1 = head(g, -1), gene2 = tail(g, -1), r = 0.3)))
#' pathways <- c(list(PWY_signal = paste0("g", 1:15)),
#'               setNames(lapply(1:19, function(i) sample(genes$GENE[-(1:15)], 15)),
#'                        paste0("PWY_", 1:19)))
#'
#' # analytic omnibus only
#' res0 <- catfish_omni2_pathways(genes, pathways = pathways, z_col = "ZSTAT",
#'                                min_genes = 10)
#'
#' # Westfall--Young calibration (small B for the example)
#' res_wy <- catfish_omni2_pathways(genes, pathways = pathways, z_col = "ZSTAT",
#'                                  chr_col = "CHR", min_genes = 10,
#'                                  perm_mode = "wy", B_perm = 200,
#'                                  magma_cor_pairs = cor_pairs, seed = 1)
#' head(res_wy[, c("pathway_id", "omni_p_analytic", "omni_p_cal",
#'                 "omni_p_fwer_wy_std")])
#'
#' \dontrun{
#' gene_results <- read.table("magma_output.genes.out", header = TRUE)
#'
#' # Analytic p-values only
#' omni_res <- catfish_omni2_pathways(
#'   gene_results = gene_results,
#'   species = "maize",
#'   gene_col = "GENE",
#'   p_raw_col = "P",
#'   z_col = "ZSTAT",
#'   perm_mode = "none"
#' )
#'
#' # Westfall--Young calibration (LD-aware, family-wise)
#' omni_wy <- catfish_omni2_pathways(
#'   gene_results = gene_results,
#'   species = "maize",
#'   chr_col = "CHR",
#'   perm_mode = "wy",
#'   B_perm = 10000,
#'   magma_cor_file = "magma_gene_cor_pairs.txt",
#'   seed = 1
#' )
#' attr(omni_wy, "wy_threshold_05")
#' }
#'
#' @seealso \code{\link{catfish_wy_pathways}}, \code{\link{catfish_acat_pathways}},
#'   \code{\link{catfish_fisher_pathways}}, \code{\link{catfish_minp_pathways}},
#'   \code{\link{catfish_stoufferZ_pathways}}
#' @export
catfish_omni2_pathways <- function(gene_results,
                                 pathways = NULL,
                                 species = NULL,
                                 pmn_gene_col = NULL,
                                 gene_col = "GENE",
                                 p_raw_col = "P",
                                 p_adj_col = NULL,
                                 z_col = NULL,
                                 weight_col = NULL,
                                 chr_col = NULL,
                                 tau_grid = "auto",
                                 tau_cap = 1,
                                 min_p = 1e-50,
                                 min_genes = 2L,
                                 max_genes = Inf,
                                 # component knobs
                                 use_tfisher = TRUE,
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
                                 perm_mode = c("none","global","wy","both"),
                                 wy_tau_grid = c(0.01, 0.05, 0.5),
                                 # GPD tail extrapolation (Knijnenburg et al. 2009)
                                 tail_mode = c("empirical", "hybrid_gpd"),
                                 tail_switch_exceed = 10L,
                                 tail_gpd_k = 250L,
                                 tail_min_B = 10000L,
                                 tail_min_tail = 50L,
                                 # back-compat / low-level
                                 B_global = 0L,
                                 B_wy = 0L,
                                 perm_pool = c("raw","obs"),
                                 magma_cor_file = NULL,
                                 magma_cor_pairs = NULL,
                                 seed = NULL,
                                 output = FALSE,
                                 out_dir = "catfish_omni2") {

  # ---- basic arg matching / validation ----
  omnibus   <- match.arg(omnibus)
  perm_mode <- match.arg(perm_mode)
  perm_pool <- match.arg(perm_pool)
  stouffer_alternative <- match.arg(stouffer_alternative)
  tail_mode <- match.arg(tail_mode)

  if (!is.null(species)) {
    species <- match.arg(species, c("maize","sorghum","arabidopsis","plant","fly",
                                    "human","human_reactome"))
  }

  min_p <- as.numeric(min_p)
  if (!is.finite(min_p) || is.na(min_p) || min_p <= 0 || min_p >= 0.5) {
    stop("min_p must be in (0, 0.5).", call. = FALSE)
  }

  min_genes <- as.integer(min_genes)
  if (!is.finite(min_genes) || is.na(min_genes) || min_genes < 1L) {
    stop("min_genes must be a positive integer.", call. = FALSE)
  }

  # ---- map B_perm + perm_mode -> B_global/B_wy (B_perm overrides) ----
  if (!is.null(B_perm)) {
    B_perm <- as.integer(B_perm)
    if (!is.finite(B_perm) || is.na(B_perm) || B_perm < 0L) stop("B_perm must be >= 0.", call. = FALSE)

    if (perm_mode == "none" || B_perm == 0L) {
      B_global <- 0L; B_wy <- 0L
    } else if (perm_mode == "global") {
      B_global <- B_perm; B_wy <- 0L
    } else if (perm_mode == "wy") {
      B_global <- 0L; B_wy <- B_perm
    } else { # both
      B_global <- B_perm; B_wy <- B_perm
    }
  } else {
    # respect perm_mode even if user sets B_global/B_wy
    B_global <- as.integer(B_global); if (is.na(B_global) || B_global < 0L) B_global <- 0L
    B_wy     <- as.integer(B_wy);     if (is.na(B_wy)     || B_wy     < 0L) B_wy <- 0L

    if (perm_mode == "none") {
      B_global <- 0L; B_wy <- 0L
    } else if (perm_mode == "global") {
      B_wy <- 0L
    } else if (perm_mode == "wy") {
      B_global <- 0L
    } else {
      # both: keep as-is
    }
  }

  # Westfall-Young requirements
  if (B_wy > 0L) {
    if (is.null(chr_col)) {
      stop("perm_mode includes 'wy' (B_wy > 0) but chr_col is NULL.\n",
           "Provide the chromosome column of gene_results (e.g. chr_col = \"CHR\").",
           call. = FALSE)
    }
    if (isTRUE(use_tfisher) && !requireNamespace("TFisher", quietly = TRUE)) {
      stop("Package 'TFisher' is required when use_tfisher = TRUE.", call. = FALSE)
    }
  }

  # Stouffer requirements (skip silently if z_col is NULL)
  if (is.null(z_col) && !is.null(weight_col)) {
    warning("weight_col provided but z_col is NULL; Stouffer is skipped and weight_col is ignored.", call. = FALSE)
    weight_col <- NULL
  }

  # ---- prepare (standardizes gene_results + pathways) ----
  prep <- catfish_omni2_prepare(
    gene_results = gene_results,
    pathways = pathways,
    species = species,
    pmn_gene_col = pmn_gene_col,
    gene_col = gene_col,
    p_raw_col = p_raw_col,
    p_adj_col = p_adj_col,
    z_col = z_col,
    weight_col = weight_col,
    chr_col = chr_col,
    tau_grid = tau_grid,
    tau_cap = tau_cap,
    min_p = min_p,
    min_genes = min_genes,
    max_genes = max_genes,
    seed = seed
  )

  # ---- run (components analytic; calibration ONLY for omnibus) ----
  catfish_omni2_run(
    prep = prep,
    omnibus = omnibus,
    use_tfisher = use_tfisher,
    do_fix = do_fix,
    stouffer_min_abs_w = stouffer_min_abs_w,
    stouffer_alternative = stouffer_alternative,
    magma_out = magma_out,
    include_magma_in_omni = include_magma_in_omni,
    include_magma_in_perm = include_magma_in_perm,
    B_global = B_global,
    B_wy = B_wy,
    perm_mode = perm_mode,
    perm_pool = perm_pool,
    wy_tau_grid = wy_tau_grid,
    magma_cor_file = magma_cor_file,
    magma_cor_pairs = magma_cor_pairs,
    tail_mode = tail_mode,
    tail_switch_exceed = tail_switch_exceed,
    tail_gpd_k = tail_gpd_k,
    tail_min_B = tail_min_B,
    tail_min_tail = tail_min_tail,
    output = output,
    out_dir = out_dir
  )
}
