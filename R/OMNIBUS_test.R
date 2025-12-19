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
#' @param species Character. If \code{pathways} is \code{NULL}, a species keyword (e.g. \code{"maize"})
#'   used to auto-load built-in pathway sets (e.g., PMN). Provide \code{species} \strong{or}
#'   \code{pathways} (not both).
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
# ----------------------------
.magcat_get_fun <- function(fname) {
  if (exists(fname, mode = "function")) return(get(fname, mode = "function"))
  if (requireNamespace("MAGCAT", quietly = TRUE) &&
      exists(fname, envir = asNamespace("MAGCAT"), mode = "function")) {
    return(getFromNamespace(fname, "MAGCAT"))
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
    fp <- NULL
    if (exists("fix_p_for_acat", mode = "function")) {
      fp <- get("fix_p_for_acat", mode = "function")
    } else if (requireNamespace("MAGCAT", quietly = TRUE) &&
               exists("fix_p_for_acat", envir = asNamespace("MAGCAT"), mode = "function")) {
      fp <- getFromNamespace("fix_p_for_acat", "MAGCAT")
    }
    if (!is.null(fp)) {
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

.magma_build_R_from_pairs <- function(genes, cor_pairs, make_PD = TRUE) {
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

  # load pathways if species (use your wrapper helper)
  if (is.null(pathways)) {
    if (exists("magcat_load_pathways", mode = "function")) {
      pathways <- magcat_load_pathways(species = species, gene_col = pmn_gene_col)
    } else if (requireNamespace("MAGCAT", quietly = TRUE) &&
               exists("magcat_load_pathways", envir = asNamespace("MAGCAT"), mode = "function")) {
      pathways <- MAGCAT::magcat_load_pathways(species = species, gene_col = pmn_gene_col)
    } else {
      stop("Could not find magcat_load_pathways(). Make sure MAGCAT is loaded.", call. = FALSE)
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

  # ACAT
  acat_tab <- acat_fun(
    gene_results = gene_df,
    pathways     = pw_list,
    gene_col     = prep$gene_col,
    p_col        = prep$p_use_col,
    min_p        = prep$min_p,
    do_fix       = isTRUE(do_fix),
    B            = 0L,
    seed         = prep$seed,
    output       = FALSE
  )
  res$acat_p <- acat_tab$acat_p[match(res$pathway_id, acat_tab$pathway_id)]

  # Fisher
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

  # Adaptive soft TFisher (analytic only; NO per-method permutations)
  tf_tab <- tfisher_fun(
    gene_results     = gene_df,
    pathways         = pw_list,
    gene_col         = prep$gene_col,
    p_col            = prep$p_use_col,
    tau_grid         = prep$tau_grid,
    min_p            = prep$min_p,
    do_fix           = isTRUE(do_fix),
    B_perm           = 0L,
    perm_mode        = "resample_global",  # irrelevant when B_perm=0
    magma_genes_out  = NULL,
    magma_cor_file   = NULL,
    make_PD          = make_PD,
    seed             = prep$seed,
    analytic_logical = TRUE,
    output           = FALSE
  )
  res$tfisher_p_analytic <- tf_tab$tfisher_p_analytic[match(res$pathway_id, tf_tab$pathway_id)]
  res$tau_hat            <- tf_tab$tau_hat[match(res$pathway_id, tf_tab$pathway_id)]
  res$tfisher_stat_hat   <- tf_tab$tfisher_stat_hat[match(res$pathway_id, tf_tab$pathway_id)]

  # minP (analytic only; NO per-method permutations)
  minp_tab <- minp_fun(
    gene_results     = gene_df,
    pathways         = pw_list,
    gene_col         = prep$gene_col,
    p_col            = prep$p_use_col,
    min_p            = prep$min_p,
    do_fix           = isTRUE(do_fix),
    B_perm           = 0L,
    seed             = prep$seed,
    analytic_logical = TRUE,
    output           = FALSE
  )
  res$minp_p_analytic <- minp_tab$minp_p_analytic[match(res$pathway_id, minp_tab$pathway_id)]

  # StoufferZ (analytic only) if z_col provided
  res$stouffer_p_analytic <- NA_real_
  if (!is.null(prep$z_col)) {
    if (is.null(stouf_fun)) stop("z_col provided but magcat_stoufferZ_pathways() not found.", call. = FALSE)

    st_tab <- stouf_fun(
      gene_results    = gene_df,
      pathways        = pw_list,
      gene_col        = prep$gene_col,
      z_col           = prep$z_col,
      weight_col      = prep$weight_col,
      min_abs_w       = stouffer_min_abs_w,
      B_perm          = 0L,
      perm_mode       = "resample_global", # irrelevant when B_perm=0
      magma_genes_out = NULL,
      magma_cor_file  = NULL,
      make_PD         = make_PD,
      seed            = prep$seed,
      alternative     = stouffer_alternative,
      output          = FALSE
    )
    res$stouffer_p_analytic <- st_tab$stouffer_p_analytic[match(res$pathway_id, st_tab$pathway_id)]
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
  # OMNIBUS calibration: MVN resampling (LD-aware; ONLY omnibus)
  # ============================================================
  res$omni_p_mvn <- NA_real_
  res$omni_p_mvn_BH <- NA_real_

  if (B_mvn > 0L) {

    cor_pairs <- .magma_read_gene_cor_pairs(magma_cor_file = magma_cor_file, magma_cor_pairs = magma_cor_pairs)
    if (is.null(cor_pairs) || !nrow(cor_pairs)) stop("MVN resampling requires magma_cor_file or magma_cor_pairs.", call. = FALSE)

    if (!is.null(prep$seed)) set.seed(prep$seed)

    for (i in seq_len(nrow(res))) {

      d <- res$n_genes[i]
      if (!is.finite(d) || d < prep$min_genes) next

      omni_obs <- res$omni_p_analytic[i]
      if (!is.finite(omni_obs) || is.na(omni_obs)) next

      genes_S <- prep$g_list_raw[[res$pathway_id[i]]]
      genes_S <- as.character(genes_S)
      genes_S <- genes_S[!is.na(genes_S) & genes_S != ""]
      if (length(genes_S) < prep$min_genes) next

      R <- .magma_build_R_from_pairs(genes = genes_S, cor_pairs = cor_pairs, make_PD = make_PD)
      if (is.null(R)) next

      Z <- .magma_simulate_Z_from_R(R, B = B_mvn)
      if (is.null(Z)) next

      if (mvn_marginal == "uniform") {
        P <- .mvn_z_to_p(Z, min_p = prep$min_p)
      } else {
        if (!is.null(pool_p)) {
          pool_this <- pool_p
        } else {
          # empirical pool from observed p-values (excluding pathway genes is ideal; minimal version here)
          pool_this <- prep$gene_p_all_use
        }
        P <- .mvn_z_to_empirical_p(Z, pool_p = pool_this, min_p = prep$min_p)
      }

      # ACAT/Fisher/minP from P
      pA <- vapply(seq_len(B_mvn), function(b) .magcat_acat_combine(P[b, ], min_p = prep$min_p, do_fix = isTRUE(do_fix)), numeric(1))
      pF <- stats::pchisq(-2 * rowSums(log(P)), df = 2 * ncol(P), lower.tail = FALSE)

      # TFisher adaptive from P
      pT <- rep(NA_real_, B_mvn)
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
        pT[b] <- best
      }

      pM <- {
        k <- ncol(P)
        pmin <- apply(P, 1, min)
        1 - (1 - pmin)^k
      }

      # Stouffer from MVN Z (use observed pathway weights if provided)
      pS <- rep(NA_real_, B_mvn)
      if (!is.null(prep$z_col) && !is.null(prep$z_list)) {
        w_i <- NULL
        if (!is.null(prep$w_list)) w_i <- prep$w_list[[res$pathway_id[i]]]

        # clean weights like wrapper
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
          pS <- stats::pnorm(Zs, lower.tail = FALSE)
        } else if (stouffer_alternative == "less") {
          pS <- stats::pnorm(Zs, lower.tail = TRUE)
        } else {
          pS <- 2 * stats::pnorm(-abs(Zs))
        }
      }

      omni_null <- vapply(seq_len(B_mvn), function(b) {
        comps <- c(pA[b], pF[b], pT[b], pM[b], pS[b])
        if (isTRUE(include_magma_eff) && is.finite(res$magma_pvalue[i])) comps <- c(comps, res$magma_pvalue[i])
        .magcat_omni_methods(comps, omnibus = omnibus, min_p = prep$min_p, do_fix = isTRUE(do_fix))
      }, numeric(1))

      res$omni_p_mvn[i] <- (1 + sum(omni_null <= omni_obs, na.rm = TRUE)) / (B_mvn + 1)
    }

    res$omni_p_mvn_BH <- .magcat_bh(res$omni_p_mvn)
  }

  # final p priority: MVN > global > analytic
  res$omni_p_final <- res$omni_p_analytic
  if (B_global > 0L) res$omni_p_final <- res$omni_p_global
  if (B_mvn > 0L)    res$omni_p_final <- res$omni_p_mvn
  res$omni_p_final_BH <- .magcat_bh(res$omni_p_final)

  res <- res[order(res$omni_p_final, decreasing = FALSE, na.last = TRUE), , drop = FALSE]

  if (isTRUE(output)) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    tag <- paste0("omni_", tolower(omnibus),
                  if (B_mvn > 0L) "_mvn" else if (B_global > 0L) "_global" else "_analytic")
    out_path <- file.path(out_dir, paste0(tag, ".csv"))
    utils::write.csv(res, out_path, row.names = FALSE)
    attr(res, "file") <- out_path
  }

  res
}

# ----------------------------
# Main user function
# ----------------------------
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
                                 # use YOUR wrapper knobs
                                 do_fix = TRUE,
                                 stouffer_min_abs_w = 1e-8,
                                 stouffer_alternative = c("greater","two.sided","less"),
                                 # magma competitive
                                 magma_out = NULL,
                                 include_magma_in_omni = TRUE,
                                 include_magma_in_perm = FALSE,
                                 omnibus = c("ACAT","minP"),
                                 # NEW:
                                 B_perm = NULL,
                                 perm_mode = c("none","global","mvn","both"),
                                 # NEW (MVN marginal control)
                                 mvn_marginal = c("uniform","empirical"),
                                 mvn_pool     = c("use","raw"),
                                 pool_p       = NULL,
                                 # BACKCOMPAT:
                                 B_global = 0L,
                                 B_mvn = 0L,
                                 perm_pool = c("raw","obs"),
                                 magma_cor_file = NULL,
                                 magma_cor_pairs = NULL,
                                 make_PD = TRUE,
                                 seed = NULL,
                                 output = FALSE,
                                 out_dir = "magcat_omni2") {

  omnibus   <- match.arg(omnibus)
  perm_pool <- match.arg(perm_pool)
  perm_mode <- match.arg(perm_mode)
  mvn_marginal <- match.arg(mvn_marginal)
  mvn_pool     <- match.arg(mvn_pool)
  stouffer_alternative <- match.arg(stouffer_alternative)

  # map B_perm + perm_mode -> B_global / B_mvn
  if (!is.null(B_perm)) {
    B_perm <- as.integer(B_perm)
    if (perm_mode == "global") {
      B_global <- B_perm; B_mvn <- 0L
    } else if (perm_mode == "mvn") {
      B_global <- 0L; B_mvn <- B_perm
    } else if (perm_mode == "both") {
      B_global <- B_perm; B_mvn <- B_perm
    } else {
      B_global <- 0L; B_mvn <- 0L
    }
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
    mvn_pool     = mvn_pool,
    pool_p       = pool_p,
    magma_cor_file = magma_cor_file,
    magma_cor_pairs = magma_cor_pairs,
    make_PD = make_PD,
    output = output,
    out_dir = out_dir
  )
}





















# magcat_omni2_pathways <- function(gene_results,
#                                   pathways = NULL,
#                                   species = NULL,
#                                   pmn_gene_col = NULL,
#                                   gene_col = "GENE",
#                                   p_raw_col = "P",
#                                   p_adj_col = NULL,
#                                   z_col = NULL,
#                                   weight_col = NULL,
#                                   tau_grid = c(0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.001),
#                                   tau_cap = 1,
#                                   min_p = 1e-15,
#                                   min_genes = 2L,
#                                   magma_out = NULL,
#                                   include_magma_in_omni = TRUE,
#                                   include_magma_in_perm = FALSE,
#                                   omnibus = c("ACAT","minP"),
#                                   ## NEW:
#                                   B_perm = NULL,
#                                   perm_mode = c("none","global","mvn","both"),
#                                   ## NEW (MVN marginal control)
#                                   mvn_marginal = c("uniform","empirical"),
#                                   mvn_pool     = c("use","raw"),
#                                   pool_p       = NULL,
#                                   ## BACKCOMPAT:
#                                   B_global = 0L,
#                                   B_mvn = 0L,
#                                   perm_pool = c("raw","obs"),
#                                   magma_cor_file = NULL,
#                                   magma_cor_pairs = NULL,
#                                   make_PD = TRUE,
#                                   seed = NULL,
#                                   output = FALSE,
#                                   out_dir = "magcat_omni2") {

#   omnibus   <- match.arg(omnibus)
#   perm_pool <- match.arg(perm_pool)
#   perm_mode <- match.arg(perm_mode)

#   # If user uses the new interface, map it
#   if (!is.null(B_perm)) {
#     B_perm <- as.integer(B_perm)
#     if (perm_mode == "global") {
#       B_global <- B_perm; B_mvn <- 0L
#     } else if (perm_mode == "mvn") {
#       B_global <- 0L; B_mvn <- B_perm
#     } else if (perm_mode == "both") {
#       B_global <- B_perm; B_mvn <- B_perm
#     } else {
#       B_global <- 0L; B_mvn <- 0L
#     }
#   }

#   prep <- magcat_omni2_prepare(
#     gene_results = gene_results,
#     pathways = pathways,
#     species = species,
#     pmn_gene_col = pmn_gene_col,
#     gene_col = gene_col,
#     p_raw_col = p_raw_col,
#     p_adj_col = p_adj_col,
#     z_col = z_col,
#     weight_col = weight_col,
#     tau_grid = tau_grid,
#     tau_cap = tau_cap,
#     min_p = min_p,
#     min_genes = min_genes,
#     seed = seed
#   )

#   magcat_omni2_run(
#     prep = prep,
#     omnibus = omnibus,
#     magma_out = magma_out,
#     include_magma_in_omni = include_magma_in_omni,
#     include_magma_in_perm = include_magma_in_perm,
#     B_global = B_global,
#     B_mvn = B_mvn,
#     mvn_marginal = mvn_marginal,
#     mvn_pool     = mvn_pool,
#     pool_p       = pool_p,
#     perm_pool = perm_pool,
#     magma_cor_file = magma_cor_file,
#     magma_cor_pairs = magma_cor_pairs,
#     make_PD = make_PD,
#     output = output,
#     out_dir = out_dir
#   )
# }


# ## ----------------------------
# ## helpers: p fixing + combining
# ## ----------------------------

# ## ----------------------------
# ## MVN: Z -> p helpers (uniform vs empirical marginals)
# ## ----------------------------

# .mvn_z_to_p <- function(Z, min_p = 1e-15) {
#   U <- stats::pnorm(Z)
#   P <- 2 * pmin(U, 1 - U)
#   P <- pmax(pmin(P, 1 - min_p), min_p)
#   P
# }

# .empirical_qfun <- function(pool_p) {
#   pool_p <- suppressWarnings(as.numeric(pool_p))
#   pool_p <- pool_p[is.finite(pool_p) & !is.na(pool_p) & pool_p > 0 & pool_p < 1]
#   pool_p <- sort(pool_p)
#   n <- length(pool_p)
#   if (n < 1000) warning("pool_p is small (n=", n, "); empirical mapping may be noisy.")

#   stats::approxfun(
#     x = (1:n) / (n + 1),
#     y = pool_p,
#     method = "linear",
#     rule = 2,
#     ties = "ordered"
#   )
# }

# .mvn_z_to_empirical_p <- function(Z, pool_p, min_p = 1e-15) {
#   q0 <- .empirical_qfun(pool_p)
#   U  <- stats::pnorm(Z)
#   P  <- matrix(q0(as.vector(U)), nrow = nrow(Z), ncol = ncol(Z))
#   P  <- pmax(pmin(P, 1 - min_p), min_p)
#   P
# }

# .magcat_fix_p <- function(p, min_p = 1e-15) {
#   p <- suppressWarnings(as.numeric(p))
#   p <- p[is.finite(p) & !is.na(p)]
#   if (!length(p)) return(numeric(0))
#   p[p <= 0] <- min_p
#   p[p >= 1] <- 1 - min_p
#   p
# }

# .magcat_bh <- function(p) {
#   p <- as.numeric(p)
#   out <- rep(NA_real_, length(p))
#   idx <- which(is.finite(p) & !is.na(p))
#   if (length(idx)) out[idx] <- stats::p.adjust(p[idx], method = "BH")
#   out
# }

# ## ACAT combine for a vector of p-values (equal weights)
# .magcat_acat_combine <- function(p, min_p = 1e-15) {
#   p <- .magcat_fix_p(p, min_p = min_p)
#   p <- p[p > 0 & p < 1]
#   if (!length(p)) return(NA_real_)
#   tstat <- mean(tan((0.5 - p) * pi))
#   out <- 0.5 - atan(tstat) / pi
#   out <- max(min(out, 1 - min_p), min_p)
#   as.numeric(out)
# }

# ## ACAT combine row-wise for a matrix P (B x d)
# .magcat_acat_combine_mat <- function(P, min_p = 1e-15) {
#   P <- suppressWarnings(matrix(as.numeric(P), nrow = nrow(P)))
#   P[!is.finite(P) | is.na(P)] <- NA_real_
#   P[P <= 0] <- min_p
#   P[P >= 1] <- 1 - min_p

#   # rowMeans(tan(.)) but safely skip NAs
#   T <- rowMeans(tan((0.5 - P) * pi), na.rm = TRUE)
#   out <- 0.5 - atan(T) / pi
#   out[out <= 0] <- min_p
#   out[out >= 1] <- 1 - min_p
#   as.numeric(out)
# }

# ## Sidak minP across genes (vector)
# .magcat_minp_sidak <- function(p, min_p = 1e-15) {
#   p <- .magcat_fix_p(p, min_p = min_p)
#   p <- p[p > 0 & p < 1]
#   k <- length(p)
#   if (k < 1L) return(NA_real_)
#   pmin <- min(p)
#   1 - (1 - pmin)^k
# }

# ## Fisher across genes (vector)
# .magcat_fisher <- function(p, min_p = 1e-15) {
#   p <- .magcat_fix_p(p, min_p = min_p)
#   p <- p[p > 0 & p < 1]
#   k <- length(p)
#   if (k < 2L) return(NA_real_)
#   stats::pchisq(-2 * sum(log(p)), df = 2 * k, lower.tail = FALSE)
# }

# ## Stouffer from Z (vector) with optional weights
# .magcat_stouffer_p <- function(z, w = NULL, min_abs_w = 1e-8) {
#   z <- suppressWarnings(as.numeric(z))
#   keep <- is.finite(z) & !is.na(z)
#   z <- z[keep]
#   if (length(z) < 2L) return(NA_real_)

#   if (is.null(w)) {
#     w <- rep(1, length(z))
#   } else {
#     w <- suppressWarnings(as.numeric(w))[keep]
#     bad <- !is.finite(w) | is.na(w) | abs(w) <= min_abs_w
#     if (any(bad)) {
#       w_ok <- w[!bad]
#       repl <- if (length(w_ok)) stats::median(abs(w_ok)) else 1
#       w[bad] <- repl
#     }
#   }

#   Zs <- sum(w * z) / sqrt(sum(w^2))
#   2 * stats::pnorm(-abs(Zs))
# }

# ## Adaptive soft TFisher analytic p-value across tau_grid
# .magcat_tfisher_adapt <- function(p, tau_grid, min_p = 1e-15) {
#   if (!requireNamespace("TFisher", quietly = TRUE)) {
#     stop("TFisher package required for adaptive soft TFisher.", call. = FALSE)
#   }
#   p <- .magcat_fix_p(p, min_p = min_p)
#   p <- p[p > 0 & p < 1]
#   n <- length(p)
#   if (n < 2L) return(list(p = NA_real_, tau = NA_real_, stat = NA_real_))

#   best_p <- Inf
#   best_tau <- tau_grid[1]
#   best_stat <- NA_real_

#   for (tau in tau_grid) {
#     st <- TFisher::stat.soft(p = p, tau1 = tau)
#     # p.soft returns CDF -> right tail = 1 - CDF
#     pv <- 1 - as.numeric(TFisher::p.soft(q = st, n = n, tau1 = tau, M = NULL))
#     if (is.finite(pv) && !is.na(pv) && pv < best_p) {
#       best_p <- pv
#       best_tau <- tau
#       best_stat <- st
#     }
#   }

#   list(p = as.numeric(best_p), tau = as.numeric(best_tau), stat = as.numeric(best_stat))
# }

# ## omnibus across METHODS (vector of method p-values)
# .magcat_omni_methods <- function(p_methods, omnibus = c("ACAT", "minP"), min_p = 1e-15) {
#   omnibus <- match.arg(omnibus)
#   pv <- .magcat_fix_p(p_methods, min_p = min_p)
#   pv <- pv[pv > 0 & pv < 1]
#   if (!length(pv)) return(NA_real_)

#   if (omnibus == "ACAT") {
#     .magcat_acat_combine(pv, min_p = min_p)
#   } else {
#     k <- length(pv)
#     pmin <- min(pv)
#     1 - (1 - pmin)^k
#   }
# }

# ## ----------------------------
# ## MAGMA correlation helpers (MVN)
# ## ----------------------------
# .magma_read_gene_cor_pairs <- function(magma_cor_file = NULL, magma_cor_pairs = NULL) {
#   if (!is.null(magma_cor_pairs)) {
#     tab <- magma_cor_pairs
#     if (!is.data.frame(tab) || ncol(tab) < 3) stop("magma_cor_pairs must be a data.frame with >=3 cols.", call. = FALSE)
#     tab <- tab[, 1:3]
#     names(tab) <- c("gene1", "gene2", "r")
#     tab$gene1 <- as.character(tab$gene1)
#     tab$gene2 <- as.character(tab$gene2)
#     tab$r <- suppressWarnings(as.numeric(tab$r))
#     tab <- tab[is.finite(tab$r) & !is.na(tab$r), , drop = FALSE]
#     return(tab)
#   }

#   if (is.null(magma_cor_file)) return(NULL)
#   if (!file.exists(magma_cor_file)) stop("magma_cor_file does not exist: ", magma_cor_file, call. = FALSE)

#   if (requireNamespace("data.table", quietly = TRUE)) {
#     tab <- data.table::fread(magma_cor_file, header = FALSE, data.table = FALSE, fill = TRUE, showProgress = FALSE)
#   } else {
#     tab <- utils::read.table(magma_cor_file, header = FALSE, stringsAsFactors = FALSE, fill = TRUE)
#   }
#   if (ncol(tab) < 3) stop("magma_cor_file must have 3 columns: gene1 gene2 r", call. = FALSE)

#   tab <- tab[, 1:3]
#   names(tab) <- c("gene1", "gene2", "r")
#   tab$gene1 <- as.character(tab$gene1)
#   tab$gene2 <- as.character(tab$gene2)
#   tab$r <- suppressWarnings(as.numeric(tab$r))
#   tab <- tab[is.finite(tab$r) & !is.na(tab$r), , drop = FALSE]
#   tab
# }

# .magma_build_R_from_pairs <- function(genes, cor_pairs, make_PD = TRUE) {
#   genes <- as.character(genes)
#   genes <- genes[!is.na(genes) & genes != ""]
#   d <- length(genes)
#   if (d < 2L) return(NULL)

#   R <- diag(1, d)
#   rownames(R) <- genes
#   colnames(R) <- genes

#   if (!is.null(cor_pairs) && nrow(cor_pairs)) {
#     idx <- seq_len(d)
#     pos <- stats::setNames(idx, genes)

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
#     if (!requireNamespace("Matrix", quietly = TRUE)) {
#       stop("Matrix package required for make_PD=TRUE.", call. = FALSE)
#     }
#     npd <- Matrix::nearPD(R, corr = TRUE)
#     R <- as.matrix(npd$mat)
#   }

#   R
# }

# .magma_simulate_Z_from_R <- function(R, B) {
#   B <- as.integer(B)
#   d <- ncol(R)
#   if (B < 1L || d < 2L) return(NULL)

#   # Cholesky for correlated normals
#   U <- tryCatch(chol(R), error = function(e) NULL)
#   if (is.null(U)) stop("chol(R) failed even after PD-fix. Check correlation input.", call. = FALSE)

#   Z0 <- matrix(stats::rnorm(B * d), nrow = B, ncol = d)
#   Z  <- Z0 %*% U
#   colnames(Z) <- colnames(R)
#   Z
# }

# ## ----------------------------
# ## prepare (fast reuse)
# ## ----------------------------
# magcat_omni2_prepare <- function(gene_results,
#                                  pathways = NULL,
#                                  species = NULL,
#                                  pmn_gene_col = NULL,
#                                  gene_col = "GENE",
#                                  p_raw_col = "P",
#                                  p_adj_col = NULL,
#                                  z_col = NULL,
#                                  weight_col = NULL,
#                                  tau_grid = c(0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.001),
#                                  tau_cap = 1,
#                                  min_p = 1e-15,
#                                  min_genes = 2L,
#                                  seed = NULL) {

#   if (!is.null(seed)) set.seed(seed)

#   if (!is.data.frame(gene_results)) stop("gene_results must be a data.frame.", call. = FALSE)
#   if (!(gene_col %in% names(gene_results))) stop("gene_results missing gene_col: ", gene_col, call. = FALSE)

#   if (!is.null(pathways) && !is.null(species)) {
#     stop("Provide either pathways OR species, not both.", call. = FALSE)
#   }
#   if (is.null(pathways) && is.null(species)) {
#     stop("Need pathways (list/df) or species.", call. = FALSE)
#   }

#   # load pathways if species
#   if (is.null(pathways)) {
#     fn <- NULL
#     if (requireNamespace("MAGCAT", quietly = TRUE)) {
#       fn <- tryCatch(getFromNamespace("magcat_load_pathways", "MAGCAT"), error = function(e) NULL)
#     }
#     if (is.null(fn)) stop("Could not find MAGCAT::magcat_load_pathways().", call. = FALSE)
#     pathways <- fn(species = species, gene_col = pmn_gene_col)
#   }

#   # choose p column for non-MAGMA methods: prefer adjusted if provided and present
#   p_use_col <- p_raw_col
#   if (!is.null(p_adj_col) && (p_adj_col %in% names(gene_results))) {
#     p_use_col <- p_adj_col
#   }
#   if (!(p_use_col %in% names(gene_results))) {
#     stop("Chosen p column not found in gene_results: ", p_use_col, call. = FALSE)
#   }
#   if (!(p_raw_col %in% names(gene_results))) {
#     # raw p might be missing; that's OK unless you want raw pool / MAGMA tidy from gene_results
#     p_raw_col <- NULL
#   }

#   # gene ids
#   g_raw  <- as.character(gene_results[[gene_col]])
#   g_norm <- tolower(trimws(g_raw))
#   okg <- !is.na(g_norm) & g_norm != "" & g_norm != "unknown"
#   g_raw <- g_raw[okg]
#   g_norm <- g_norm[okg]

#     # p (use)
#   p_use <- suppressWarnings(as.numeric(gene_results[[p_use_col]]))[okg]
#   okp <- is.finite(p_use) & !is.na(p_use)
#   g_raw  <- g_raw[okp]
#   g_norm <- g_norm[okp]
#   p_use  <- p_use[okp]

#   # ALSO capture raw p aligned to the same kept genes (optional)
#   p_raw_vec <- NULL
#   if (!is.null(p_raw_col) && (p_raw_col %in% names(gene_results))) {
#     p_raw_vec <- suppressWarnings(as.numeric(gene_results[[p_raw_col]]))[okg]
#     p_raw_vec <- p_raw_vec[okp]
#   }

#   # z
#   z_use <- NULL
#   if (!is.null(z_col)) {
#     if (!(z_col %in% names(gene_results))) stop("z_col not found: ", z_col, call. = FALSE)
#     z_use <- suppressWarnings(as.numeric(gene_results[[z_col]]))[okg][okp]
#   }

#   # weights
#   w_use <- NULL
#   if (!is.null(weight_col)) {
#     if (!(weight_col %in% names(gene_results))) stop("weight_col not found: ", weight_col, call. = FALSE)
#     w_use <- suppressWarnings(as.numeric(gene_results[[weight_col]]))[okg][okp]
#   }

#   # dedup by normalized gene id (keep first)
#   keep1 <- !duplicated(g_norm)
#   g_raw  <- g_raw[keep1]
#   g_norm <- g_norm[keep1]
#   p_use  <- p_use[keep1]
#   if (!is.null(p_raw_vec)) p_raw_vec <- p_raw_vec[keep1]
#   if (!is.null(z_use)) z_use <- z_use[keep1]
#   if (!is.null(w_use)) w_use <- w_use[keep1]

#   gene_map <- stats::setNames(g_raw, g_norm)
  
#   gene_p_use_by_norm <- stats::setNames(p_use, g_norm)
#   gene_p_raw_by_norm <- if (!is.null(p_raw_vec)) stats::setNames(p_raw_vec, g_norm) else NULL

#   gene_p   <- stats::setNames(p_use, g_norm)
#   gene_z   <- if (!is.null(z_use)) stats::setNames(z_use, g_norm) else NULL
#   gene_w   <- if (!is.null(w_use)) stats::setNames(w_use, g_norm) else NULL

#   # coerce pathways to long
#   if (is.data.frame(pathways)) {
#     if (!all(c("pathway_id","gene_id") %in% names(pathways))) {
#       stop("pathways df must have pathway_id and gene_id.", call. = FALSE)
#     }
#     if (!("pathway_name" %in% names(pathways))) pathways$pathway_name <- pathways$pathway_id
#     pw_long <- pathways[, c("pathway_id","pathway_name","gene_id"), drop = FALSE]
#     pw_long$gene_id <- as.character(pw_long$gene_id)
#   } else if (is.list(pathways)) {
#     pid <- names(pathways); if (is.null(pid)) pid <- paste0("PWY_", seq_along(pathways))
#     pw_long <- data.frame(
#       pathway_id   = rep(pid, lengths(pathways)),
#       pathway_name = rep(pid, lengths(pathways)),
#       gene_id      = unlist(pathways, use.names = FALSE),
#       stringsAsFactors = FALSE
#     )
#   } else {
#     stop("pathways must be a data.frame or list.", call. = FALSE)
#   }

#   pw_long$gene_id   <- gsub("^gene:", "", pw_long$gene_id, ignore.case = TRUE)
#   pw_long$gene_norm <- tolower(trimws(pw_long$gene_id))
#   pw_long <- pw_long[!is.na(pw_long$gene_norm) & pw_long$gene_norm != "" & pw_long$gene_norm != "unknown", , drop = FALSE]

#   g_list_norm <- split(pw_long$gene_norm, pw_long$pathway_id)
#   pw_name     <- tapply(pw_long$pathway_name, pw_long$pathway_id, function(x) x[1])

#   # intersect with available genes
#   g_list_norm <- lapply(g_list_norm, function(g) unique(g[g %in% names(gene_p)]))

#   # drop small pathways early (optional)
#   min_genes <- as.integer(min_genes)
#   if (!is.na(min_genes) && min_genes >= 1L) {
#     keep_pw <- lengths(g_list_norm) >= min_genes
#     g_list_norm <- g_list_norm[keep_pw]
#   }

#   pid <- names(g_list_norm)

#   # per pathway vectors
#   p_list <- lapply(g_list_norm, function(g) unname(gene_p[g]))
#   z_list <- if (!is.null(gene_z)) lapply(g_list_norm, function(g) unname(gene_z[g])) else NULL
#   w_list <- if (!is.null(gene_w)) lapply(g_list_norm, function(g) unname(gene_w[g])) else NULL
#   g_list_raw <- lapply(g_list_norm, function(g) unname(gene_map[g]))

#   # tau sanitize (allow tau=1)
#   tau_grid <- suppressWarnings(as.numeric(tau_grid))
#   tau_grid <- tau_grid[is.finite(tau_grid) & !is.na(tau_grid) & tau_grid > 0 & tau_grid <= tau_cap]
#   tau_grid <- sort(unique(tau_grid), decreasing = TRUE)
#   if (!length(tau_grid)) tau_grid <- 0.05

#   structure(
#     list(
#       pathway_id   = pid,
#       pathway_name = as.character(pw_name[pid]),
#       g_list_norm  = g_list_norm,
#       g_list_raw   = g_list_raw,
#       p_list       = p_list,
#       z_list       = z_list,
#       w_list       = w_list,
#       gene_p_all_use = p_use,              # for obs pool if needed
#       gene_p_all_raw = if (!is.null(p_raw_col)) suppressWarnings(as.numeric(gene_results[[p_raw_col]])) else NULL,
#       gene_z_all     = if (!is.null(z_col)) suppressWarnings(as.numeric(gene_results[[z_col]])) else NULL,
#       gene_w_all     = if (!is.null(weight_col)) suppressWarnings(as.numeric(gene_results[[weight_col]])) else NULL,
#       gene_col = gene_col,
#       p_use_col = p_use_col,
#       p_raw_col = p_raw_col,
#       p_adj_col = p_adj_col,
#       z_col = z_col,
#       weight_col = weight_col,
#       tau_grid = tau_grid,
#       tau_cap = tau_cap,
#       min_p = min_p,
#       min_genes = min_genes,
#       seed = seed,
#       gene_p_by_norm_use = gene_p_use_by_norm,
#       gene_p_by_norm_raw = gene_p_raw_by_norm
#     ),
#     class = "magcat_omni2_prep"
#   )
# }

# ## ============================================================
# ## PATCH: fix global sampler + add perm_mode/B_perm interface
# ## Paste this AFTER your previous definitions to override them.
# ## ============================================================

# .magcat_sample_idx_mat <- function(n_pool, d, B) {
#   n_pool <- as.integer(n_pool)
#   d      <- as.integer(d)
#   B      <- as.integer(B)
#   if (n_pool < 1L || d < 1L || B < 1L) return(NULL)

#   if (n_pool >= d) {
#     # IMPORTANT: sample WITHOUT replacement *within each permutation*
#     idx <- replicate(B, sample.int(n_pool, size = d, replace = FALSE))
#     t(idx) # B x d
#   } else {
#     # fallback when pathway is larger than pool: sample WITH replacement
#     matrix(sample.int(n_pool, size = B * d, replace = TRUE), nrow = B, ncol = d)
#   }
# }

# magcat_omni2_run <- function(prep,
#                              omnibus = c("ACAT","minP"),
#                              magma_out = NULL,
#                              include_magma_in_omni = TRUE,
#                              include_magma_in_perm = FALSE,
#                              B_global = 0L,
#                              B_mvn = 0L,
#                              perm_pool = c("raw","obs"),
#                              mvn_marginal = c("uniform","empirical"),
#                              mvn_pool     = c("use","raw"),
#                              pool_p       = NULL,
#                              magma_cor_file = NULL,
#                              magma_cor_pairs = NULL,
#                              make_PD = TRUE,
#                              output = FALSE,
#                              out_dir = "magcat_omni2") {


#   if (!inherits(prep, "magcat_omni2_prep")) stop("prep must come from magcat_omni2_prepare().", call. = FALSE)
#   omnibus   <- match.arg(omnibus)
#   perm_pool <- match.arg(perm_pool)
#   mvn_marginal <- match.arg(mvn_marginal)
#   mvn_pool     <- match.arg(mvn_pool)

#   pid <- prep$pathway_id
#   n_genes <- lengths(prep$g_list_norm)

#   res <- data.frame(
#     pathway_id   = pid,
#     pathway_name = prep$pathway_name,
#     n_genes      = as.integer(n_genes),
#     stringsAsFactors = FALSE
#   )

#   ## observed component p-values (analytic)
#   res$acat_p <- vapply(pid, function(k) .magcat_acat_combine(prep$p_list[[k]], min_p = prep$min_p), numeric(1))
#   res$fisher_p <- vapply(pid, function(k) .magcat_fisher(prep$p_list[[k]], min_p = prep$min_p), numeric(1))

#   tf_obs <- lapply(pid, function(k) .magcat_tfisher_adapt(prep$p_list[[k]], prep$tau_grid, min_p = prep$min_p))
#   res$tfisher_p_analytic <- vapply(tf_obs, `[[`, numeric(1), "p")
#   res$tau_hat            <- vapply(tf_obs, `[[`, numeric(1), "tau")
#   res$tfisher_stat_hat   <- vapply(tf_obs, `[[`, numeric(1), "stat")

#   res$minp_p_analytic <- vapply(pid, function(k) .magcat_minp_sidak(prep$p_list[[k]], min_p = prep$min_p), numeric(1))

#   if (!is.null(prep$z_col) && !is.null(prep$z_list)) {
#     res$stouffer_p_analytic <- vapply(pid, function(k) .magcat_stouffer_p(prep$z_list[[k]], prep$w_list[[k]]), numeric(1))
#   } else {
#     res$stouffer_p_analytic <- NA_real_
#   }

#   ## MAGMA competitive (you pass magma_out data.frame "out")
#   res$magma_pvalue <- NA_real_
#   if (!is.null(magma_out)) {
#     if (!is.data.frame(magma_out)) stop("magma_out must be a data.frame (MAGMA competitive output).", call. = FALSE)
#     if (!all(c("pathway_id","magma_pvalue") %in% names(magma_out))) {
#       stop("magma_out must have columns pathway_id and magma_pvalue.", call. = FALSE)
#     }
#     res$magma_pvalue <- magma_out$magma_pvalue[match(res$pathway_id, magma_out$pathway_id)]
#   }

#   ## omnibus analytic (across METHODS)
#   res$omni_p_analytic <- vapply(seq_len(nrow(res)), function(i) {
#     comps <- c(res$acat_p[i],
#                res$fisher_p[i],
#                res$tfisher_p_analytic[i],
#                res$minp_p_analytic[i],
#                res$stouffer_p_analytic[i])

#     if (isTRUE(include_magma_in_omni) && is.finite(res$magma_pvalue[i])) {
#       comps <- c(comps, res$magma_pvalue[i])
#     }
#     .magcat_omni_methods(comps, omnibus = omnibus, min_p = prep$min_p)
#   }, numeric(1))

#   res$omni_p_analytic_BH <- .magcat_bh(res$omni_p_analytic)

#   ## ----------------------------
#   ## global resampling (omnibus only)  [FIXED]
#   ## ----------------------------
#   res$omni_p_global <- NA_real_
#   res$omni_p_global_BH <- NA_real_

#   B_global <- as.integer(B_global)
#   if (B_global > 0L) {

#     if (!is.null(prep$seed)) set.seed(prep$seed)

#     pool_p <- NULL
#     pool_z <- NULL

#     if (perm_pool == "raw" && !is.null(prep$p_raw_col) && !is.null(prep$gene_p_all_raw)) {
#       pool_p <- suppressWarnings(as.numeric(prep$gene_p_all_raw))
#     } else {
#       pool_p <- suppressWarnings(as.numeric(prep$gene_p_all_use))
#     }

#     okp <- is.finite(pool_p) & !is.na(pool_p) & pool_p > 0 & pool_p < 1
#     pool_p <- pool_p[okp]

#     if (!is.null(prep$z_col) && !is.null(prep$gene_z_all)) {
#       pool_z_full <- suppressWarnings(as.numeric(prep$gene_z_all))
#       pool_z_full <- pool_z_full[okp]
#       okz <- is.finite(pool_z_full) & !is.na(pool_z_full)
#       pool_p <- pool_p[okz]
#       pool_z <- pool_z_full[okz]
#     }

#     pool_p <- .magcat_fix_p(pool_p, min_p = prep$min_p)
#     if (length(pool_p) < 2L) stop("Global resampling pool has <2 valid p-values.", call. = FALSE)

#     n_pool <- length(pool_p)

#     for (i in seq_len(nrow(res))) {

#       d <- res$n_genes[i]
#       if (!is.finite(d) || d < prep$min_genes) next

#       omni_obs <- res$omni_p_analytic[i]
#       if (!is.finite(omni_obs) || is.na(omni_obs)) next

#       idx_mat <- .magcat_sample_idx_mat(n_pool = n_pool, d = d, B = B_global)
#       if (is.null(idx_mat)) next

#       P <- matrix(pool_p[idx_mat], nrow = B_global, ncol = d)

#       pA <- .magcat_acat_combine_mat(P, min_p = prep$min_p)
#       pF <- stats::pchisq(-2 * rowSums(log(P)), df = 2 * d, lower.tail = FALSE)

#       pT <- rep(NA_real_, B_global)
#       for (b in seq_len(B_global)) {
#         pT[b] <- .magcat_tfisher_adapt(P[b, ], prep$tau_grid, min_p = prep$min_p)$p
#       }

#       pM <- {
#         pmin <- apply(P, 1, min)
#         1 - (1 - pmin)^d
#       }

#       pS <- rep(NA_real_, B_global)
#       if (!is.null(pool_z)) {
#         Z <- matrix(pool_z[idx_mat], nrow = B_global, ncol = d)
#         Zs <- rowSums(Z) / sqrt(d)
#         pS <- 2 * stats::pnorm(-abs(Zs))
#       }

#       omni_null <- vapply(seq_len(B_global), function(b) {
#         comps <- c(pA[b], pF[b], pT[b], pM[b], pS[b])

#         # Not recommended: including fixed MAGMA p in null
#         if (isTRUE(include_magma_in_perm) && isTRUE(include_magma_in_omni) && is.finite(res$magma_pvalue[i])) {
#           comps <- c(comps, res$magma_pvalue[i])
#         }

#         .magcat_omni_methods(comps, omnibus = omnibus, min_p = prep$min_p)
#       }, numeric(1))

#       res$omni_p_global[i] <- (1 + sum(omni_null <= omni_obs, na.rm = TRUE)) / (B_global + 1)
#     }

#     res$omni_p_global_BH <- .magcat_bh(res$omni_p_global)
#   }

#   ## ----------------------------
#   ## MVN resampling (LD-aware)
#   ## ----------------------------
#   res$omni_p_mvn <- NA_real_
#   res$omni_p_mvn_BH <- NA_real_

#   B_mvn <- as.integer(B_mvn)
#   if (B_mvn > 0L) {

#     cor_pairs <- .magma_read_gene_cor_pairs(magma_cor_file = magma_cor_file, magma_cor_pairs = magma_cor_pairs)
#     if (is.null(cor_pairs) || !nrow(cor_pairs)) {
#       stop("MVN resampling requires magma_cor_file or magma_cor_pairs (gene1, gene2, r).", call. = FALSE)
#     }

#     if (!is.null(prep$seed)) set.seed(prep$seed)

#     for (i in seq_len(nrow(res))) {

#       d <- res$n_genes[i]
#       if (!is.finite(d) || d < prep$min_genes) next

#       omni_obs <- .magcat_omni_methods(
#         c(res$acat_p[i], res$fisher_p[i], res$tfisher_p_analytic[i], res$minp_p_analytic[i], res$stouffer_p_analytic[i]),
#         omnibus = omnibus, min_p = prep$min_p
#       )
#       if (!is.finite(omni_obs) || is.na(omni_obs)) next

#       genes_S <- prep$g_list_raw[[res$pathway_id[i]]]
#       genes_S <- as.character(genes_S)
#       genes_S <- genes_S[!is.na(genes_S) & genes_S != ""]
#       if (length(genes_S) < prep$min_genes) next

#       R <- .magma_build_R_from_pairs(genes = genes_S, cor_pairs = cor_pairs, make_PD = make_PD)
#       if (is.null(R)) next

#       Z <- .magma_simulate_Z_from_R(R, B = B_mvn)
#       if (is.null(Z)) next

#       if (mvn_marginal == "uniform") {
#         P <- .mvn_z_to_p(Z, min_p = prep$min_p)
#       } else {
#         # Build a null pool (exclude pathway genes to reduce leakage)
#         if (!is.null(pool_p)) {
#           pool_this <- pool_p
#         } else {
#           g_norm_S <- prep$g_list_norm[[res$pathway_id[i]]]
#           if (mvn_pool == "raw" && !is.null(prep$gene_p_by_norm_raw)) {
#             pool_all <- prep$gene_p_by_norm_raw
#           } else {
#             pool_all <- prep$gene_p_by_norm_use
#           }
#           # drop pathway genes
#           pool_this <- pool_all[!(names(pool_all) %in% g_norm_S)]
#         }

#         P <- .mvn_z_to_empirical_p(Z, pool_p = pool_this, min_p = prep$min_p)
#       }

#       pA <- .magcat_acat_combine_mat(P, min_p = prep$min_p)
#       pF <- stats::pchisq(-2 * rowSums(log(P)), df = 2 * ncol(P), lower.tail = FALSE)

#       pT <- rep(NA_real_, B_mvn)
#       for (b in seq_len(B_mvn)) {
#         pT[b] <- .magcat_tfisher_adapt(P[b, ], prep$tau_grid, min_p = prep$min_p)$p
#       }

#       pM <- {
#         k <- ncol(P)
#         pmin <- apply(P, 1, min)
#         1 - (1 - pmin)^k
#       }

#       pS <- {
#         k <- ncol(Z)
#         Zs <- rowSums(Z) / sqrt(k)
#         2 * stats::pnorm(-abs(Zs))
#       }

#       omni_null <- vapply(seq_len(B_mvn), function(b) {
#         comps <- c(pA[b], pF[b], pT[b], pM[b], pS[b])

#         if (isTRUE(include_magma_in_perm) && isTRUE(include_magma_in_omni) && is.finite(res$magma_pvalue[i])) {
#           comps <- c(comps, res$magma_pvalue[i])
#         }

#         .magcat_omni_methods(comps, omnibus = omnibus, min_p = prep$min_p)
#       }, numeric(1))

#       res$omni_p_mvn[i] <- (1 + sum(omni_null <= omni_obs, na.rm = TRUE)) / (B_mvn + 1)
#     }

#     res$omni_p_mvn_BH <- .magcat_bh(res$omni_p_mvn)
#   }

#   ## final p priority: MVN > global > analytic
#   res$omni_p_final <- res$omni_p_analytic
#   if (B_global > 0L) res$omni_p_final <- res$omni_p_global
#   if (B_mvn > 0L)    res$omni_p_final <- res$omni_p_mvn
#   res$omni_p_final_BH <- .magcat_bh(res$omni_p_final)

#   res <- res[order(res$omni_p_final, decreasing = FALSE, na.last = TRUE), , drop = FALSE]

#   if (isTRUE(output)) {
#     if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#     tag <- paste0("omni_", tolower(omnibus),
#                   if (B_mvn > 0L) "_mvn" else if (B_global > 0L) "_global" else "_analytic")
#     out_path <- file.path(out_dir, paste0(tag, ".csv"))
#     utils::write.csv(res, out_path, row.names = FALSE)
#     attr(res, "file") <- out_path
#   }

#   res
# }



## ============================================================
## HOW TO CALL (what you expected)
## ============================================================
## 1) Global resampling only:
# omni <- magcat_omni2_pathways(
#   gene_results = genes_adj, species="maize", gene_col="GENE",
#   p_raw_col="P", p_adj_col="P_adj", z_col="Z_adj",
#   magma_out = out,
#   omnibus="ACAT",
#   B_perm = 1000L, perm_mode="global",
#   perm_pool="raw",
#   seed=123, output=TRUE, out_dir="magcat_omni2_results"
# )

## 2) MVN only (LD-aware) — REQUIRES magma_cor_file (or magma_cor_pairs):
# omni <- magcat_omni2_pathways(
#   gene_results = genes_adj, species="maize", gene_col="GENE",
#   p_raw_col="P", p_adj_col="P_adj", z_col="Z_adj",
#   magma_out = out,
#   omnibus="ACAT",
#   B_perm = 1000L, perm_mode="mvn",
#   magma_cor_file="magma_gene_cor_pairs.txt",
#   make_PD=TRUE,
#   seed=123, output=TRUE, out_dir="magcat_omni2_results"
# )

## 3) BOTH global + MVN:
# omni <- magcat_omni2_pathways(
#   gene_results = genes_adj, species="maize", gene_col="GENE",
#   p_raw_col="P", p_adj_col="P_adj", z_col="Z_adj",
#   magma_out = out,
#   omnibus="ACAT",
#   B_perm = 1000L, perm_mode="both",
#   perm_pool="raw",
#   magma_cor_file="magma_gene_cor_pairs.txt",
#   make_PD=TRUE,
#   seed=123, output=TRUE, out_dir="magcat_omni2_results"
# )






# ## ============================================================
# ##  MAGCAT omnibus (fresh, fast, 6-method, with global + MVN resampling)
# ##
# ##  Methods (component p-values per pathway):
# ##    1) ACAT (gene-level)
# ##    2) Fisher (gene-level)
# ##    3) Adaptive soft TFisher (analytic p.soft across tau_grid)
# ##    4) minP (Sidak across genes)
# ##    5) Stouffer Z (from gene-level Z)
# ##    6) MAGMA competitive (optional; you pass the data.frame output "out")
# ##
# ##  Omnibus across METHODS (per pathway): ACAT or minP (Sidak across methods)
# ##
# ##  Permutation / resampling (omnibus only):
# ##    - global_resampling: resample genes from a global pool (fast heuristic)
# ##    - mvn_resampling: simulate correlated Z ~ MVN(0, R_S) using gene-gene correlations
# ##      (requires magma_cor_file or magma_cor_pairs)
# ##
# ##  Key design choices:
# ##    - Non-MAGMA methods use adjusted p-values if you provide p_adj_col; else raw.
# ##    - MAGMA competitive uses your provided magma_out and your raw/unadjusted genes_all.
# ##    - Permutations NEVER include MAGMA competitive by default (no cheap null for MAGMA).
# ## ============================================================

# ## ----------------------------
# ## helpers: p fixing + combining
# ## ----------------------------
# .magcat_fix_p <- function(p, min_p = 1e-15) {
#   p <- suppressWarnings(as.numeric(p))
#   p <- p[is.finite(p) & !is.na(p)]
#   if (!length(p)) return(numeric(0))
#   p[p <= 0] <- min_p
#   p[p >= 1] <- 1 - min_p
#   p
# }

# .magcat_bh <- function(p) {
#   p <- as.numeric(p)
#   out <- rep(NA_real_, length(p))
#   idx <- which(is.finite(p) & !is.na(p))
#   if (length(idx)) out[idx] <- stats::p.adjust(p[idx], method = "BH")
#   out
# }

# ## ACAT combine for a vector of p-values (equal weights)
# .magcat_acat_combine <- function(p, min_p = 1e-15) {
#   p <- .magcat_fix_p(p, min_p = min_p)
#   p <- p[p > 0 & p < 1]
#   if (!length(p)) return(NA_real_)
#   tstat <- mean(tan((0.5 - p) * pi))
#   out <- 0.5 - atan(tstat) / pi
#   out <- max(min(out, 1 - min_p), min_p)
#   as.numeric(out)
# }

# ## ACAT combine row-wise for a matrix P (B x d)
# .magcat_acat_combine_mat <- function(P, min_p = 1e-15) {
#   P <- suppressWarnings(matrix(as.numeric(P), nrow = nrow(P)))
#   P[!is.finite(P) | is.na(P)] <- NA_real_
#   P[P <= 0] <- min_p
#   P[P >= 1] <- 1 - min_p

#   # rowMeans(tan(.)) but safely skip NAs
#   T <- rowMeans(tan((0.5 - P) * pi), na.rm = TRUE)
#   out <- 0.5 - atan(T) / pi
#   out[out <= 0] <- min_p
#   out[out >= 1] <- 1 - min_p
#   as.numeric(out)
# }

# ## Sidak minP across genes (vector)
# .magcat_minp_sidak <- function(p, min_p = 1e-15) {
#   p <- .magcat_fix_p(p, min_p = min_p)
#   p <- p[p > 0 & p < 1]
#   k <- length(p)
#   if (k < 1L) return(NA_real_)
#   pmin <- min(p)
#   1 - (1 - pmin)^k
# }

# ## Fisher across genes (vector)
# .magcat_fisher <- function(p, min_p = 1e-15) {
#   p <- .magcat_fix_p(p, min_p = min_p)
#   p <- p[p > 0 & p < 1]
#   k <- length(p)
#   if (k < 2L) return(NA_real_)
#   stats::pchisq(-2 * sum(log(p)), df = 2 * k, lower.tail = FALSE)
# }

# ## Stouffer from Z (vector) with optional weights
# .magcat_stouffer_p <- function(z, w = NULL, min_abs_w = 1e-8) {
#   z <- suppressWarnings(as.numeric(z))
#   keep <- is.finite(z) & !is.na(z)
#   z <- z[keep]
#   if (length(z) < 2L) return(NA_real_)

#   if (is.null(w)) {
#     w <- rep(1, length(z))
#   } else {
#     w <- suppressWarnings(as.numeric(w))[keep]
#     bad <- !is.finite(w) | is.na(w) | abs(w) <= min_abs_w
#     if (any(bad)) {
#       w_ok <- w[!bad]
#       repl <- if (length(w_ok)) stats::median(abs(w_ok)) else 1
#       w[bad] <- repl
#     }
#   }

#   Zs <- sum(w * z) / sqrt(sum(w^2))
#   2 * stats::pnorm(-abs(Zs))
# }

# ## Adaptive soft TFisher analytic p-value across tau_grid
# .magcat_tfisher_adapt <- function(p, tau_grid, min_p = 1e-15) {
#   if (!requireNamespace("TFisher", quietly = TRUE)) {
#     stop("TFisher package required for adaptive soft TFisher.", call. = FALSE)
#   }
#   p <- .magcat_fix_p(p, min_p = min_p)
#   p <- p[p > 0 & p < 1]
#   n <- length(p)
#   if (n < 2L) return(list(p = NA_real_, tau = NA_real_, stat = NA_real_))

#   best_p <- Inf
#   best_tau <- tau_grid[1]
#   best_stat <- NA_real_

#   for (tau in tau_grid) {
#     st <- TFisher::stat.soft(p = p, tau1 = tau)
#     # p.soft returns CDF -> right tail = 1 - CDF
#     pv <- 1 - as.numeric(TFisher::p.soft(q = st, n = n, tau1 = tau, M = NULL))
#     if (is.finite(pv) && !is.na(pv) && pv < best_p) {
#       best_p <- pv
#       best_tau <- tau
#       best_stat <- st
#     }
#   }

#   list(p = as.numeric(best_p), tau = as.numeric(best_tau), stat = as.numeric(best_stat))
# }

# ## omnibus across METHODS (vector of method p-values)
# .magcat_omni_methods <- function(p_methods, omnibus = c("ACAT", "minP"), min_p = 1e-15) {
#   omnibus <- match.arg(omnibus)
#   pv <- .magcat_fix_p(p_methods, min_p = min_p)
#   pv <- pv[pv > 0 & pv < 1]
#   if (!length(pv)) return(NA_real_)

#   if (omnibus == "ACAT") {
#     .magcat_acat_combine(pv, min_p = min_p)
#   } else {
#     k <- length(pv)
#     pmin <- min(pv)
#     1 - (1 - pmin)^k
#   }
# }

# ## ----------------------------
# ## MAGMA correlation helpers (MVN)
# ## ----------------------------
# .magma_read_gene_cor_pairs <- function(magma_cor_file = NULL, magma_cor_pairs = NULL) {
#   if (!is.null(magma_cor_pairs)) {
#     tab <- magma_cor_pairs
#     if (!is.data.frame(tab) || ncol(tab) < 3) stop("magma_cor_pairs must be a data.frame with >=3 cols.", call. = FALSE)
#     tab <- tab[, 1:3]
#     names(tab) <- c("gene1", "gene2", "r")
#     tab$gene1 <- as.character(tab$gene1)
#     tab$gene2 <- as.character(tab$gene2)
#     tab$r <- suppressWarnings(as.numeric(tab$r))
#     tab <- tab[is.finite(tab$r) & !is.na(tab$r), , drop = FALSE]
#     return(tab)
#   }

#   if (is.null(magma_cor_file)) return(NULL)
#   if (!file.exists(magma_cor_file)) stop("magma_cor_file does not exist: ", magma_cor_file, call. = FALSE)

#   if (requireNamespace("data.table", quietly = TRUE)) {
#     tab <- data.table::fread(magma_cor_file, header = FALSE, data.table = FALSE, fill = TRUE, showProgress = FALSE)
#   } else {
#     tab <- utils::read.table(magma_cor_file, header = FALSE, stringsAsFactors = FALSE, fill = TRUE)
#   }
#   if (ncol(tab) < 3) stop("magma_cor_file must have 3 columns: gene1 gene2 r", call. = FALSE)

#   tab <- tab[, 1:3]
#   names(tab) <- c("gene1", "gene2", "r")
#   tab$gene1 <- as.character(tab$gene1)
#   tab$gene2 <- as.character(tab$gene2)
#   tab$r <- suppressWarnings(as.numeric(tab$r))
#   tab <- tab[is.finite(tab$r) & !is.na(tab$r), , drop = FALSE]
#   tab
# }

# .magma_build_R_from_pairs <- function(genes, cor_pairs, make_PD = TRUE) {
#   genes <- as.character(genes)
#   genes <- genes[!is.na(genes) & genes != ""]
#   d <- length(genes)
#   if (d < 2L) return(NULL)

#   R <- diag(1, d)
#   rownames(R) <- genes
#   colnames(R) <- genes

#   if (!is.null(cor_pairs) && nrow(cor_pairs)) {
#     idx <- seq_len(d)
#     pos <- stats::setNames(idx, genes)

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
#     if (!requireNamespace("Matrix", quietly = TRUE)) {
#       stop("Matrix package required for make_PD=TRUE.", call. = FALSE)
#     }
#     npd <- Matrix::nearPD(R, corr = TRUE)
#     R <- as.matrix(npd$mat)
#   }

#   R
# }

# .magma_simulate_Z_from_R <- function(R, B) {
#   B <- as.integer(B)
#   d <- ncol(R)
#   if (B < 1L || d < 2L) return(NULL)

#   # Cholesky for correlated normals
#   U <- tryCatch(chol(R), error = function(e) NULL)
#   if (is.null(U)) stop("chol(R) failed even after PD-fix. Check correlation input.", call. = FALSE)

#   Z0 <- matrix(stats::rnorm(B * d), nrow = B, ncol = d)
#   Z  <- Z0 %*% U
#   colnames(Z) <- colnames(R)
#   Z
# }

# ## ----------------------------
# ## prepare (fast reuse)
# ## ----------------------------
# magcat_omni2_prepare <- function(gene_results,
#                                  pathways = NULL,
#                                  species = NULL,
#                                  pmn_gene_col = NULL,
#                                  gene_col = "GENE",
#                                  p_raw_col = "P",
#                                  p_adj_col = NULL,
#                                  z_col = NULL,
#                                  weight_col = NULL,
#                                  tau_grid = c(0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.001),
#                                  tau_cap = 1,
#                                  min_p = 1e-15,
#                                  min_genes = 2L,
#                                  seed = NULL) {

#   if (!is.null(seed)) set.seed(seed)

#   if (!is.data.frame(gene_results)) stop("gene_results must be a data.frame.", call. = FALSE)
#   if (!(gene_col %in% names(gene_results))) stop("gene_results missing gene_col: ", gene_col, call. = FALSE)

#   if (!is.null(pathways) && !is.null(species)) {
#     stop("Provide either pathways OR species, not both.", call. = FALSE)
#   }
#   if (is.null(pathways) && is.null(species)) {
#     stop("Need pathways (list/df) or species.", call. = FALSE)
#   }

#   # load pathways if species
#   if (is.null(pathways)) {
#     fn <- NULL
#     if (requireNamespace("MAGCAT", quietly = TRUE)) {
#       fn <- tryCatch(getFromNamespace("magcat_load_pathways", "MAGCAT"), error = function(e) NULL)
#     }
#     if (is.null(fn)) stop("Could not find MAGCAT::magcat_load_pathways().", call. = FALSE)
#     pathways <- fn(species = species, gene_col = pmn_gene_col)
#   }

#   # choose p column for non-MAGMA methods: prefer adjusted if provided and present
#   p_use_col <- p_raw_col
#   if (!is.null(p_adj_col) && (p_adj_col %in% names(gene_results))) {
#     p_use_col <- p_adj_col
#   }
#   if (!(p_use_col %in% names(gene_results))) {
#     stop("Chosen p column not found in gene_results: ", p_use_col, call. = FALSE)
#   }
#   if (!(p_raw_col %in% names(gene_results))) {
#     # raw p might be missing; that's OK unless you want raw pool / MAGMA tidy from gene_results
#     p_raw_col <- NULL
#   }

#   # gene ids
#   g_raw  <- as.character(gene_results[[gene_col]])
#   g_norm <- tolower(trimws(g_raw))
#   okg <- !is.na(g_norm) & g_norm != "" & g_norm != "unknown"
#   g_raw <- g_raw[okg]
#   g_norm <- g_norm[okg]

#   # p (use)
#   p_use <- suppressWarnings(as.numeric(gene_results[[p_use_col]]))[okg]
#   okp <- is.finite(p_use) & !is.na(p_use)
#   g_raw <- g_raw[okp]; g_norm <- g_norm[okp]; p_use <- p_use[okp]

#   # z
#   z_use <- NULL
#   if (!is.null(z_col)) {
#     if (!(z_col %in% names(gene_results))) stop("z_col not found: ", z_col, call. = FALSE)
#     z_use <- suppressWarnings(as.numeric(gene_results[[z_col]]))[okg][okp]
#   }

#   # weights
#   w_use <- NULL
#   if (!is.null(weight_col)) {
#     if (!(weight_col %in% names(gene_results))) stop("weight_col not found: ", weight_col, call. = FALSE)
#     w_use <- suppressWarnings(as.numeric(gene_results[[weight_col]]))[okg][okp]
#   }

#   # dedup by normalized gene id (keep first)
#   keep1 <- !duplicated(g_norm)
#   g_raw <- g_raw[keep1]; g_norm <- g_norm[keep1]; p_use <- p_use[keep1]
#   if (!is.null(z_use)) z_use <- z_use[keep1]
#   if (!is.null(w_use)) w_use <- w_use[keep1]

#   gene_map <- stats::setNames(g_raw, g_norm)
#   gene_p   <- stats::setNames(p_use, g_norm)
#   gene_z   <- if (!is.null(z_use)) stats::setNames(z_use, g_norm) else NULL
#   gene_w   <- if (!is.null(w_use)) stats::setNames(w_use, g_norm) else NULL

#   # coerce pathways to long
#   if (is.data.frame(pathways)) {
#     if (!all(c("pathway_id","gene_id") %in% names(pathways))) {
#       stop("pathways df must have pathway_id and gene_id.", call. = FALSE)
#     }
#     if (!("pathway_name" %in% names(pathways))) pathways$pathway_name <- pathways$pathway_id
#     pw_long <- pathways[, c("pathway_id","pathway_name","gene_id"), drop = FALSE]
#     pw_long$gene_id <- as.character(pw_long$gene_id)
#   } else if (is.list(pathways)) {
#     pid <- names(pathways); if (is.null(pid)) pid <- paste0("PWY_", seq_along(pathways))
#     pw_long <- data.frame(
#       pathway_id   = rep(pid, lengths(pathways)),
#       pathway_name = rep(pid, lengths(pathways)),
#       gene_id      = unlist(pathways, use.names = FALSE),
#       stringsAsFactors = FALSE
#     )
#   } else {
#     stop("pathways must be a data.frame or list.", call. = FALSE)
#   }

#   pw_long$gene_id   <- gsub("^gene:", "", pw_long$gene_id, ignore.case = TRUE)
#   pw_long$gene_norm <- tolower(trimws(pw_long$gene_id))
#   pw_long <- pw_long[!is.na(pw_long$gene_norm) & pw_long$gene_norm != "" & pw_long$gene_norm != "unknown", , drop = FALSE]

#   g_list_norm <- split(pw_long$gene_norm, pw_long$pathway_id)
#   pw_name     <- tapply(pw_long$pathway_name, pw_long$pathway_id, function(x) x[1])

#   # intersect with available genes
#   g_list_norm <- lapply(g_list_norm, function(g) unique(g[g %in% names(gene_p)]))

#   # drop small pathways early (optional)
#   min_genes <- as.integer(min_genes)
#   if (!is.na(min_genes) && min_genes >= 1L) {
#     keep_pw <- lengths(g_list_norm) >= min_genes
#     g_list_norm <- g_list_norm[keep_pw]
#   }

#   pid <- names(g_list_norm)

#   # per pathway vectors
#   p_list <- lapply(g_list_norm, function(g) unname(gene_p[g]))
#   z_list <- if (!is.null(gene_z)) lapply(g_list_norm, function(g) unname(gene_z[g])) else NULL
#   w_list <- if (!is.null(gene_w)) lapply(g_list_norm, function(g) unname(gene_w[g])) else NULL
#   g_list_raw <- lapply(g_list_norm, function(g) unname(gene_map[g]))

#   # tau sanitize (allow tau=1)
#   tau_grid <- suppressWarnings(as.numeric(tau_grid))
#   tau_grid <- tau_grid[is.finite(tau_grid) & !is.na(tau_grid) & tau_grid > 0 & tau_grid <= tau_cap]
#   tau_grid <- sort(unique(tau_grid), decreasing = TRUE)
#   if (!length(tau_grid)) tau_grid <- 0.05

#   structure(
#     list(
#       pathway_id   = pid,
#       pathway_name = as.character(pw_name[pid]),
#       g_list_norm  = g_list_norm,
#       g_list_raw   = g_list_raw,
#       p_list       = p_list,
#       z_list       = z_list,
#       w_list       = w_list,
#       gene_p_all_use = p_use,              # for obs pool if needed
#       gene_p_all_raw = if (!is.null(p_raw_col)) suppressWarnings(as.numeric(gene_results[[p_raw_col]])) else NULL,
#       gene_z_all     = if (!is.null(z_col)) suppressWarnings(as.numeric(gene_results[[z_col]])) else NULL,
#       gene_w_all     = if (!is.null(weight_col)) suppressWarnings(as.numeric(gene_results[[weight_col]])) else NULL,
#       gene_col = gene_col,
#       p_use_col = p_use_col,
#       p_raw_col = p_raw_col,
#       p_adj_col = p_adj_col,
#       z_col = z_col,
#       weight_col = weight_col,
#       tau_grid = tau_grid,
#       tau_cap = tau_cap,
#       min_p = min_p,
#       min_genes = min_genes,
#       seed = seed
#     ),
#     class = "magcat_omni2_prep"
#   )
# }

# ## ============================================================
# ## PATCH: fix global sampler + add perm_mode/B_perm interface
# ## Paste this AFTER your previous definitions to override them.
# ## ============================================================

# .magcat_sample_idx_mat <- function(n_pool, d, B) {
#   n_pool <- as.integer(n_pool)
#   d      <- as.integer(d)
#   B      <- as.integer(B)
#   if (n_pool < 1L || d < 1L || B < 1L) return(NULL)

#   if (n_pool >= d) {
#     # IMPORTANT: sample WITHOUT replacement *within each permutation*
#     idx <- replicate(B, sample.int(n_pool, size = d, replace = FALSE))
#     t(idx) # B x d
#   } else {
#     # fallback when pathway is larger than pool: sample WITH replacement
#     matrix(sample.int(n_pool, size = B * d, replace = TRUE), nrow = B, ncol = d)
#   }
# }

# magcat_omni2_run <- function(prep,
#                              omnibus = c("ACAT","minP"),
#                              magma_out = NULL,
#                              include_magma_in_omni = TRUE,
#                              include_magma_in_perm = FALSE,
#                              B_global = 0L,
#                              B_mvn = 0L,
#                              perm_pool = c("raw","obs"),
#                              magma_cor_file = NULL,
#                              magma_cor_pairs = NULL,
#                              make_PD = TRUE,
#                              output = FALSE,
#                              out_dir = "magcat_omni2") {

#   if (!inherits(prep, "magcat_omni2_prep")) stop("prep must come from magcat_omni2_prepare().", call. = FALSE)
#   omnibus   <- match.arg(omnibus)
#   perm_pool <- match.arg(perm_pool)

#   pid <- prep$pathway_id
#   n_genes <- lengths(prep$g_list_norm)

#   res <- data.frame(
#     pathway_id   = pid,
#     pathway_name = prep$pathway_name,
#     n_genes      = as.integer(n_genes),
#     stringsAsFactors = FALSE
#   )

#   ## observed component p-values (analytic)
#   res$acat_p <- vapply(pid, function(k) .magcat_acat_combine(prep$p_list[[k]], min_p = prep$min_p), numeric(1))
#   res$fisher_p <- vapply(pid, function(k) .magcat_fisher(prep$p_list[[k]], min_p = prep$min_p), numeric(1))

#   tf_obs <- lapply(pid, function(k) .magcat_tfisher_adapt(prep$p_list[[k]], prep$tau_grid, min_p = prep$min_p))
#   res$tfisher_p_analytic <- vapply(tf_obs, `[[`, numeric(1), "p")
#   res$tau_hat            <- vapply(tf_obs, `[[`, numeric(1), "tau")
#   res$tfisher_stat_hat   <- vapply(tf_obs, `[[`, numeric(1), "stat")

#   res$minp_p_analytic <- vapply(pid, function(k) .magcat_minp_sidak(prep$p_list[[k]], min_p = prep$min_p), numeric(1))

#   if (!is.null(prep$z_col) && !is.null(prep$z_list)) {
#     res$stouffer_p_analytic <- vapply(pid, function(k) .magcat_stouffer_p(prep$z_list[[k]], prep$w_list[[k]]), numeric(1))
#   } else {
#     res$stouffer_p_analytic <- NA_real_
#   }

#   ## MAGMA competitive (you pass magma_out data.frame "out")
#   res$magma_pvalue <- NA_real_
#   if (!is.null(magma_out)) {
#     if (!is.data.frame(magma_out)) stop("magma_out must be a data.frame (MAGMA competitive output).", call. = FALSE)
#     if (!all(c("pathway_id","magma_pvalue") %in% names(magma_out))) {
#       stop("magma_out must have columns pathway_id and magma_pvalue.", call. = FALSE)
#     }
#     res$magma_pvalue <- magma_out$magma_pvalue[match(res$pathway_id, magma_out$pathway_id)]
#   }

#   ## omnibus analytic (across METHODS)
#   res$omni_p_analytic <- vapply(seq_len(nrow(res)), function(i) {
#     comps <- c(res$acat_p[i],
#                res$fisher_p[i],
#                res$tfisher_p_analytic[i],
#                res$minp_p_analytic[i],
#                res$stouffer_p_analytic[i])

#     if (isTRUE(include_magma_in_omni) && is.finite(res$magma_pvalue[i])) {
#       comps <- c(comps, res$magma_pvalue[i])
#     }
#     .magcat_omni_methods(comps, omnibus = omnibus, min_p = prep$min_p)
#   }, numeric(1))

#   res$omni_p_analytic_BH <- .magcat_bh(res$omni_p_analytic)

#   ## ----------------------------
#   ## global resampling (omnibus only)  [FIXED]
#   ## ----------------------------
#   res$omni_p_global <- NA_real_
#   res$omni_p_global_BH <- NA_real_

#   B_global <- as.integer(B_global)
#   if (B_global > 0L) {

#     if (!is.null(prep$seed)) set.seed(prep$seed)

#     pool_p <- NULL
#     pool_z <- NULL

#     if (perm_pool == "raw" && !is.null(prep$p_raw_col) && !is.null(prep$gene_p_all_raw)) {
#       pool_p <- suppressWarnings(as.numeric(prep$gene_p_all_raw))
#     } else {
#       pool_p <- suppressWarnings(as.numeric(prep$gene_p_all_use))
#     }

#     okp <- is.finite(pool_p) & !is.na(pool_p) & pool_p > 0 & pool_p < 1
#     pool_p <- pool_p[okp]

#     if (!is.null(prep$z_col) && !is.null(prep$gene_z_all)) {
#       pool_z_full <- suppressWarnings(as.numeric(prep$gene_z_all))
#       pool_z_full <- pool_z_full[okp]
#       okz <- is.finite(pool_z_full) & !is.na(pool_z_full)
#       pool_p <- pool_p[okz]
#       pool_z <- pool_z_full[okz]
#     }

#     pool_p <- .magcat_fix_p(pool_p, min_p = prep$min_p)
#     if (length(pool_p) < 2L) stop("Global resampling pool has <2 valid p-values.", call. = FALSE)

#     n_pool <- length(pool_p)

#     for (i in seq_len(nrow(res))) {

#       d <- res$n_genes[i]
#       if (!is.finite(d) || d < prep$min_genes) next

#       omni_obs <- res$omni_p_analytic[i]
#       if (!is.finite(omni_obs) || is.na(omni_obs)) next

#       idx_mat <- .magcat_sample_idx_mat(n_pool = n_pool, d = d, B = B_global)
#       if (is.null(idx_mat)) next

#       P <- matrix(pool_p[idx_mat], nrow = B_global, ncol = d)

#       pA <- .magcat_acat_combine_mat(P, min_p = prep$min_p)
#       pF <- stats::pchisq(-2 * rowSums(log(P)), df = 2 * d, lower.tail = FALSE)

#       pT <- rep(NA_real_, B_global)
#       for (b in seq_len(B_global)) {
#         pT[b] <- .magcat_tfisher_adapt(P[b, ], prep$tau_grid, min_p = prep$min_p)$p
#       }

#       pM <- {
#         pmin <- apply(P, 1, min)
#         1 - (1 - pmin)^d
#       }

#       pS <- rep(NA_real_, B_global)
#       if (!is.null(pool_z)) {
#         Z <- matrix(pool_z[idx_mat], nrow = B_global, ncol = d)
#         Zs <- rowSums(Z) / sqrt(d)
#         pS <- 2 * stats::pnorm(-abs(Zs))
#       }

#       omni_null <- vapply(seq_len(B_global), function(b) {
#         comps <- c(pA[b], pF[b], pT[b], pM[b], pS[b])

#         # Not recommended: including fixed MAGMA p in null
#         if (isTRUE(include_magma_in_perm) && isTRUE(include_magma_in_omni) && is.finite(res$magma_pvalue[i])) {
#           comps <- c(comps, res$magma_pvalue[i])
#         }

#         .magcat_omni_methods(comps, omnibus = omnibus, min_p = prep$min_p)
#       }, numeric(1))

#       res$omni_p_global[i] <- (1 + sum(omni_null <= omni_obs, na.rm = TRUE)) / (B_global + 1)
#     }

#     res$omni_p_global_BH <- .magcat_bh(res$omni_p_global)
#   }

#   ## ----------------------------
#   ## MVN resampling (LD-aware)
#   ## ----------------------------
#   res$omni_p_mvn <- NA_real_
#   res$omni_p_mvn_BH <- NA_real_

#   B_mvn <- as.integer(B_mvn)
#   if (B_mvn > 0L) {

#     cor_pairs <- .magma_read_gene_cor_pairs(magma_cor_file = magma_cor_file, magma_cor_pairs = magma_cor_pairs)
#     if (is.null(cor_pairs) || !nrow(cor_pairs)) {
#       stop("MVN resampling requires magma_cor_file or magma_cor_pairs (gene1, gene2, r).", call. = FALSE)
#     }

#     if (!is.null(prep$seed)) set.seed(prep$seed)

#     for (i in seq_len(nrow(res))) {

#       d <- res$n_genes[i]
#       if (!is.finite(d) || d < prep$min_genes) next

#       omni_obs <- .magcat_omni_methods(
#         c(res$acat_p[i], res$fisher_p[i], res$tfisher_p_analytic[i], res$minp_p_analytic[i], res$stouffer_p_analytic[i]),
#         omnibus = omnibus, min_p = prep$min_p
#       )
#       if (!is.finite(omni_obs) || is.na(omni_obs)) next

#       genes_S <- prep$g_list_raw[[res$pathway_id[i]]]
#       genes_S <- as.character(genes_S)
#       genes_S <- genes_S[!is.na(genes_S) & genes_S != ""]
#       if (length(genes_S) < prep$min_genes) next

#       R <- .magma_build_R_from_pairs(genes = genes_S, cor_pairs = cor_pairs, make_PD = make_PD)
#       if (is.null(R)) next

#       Z <- .magma_simulate_Z_from_R(R, B = B_mvn)
#       if (is.null(Z)) next

#       P <- 2 * stats::pnorm(-abs(Z))  # two-sided p

#       pA <- .magcat_acat_combine_mat(P, min_p = prep$min_p)
#       pF <- stats::pchisq(-2 * rowSums(log(P)), df = 2 * ncol(P), lower.tail = FALSE)

#       pT <- rep(NA_real_, B_mvn)
#       for (b in seq_len(B_mvn)) {
#         pT[b] <- .magcat_tfisher_adapt(P[b, ], prep$tau_grid, min_p = prep$min_p)$p
#       }

#       pM <- {
#         k <- ncol(P)
#         pmin <- apply(P, 1, min)
#         1 - (1 - pmin)^k
#       }

#       pS <- {
#         k <- ncol(Z)
#         Zs <- rowSums(Z) / sqrt(k)
#         2 * stats::pnorm(-abs(Zs))
#       }

#       omni_null <- vapply(seq_len(B_mvn), function(b) {
#         comps <- c(pA[b], pF[b], pT[b], pM[b], pS[b])

#         if (isTRUE(include_magma_in_perm) && isTRUE(include_magma_in_omni) && is.finite(res$magma_pvalue[i])) {
#           comps <- c(comps, res$magma_pvalue[i])
#         }

#         .magcat_omni_methods(comps, omnibus = omnibus, min_p = prep$min_p)
#       }, numeric(1))

#       res$omni_p_mvn[i] <- (1 + sum(omni_null <= omni_obs, na.rm = TRUE)) / (B_mvn + 1)
#     }

#     res$omni_p_mvn_BH <- .magcat_bh(res$omni_p_mvn)
#   }

#   ## final p priority: MVN > global > analytic
#   res$omni_p_final <- res$omni_p_analytic
#   if (B_global > 0L) res$omni_p_final <- res$omni_p_global
#   if (B_mvn > 0L)    res$omni_p_final <- res$omni_p_mvn
#   res$omni_p_final_BH <- .magcat_bh(res$omni_p_final)

#   res <- res[order(res$omni_p_final, decreasing = FALSE, na.last = TRUE), , drop = FALSE]

#   if (isTRUE(output)) {
#     if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#     tag <- paste0("omni_", tolower(omnibus),
#                   if (B_mvn > 0L) "_mvn" else if (B_global > 0L) "_global" else "_analytic")
#     out_path <- file.path(out_dir, paste0(tag, ".csv"))
#     utils::write.csv(res, out_path, row.names = FALSE)
#     attr(res, "file") <- out_path
#   }

#   res
# }

# magcat_omni2_pathways <- function(gene_results,
#                                   pathways = NULL,
#                                   species = NULL,
#                                   pmn_gene_col = NULL,
#                                   gene_col = "GENE",
#                                   p_raw_col = "P",
#                                   p_adj_col = NULL,
#                                   z_col = NULL,
#                                   weight_col = NULL,
#                                   tau_grid = c(0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.001),
#                                   tau_cap = 1,
#                                   min_p = 1e-15,
#                                   min_genes = 2L,
#                                   magma_out = NULL,
#                                   include_magma_in_omni = TRUE,
#                                   include_magma_in_perm = FALSE,
#                                   omnibus = c("ACAT","minP"),
#                                   ## NEW:
#                                   B_perm = NULL,
#                                   perm_mode = c("none","global","mvn","both"),
#                                   ## BACKCOMPAT:
#                                   B_global = 0L,
#                                   B_mvn = 0L,
#                                   perm_pool = c("raw","obs"),
#                                   magma_cor_file = NULL,
#                                   magma_cor_pairs = NULL,
#                                   make_PD = TRUE,
#                                   seed = NULL,
#                                   output = FALSE,
#                                   out_dir = "magcat_omni2") {

#   omnibus   <- match.arg(omnibus)
#   perm_pool <- match.arg(perm_pool)
#   perm_mode <- match.arg(perm_mode)

#   # If user uses the new interface, map it
#   if (!is.null(B_perm)) {
#     B_perm <- as.integer(B_perm)
#     if (perm_mode == "global") {
#       B_global <- B_perm; B_mvn <- 0L
#     } else if (perm_mode == "mvn") {
#       B_global <- 0L; B_mvn <- B_perm
#     } else if (perm_mode == "both") {
#       B_global <- B_perm; B_mvn <- B_perm
#     } else {
#       B_global <- 0L; B_mvn <- 0L
#     }
#   }

#   prep <- magcat_omni2_prepare(
#     gene_results = gene_results,
#     pathways = pathways,
#     species = species,
#     pmn_gene_col = pmn_gene_col,
#     gene_col = gene_col,
#     p_raw_col = p_raw_col,
#     p_adj_col = p_adj_col,
#     z_col = z_col,
#     weight_col = weight_col,
#     tau_grid = tau_grid,
#     tau_cap = tau_cap,
#     min_p = min_p,
#     min_genes = min_genes,
#     seed = seed
#   )

#   magcat_omni2_run(
#     prep = prep,
#     omnibus = omnibus,
#     magma_out = magma_out,
#     include_magma_in_omni = include_magma_in_omni,
#     include_magma_in_perm = include_magma_in_perm,
#     B_global = B_global,
#     B_mvn = B_mvn,
#     perm_pool = perm_pool,
#     magma_cor_file = magma_cor_file,
#     magma_cor_pairs = magma_cor_pairs,
#     make_PD = make_PD,
#     output = output,
#     out_dir = out_dir
#   )
# }

# ## ============================================================
# ## HOW TO CALL (what you expected)
# ## ============================================================
# ## 1) Global resampling only:
# # omni <- magcat_omni2_pathways(
# #   gene_results = genes_adj, species="maize", gene_col="GENE",
# #   p_raw_col="P", p_adj_col="P_adj", z_col="Z_adj",
# #   magma_out = out,
# #   omnibus="ACAT",
# #   B_perm = 1000L, perm_mode="global",
# #   perm_pool="raw",
# #   seed=123, output=TRUE, out_dir="magcat_omni2_results"
# # )

# ## 2) MVN only (LD-aware) — REQUIRES magma_cor_file (or magma_cor_pairs):
# # omni <- magcat_omni2_pathways(
# #   gene_results = genes_adj, species="maize", gene_col="GENE",
# #   p_raw_col="P", p_adj_col="P_adj", z_col="Z_adj",
# #   magma_out = out,
# #   omnibus="ACAT",
# #   B_perm = 1000L, perm_mode="mvn",
# #   magma_cor_file="magma_gene_cor_pairs.txt",
# #   make_PD=TRUE,
# #   seed=123, output=TRUE, out_dir="magcat_omni2_results"
# # )

# ## 3) BOTH global + MVN:
# # omni <- magcat_omni2_pathways(
# #   gene_results = genes_adj, species="maize", gene_col="GENE",
# #   p_raw_col="P", p_adj_col="P_adj", z_col="Z_adj",
# #   magma_out = out,
# #   omnibus="ACAT",
# #   B_perm = 1000L, perm_mode="both",
# #   perm_pool="raw",
# #   magma_cor_file="magma_gene_cor_pairs.txt",
# #   make_PD=TRUE,
# #   seed=123, output=TRUE, out_dir="magcat_omni2_results"
# # )
