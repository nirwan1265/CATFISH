# =============================================================================
# Westfall--Young family-wise calibration for CATFISH
#
# One genome-wide correlated null
# draw is scored through EVERY pathway, so pathways that share genes are
# perturbed jointly. The resulting B x K matrix of null omnibus statistics is
# read two ways:
#   * down each column  -> per-pathway calibrated p-value
#   * across each row   -> minimum over pathways -> Westfall--Young FWER p-value
#
# All component statistics are computed in log space (the GLGC lipid GWAS has
# gene-level p < 1e-300), with the five numerical safeguards of the CATFISH-WY
# paper (Section 2.4) and the soft TFisher definition of Section 2.1:
#   * Wo = min over a FIXED tau grid of the right-tail p (stat.soft.omni$omni),
#     not TFisher::p.soft.omni, which is miscalibrated in both tails;
#   * an empty soft-thresholding statistic is undefined (NA) and the pathway is
#     combined from the components that are defined;
#   * Wo is capped at 1e-16, the resolution of 1 - CDF in double precision.
# The tau grid must be fixed in advance: a grid estimated from all gene-level
# p-values makes every pathway's statistic depend on genes outside it and
# breaks the subset pivotality that gives strong FWER control.
# =============================================================================

## ---------------------------------------------------------------- helpers ---

.wy_logsumexp <- function(lx) {
  m <- max(lx)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(lx - m)))
}

## ACAT on log p-values, returning log(p_ACAT).
##   T = mean tan((1/2 - p) pi). For p < 1e-15, tan(.) ~ 1/(p pi) and those
##   terms are accumulated by log-sum-exp; for T > 1e6, p ~ 1/(T pi), avoiding
##   cancellation near 1/2; and T may be NEGATIVE when most inputs exceed 1/2,
##   which is handled directly rather than clamped (clamping sets p = 0.5 and
##   destroys the ordering of the upper half).
.wy_acat_log <- function(lp) {
  lp <- lp[is.finite(lp)]
  m  <- length(lp)
  if (!m) return(NA_real_)
  lp <- pmin(lp, log1p(-1e-15))          # keep p < 1 so tan() stays finite
  small <- lp < log(1e-15)
  S_big <- if (any(!small)) sum(tan((0.5 - exp(lp[!small])) * pi)) else 0
  if (any(small)) {
    logS_small <- .wy_logsumexp(-lp[small] - log(pi))
    ## |S_big| <= m / (1e-15 pi); once S_small dwarfs that, S is S_small
    if (logS_small > log(m) + log(1 / (1e-15 * pi)) + 40)
      return(-(logS_small - log(m)) - log(pi))
    S <- exp(logS_small) + S_big
  } else {
    S <- S_big
  }
  Tm <- S / m
  if (Tm > 1e6) return(-log(Tm) - log(pi))
  log(0.5 - atan(Tm) / pi)
}

.wy_fisher_log <- function(lp) {
  stats::pchisq(-2 * sum(lp), df = 2 * length(lp),
                lower.tail = FALSE, log.p = TRUE)
}

.wy_stouffer_log <- function(lp) {
  z <- stats::qnorm(lp, lower.tail = FALSE, log.p = TRUE)
  z[!is.finite(z)] <- 0
  stats::pnorm(sum(z) / sqrt(length(z)), lower.tail = FALSE, log.p = TRUE)
}

## minP with Sidak inversion 1 - (1 - q)^m, which underflows once q is below
## about 1e-16 / m and would silently drop minP from the omnibus; the m * q
## asymptote is used whenever it is below 1e-8.
.wy_minp_log <- function(lp) {
  m <- length(lp); lmin <- min(lp)
  approx <- lmin + log(m)
  if (approx < log(1e-8)) return(approx)
  r <- log1p(-exp(m * log1p(-exp(lmin))))
  if (is.finite(r)) r else approx
}

## Adaptive soft TFisher as Wo, with the empty case undefined and a precision cap.
.wy_tfisher_log <- function(lp, tau_grid, tf_min_p = 1e-16) {
  p  <- pmin(pmax(exp(lp), 1e-300), 1 - 1e-15)
  Wo <- tryCatch(as.numeric(TFisher::stat.soft.omni(p = p, TAU1 = tau_grid,
                                                    M = NULL)$omni),
                 error = function(e) NA_real_)
  if (!is.finite(Wo) || Wo >= 1 - 1e-12) return(NA_real_)   # empty: undefined
  log(max(Wo, tf_min_p))                                      # precision cap
}

## Omnibus = ACAT across the component log p-values that are defined.
.wy_omni_log <- function(lp, use_tfisher = TRUE, tau_grid = c(0.01, 0.05, 0.5)) {
  comps <- c(.wy_acat_log(lp), .wy_fisher_log(lp),
             .wy_stouffer_log(lp), .wy_minp_log(lp))
  if (isTRUE(use_tfisher)) comps <- c(comps, .wy_tfisher_log(lp, tau_grid))
  comps <- comps[is.finite(comps)]
  if (!length(comps)) return(0)
  .wy_acat_log(comps)
}

## Per-chromosome correlation matrices and Cholesky factors from gene pairs.
## Unreported pairs are zero; genes on different chromosomes are independent.
## A matrix that is not positive definite has its eigenvalues floored at 1e-6
## and its diagonal renormalised.
.wy_chol_by_chr <- function(gene_id, gene_chr, cor_pairs, verbose = TRUE) {
  chr_of <- stats::setNames(gene_chr, gene_id)
  cp_by_chr <- list()
  if (!is.null(cor_pairs) && nrow(cor_pairs)) {
    cp <- cor_pairs
    cp$chr <- chr_of[cp$g1]
    cp <- cp[!is.na(cp$chr) & !is.na(chr_of[cp$g2]), , drop = FALSE]
    if (nrow(cp)) cp_by_chr <- split(cp[, c("g1", "g2", "r")], cp$chr)
  }
  chroms <- unique(gene_chr)
  Ut <- list(); idx <- list(); repaired <- character(0)
  for (cc in chroms) {
    ii  <- which(gene_chr == cc)
    ids <- gene_id[ii]
    loc <- stats::setNames(seq_along(ids), ids)
    m   <- length(ids)
    R   <- diag(1, m)
    sub <- cp_by_chr[[cc]]
    npair <- if (is.null(sub)) 0L else nrow(sub)
    if (npair) {
      a <- loc[sub$g1]; b <- loc[sub$g2]; rr <- as.numeric(sub$r)
      ok <- !is.na(a) & !is.na(b) & is.finite(rr)
      R[cbind(a[ok], b[ok])] <- rr[ok]
      R[cbind(b[ok], a[ok])] <- rr[ok]
    }
    U <- tryCatch(chol(R), error = function(e) NULL)
    if (is.null(U)) {
      ev <- eigen(R, symmetric = TRUE)
      ev$values[ev$values < 1e-6] <- 1e-6
      R <- ev$vectors %*% diag(ev$values, nrow = m) %*% t(ev$vectors)
      d <- sqrt(diag(R)); R <- R / outer(d, d)
      U <- chol(R)
      repaired <- c(repaired, cc)
    }
    if (isTRUE(verbose))
      message(sprintf("    chr %-3s genes=%-5d pairs=%-7d%s", cc, m, npair,
                      if (cc %in% repaired) " (eigen-repaired)" else ""))
    Ut[[cc]]  <- t(U)            # draw as t(U) %*% eps
    idx[[cc]] <- ii
  }
  covered <- chroms[vapply(chroms, function(cc) {
    L <- Ut[[cc]]; any(L[lower.tri(L)] != 0)
  }, logical(1))]
  list(Ut = Ut, idx = idx, chroms = chroms, repaired = repaired, covered = covered)
}

## The engine: observed statistics, the genome-wide null, and both readings.
.wy_engine <- function(gene_id, gene_lp, gene_chr, paths, cor_pairs,
                       B, use_tfisher, tau_grid, seed, verbose) {
  if (isTRUE(verbose)) message("[WY] per-chromosome correlation and Cholesky factors")
  ch <- .wy_chol_by_chr(gene_id, gene_chr, cor_pairs, verbose)
  if (length(ch$covered) < length(ch$chroms)) {
    warning(sprintf(paste0(
      "gene correlations cover only %d of %d chromosomes (%s); genes on the ",
      "rest are treated as INDEPENDENT, so the LD-aware null is incomplete."),
      length(ch$covered), length(ch$chroms), paste(ch$covered, collapse = ",")),
      call. = FALSE)
  }

  if (isTRUE(verbose)) message("[WY] observed statistics")
  obs_log <- vapply(paths, function(ii)
    .wy_omni_log(gene_lp[ii], use_tfisher, tau_grid), numeric(1))

  if (isTRUE(verbose))
    message(sprintf("[WY] null: B = %d genome-wide draws, %d pathways, TFisher = %s",
                    B, length(paths), use_tfisher))
  if (!is.null(seed)) set.seed(seed)
  NG <- length(gene_id); NP <- length(paths)
  minp_null <- numeric(B)
  null_mat  <- matrix(NA_real_, nrow = B, ncol = NP)
  at <- max(1L, B %/% 20L)
  for (b in seq_len(B)) {
    znull <- numeric(NG)
    for (cc in ch$chroms) {
      L <- ch$Ut[[cc]]
      znull[ch$idx[[cc]]] <- as.vector(L %*% stats::rnorm(ncol(L)))
    }
    lpnull <- stats::pnorm(znull, lower.tail = FALSE, log.p = TRUE)
    st <- vapply(paths, function(ii)
      .wy_omni_log(lpnull[ii], use_tfisher, tau_grid), numeric(1))
    minp_null[b] <- min(st, na.rm = TRUE)
    null_mat[b, ] <- st
    if (isTRUE(verbose) && b %% at == 0L)
      message(sprintf("    %d/%d (%.0f%%)", b, B, 100 * b / B))
  }

  rd <- .wy_read_null(obs_log, null_mat, minp_null)
  list(obs_log = obs_log, z_std = rd$z_std, p_cal = rd$p_cal,
       p_fwer_wy_std = rd$p_fwer_std, p_fwer_wy = rd$p_fwer,
       minp_null = minp_null, null_mat = null_mat, chol = ch)
}

## Both readings of a B x K null matrix of pathway log p-values: each pathway's
## own column (calibrated p) and the family minimum in each replicate
## (Westfall--Young), on the standardised and on the raw scale.
.wy_read_null <- function(obs_log, null_mat, minp_null) {
  B <- nrow(null_mat); NP <- ncol(null_mat)
  ## per-pathway standardisation against each pathway's own null, so the
  ## family minimum compares pathways of different size on a common scale
  mu  <- colMeans(null_mat)
  sdv <- apply(null_mat, 2, stats::sd)
  sdv[!is.finite(sdv) | sdv <= 0] <- NA_real_
  z_obs  <- (obs_log - mu) / sdv
  z_null <- sweep(sweep(null_mat, 2, mu, "-"), 2, sdv, "/")
  ok <- is.finite(z_obs)
  minz_null <- apply(z_null[, ok, drop = FALSE], 1, min, na.rm = TRUE)

  ## counting form: null minima equal to the observed value count as extreme
  p_fwer_std <- rep(NA_real_, NP)
  p_fwer_std[ok] <- vapply(z_obs[ok], function(o)
    (1 + sum(minz_null <= o)) / (B + 1), numeric(1))
  p_fwer <- vapply(obs_log, function(o)
    (1 + sum(minp_null <= o)) / (B + 1), numeric(1))
  p_cal <- vapply(seq_len(NP), function(k) {
    v <- null_mat[, k]; v <- v[is.finite(v)]
    (1 + sum(v <= obs_log[k])) / (length(v) + 1)
  }, numeric(1))
  list(z_std = z_obs, p_fwer_std = p_fwer_std, p_fwer = p_fwer, p_cal = p_cal)
}

## Competitive null by circular genomic permutation. Genes are laid out in
## genome order and, in each replicate, every gene takes the observed score of
## the gene s places further along (wrapping at the end), for one random shift s
## shared by all pathways. A pathway's genes therefore keep their spacing, so
## neighbouring genes in LD stay together, while their scores come from a random
## place in the genome. The question is whether the pathway's genes score higher
## than genes elsewhere in the same GWAS: genome-wide polygenic signal is in the
## null as well as in the observed statistic. Shifts are drawn without
## replacement, so at most (number of genes - 1) replicates exist.
.wy_competitive <- function(gene_lp, ord, paths, obs_log, B, use_tfisher, tau_grid,
                            seed, verbose) {
  NG <- length(gene_lp); NP <- length(paths)
  if (NG < 3L) stop("the competitive test needs at least 3 genes.", call. = FALSE)
  B_eff <- min(as.integer(B), NG - 1L)
  if (isTRUE(verbose)) {
    message(sprintf("[WY competitive] circular genomic shift: %d replicates over %d genes",
                    B_eff, NG))
    if (B_eff < B)
      message(sprintf("[WY competitive] only %d distinct shifts exist; B reduced from %d",
                      B_eff, B))
  }
  lp_ord <- gene_lp[ord]                       # scores in genome order
  pos <- integer(NG); pos[ord] <- seq_len(NG)  # genome position of each gene
  if (!is.null(seed)) set.seed(seed + 1L)      # a different stream from the MVN null
  shifts <- sample.int(NG - 1L, B_eff)
  minp_null <- numeric(B_eff)
  null_mat  <- matrix(NA_real_, nrow = B_eff, ncol = NP)
  at <- max(1L, B_eff %/% 20L)
  for (b in seq_len(B_eff)) {
    lpn <- lp_ord[((pos - 1L + shifts[b]) %% NG) + 1L]
    st <- vapply(paths, function(ii)
      .wy_omni_log(lpn[ii], use_tfisher, tau_grid), numeric(1))
    minp_null[b] <- min(st, na.rm = TRUE)
    null_mat[b, ] <- st
    if (isTRUE(verbose) && b %% at == 0L)
      message(sprintf("    %d/%d (%.0f%%)", b, B_eff, 100 * b / B_eff))
  }
  rd <- .wy_read_null(obs_log, null_mat, minp_null)
  c(rd, list(B = B_eff, minp_null = minp_null))
}

## Read a gene-pair correlation table (gene1 gene2 r) from a file or data.frame.
.wy_read_cor_pairs <- function(cor_pairs = NULL, cor_file = NULL) {
  if (is.null(cor_pairs)) {
    if (is.null(cor_file)) return(NULL)
    if (!file.exists(cor_file)) stop("cor_file does not exist: ", cor_file, call. = FALSE)
    cor_pairs <- if (requireNamespace("data.table", quietly = TRUE)) {
      as.data.frame(data.table::fread(cor_file, showProgress = FALSE))
    } else {
      utils::read.delim(cor_file, header = TRUE, stringsAsFactors = FALSE)
    }
  }
  if (!is.data.frame(cor_pairs) || ncol(cor_pairs) < 3)
    stop("gene correlations must have three columns: gene1, gene2, r.", call. = FALSE)
  cp <- cor_pairs[, 1:3]
  names(cp) <- c("g1", "g2", "r")
  cp$g1 <- tolower(trimws(as.character(cp$g1)))
  cp$g2 <- tolower(trimws(as.character(cp$g2)))
  cp$r  <- suppressWarnings(as.numeric(cp$r))
  cp
}

## Gene-level genomic control: centre the gene z-scores of one GWAS at their
## median. A GWAS can carry a genome-wide shift in its gene z-scores that differs
## from one GWAS to another and is shared across chromosomes, so the correlated
## null, which centres every gene at zero, does not contain it; a self-contained
## test would read it as pathway signal. z is taken from log p, so p-values below
## 1e-300 keep their order, and the centred p-value is returned on the log scale.
## The order of genes is unchanged.
.catfish_gc <- function(lp, apply = TRUE) {
  z <- stats::qnorm(lp, lower.tail = FALSE, log.p = TRUE)
  med <- stats::median(z)
  if (!isTRUE(apply)) return(list(lp = lp, shift = 0, median_z = med))
  list(lp = stats::pnorm(z - med, lower.tail = FALSE, log.p = TRUE), shift = med,
       median_z = med)
}

## What to tell the user about the centring. A large positive median z usually
## means the median gene is itself associated (a saturated, highly polygenic
## GWAS), and centring then turns the self-contained question part of the way
## into a competitive one; a large negative one means deflated gene p-values.
.catfish_gc_message <- function(med, gc, tag) {
  if (!isTRUE(gc))
    return(sprintf("[%s] genomic control off (gc = FALSE): median gene z %.3f left in", tag, med))
  msg <- sprintf("[%s] genomic control: median gene z %.3f removed", tag, med)
  if (med > 0.5)
    msg <- paste0(msg, sprintf(paste0(
      "\n[%s] note: a median gene z of %.2f suggests the median gene is itself associated. ",
      "Centring compares each pathway with the median gene, which is closer to a ",
      "competitive test; use gc = FALSE for the strict self-contained test."), tag, med))
  if (med < -0.5)
    msg <- paste0(msg, sprintf(paste0(
      "\n[%s] note: a median gene z of %.2f means the gene p-values are deflated; ",
      "check the calibration of the GWAS."), tag, med))
  msg
}

## --------------------------------------------------------------- exported ---

#' Pathway testing with Westfall--Young family-wise calibration
#'
#' Scores every pathway with the CATFISH omnibus (ACAT across ACAT, Fisher,
#' Stouffer, minP and, optionally, adaptive soft TFisher) and calibrates it
#' against a genome-wide correlated null. Each of \code{B} null replicates
#' draws one vector of gene-level \eqn{z}-scores from a multivariate normal
#' model of gene--gene correlation, chromosome by chromosome, and scores every
#' pathway on that same draw, so pathways that share genes are perturbed
#' jointly. The resulting null matrix gives each pathway a calibrated
#' marginal \eqn{p}-value (its own column) and a Westfall--Young family-wise
#' \eqn{p}-value (the minimum across pathways in each replicate).
#'
#' All statistics are computed in log space, so gene-level \eqn{p}-values far
#' below \eqn{10^{-300}} keep their ordering. Soft TFisher contributes
#' \eqn{W_o}, the minimum over \code{tau_grid} of the single-threshold
#' right-tail \eqn{p}; it is treated as undefined, and omitted from that
#' pathway's omnibus, when no gene falls below any threshold, and is capped at
#' \eqn{10^{-16}}, the resolution of its tail computation.
#'
#' Family-wise control holds in the strong sense under the multivariate normal
#' model, because each pathway's statistic depends only on its own genes. This
#' requires \code{tau_grid} to be fixed in advance rather than estimated from
#' the gene-level \eqn{p}-values.
#'
#' By default (\code{gc = TRUE}) gene-level genomic control is applied. Each gene
#' \eqn{p}-value is converted to \eqn{z_g = \Phi^{-1}(1 - p_g)} (from
#' \eqn{\log p_g}), the median over all genes analysed is subtracted, and the
#' result is converted back. A GWAS can carry a genome-wide shift in its gene
#' \eqn{z}-scores that is shared across chromosomes; the correlated null does not
#' contain it, and a self-contained test would read it as pathway signal.
#' Centring is one constant per GWAS: it leaves the order of genes and the
#' gene--gene correlations unchanged and moves the median gene \eqn{p}-value to
#' 0.5. Pass the results of the whole genome, not a subset, because the median is
#' taken over the genes supplied.
#'
#' Centring cannot tell a technical shift from real polygenic signal. In a highly
#' polygenic, well-powered GWAS the median gene is itself associated, and
#' centring then asks whether a pathway is more associated than the median gene,
#' which is closer to a competitive test; a note is printed when the median gene
#' \eqn{z} exceeds 0.5. Set \code{gc = FALSE} for the strict self-contained test,
#' for gene results that do not cover the whole genome, for gene \eqn{p}-values
#' already corrected in another way, or to reproduce uncentred results.
#'
#' \strong{Competitive test.} The test above is self-contained: it asks whether
#' a pathway carries any association. With \code{competitive = TRUE} a second,
#' competitive test asks whether the pathway's genes are more associated than
#' genes elsewhere in the same GWAS. It uses the same omnibus statistic and the
#' same two readings (calibrated and Westfall--Young family-wise), with a
#' different null: genes are put in genome order (\code{chr_col}, then
#' \code{pos_col}) and in each replicate every gene takes the observed score of
#' the gene a random number of places further along the genome, wrapping at the
#' end, with one shift shared by all pathways (circular genomic permutation). A
#' pathway's genes keep their spacing, so genes in LD with each other stay
#' together, and genome-wide polygenic signal is present in the null as well
#' as in the observed statistic. It uses the gene \eqn{p}-values the
#' self-contained test uses (after genomic control when \code{gc = TRUE}). There
#' are at most (number of genes \eqn{-} 1) distinct shifts, which caps the
#' number of replicates and the resolution of the competitive \eqn{p}-values.
#' The null does not match genes on size or SNP count.
#'
#' @param gene_results data.frame of gene-level results (e.g. MAGMA
#'   \code{.genes.out}), with columns for gene ID, \eqn{p}-value and chromosome.
#' @param pathways Named list of gene-ID vectors, or a data.frame with
#'   \code{pathway_id} and \code{gene_id}. Provide this or \code{species}.
#' @param species Optional; built-in PMN pathways via
#'   \code{catfish_load_pathways()}.
#' @param pmn_gene_col Passed to \code{catfish_load_pathways()} when
#'   \code{species} is used.
#' @param gene_col,p_col,chr_col Columns of \code{gene_results} holding the
#'   gene ID, gene-level \eqn{p}-value and chromosome.
#' @param z_col Optional column (e.g. MAGMA \code{ZSTAT}); if given, genes with
#'   a non-finite value are dropped, as MAGMA output is usually filtered.
#' @param cor_pairs data.frame of gene pairs and correlations
#'   (gene1, gene2, r), e.g. from \code{magma_genesraw_to_cor_pairs_banded()}.
#' @param cor_file Path to the same table on disk; used if \code{cor_pairs} is
#'   \code{NULL}.
#' @param B Number of genome-wide null replicates.
#' @param use_tfisher Include adaptive soft TFisher as a fifth component.
#' @param tau_grid Fixed soft TFisher truncation grid. The default is the grid
#'   of Zhang et al. (2020) without its \eqn{\tau=1} corner, which is Fisher.
#' @param min_genes,max_genes Pathway size limits, counted over genes present
#'   in \code{gene_results}.
#' @param seed Random seed for the null draws (\code{NULL} leaves it unset).
#' @param gc Gene-level genomic control: centre the gene \eqn{z}-scores at their
#'   median before testing (default \code{TRUE}). See Details.
#' @param competitive Also run the competitive test (default \code{FALSE}).
#'   See Details.
#' @param pos_col Column of \code{gene_results} with the gene position (e.g.
#'   MAGMA \code{START}), used to put genes in genome order for the competitive
#'   test. If \code{NULL}, genes are taken in the order supplied within each
#'   chromosome, and a message says so.
#' @param B_comp Replicates of the competitive null (default \code{B}); at most
#'   the number of genes minus one.
#' @param verbose Print progress.
#' @param output If \code{TRUE}, write a CSV to \code{out_dir}.
#' @param out_dir Output directory when \code{output = TRUE}.
#'
#' @return A data.frame, one row per pathway, sorted by \code{z_std}:
#' \describe{
#'   \item{pathway_id}{Pathway identifier}
#'   \item{pathway_name}{Pathway name (the \code{pathway_name} column of a
#'     data.frame \code{pathways}, or the built-in set; the ID otherwise)}
#'   \item{z_std}{Observed omnibus standardised against the pathway's own null
#'     (more negative = more extreme)}
#'   \item{p_fwer_std}{Westfall--Young family-wise \eqn{p}-value on the
#'     standardised scale; the headline family-wise result}
#'   \item{n_genes}{Genes in the pathway present in \code{gene_results}}
#'   \item{log10_p_omnibus, p_omnibus}{Uncalibrated analytic omnibus}
#'   \item{p_fwer_wy}{Westfall--Young family-wise \eqn{p}-value on the raw
#'     scale (valid when pathways share one marginal null)}
#'   \item{p_bonferroni}{Bonferroni applied to the analytic omnibus, for
#'     comparison}
#'   \item{p_calibrated}{Per-pathway calibrated \eqn{p}-value}
#'   \item{gc_shift}{Median gene \eqn{z}-score removed by genomic control
#'     (the same for every row; 0 when \code{gc = FALSE})}
#'   \item{z_std_comp, p_fwer_comp, p_calibrated_comp}{With
#'     \code{competitive = TRUE}: the standardised statistic, family-wise
#'     \eqn{p}-value and calibrated \eqn{p}-value of the competitive test}
#' }
#' Attributes \code{wy_threshold_05}, \code{B}, \code{tau_grid},
#' \code{gc_shift} and \code{median_gene_z} (the median gene \eqn{z} before
#' any centring) are attached, and \code{B_comp} with \code{competitive = TRUE}.
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
#' res <- catfish_wy_pathways(genes, pathways = pathways, cor_pairs = cor_pairs,
#'                            B = 200, min_genes = 10, seed = 1, verbose = FALSE)
#' head(res[, c("pathway_id", "n_genes", "p_calibrated", "p_fwer_std")])
#' attr(res, "wy_threshold_05")
#'
#' \dontrun{
#' # MAGMA output (gene-level p-values with chromosome) and MAGMA gene-gene
#' # correlations; use B = 10000 for final results
#' genes <- read.table("trait.genes.out", header = TRUE)
#' res <- catfish_wy_pathways(genes, pathways = my_pathways,
#'                            gene_col = "GENE", p_col = "P", chr_col = "CHR",
#'                            z_col = "ZSTAT", cor_file = "trait.cor_pairs.txt",
#'                            B = 10000, min_genes = 10, max_genes = 200)
#' }
#' @export
catfish_wy_pathways <- function(gene_results,
                                pathways     = NULL,
                                species      = NULL,
                                pmn_gene_col = NULL,
                                gene_col     = "GENE",
                                p_col        = "P",
                                chr_col      = "CHR",
                                z_col        = NULL,
                                cor_pairs    = NULL,
                                cor_file     = NULL,
                                B            = 1000L,
                                use_tfisher  = TRUE,
                                tau_grid     = c(0.01, 0.05, 0.5),
                                min_genes    = 10L,
                                max_genes    = 200L,
                                seed         = 20260706L,
                                gc           = TRUE,
                                competitive  = FALSE,
                                pos_col      = NULL,
                                B_comp       = B,
                                verbose      = TRUE,
                                output       = FALSE,
                                out_dir      = "catfish_wy") {
  if (!is.data.frame(gene_results)) stop("gene_results must be a data.frame.", call. = FALSE)
  for (cl in c(gene_col, p_col, chr_col, pos_col))
    if (!cl %in% names(gene_results)) stop("gene_results is missing column: ", cl, call. = FALSE)
  force(B_comp)
  if (!is.null(pathways) && !is.null(species))
    stop("Provide either 'pathways' OR 'species', not both.", call. = FALSE)
  if (is.null(pathways) && is.null(species))
    stop("Provide 'pathways' (list or data.frame) or 'species'.", call. = FALSE)
  if (isTRUE(use_tfisher) && !requireNamespace("TFisher", quietly = TRUE))
    stop("Package 'TFisher' is required when use_tfisher = TRUE.", call. = FALSE)
  if (is.character(tau_grid))
    stop("tau_grid must be a fixed numeric grid; a data-adaptive grid breaks ",
         "the strong family-wise control of the Westfall-Young procedure.", call. = FALSE)
  tau_grid <- sort(unique(as.numeric(tau_grid)))
  if (!length(tau_grid) || any(!is.finite(tau_grid)) || any(tau_grid <= 0 | tau_grid >= 1))
    stop("tau_grid values must lie in (0, 1).", call. = FALSE)
  B <- as.integer(B)
  if (!is.finite(B) || B < 1L) stop("B must be a positive integer.", call. = FALSE)
  if (is.null(cor_pairs) && is.null(cor_file))
    warning("No gene correlations supplied: every gene is treated as independent.",
            call. = FALSE)

  ## ---- genes (finite p, optional finite z, first occurrence of each ID)
  g <- gene_results
  P <- suppressWarnings(as.numeric(g[[p_col]]))
  keep <- is.finite(P)
  if (!is.null(z_col)) {
    if (!z_col %in% names(g)) stop("z_col not found: ", z_col, call. = FALSE)
    keep <- keep & is.finite(suppressWarnings(as.numeric(g[[z_col]])))
  }
  gid <- tolower(trimws(as.character(g[[gene_col]])))[keep]
  chr <- as.character(g[[chr_col]])[keep]
  pos <- if (is.null(pos_col)) seq_along(gid) else
    suppressWarnings(as.numeric(g[[pos_col]]))[keep]
  P   <- P[keep]
  first <- !duplicated(gid)
  gid <- gid[first]; chr <- chr[first]; P <- P[first]; pos <- pos[first]
  gcr     <- .catfish_gc(log(pmin(pmax(P, 1e-320), 1)), apply = gc)
  gene_lp <- gcr$lp
  idx_of  <- stats::setNames(seq_along(gid), gid)
  if (isTRUE(verbose)) {
    message(sprintf("[WY] genes: %d across %d chromosomes", length(gid),
                    length(unique(chr))))
    message(.catfish_gc_message(gcr$median_z, gc, "WY"))
  }

  ## ---- pathways
  if (is.null(pathways)) {
    pathways <- if (is.null(pmn_gene_col)) catfish_load_pathways(species = species)
                else catfish_load_pathways(species = species, gene_col = pmn_gene_col)
  }
  pw <- .catfish_pathways_to_list(pathways)
  paths <- list()
  for (nm in names(pw$p_list)) {
    gs <- unique(tolower(trimws(as.character(pw$p_list[[nm]]))))
    gs <- gs[nzchar(gs)]
    ii <- idx_of[gs]; ii <- ii[!is.na(ii)]
    if (length(ii) >= min_genes && length(ii) <= max_genes) paths[[nm]] <- as.integer(ii)
  }
  NP <- length(paths)
  if (!NP) stop("No pathway has between min_genes and max_genes genes present.", call. = FALSE)
  if (isTRUE(verbose))
    message(sprintf("[WY] pathways kept (%d-%s genes present): %d", min_genes,
                    format(max_genes), NP))

  ## ---- engine
  cp  <- .wy_read_cor_pairs(cor_pairs, cor_file)
  eng <- .wy_engine(gid, gene_lp, chr, paths, cp, B, use_tfisher, tau_grid,
                    seed, verbose)

  ## ---- competitive test (circular genomic permutation)
  comp <- NULL
  if (isTRUE(competitive)) {
    if (is.null(pos_col) && isTRUE(verbose))
      message("[WY competitive] pos_col not given: genes taken in the order supplied within each chromosome")
    if (!is.null(pos_col) && any(!is.finite(pos)))
      warning(sprintf("%d genes have no position in '%s' and are placed at the end of their chromosome.",
                      sum(!is.finite(pos)), pos_col), call. = FALSE)
    chr_num <- suppressWarnings(as.numeric(chr))
    ord <- order(is.na(chr_num), chr_num, chr, pos, seq_along(gid), na.last = TRUE)
    comp <- .wy_competitive(gene_lp, ord, paths, eng$obs_log, B_comp, use_tfisher,
                            tau_grid, seed, verbose)
  }

  res <- data.frame(
    pathway_id      = names(paths),
    pathway_name    = unname(as.character(pw$p_names[names(paths)])),
    z_std           = eng$z_std,
    p_fwer_std      = eng$p_fwer_wy_std,
    n_genes         = vapply(paths, length, integer(1)),
    log10_p_omnibus = eng$obs_log / log(10),
    p_omnibus       = exp(eng$obs_log),
    p_fwer_wy       = eng$p_fwer_wy,
    p_bonferroni    = pmin(exp(eng$obs_log) * NP, 1),
    p_calibrated    = eng$p_cal,
    gc_shift        = gcr$shift,
    stringsAsFactors = FALSE
  )
  if (!is.null(comp)) {
    res$z_std_comp        <- comp$z_std
    res$p_fwer_comp       <- comp$p_fwer_std
    res$p_calibrated_comp <- comp$p_cal
  }
  res <- res[order(res$z_std), , drop = FALSE]
  rownames(res) <- NULL
  attr(res, "wy_threshold_05") <- exp(stats::quantile(eng$minp_null, 0.05, names = FALSE))
  attr(res, "B") <- B
  attr(res, "tau_grid") <- tau_grid
  attr(res, "gc_shift") <- gcr$shift
  attr(res, "median_gene_z") <- gcr$median_z
  if (!is.null(comp)) attr(res, "B_comp") <- comp$B

  if (isTRUE(output)) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    f <- file.path(out_dir, "catfish_wy_pathways.csv")
    utils::write.csv(res, f, row.names = FALSE)
    attr(res, "file") <- f
  }
  res
}
