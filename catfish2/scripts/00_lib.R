## =============================================================================
## 00_lib.R -- shared statistics for the CATFISH-WY paper
##
## Everything is computed as log(p). On pathways of 100-200 genes with dense
## signal the Fisher and Stouffer p-values fall below ~1e-308 and underflow to
## zero in linear space, which ties large blocks of pathways at the floor and
## destroys any ranking. R's *_log.p arguments give exact logs instead.
## =============================================================================

TAU_GRID_DEFAULT <- c(0.1, 0.05, 0.01, 0.005, 0.001)

.logsumexp <- function(lx) {
  m <- max(lx)
  if (!is.finite(m)) return(-Inf)
  m + log(sum(exp(lx - m)))
}

## ACAT on log p-values, returning log(p_ACAT).
##   T = mean_i tan((0.5 - p_i) * pi);  tan(.) -> 1/(p*pi) as p -> 0
##   large T  =>  p_ACAT -> 1/(T*pi)
acat_log <- function(lp) {
  lp <- pmin(lp, -1e-300)
  m  <- length(lp)
  small <- lp < log(1e-15)
  S_big <- if (any(!small)) sum(tan((0.5 - exp(lp[!small])) * pi)) else 0

  ## No extreme p-value: T = S_big/m is finite and may be of EITHER sign. The
  ## log-space branch below cannot represent a negative T, and clamping it (as
  ## an earlier version did) pins every such pathway at p = 0.5 exactly, which
  ## destroys the ordering over the whole upper half of the distribution.
  if (!any(small)) {
    pv <- 0.5 - atan(S_big / m) / pi
    return(log(max(pv, 1e-323)))
  }

  ## At least one extreme p: S is dominated by the 1/(p*pi) terms, so work in
  ## logs to avoid overflow.
  logS_small <- .logsumexp(-lp[small] - log(pi))
  logS <- if (logS_small > log(abs(S_big) + 1e-300) + 50) {
    logS_small
  } else {
    log(max(exp(logS_small) + S_big, 1e-300))
  }
  logT <- logS - log(m)
  ## switch to the 1/(T*pi) asymptote well before 1e15: beyond T ~ 1e5 the
  ## 0.5 - atan(T)/pi form loses more to cancellation than the asymptote loses
  ## to truncation.
  if (logT > log(1e6)) return(-logT - log(pi))
  pv <- 0.5 - atan(exp(logT)) / pi
  if (pv <= 0) -745 else log(pv)
}

fisher_log <- function(lp) {
  stats::pchisq(-2 * sum(lp), df = 2 * length(lp), lower.tail = FALSE, log.p = TRUE)
}

stouffer_log <- function(lp) {
  z <- stats::qnorm(lp, lower.tail = FALSE, log.p = TRUE)
  z[!is.finite(z)] <- 0
  stats::pnorm(sum(z) / sqrt(length(z)), lower.tail = FALSE, log.p = TRUE)
}

## p_minP = 1 - (1 - q)^m with q = exp(lmin).
## For small m*q this equals m*q to full precision, and the 1-(1-q)^m form
## CANNOT be used there: once q < ~1e-16/m, exp(m*log1p(-q)) rounds to exactly
## 1, log1p(-1) is -Inf, and the component is silently dropped from the omnibus
## by the is.finite() filter -- i.e. minP disappears exactly when a pathway has
## a very strong gene, which is the case it exists to detect.
minp_log <- function(lp) {
  m <- length(lp); lmin <- min(lp)
  approx <- lmin + log(m)
  if (approx < log(1e-8)) return(approx)
  v <- m * log1p(-exp(lmin))
  if (!is.finite(v)) return(approx)
  r <- log1p(-exp(v))
  if (!is.finite(r)) approx else r
}

tfisher_log <- function(lp, tau_grid = TAU_GRID_DEFAULT) {
  if (!requireNamespace("TFisher", quietly = TRUE)) return(NA_real_)
  p <- exp(pmax(lp, -700)); m <- length(p); best <- 0
  for (tau in tau_grid) {
    st <- TFisher::stat.soft(p = p, tau1 = tau)
    pv <- TFisher::p.soft(q = st, n = m, tau1 = tau, M = NULL)
    if (is.finite(pv) && pv > 0 && log(pv) < best) best <- log(pv)
  }
  if (best == 0) NA_real_ else best
}

## Component vector (log scale). `drop` lets the breakage tests remove or
## corrupt individual components.
components_log <- function(lp, use_tfisher = TRUE, tau_grid = TAU_GRID_DEFAULT,
                           drop = character(0)) {
  out <- c(ACAT = acat_log(lp), Fisher = fisher_log(lp),
           Stouffer = stouffer_log(lp), minP = minp_log(lp))
  if (isTRUE(use_tfisher)) out <- c(out, TFisher = tfisher_log(lp, tau_grid))
  out[setdiff(names(out), drop)]
}

## Omnibus = ACAT across component log-p values.
omnibus_log <- function(lp, use_tfisher = TRUE, tau_grid = TAU_GRID_DEFAULT,
                        drop = character(0), gamma = NULL) {
  comps <- components_log(lp, use_tfisher, tau_grid, drop)
  ## gamma: multiply selected component log-p by gamma < 1 to make them
  ## anti-conservative (p^gamma). Used by the breakage stress test.
  if (!is.null(gamma)) {
    for (nm in names(gamma)) if (nm %in% names(comps)) comps[nm] <- comps[nm] * gamma[[nm]]
  }
  comps <- comps[is.finite(comps)]
  if (!length(comps)) return(0)
  acat_log(comps)
}

## -----------------------------------------------------------------------------
## Null generators
## -----------------------------------------------------------------------------

## Exchangeable correlation matrix: 1 on the diagonal, rho off-diagonal.
R_exch <- function(m, rho) {
  R <- matrix(rho, m, m); diag(R) <- 1; R
}

## Cholesky with an eigen repair fallback (MAGMA's banded correlation estimates
## are frequently not positive definite).
safe_chol <- function(R) {
  U <- tryCatch(chol(R), error = function(e) NULL)
  if (!is.null(U)) return(U)
  ev <- eigen(R, symmetric = TRUE)
  ev$values[ev$values < 1e-6] <- 1e-6
  R2 <- ev$vectors %*% diag(ev$values) %*% t(ev$vectors)
  d <- sqrt(diag(R2)); R2 <- R2 / outer(d, d)
  chol(R2)
}

## One draw of correlated log p-values (one-sided, uniform marginals).
draw_lp <- function(U) {
  z <- as.vector(t(U) %*% stats::rnorm(ncol(U)))
  stats::pnorm(z, lower.tail = FALSE, log.p = TRUE)
}

## -----------------------------------------------------------------------------
## Calibration
## -----------------------------------------------------------------------------

## Per-pathway MVN calibration: empirical tail probability of the observed
## statistic against B draws from this pathway's own null.
mvn_calibrate <- function(obs_log, null_log) {
  ## Drop non-finite null draws before counting. A component may legitimately
  ## be undefined on some draws -- TFisher returns NA when no gene falls below
  ## the smallest tau -- and without this filter one such draw turns the whole
  ## calibrated p-value into NA.
  null_log <- null_log[is.finite(null_log)]
  if (!length(null_log)) return(NA_real_)
  (1 + sum(null_log <= obs_log)) / (length(null_log) + 1)
}

## Genomic control lambda from a vector of p-values (median chi-square ratio).
lambda_gc <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  if (!length(p)) return(NA_real_)
  stats::median(stats::qchisq(p, df = 1, lower.tail = FALSE)) / stats::qchisq(0.5, df = 1, lower.tail = FALSE)
}

type1 <- function(p, alpha) mean(p < alpha, na.rm = TRUE)

## =============================================================================
## Shared gene universe with overlapping pathways
##
## Scripts that draw each pathway from its own R_S contain NO between-pathway
## dependence, so a Westfall-Young threshold computed on them is just a
## Bonferroni threshold with extra steps. These helpers build one gene universe
## with block LD and define pathways as OVERLAPPING subsets of it, so a single
## draw perturbs every pathway containing a given gene at once -- the structure
## a real pathway database has.
## =============================================================================

## G genes in independent blocks of `block` genes; exchangeable rho within a
## block. Blocks stand in for chromosomes / LD blocks.
make_universe <- function(G = 3000L, block = 50L, rho = 0.3) {
  stopifnot(G %% block == 0L)
  list(U = safe_chol(R_exch(block, rho)), nblock = G %/% block,
       G = G, block = block, rho = rho)
}

## One genome-wide draw -> vector of G log p-values.
draw_universe_lp <- function(uni) {
  z <- as.vector(vapply(seq_len(uni$nblock),
         function(i) as.vector(t(uni$U) %*% stats::rnorm(uni$block)),
         numeric(uni$block)))
  stats::pnorm(z, lower.tail = FALSE, log.p = TRUE)
}

## K pathways of m genes each, drawn from the first `pool` genes of the
## universe. Smaller pool => each gene sits in more pathways => more overlap.
## Expected pathways per gene is K * m / pool.
make_pathways <- function(K, m, pool, G) {
  stopifnot(pool <= G, m <= pool)
  p <- lapply(seq_len(K), function(i) sample.int(pool, m))
  names(p) <- sprintf("PW%04d", seq_len(K))
  p
}

mean_memberships <- function(paths, G) {
  tab <- tabulate(unlist(paths), nbins = G)
  mean(tab[tab > 0])
}

## B x K x C array of component log p-values: every pathway scored on every
## draw. This is the expensive step, and it is deliberately separated from the
## omnibus so that breakage and leave-one-out variants can be evaluated from
## ONE set of component evaluations instead of recomputing them per variant.
components_array <- function(uni, paths, B, use_tfisher = TRUE,
                             tau_grid = TAU_GRID_DEFAULT, verbose = TRUE) {
  K  <- length(paths)
  cn <- names(components_log(rep(-1, 5), use_tfisher, tau_grid))
  out <- array(NA_real_, c(B, K, length(cn)),
               dimnames = list(NULL, names(paths), cn))
  at <- max(1L, B %/% 10L)
  for (b in seq_len(B)) {
    lp <- draw_universe_lp(uni)
    for (k in seq_len(K)) out[b, k, ] <- components_log(lp[paths[[k]]],
                                                        use_tfisher, tau_grid)
    if (verbose && b %% at == 0L) cat(sprintf("      draws %d/%d\n", b, B))
  }
  out
}

## Omnibus from a precomputed component array. `gamma` corrupts components
## (log p * gamma, i.e. p^gamma); `drop` removes them entirely.
omnibus_from_components <- function(arr, gamma = NULL, drop = character(0)) {
  keep <- setdiff(dimnames(arr)[[3]], drop)
  sub  <- arr[, , keep, drop = FALSE]
  if (!is.null(gamma)) {
    for (nm in names(gamma)) {
      if (nm %in% keep) sub[, , nm] <- sub[, , nm] * gamma[[nm]]
    }
  }
  B <- dim(sub)[1]; K <- dim(sub)[2]
  out <- matrix(NA_real_, B, K)
  for (b in seq_len(B)) {
    for (k in seq_len(K)) {
      v <- sub[b, k, ]; v <- v[is.finite(v)]
      out[b, k] <- if (length(v)) acat_log(v) else 0
    }
  }
  out
}

## Column-wise (marginal) calibration: each column of `obs` against the same
## column of `cal`. This is the per-pathway MVN p-value.
calibrate_cols <- function(obs, cal) {
  res <- obs
  for (k in seq_len(ncol(cal))) {
    res[, k] <- vapply(obs[, k], function(o) mvn_calibrate(o, cal[, k]), numeric(1))
  }
  res
}

## Row-wise (family-wise) threshold: the alpha-quantile of the family minimum.
## Computed on the RAW omnibus scale, which is comparable across pathways here
## because every pathway in a given setting has the same size.
wy_threshold <- function(cal_om, alpha) {
  stats::quantile(apply(cal_om, 1, min, na.rm = TRUE), alpha, names = FALSE)
}

## Family-wise error rates on held-out draws.
fwer_rates <- function(ev_om, thr_wy, K, alpha) {
  mn <- apply(ev_om, 1, min, na.rm = TRUE)
  c(wy = mean(mn <= thr_wy), bonferroni = mean(exp(mn) < alpha / K))
}
