## =============================================================================
## 00_lib.R -- simulation infrastructure for the CATFISH-WY paper.
##
## The component tests and the omnibus are NOT implemented here. They are taken
## from the CATFISH package, so these simulations exercise the same code the
## real-data analyses run. What remains in this file is the simulation scaffold
## only: the correlated gene universe, the pathway sampler, the per-pathway
## (marginal) calibration, the Westfall-Young family-wise threshold, and the
## error-rate summaries.
##
## Statistics are carried as log(p) inside the scaffold because the family
## minimum and the calibration comparisons are done on that scale; the package
## returns linear p-values floored at min_p and they are logged on arrival.
## =============================================================================

## -----------------------------------------------------------------------------
## Component tests -- ALL of these come from the CATFISH package.
##
## Nothing in this file reimplements a component test any more. Every component
## p-value and the ACAT-O omnibus are produced by
## CATFISH::catfish_omni2_pathways(), which is the same code path the real-data
## analyses use, so the simulations in this paper validate the shipped package
## rather than a parallel implementation of it.
##
## Columns returned by the package and the names used here:
##   acat_p               -> ACAT
##   fisher_p             -> Fisher
##   stouffer_p_analytic  -> Stouffer
##   minp_p_analytic      -> minP
##   tfisher_p_analytic   -> TFisher
##   omni_p_analytic      -> Omnibus
##
## The package floors p-values at min_p (1e-50 in omni2). Every threshold used
## in this paper -- alpha = 0.05, alpha/K = 3.3e-4, and the Westfall-Young
## family-wise threshold near 1e-3 -- is many orders of magnitude above that
## floor, so the cap never binds in these simulations.
## -----------------------------------------------------------------------------

## Truncation thresholds for adaptive soft TFisher. This is the grid of
## Zhang, Tong, Landers & Wu (2020), TAU1 = c(0.01, 0.05, 0.5, 1), with the
## tau = 1 corner dropped because tau = 1 IS Fisher and Fisher is already a
## separate component. Spanning small, intermediate and large tau is their
## stated recommendation; a grid confined to small tau leaves m*tau << 1 at
## realistic pathway sizes and the statistic is empty on most pathways.
TAU_GRID_DEFAULT <- c(0.01, 0.05, 0.5)
## Lower floor for component and omnibus p-values. Raised from 1e-50 because
## ACAT has a Cauchy tail: with every gene in a pathway correlated at rho=0.6
## the analytic omnibus p underflows past 1e-50 on roughly a quarter of null
## draws, those draws tie at the floor, and the Westfall-Young threshold then
## lands on the tie. 1e-300 is near the limit of double precision and is not
## reached in any configuration tested.
MIN_P_DEFAULT    <- 1e-300

## TFisher reports Wo = 1 - p.soft(...). In double precision 1 - x is exactly 0
## once x > 1 - 1e-16, so a strongly enriched pathway returns Wo = 0 and would
## be clamped to MIN_P_DEFAULT -- reporting 1e-300 for a statistic that cannot
## resolve past 1e-16. Cap TFisher at its own precision limit instead. Fisher
## and Stouffer are computed in closed form and genuinely reach 1e-100, so they
## keep MIN_P_DEFAULT.
TF_MIN_P         <- 1e-16

## Component p-value column returned by each CATFISH wrapper.
COMPONENT_FUN <- list(
  ACAT     = list(fn = "catfish_acat_pathways",                   col = "acat_p"),
  Fisher   = list(fn = "catfish_fisher_pathways",                 col = "fisher_p"),
  Stouffer = list(fn = "catfish_stoufferZ_pathways",              col = "stouffer_p"),
  minP     = list(fn = "catfish_minp_pathways",                   col = "minp_p"),
  ## TFisher contributes Wo = min over the tau grid of the single-tau right-tail
  ## p (the package column tfisher_p_min), NOT the adaptive p-value
  ## tfisher_p_omni. TFisher::p.soft.omni, which converts Wo into a p-value, is
  ## miscalibrated in both directions: lambda ~ 0.12 in the bulk, and it
  ## underflows in the tail so that strongly enriched pathways clamp to min_p
  ## and tie. Wo is a perfectly good monotone statistic (smaller = more
  ## extreme) and the resampled null bank calibrates it, which is what every
  ## other component here gets anyway.
  ##
  ## CONSEQUENCE: Wo is not a p-value -- under the null it is the minimum of
  ## five right-tail p-values and so is stochastically smaller than uniform.
  ## Any column that treats a component value AS a p-value is therefore not
  ## interpretable for TFisher: t1_analytic, fwer_bonferroni and
  ## lambda_analytic. Judge TFisher on the calibrated and family-wise columns.
  TFisher  = list(fn = "catfish_soft_tfisher_adaptive_pathways",  col = "tfisher_p_min")
)
STAT_NAMES <- c(names(COMPONENT_FUN), "Omnibus")

.require_catfish <- function() {
  if (!requireNamespace("CATFISH", quietly = TRUE))
    stop("the CATFISH package is required; install/load it first", call. = FALSE)
}

## Gene identifiers for a simulated universe of G genes.
universe_gene_ids <- function(G) sprintf("g%05d", seq_len(G))

## Convert integer-index pathways (what make_pathways returns) to the character
## gene-id lists the package expects.
paths_to_chr <- function(paths, G) {
  ids <- universe_gene_ids(G)
  lapply(paths, function(ix) ids[ix])
}

## One draw -> a K x 6 matrix of log p (5 components + Omnibus), pathway order
## preserved.
##
## The five components come from the CATFISH wrappers and the omnibus from
## CATFISH::catfish_acat_pathways applied to those five. The wrappers are used
## rather than catfish_omni2_pathways because omni2 additionally runs IHW,
## qvalue, BKY and BH on every call; that multiple-testing layer is what the
## real-data report needs and is dead weight inside a resampling loop. The
## component math is identical either way.
##
## Only genes that appear in at least one pathway are passed in: a gene in no
## pathway cannot affect any component, and dropping them keeps the per-call
## gene table small.
components_from_draw <- function(lp, z, paths_chr, G,
                                 tau_grid = TAU_GRID_DEFAULT,
                                 min_p = MIN_P_DEFAULT) {
  .require_catfish()
  ids  <- universe_gene_ids(G)
  used <- sort(unique(unlist(paths_chr, use.names = FALSE)))
  ix   <- match(used, ids)
  gr <- data.frame(GENE  = used,
                   P     = pmin(pmax(exp(lp[ix]), min_p), 1 - 1e-15),
                   ZSTAT = z[ix],
                   stringsAsFactors = FALSE)
  K  <- length(paths_chr)
  pm <- matrix(NA_real_, K, length(COMPONENT_FUN),
               dimnames = list(names(paths_chr), names(COMPONENT_FUN)))
  for (nm in names(COMPONENT_FUN)) {
    spec <- COMPONENT_FUN[[nm]]
    fn   <- getExportedValue("CATFISH", spec$fn)
    args <- list(gene_results = gr, pathways = paths_chr,
                 gene_col = "GENE", min_p = min_p, output = FALSE)
    if (nm == "Stouffer")      args$z_col    <- "ZSTAT"
    else                       args$p_col    <- "P"
    if (nm == "TFisher") { args$p_col <- "P"; args$tau_grid <- tau_grid }
    out <- do.call(fn, args)
    v <- as.numeric(out[[spec$col]])[match(names(paths_chr), out$pathway_id)]
    ## Soft TFisher is UNDEFINED, not non-significant, when no gene in the
    ## pathway beats any tau: stat.soft is 0 at every tau and the right-tail p
    ## returns exactly 1. Handing that 1 to ACAT contributes tan(-0.3*pi) ~
    ## -1.38, so a test that could not be computed votes AGAINST the pathway.
    ## NA instead, and the pathway is combined from whichever components ARE
    ## defined. This is decided per pathway, not per size: a large pathway in
    ## which every gene is weak is empty too.
    ##
    ## This conditions TFisher on having at least one gene below tau, so where
    ## it does contribute its p-values are stochastically small. The null bank
    ## is built under the same rule, so the calibration and the WY threshold
    ## absorb that bias -- an ANALYTIC omnibus p under this rule would be
    ## anti-conservative.
    if (nm == "TFisher") {
      v[is.finite(v) & v >= 1 - 1e-12] <- NA_real_          # empty statistic
      v[is.finite(v)] <- pmax(v[is.finite(v)], TF_MIN_P)    # precision limit
    }
    pm[, nm] <- v
  }
  cbind(log(pmax(pm, min_p)),
        Omnibus = acat_combine_log(log(pmax(pm, min_p)), min_p))
}

## Many INDEPENDENT pathways scored in one go. `zmat` is n x m: one row per
## pathway, each row an independent draw of m gene-level z values. Each pathway
## is given its own private genes, so the pathways stay independent -- this is
## purely a batching device for scripts 01 and 04, which otherwise call the
## package once per pathway and pay its per-call setup tens of thousands of
## times. Processed in chunks to keep the gene table a sane size.
components_from_zmat <- function(zmat, tau_grid = TAU_GRID_DEFAULT,
                                 min_p = MIN_P_DEFAULT, chunk = 2000L) {
  if (is.null(dim(zmat))) zmat <- matrix(zmat, nrow = 1L)
  n <- nrow(zmat); m <- ncol(zmat)
  out <- matrix(NA_real_, n, length(STAT_NAMES),
                dimnames = list(NULL, STAT_NAMES))
  for (lo in seq(1L, n, by = chunk)) {
    hi  <- min(lo + chunk - 1L, n)
    nn  <- hi - lo + 1L
    G   <- nn * m
    ids <- universe_gene_ids(G)
    pw  <- split(ids, rep(seq_len(nn), each = m))
    names(pw) <- sprintf("PW%06d", seq_len(nn))
    z   <- as.vector(t(zmat[lo:hi, , drop = FALSE]))
    lp  <- stats::pnorm(z, lower.tail = FALSE, log.p = TRUE)
    out[lo:hi, ] <- components_from_draw(lp, z, pw, G, tau_grid, min_p)
  }
  out
}

## n independent pathway draws of m genes from a Cholesky factor, as a matrix.
draw_zmat <- function(U, n, mu = 0) {
  m <- ncol(U)
  t(matrix(as.vector(t(U) %*% matrix(stats::rnorm(m * n), m, n)), m, n)) + mu
}

## Single pathway of m genes scored through the package. Used by the scripts
## that test one pathway at a time (01, 04) rather than a whole family.
## `z` is the gene-level z vector; one-sided p-values are taken from it.
components_single <- function(z, tau_grid = TAU_GRID_DEFAULT,
                              min_p = MIN_P_DEFAULT) {
  m  <- length(z)
  lp <- stats::pnorm(z, lower.tail = FALSE, log.p = TRUE)
  pw <- list(PW0001 = universe_gene_ids(m))
  components_from_draw(lp, z, pw, m, tau_grid, min_p)[1, ]
}

## Omnibus log p for a single pathway, from the package.
omnibus_log_z <- function(z, ...) unname(components_single(z, ...)[["Omnibus"]])

## One correlated draw of m gene-level z values from a Cholesky factor.
draw_z <- function(U, mu = 0) as.vector(t(U) %*% stats::rnorm(ncol(U))) + mu

acat_combine_log <- function(log_comp_mat, min_p = MIN_P_DEFAULT) {
  .require_catfish()
  if (is.null(dim(log_comp_mat)))
    log_comp_mat <- matrix(log_comp_mat, nrow = 1L,
                           dimnames = list(NULL, names(log_comp_mat)))
  K  <- nrow(log_comp_mat); cn <- colnames(log_comp_mat)
  gid <- character(0); pv <- numeric(0)
  pw  <- vector("list", K); names(pw) <- sprintf("PW%06d", seq_len(K))
  for (k in seq_len(K)) {
    ok  <- which(is.finite(log_comp_mat[k, ]))   ## undefined components drop out
    ids <- sprintf("%s_%06d", cn[ok], k)
    gid <- c(gid, ids); pv <- c(pv, exp(log_comp_mat[k, ok]))
    pw[[k]] <- ids
  }
  res  <- rep(NA_real_, K)
  keep <- lengths(pw) > 0L
  if (!any(keep)) return(res)
  gr <- data.frame(GENE = gid, P = pmin(pmax(pv, min_p), 1 - 1e-15),
                   stringsAsFactors = FALSE)
  out <- CATFISH::catfish_acat_pathways(gene_results = gr, pathways = pw[keep],
                                        gene_col = "GENE", p_col = "P",
                                        min_p = min_p, output = FALSE)
  out <- out[match(names(pw)[keep], out$pathway_id), , drop = FALSE]
  res[keep] <- log(pmax(as.numeric(out$acat_p), min_p))
  res
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

## One genome-wide draw on the z scale, so a mean shift can be added before the
## one-sided p-values are taken.
draw_universe_z <- function(uni) {
  as.vector(vapply(seq_len(uni$nblock),
    function(i) as.vector(t(uni$U) %*% stats::rnorm(uni$block)),
    numeric(uni$block)))
}

## One genome-wide draw -> vector of G log p-values.
draw_universe_lp <- function(uni) {
  stats::pnorm(draw_universe_z(uni), lower.tail = FALSE, log.p = TRUE)
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
  paths_chr <- paths_to_chr(paths, uni$G)
  cn  <- STAT_NAMES
  out <- array(NA_real_, c(B, length(paths), length(cn)),
               dimnames = list(NULL, names(paths), cn))
  at  <- max(1L, B %/% 10L)
  for (b in seq_len(B)) {
    z  <- draw_universe_z(uni)
    lp <- stats::pnorm(z, lower.tail = FALSE, log.p = TRUE)
    out[b, , ] <- components_from_draw(lp, z, paths_chr, uni$G, tau_grid)
    if (verbose && b %% at == 0L) cat(sprintf("      draws %d/%d\n", b, B))
  }
  if (!isTRUE(use_tfisher)) out <- out[, , setdiff(cn, "TFisher"), drop = FALSE]
  out
}

## Omnibus from a precomputed component array. `gamma` corrupts components
## (log p * gamma, i.e. p^gamma); `drop` removes them entirely.
omnibus_from_components <- function(arr, gamma = NULL, drop = character(0)) {
  comps <- setdiff(dimnames(arr)[[3]], c("Omnibus", drop))
  ## untouched: the package already returned its own ACAT-O omnibus
  if (is.null(gamma) && !length(drop) && "Omnibus" %in% dimnames(arr)[[3]])
    return(arr[, , "Omnibus"])
  sub <- arr[, , comps, drop = FALSE]
  if (!is.null(gamma)) {
    for (nm in names(gamma)) if (nm %in% comps) sub[, , nm] <- sub[, , nm] * gamma[[nm]]
  }
  B <- dim(sub)[1]
  out <- matrix(NA_real_, B, dim(sub)[2])
  for (b in seq_len(B)) out[b, ] <- acat_combine_log(sub[b, , , drop = TRUE])
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
