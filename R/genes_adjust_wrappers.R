#' Adjust MAGMA gene p-values using gene length and NSNPS
#'
#' Regresses out the effect of gene length and number of SNPs (NSNPS) from
#' gene-level Z-statistics. This can reduce bias from longer genes having
#' more SNPs and therefore more significant p-values.
#'
#' @param gene_results data.frame. Gene-level results from MAGMA (typically
#'   from reading a `.genes.out` file). Must contain gene ID, p-value, and
#'   NSNPS columns.
#' @param gene_lengths data.frame. Gene lengths, typically from
#'   \code{\link{get_gene_lengths}}. Must contain gene ID and length columns.
#' @param gene_col character. Column name in \code{gene_results} for gene IDs.
#'   Default "GENE".
#' @param nsnp_col character. Column name in \code{gene_results} for the number
#'   of SNPs per gene. Default "NSNPS".
#' @param p_col character. Column name in \code{gene_results} for raw p-values.
#'   Default "P".
#' @param z_col character or NULL. Column name for Z-scores in \code{gene_results}.
#'   If NULL or not present, Z-scores are reconstructed from p-values (magnitude
#'   only). Default "ZSTAT".
#' @param len_gene_col character. Column name in \code{gene_lengths} for gene IDs.
#'   Default "gene_id".
#' @param len_col character. Column name in \code{gene_lengths} for gene length.
#'   Default "length".
#' @param log1p_covars logical. If TRUE (default), apply log(1+x) transformation
#'   to NSNPS and gene length covariates before regression.
#'
#' @return data.frame with columns:
#' \describe{
#'   \item{gene_id}{Gene identifier}
#'   \item{z_adj}{Adjusted Z-score (residuals from regression)}
#'   \item{p_adj}{Adjusted p-value derived from \code{z_adj}}
#'   \item{fit}{List column containing the lm fit object for each gene}
#' }
#' An attribute "lm_fit" is also attached with the linear model object.
#'
#' @examples
#' \dontrun{
#' # Load gene results from MAGMA
#' gene_results <- read.delim("magma_output.genes.out")
#'
#' # Get gene lengths from GFF3
#' gene_lengths <- get_gene_lengths("reference.gff3")
#'
#' # Adjust gene p-values
#' adjusted <- magcat_adjust_gene_p(
#'   gene_results = gene_results,
#'   gene_lengths = gene_lengths,
#'   gene_col = "GENE",
#'   nsnp_col = "NSNPS",
#'   p_col = "P"
#' )
#'
#' # Use adjusted p-values for pathway analysis
#' gene_results$P_adj <- adjusted$p_adj[match(gene_results$GENE, adjusted$gene_id)]
#' }
#'
#' @seealso \code{\link{magma_gene}}, \code{\link{get_gene_lengths}}
#' @export
magcat_adjust_gene_p <- function(gene_results,
                                 gene_lengths,
                                 gene_col     = "GENE",
                                 nsnp_col     = "NSNPS",
                                 p_col        = "P",
                                 z_col        = "ZSTAT",
                                 len_gene_col = "gene_id",
                                 len_col      = "length",
                                 log1p_covars = TRUE) {

  if (!is.data.frame(gene_results)) stop("magcat_adjust_gene_p(): gene_results must be a data.frame.", call. = FALSE)
  if (!is.data.frame(gene_lengths)) stop("magcat_adjust_gene_p(): gene_lengths must be a data.frame.", call. = FALSE)

  for (nm in c(gene_col, nsnp_col, p_col)) {
    if (!(nm %in% names(gene_results))) {
      stop("magcat_adjust_gene_p(): missing column in gene_results: ", nm, call. = FALSE)
    }
  }
  for (nm in c(len_gene_col, len_col)) {
    if (!(nm %in% names(gene_lengths))) {
      stop("magcat_adjust_gene_p(): missing column in gene_lengths: ", nm, call. = FALSE)
    }
  }

  gr <- gene_results
  gl <- gene_lengths

  gr[[gene_col]] <- as.character(gr[[gene_col]])
  gl[[len_gene_col]] <- as.character(gl[[len_gene_col]])

  # fix "gene:Zm..." -> "Zm..."
  gl[[len_gene_col]] <- sub("^gene:", "", gl[[len_gene_col]])

  # merge length onto results
  m <- merge(
    gr,
    gl[, c(len_gene_col, len_col), drop = FALSE],
    by.x = gene_col,
    by.y = len_gene_col,
    all.x = TRUE,
    sort = FALSE
  )

  # get z_raw (prefer supplied ZSTAT; else reconstruct from P)
  have_z <- (!is.null(z_col)) && (z_col %in% names(m))
  if (have_z) {
    z_raw <- suppressWarnings(as.numeric(m[[z_col]]))
  } else {
    p <- suppressWarnings(as.numeric(m[[p_col]]))
    p[p <= 0] <- 1e-300
    p[p >= 1] <- 1 - 1e-16
    z_raw <- stats::qnorm(1 - p/2)  # magnitude only
  }

  # residualize magnitude
  y <- abs(z_raw)

  nsnps <- suppressWarnings(as.numeric(m[[nsnp_col]]))
  glen  <- suppressWarnings(as.numeric(m[[len_col]]))

  if (log1p_covars) {
    nsnps <- log1p(nsnps)
    glen  <- log1p(glen)
  }

  ok <- is.finite(y) & is.finite(nsnps) & is.finite(glen)

  if (!any(ok)) {
    stop("magcat_adjust_gene_p(): no complete cases after merging covariates.", call. = FALSE)
  }

  df <- data.frame(y = y[ok], nsnps = nsnps[ok], len = glen[ok])

  fit <- stats::lm(y ~ nsnps + len, data = df)

  z_adj <- rep(NA_real_, nrow(m))
  z_adj[ok] <- stats::residuals(fit)

  p_adj <- rep(NA_real_, nrow(m))
  p_adj[ok] <- 2 * stats::pnorm(-abs(z_adj[ok]))

  out <- data.frame(
    gene_id = m[[gene_col]],
    z_adj   = z_adj,
    p_adj   = p_adj,
    stringsAsFactors = FALSE
  )

  # ✅ save model both ways
  attr(out, "lm_fit") <- fit
  out$fit <- list(fit)

  out
}
