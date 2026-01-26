#' Pathway-level (weighted) Stouffer Z test on gene-level Z statistics
#'
#' Computes a (possibly weighted) Stouffer Z within each pathway:
#' \deqn{Z_S = \frac{\sum_{g \in S} w_g Z_g}{\sqrt{\sum_{g \in S} w_g^2}}}
#'
#' and returns a p-value using the requested \code{alternative}:
#' \itemize{
#'   \item \code{"two.sided"}: \eqn{p = 2 \Phi(-|Z_S|)}
#'   \item \code{"greater"}:  \eqn{p = \Pr(Z \ge Z_S)}
#'   \item \code{"less"}:     \eqn{p = \Pr(Z \le Z_S)}
#' }
#'
#' Gene IDs are matched case-insensitively: both \code{gene_results[[gene_col]]} and
#' pathway gene IDs are normalized to lower-case before matching. Genes not present
#' in \code{gene_results} are dropped from that pathway.
#'
#' NOTE: The analytic p-value assumes (approximately) independent \eqn{Z_g} within a pathway.
#' If genes are correlated (e.g., LD-induced correlation in MAGMA gene statistics), then this
#' analytic p-value may be anti-conservative. In MAGCAT/CATFISH, correlation-aware calibration
#' should be handled in your MVN/Omni layer (not inside this wrapper).
#'
#' @param gene_results data.frame with at least a gene ID column and a Z-statistic column.
#'   For MAGMA `.genes.out`, you may have \code{GENE} and \code{ZSTAT} (or similar).
#' @param pathways Either:
#' \itemize{
#'   \item A named list: each element is a character vector of gene IDs in that pathway, OR
#'   \item A data.frame with columns \code{pathway_id}, \code{gene_id} (and optional \code{pathway_name}).
#' }
#' @param species Optional; one of "maize", "sorghum", "arabidopsis", "plant".
#'   If provided, built-in PMN pathways are loaded via \code{magcat_load_pathways()}.
#'   Provide either \code{pathways} OR \code{species}, but not both.
#' @param pmn_gene_col Optional; passed to \code{magcat_load_pathways(gene_col=...)}
#'   when \code{species} is used. If NULL, the PMN loader prefers "Gene-name"
#'   then "Gene-id".
#' @param gene_col Column name in \code{gene_results} containing gene IDs (default "GENE").
#' @param z_col Column name in \code{gene_results} containing gene-level Z statistics (default "ZSTAT").
#' @param weight_col Optional column name in \code{gene_results} for per-gene weights
#'   (default NULL = equal weights). Examples: \code{sqrt(NSNPS)}, \code{1/SE}, etc.
#' @param min_abs_w Non-finite / NA / <=0 weights are replaced by this small positive value
#'   (default 1e-8).
#' @param alternative One of \code{"greater"}, \code{"two.sided"}, \code{"less"}.
#'   \code{"greater"} tests for positive enrichment (default "greater").
#' @param output If TRUE, write a CSV of results to \code{out_dir}.
#' @param out_dir Directory to write CSV when \code{output = TRUE} (default "magcat_stouffer_z").
#'
#' @return A data.frame with columns:
#' \describe{
#'   \item{pathway_id}{Pathway identifier}
#'   \item{pathway_name}{Pathway name (if available; otherwise equals \code{pathway_id})}
#'   \item{n_genes}{Number of genes used in the pathway}
#'   \item{gene_names}{Semicolon-separated gene IDs actually used (from \code{gene_results})}
#'   \item{gene_zvals}{Semicolon-separated Z statistics actually used}
#'   \item{stouffer_Z}{Combined Stouffer Z-score}
#'   \item{stouffer_p}{Analytic p-value under the chosen \code{alternative}}
#' }
#' Rows are sorted by \code{stouffer_p} in increasing order (smallest p first; NA at bottom).
#' If \code{output = TRUE}, an attribute \code{"file"} is attached with the CSV path.
#'
#' @examples
#' \dontrun{
#' # Load gene results from MAGMA output (must include ZSTAT column)
#' gene_results <- read.delim("magma_output.genes.out")
#'
#' # Run Stouffer's Z method on maize PMN pathways
#' stouffer_res <- magcat_stoufferZ_pathways(
#'   gene_results = gene_results,
#'   species = "maize",
#'   gene_col = "GENE",
#'   z_col = "ZSTAT"
#' )
#' head(stouffer_res)
#'
#' # With weighted Stouffer (e.g., by number of SNPs)
#' stouffer_weighted <- magcat_stoufferZ_pathways(
#'   gene_results = gene_results,
#'   species = "maize",
#'   weight_col = "NSNPS"
#' )
#' }
#'
#' @seealso \code{\link{magcat_acat_pathways}}, \code{\link{magcat_fisher_pathways}}
#' @export
magcat_stoufferZ_pathways <- function(gene_results,
                                      pathways     = NULL,
                                      species      = NULL,
                                      pmn_gene_col = NULL,
                                      gene_col     = "GENE",
                                      z_col        = "ZSTAT",
                                      weight_col   = NULL,
                                      min_abs_w    = 1e-8,
                                      alternative  = c("greater", "two.sided", "less"),
                                      output       = FALSE,
                                      out_dir      = "magcat_stouffer_z") {

  alternative <- match.arg(alternative)

  ## -------- decide pathway source: pathways vs species ----------
  if (!is.null(pathways) && !is.null(species)) {
    stop("Provide either 'pathways' OR 'species', not both.", call. = FALSE)
  }
  if (is.null(pathways) && is.null(species)) {
    stop("You must provide either:\n",
         "  * 'pathways' (list/data.frame), OR\n",
         "  * 'species' = 'maize' | 'sorghum' | 'arabidopsis' | 'plant'.",
         call. = FALSE)
  }

  # If user gave species, load PMN pathways
  if (is.null(pathways) && !is.null(species)) {
    pathways <- if (is.null(pmn_gene_col)) {
      magcat_load_pathways(species = species)
    } else {
      magcat_load_pathways(species = species, gene_col = pmn_gene_col)
    }
  }

  ## -------- standardize gene_results ----------
  if (!all(c(gene_col, z_col) %in% names(gene_results))) {
    stop("gene_results must contain columns '", gene_col, "' and '", z_col, "'.",
         call. = FALSE)
  }
  if (!is.null(weight_col) && !(weight_col %in% names(gene_results))) {
    stop("weight_col='", weight_col, "' not found in gene_results.", call. = FALSE)
  }

  genes_all  <- as.character(gene_results[[gene_col]])
  z_all      <- as.numeric(gene_results[[z_col]])
  genes_norm <- tolower(genes_all)

  ok <- !is.na(genes_norm) & genes_norm != "" & is.finite(z_all) & !is.na(z_all)
  if (!any(ok)) {
    stop("No valid (non-NA) gene IDs and Z statistics found in gene_results.", call. = FALSE)
  }

  # Robust to duplicates: keep first Z per gene (MAGMA genes.out should be unique anyway)
  z_vec <- tapply(z_all[ok], genes_norm[ok], function(x) x[1])

  # Map normalized -> canonical gene ID (first occurrence)
  gene_map <- tapply(genes_all[ok], genes_norm[ok], function(x) x[1])

  # Optional weights
  w_vec <- NULL
  if (!is.null(weight_col)) {
    w_all <- as.numeric(gene_results[[weight_col]])
    ok_w  <- !is.na(genes_norm) & genes_norm != "" & is.finite(w_all) & !is.na(w_all)
    if (any(ok_w)) {
      w_vec <- tapply(w_all[ok_w], genes_norm[ok_w], function(x) x[1])
    } else {
      # weight_col provided but all weights invalid -> fallback to equal weights
      w_vec <- NULL
    }
  }

  ## -------- pathways -> named list + names ----------
  if (is.data.frame(pathways)) {
    if (!("pathway_id" %in% names(pathways)) || !("gene_id" %in% names(pathways))) {
      stop("If 'pathways' is a data.frame it must have columns 'pathway_id' and 'gene_id'.",
           call. = FALSE)
    }
    if (!("pathway_name" %in% names(pathways))) {
      pathways$pathway_name <- pathways$pathway_id
    }

    p_list  <- split(as.character(pathways$gene_id), pathways$pathway_id)
    p_names <- tapply(pathways$pathway_name, pathways$pathway_id, function(x) x[1])

  } else if (is.list(pathways)) {
    p_list <- pathways
    nm <- names(p_list)
    if (is.null(nm)) {
      nm <- paste0("PWY_", seq_along(p_list))
      names(p_list) <- nm
    }
    p_names <- stats::setNames(nm, nm)

  } else {
    stop("'pathways' must be either a list or a data.frame.", call. = FALSE)
  }

  # normalize IDs inside pathways
  p_list <- lapply(p_list, function(g) tolower(as.character(g)))

  ## -------- helpers ----------
  .fix_w <- function(w) {
    w <- as.numeric(w)
    bad <- !is.finite(w) | is.na(w) | w <= 0
    if (any(bad)) w[bad] <- min_abs_w
    w
  }

  .p_from_Z <- function(Zs) {
    if (!is.finite(Zs) || is.na(Zs)) return(NA_real_)
    switch(alternative,
      "two.sided" = 2 * stats::pnorm(-abs(Zs)),
      "greater"   = stats::pnorm(Zs, lower.tail = FALSE),
      "less"      = stats::pnorm(Zs, lower.tail = TRUE)
    )
  }

  .stouffer_from_z <- function(z_i, w_i = NULL) {
    z_i <- as.numeric(z_i)
    z_i <- z_i[is.finite(z_i) & !is.na(z_i)]
    if (!length(z_i)) return(list(Z = NA_real_, p = NA_real_))

    if (is.null(w_i)) {
      w_i <- rep(1, length(z_i))
    } else {
      w_i <- .fix_w(w_i)
      if (length(w_i) != length(z_i)) w_i <- rep(1, length(z_i))
    }

    denom <- sqrt(sum(w_i^2))
    if (!is.finite(denom) || denom <= 0) return(list(Z = NA_real_, p = NA_real_))

    Zs <- sum(w_i * z_i) / denom
    list(Z = Zs, p = .p_from_Z(Zs))
  }

  ## -------- result container ----------
  pw_ids <- names(p_list)
  pw_names <- unname(p_names[pw_ids])
  pw_names[is.na(pw_names) | pw_names == ""] <- pw_ids[is.na(pw_names) | pw_names == ""]

  res <- data.frame(
    pathway_id   = pw_ids,
    pathway_name = pw_names,
    n_genes      = NA_integer_,
    gene_names   = NA_character_,
    gene_zvals   = NA_character_,
    stouffer_Z   = NA_real_,
    stouffer_p   = NA_real_,
    stringsAsFactors = FALSE
  )

  ## -------- main loop ----------
  for (i in seq_along(p_list)) {
    genes_i <- p_list[[i]]

    z_i <- z_vec[genes_i]
    z_i <- as.numeric(z_i[!is.na(z_i) & is.finite(z_i)])
    genes_used_norm <- names(z_i)

    d <- length(z_i)
    res$n_genes[i] <- d
    if (d == 0L) next

    canon_ids <- unname(gene_map[genes_used_norm])
    canon_ids <- unique(canon_ids[!is.na(canon_ids) & canon_ids != ""])
    res$gene_names[i] <- paste(canon_ids, collapse = ";")
    res$gene_zvals[i] <- paste(z_i, collapse = ";")

    w_i <- NULL
    if (!is.null(w_vec)) {
      # align weights to the genes actually used
      w_i <- as.numeric(w_vec[genes_used_norm])
    }

    st <- .stouffer_from_z(z_i, w_i)
    res$stouffer_Z[i] <- st$Z
    res$stouffer_p[i] <- st$p
  }

  ## -------- sort by p (smallest first) ----------
  ord <- order(res$stouffer_p, decreasing = FALSE, na.last = TRUE)
  res <- res[ord, , drop = FALSE]

  ## -------- optional CSV output ----------
  if (isTRUE(output)) {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    species_tag <- if (is.null(species)) "custom" else species
    out_path <- file.path(out_dir, paste0("magcat_stoufferZ_pathways_", species_tag, ".csv"))
    utils::write.csv(res, out_path, row.names = FALSE)
    attr(res, "file") <- out_path
  }

  res
}
