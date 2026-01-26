#' Pathway-level *adaptive* soft TFisher (omnibus) on gene-level p-values
#'
#' For each pathway, this function:
#' \itemize{
#'   \item extracts the pathway's gene-level p-values (e.g., from MAGMA `.genes.out`)
#'   \item computes soft-thresholding Fisher statistics across a grid of \code{tau} values
#'   \item computes per-\code{tau} right-tail p-values under the null (via \code{TFisher::p.soft})
#'   \item defines the adaptive score as \eqn{W_o = \min_{\tau \in \text{grid}} p_\tau}
#'   \item computes a calibrated adaptive p-value under independence using
#'         \code{TFisher::p.soft.omni(q = W_o, ...)}
#' }
#'
#' Gene IDs are matched case-insensitively: both \code{gene_results[[gene_col]]} and
#' pathway gene IDs are normalized to lower-case before matching. Genes not present
#' in \code{gene_results} are dropped from that pathway.
#'
#' @param gene_results data.frame with at least a gene ID column and a p-value column.
#'   For MAGMA `.genes.out`, these are commonly \code{GENE} and \code{P}.
#' @param pathways Either:
#' \itemize{
#'   \item A named list: each element is a character vector of gene IDs in that pathway, OR
#'   \item A data.frame with columns \code{pathway_id}, \code{gene_id} (and optional \code{pathway_name}).
#' }
#' @param species Optional; one of "maize", "sorghum", "arabidopsis", "plant", "fly".
#'   If provided, built-in PMN pathways are loaded via \code{magcat_load_pathways()}.
#'   Provide either \code{pathways} OR \code{species}, but not both.
#' @param pmn_gene_col Optional; passed to \code{magcat_load_pathways(gene_col=...)}
#'   when \code{species} is used. If NULL, the PMN loader prefers "Gene-name" then "Gene-id".
#' @param gene_col Column name in \code{gene_results} containing gene IDs (default "GENE").
#' @param p_col Column name in \code{gene_results} containing gene-level p-values (default "P").
#' @param tau_grid Numeric vector of \code{tau} values to scan. Internally sorted increasing,
#'   as required by \code{TFisher::stat.soft.omni()}.
#' @param min_p Lower cap for extremely small p-values before TFisher math (default 1e-15).
#' @param do_fix If TRUE (default), clean/cap p-values using \code{fix_p_for_acat()}
#'   to ensure values lie strictly in (0,1).
#' @param output If TRUE, write a CSV of results to \code{out_dir}.
#' @param out_dir Directory to write CSV when \code{output = TRUE}
#'   (default "magcat_tfisher_soft_adaptive").
#'
#' @return A data.frame with columns:
#' \describe{
#'   \item{pathway_id}{Pathway identifier}
#'   \item{pathway_name}{Pathway name (if available; otherwise equals \code{pathway_id})}
#'   \item{n_genes}{Number of genes used in the pathway}
#'   \item{gene_names}{Semicolon-separated gene IDs actually used (from \code{gene_results})}
#'   \item{gene_pvals}{Semicolon-separated p-values actually used}
#'   \item{tau_hat}{Tau achieving the minimum single-\code{tau} p-value}
#'   \item{tfisher_stat_hat}{Soft TFisher statistic at \code{tau_hat}}
#'   \item{tfisher_p_hat}{Single-\code{tau} right-tail p-value at \code{tau_hat}}
#'   \item{tfisher_p_min}{Minimum single-\code{tau} p-value across \code{tau_grid} (= \eqn{W_o})}
#'   \item{tfisher_p_omni}{Adaptive omnibus p-value for \eqn{W_o} under the null (independence)}
#' }
#'
#' Rows are sorted by \code{tfisher_p_omni} in increasing order (smallest p first; NA at bottom).
#' If \code{output = TRUE}, an attribute \code{"file"} is attached with the CSV path.
#'
#' @examples
#' \dontrun{
#' # Load gene results from MAGMA output
#' gene_results <- read.delim("magma_output.genes.out")
#'
#' # Run adaptive soft TFisher on maize PMN pathways
#' tfisher_adapt_res <- magcat_soft_tfisher_adaptive_pathways(
#'   gene_results = gene_results,
#'   species = "maize",
#'   gene_col = "GENE",
#'   p_col = "P"
#' )
#' head(tfisher_adapt_res)
#'
#' # With custom tau grid
#' tfisher_custom <- magcat_soft_tfisher_adaptive_pathways(
#'   gene_results = gene_results,
#'   species = "maize",
#'   tau_grid = c(0.1, 0.05, 0.01)
#' )
#' }
#'
#' @seealso \code{\link{magcat_tfisher_pathways}},
#'   \code{\link{magcat_soft_tfisher_pathways}}
#' @export
magcat_soft_tfisher_adaptive_pathways <- function(gene_results,
                                                 pathways     = NULL,
                                                 species      = NULL,
                                                 pmn_gene_col = NULL,
                                                 gene_col     = "GENE",
                                                 p_col        = "P",
                                                 tau_grid     = c(0.20, 0.10, 0.05, 0.02, 0.01, 0.005, 0.001),
                                                 min_p        = 1e-15,
                                                 do_fix       = TRUE,
                                                 output       = FALSE,
                                                 out_dir      = "magcat_tfisher_soft_adaptive") {
  if (!requireNamespace("TFisher", quietly = TRUE)) {
    stop("Package 'TFisher' is required for magcat_soft_tfisher_adaptive_pathways().",
         call. = FALSE)
  }

  ## -------- decide pathway source: pathways vs species ----------
  if (!is.null(pathways) && !is.null(species)) {
    stop("Provide either 'pathways' OR 'species', not both.", call. = FALSE)
  }
  if (is.null(pathways) && is.null(species)) {
    stop("You must provide either:\n",
         "  * 'pathways' (list/data.frame), OR\n",
         "  * 'species' = 'maize' | 'sorghum' | 'arabidopsis' | 'plant' | 'fly'.",
         call. = FALSE)
  }

  if (is.null(pathways) && !is.null(species)) {
    pathways <- if (is.null(pmn_gene_col)) {
      magcat_load_pathways(species = species)
    } else {
      magcat_load_pathways(species = species, gene_col = pmn_gene_col)
    }
  }

  ## -------- standardize gene_results ----------
  if (!all(c(gene_col, p_col) %in% names(gene_results))) {
    stop("gene_results must contain columns '", gene_col, "' and '", p_col, "'.",
         call. = FALSE)
  }

  genes_all <- as.character(gene_results[[gene_col]])
  p_all     <- as.numeric(gene_results[[p_col]])
  genes_norm <- tolower(genes_all)

  ok <- !is.na(genes_norm) & genes_norm != "" & is.finite(p_all) & !is.na(p_all)
  if (!any(ok)) {
    stop("No valid (non-NA) gene IDs and p-values found in gene_results.", call. = FALSE)
  }

  # Robust to duplicates: keep min p-value per gene
  gene_p_vec <- tapply(p_all[ok], genes_norm[ok], min, na.rm = TRUE)

  # Map normalized -> canonical ID (first occurrence)
  gene_map <- tapply(genes_all[ok], genes_norm[ok], function(x) x[1])

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

  # normalize gene IDs inside pathways
  p_list <- lapply(p_list, function(g) tolower(as.character(g)))

  ## -------- tau grid checks (TFisher requires non-descending TAU1) ----------
  tau_grid <- unique(as.numeric(tau_grid))
  tau_grid <- tau_grid[is.finite(tau_grid) & !is.na(tau_grid)]
  if (!length(tau_grid)) stop("tau_grid has no valid numeric values.", call. = FALSE)
  if (any(tau_grid <= 0)) stop("All tau values must be > 0.", call. = FALSE)
  tau_grid <- sort(tau_grid, decreasing = FALSE)

  ## -------- result container ----------
  pw_ids <- names(p_list)
  pw_names <- unname(p_names[pw_ids])
  pw_names[is.na(pw_names) | pw_names == ""] <- pw_ids[is.na(pw_names) | pw_names == ""]

  res <- data.frame(
    pathway_id       = pw_ids,
    pathway_name     = pw_names,
    n_genes          = NA_integer_,
    gene_names       = NA_character_,
    gene_pvals       = NA_character_,
    tau_hat          = NA_real_,
    tfisher_stat_hat = NA_real_,
    tfisher_p_hat    = NA_real_,
    tfisher_p_min    = NA_real_,
    tfisher_p_omni   = NA_real_,
    stringsAsFactors = FALSE
  )

  ## -------- main loop ----------
  for (i in seq_along(p_list)) {
    genes_i_norm <- p_list[[i]]

    p_i <- gene_p_vec[genes_i_norm]
    keep <- !is.na(p_i)
    p_i <- as.numeric(p_i[keep])
    genes_used_norm <- genes_i_norm[keep]

    d <- length(p_i)
    res$n_genes[i] <- d
    if (d == 0L) next

    # canonical gene IDs actually used
    canon_ids <- unname(gene_map[genes_used_norm])
    canon_ids <- unique(canon_ids[!is.na(canon_ids) & canon_ids != ""])
    res$gene_names[i] <- paste(canon_ids, collapse = ";")

    if (isTRUE(do_fix)) {
      p_i <- fix_p_for_acat(p_i, min_p = min_p)
    }
    res$gene_pvals[i] <- paste(p_i, collapse = ";")

    ## Per-tau stats + right-tail p-values (single-test)
    stat_vec <- vapply(tau_grid, function(tau) TFisher::stat.soft(p = p_i, tau1 = tau), numeric(1))
    cdf_vec  <- vapply(tau_grid, function(tau) TFisher::p.soft(q = TFisher::stat.soft(p = p_i, tau1 = tau),
                                                              n = d, tau1 = tau, M = NULL),
                       numeric(1))
    p_rt_vec <- 1 - as.numeric(cdf_vec)

    jhat <- which.min(p_rt_vec)
    res$tau_hat[i]          <- tau_grid[jhat]
    res$tfisher_stat_hat[i] <- stat_vec[jhat]
    res$tfisher_p_hat[i]    <- p_rt_vec[jhat]

    ## Omnibus adaptive p-value (calibrates min-over-tau under the null)
    omni_out <- TFisher::stat.soft.omni(p = p_i, TAU1 = tau_grid, M = NULL)
    Wo <- as.numeric(omni_out$omni)  # Wo = min_j p_j
    res$tfisher_p_min[i] <- Wo

    # For Wo (a minimum p-value), smaller is more extreme -> p-value is CDF, not 1-CDF
    res$tfisher_p_omni[i] <- as.numeric(TFisher::p.soft.omni(q = Wo, n = d, TAU1 = tau_grid, M = NULL))
  }

  ## -------- sort by omnibus p (smallest first) ----------
  ord <- order(res$tfisher_p_omni, decreasing = FALSE, na.last = TRUE)
  res <- res[ord, , drop = FALSE]

  ## -------- optional CSV output ----------
  if (isTRUE(output)) {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    species_tag <- if (is.null(species)) "custom" else species
    out_path <- file.path(out_dir, paste0("magcat_tfisher_soft_adaptive_pathways_", species_tag, ".csv"))
    utils::write.csv(res, out_path, row.names = FALSE)
    attr(res, "file") <- out_path
  }

  res
}
