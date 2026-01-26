#' Pathway-level Fisher test on gene-level p-values
#'
#' Computes classic Fisher's combination test within each pathway:
#' \deqn{T = -2 \sum_{g \in S} \log(p_g), \qquad T \sim \chi^2_{2m}}
#' where \eqn{m} is the number of genes with available p-values in the pathway.
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
#'   If provided, built-in pathways are loaded via \code{magcat_load_pathways()}.
#'   Provide either \code{pathways} OR \code{species}, but not both.
#' @param pmn_gene_col Optional; passed to \code{magcat_load_pathways(gene_col=...)}
#'   when \code{species} is used. If NULL, the PMN loader prefers "Gene-name"
#'   then "Gene-id".
#' @param gene_col Column name in \code{gene_results} containing gene IDs (default "GENE").
#' @param p_col Column name in \code{gene_results} containing gene-level p-values (default "P").
#' @param min_p Lower cap for very small p-values (default 1e-15).
#' @param do_fix If TRUE (default), clean/cap p-values using \code{fix_p_for_acat()}
#'   to ensure values lie strictly in (0,1).
#' @param output If TRUE, write a CSV of results to \code{out_dir}.
#' @param out_dir Directory to write CSV when \code{output = TRUE} (default "magcat_fisher").
#'
#' @return A data.frame with columns:
#' \describe{
#'   \item{pathway_id}{Pathway identifier}
#'   \item{pathway_name}{Pathway name (if available; otherwise equals \code{pathway_id})}
#'   \item{n_genes}{Number of genes used in Fisher's test for this pathway}
#'   \item{gene_names}{Semicolon-separated gene IDs actually used (from \code{gene_results})}
#'   \item{fisher_p}{Fisher combined p-value}
#' }
#' Rows are sorted by \code{fisher_p} in increasing order (smallest p first; NA at bottom).
#' If \code{output = TRUE}, an attribute \code{"file"} is attached with the CSV path.
#'
#' @examples
#' \dontrun{
#' # Load gene results from MAGMA output
#' gene_results <- read.delim("magma_output.genes.out")
#'
#' # Run Fisher's method on maize PMN pathways
#' fisher_res <- magcat_fisher_pathways(
#'   gene_results = gene_results,
#'   species = "maize",
#'   gene_col = "GENE",
#'   p_col = "P"
#' )
#' head(fisher_res)
#' }
#'
#' @seealso \code{\link{magcat_acat_pathways}}, \code{\link{magcat_minp_pathways}},
#'   \code{\link{magcat_wfisher_pathways}}
#' @export
magcat_fisher_pathways <- function(gene_results,
                                   pathways     = NULL,
                                   species      = NULL,
                                   pmn_gene_col = NULL,
                                   gene_col     = "GENE",
                                   p_col        = "P",
                                   min_p        = 1e-15,
                                   do_fix       = TRUE,
                                   output       = FALSE,
                                   out_dir      = "magcat_fisher") {

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

  # If user gave species, load PMN pathways
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

  genes_all  <- as.character(gene_results[[gene_col]])
  p_all      <- as.numeric(gene_results[[p_col]])
  genes_norm <- tolower(genes_all)

  ok <- !is.na(genes_norm) & genes_norm != "" & is.finite(p_all) & !is.na(p_all)
  if (!any(ok)) {
    stop("No valid (non-NA) gene IDs and p-values found in gene_results.", call. = FALSE)
  }

  # Robust to duplicates: keep minimum p-value per gene
  gene_p_vec <- tapply(p_all[ok], genes_norm[ok], min, na.rm = TRUE)

  # Map normalized -> canonical gene ID (first occurrence)
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

  # normalize IDs inside pathways
  p_list <- lapply(p_list, function(g) tolower(as.character(g)))

  ## -------- result container ----------
  pw_ids <- names(p_list)
  pw_names <- unname(p_names[pw_ids])
  pw_names[is.na(pw_names) | pw_names == ""] <- pw_ids[is.na(pw_names) | pw_names == ""]

  res <- data.frame(
    pathway_id   = pw_ids,
    pathway_name = pw_names,
    n_genes      = NA_integer_,
    gene_names   = NA_character_,
    fisher_p     = NA_real_,
    stringsAsFactors = FALSE
  )

  ## -------- loop over pathways ----------
  for (i in seq_along(p_list)) {
    genes_i_norm <- p_list[[i]]

    p_i <- gene_p_vec[genes_i_norm]
    p_i <- as.numeric(p_i[!is.na(p_i)])
    genes_used_norm <- names(p_i)

    d <- length(p_i)
    res$n_genes[i] <- d
    if (d == 0L) next

    # canonical gene IDs (original case) for genes actually used
    canon_ids <- unname(gene_map[genes_used_norm])
    canon_ids <- unique(canon_ids[!is.na(canon_ids) & canon_ids != ""])
    res$gene_names[i] <- paste(canon_ids, collapse = ";")

    if (isTRUE(do_fix)) {
      p_i <- fix_p_for_acat(p_i, min_p = min_p)
    }

    # Fisher's method is defined for m >= 1; keep it valid for d=1 too.
    stat <- -2 * sum(log(p_i))
    res$fisher_p[i] <- stats::pchisq(stat, df = 2 * d, lower.tail = FALSE)
  }

  ## -------- sort by fisher_p (smallest first) ----------
  ord <- order(res$fisher_p, decreasing = FALSE, na.last = TRUE)
  res <- res[ord, , drop = FALSE]

  ## -------- optional CSV output ----------
  if (isTRUE(output)) {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    species_tag <- if (is.null(species)) "custom" else species
    out_path <- file.path(out_dir, paste0("magcat_fisher_pathways_", species_tag, ".csv"))
    utils::write.csv(res, out_path, row.names = FALSE)
    attr(res, "file") <- out_path
  }

  res
}
