#' Pathway-level wFisher test on MAGMA gene p-values
#'
#' @param gene_results data.frame with at least gene + p-value columns.
#'   Typically a MAGMA `.genes.out` file, with "GENE" and "P" (and "ZSTAT").
#' @param pathways either:
#'   * a named list: each element is a character vector of gene IDs
#'     in that pathway, OR
#'   * a data.frame with columns `pathway_id`, `gene_id`
#'     (and optional `pathway_name`).
#' @param species optional; one of "maize", "sorghum", "arabidopsis", "plant", "fly".
#'   If provided, built-in PMN pathways are loaded via `magcat_load_pathways()`.
#'   You must provide either `pathways` OR `species`, but not both.
#' @param pmn_gene_col optional; passed to `magcat_load_pathways(gene_col=...)`
#'   when `species` is used.
#' @param gene_col column name in `gene_results` containing gene IDs
#'   (default "GENE").
#' @param p_col column name in `gene_results` containing gene-level p-values
#'   (default "P").
#' @param effect_col character. Column name containing effect sizes (beta/log-odds
#'   or Z-statistic) to derive effect directions. Default "ZSTAT". Required when
#'   `is_onetail = FALSE` for directional weighted Fisher combination.
#' @param weight_col optional column in `gene_results` to use as weight
#'   (e.g. "NSNPS"). If NULL, all weights = 1 (classic Fisher).
#' @param is_onetail logical; passed to `wFisher(is.onetail=...)`.
#'   Default FALSE (two-sided p-values + directions).
#' @param min_p lower cap for very small p-values (default 1e-15).
#' @param do_fix logical; if TRUE, clean p-values with `fix_p_for_acat()`.
#' @param output logical; if TRUE, write a CSV of results to `out_dir`.
#' @param out_dir directory to write CSV when `output = TRUE`
#'   (default "magcat_wfisher").
#'
#' @return data.frame with columns:
#'   * pathway_id
#'   * pathway_name
#'   * n_genes
#'   * gene_names
#'   * wfisher_p
#'   * overall_eff_direction
#'   Sorted by `wfisher_p` ascending (most significant first).
#'   If `output = TRUE`, an attribute `"file"` is attached with the CSV path.
#'
#' @examples
#' \dontrun{
#' # Load gene results from MAGMA output
#' gene_results <- read.delim("magma_output.genes.out")
#'
#' # Run weighted Fisher on maize PMN pathways
#' wfisher_res <- magcat_wfisher_pathways(
#'   gene_results = gene_results,
#'   species = "maize",
#'   gene_col = "GENE",
#'   p_col = "P",
#'   effect_col = "ZSTAT",
#'   is_onetail = FALSE
#' )
#' head(wfisher_res)
#'
#' # With weights from NSNPS column
#' wfisher_weighted <- magcat_wfisher_pathways(
#'   gene_results = gene_results,
#'   species = "maize",
#'   weight_col = "NSNPS"
#' )
#' }
#'
#' @seealso \code{\link{magcat_fisher_pathways}}
#' @export
magcat_wfisher_pathways <- function(gene_results,
                                    pathways     = NULL,
                                    species      = NULL,
                                    pmn_gene_col = NULL,
                                    gene_col     = "GENE",
                                    p_col        = "P",
                                    effect_col   = "ZSTAT",
                                    weight_col   = NULL,
                                    is_onetail   = FALSE,
                                    min_p        = 1e-15,
                                    do_fix       = TRUE,
                                    output       = FALSE,
                                    out_dir      = "magcat_wfisher") {
  if (!exists("wFisher")) {
    stop("Function 'wFisher' not found. Please load the package/code that defines wFisher().",
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
    if (is.null(pmn_gene_col)) {
      pathways <- magcat_load_pathways(species = species)
    } else {
      pathways <- magcat_load_pathways(
        species  = species,
        gene_col = pmn_gene_col
      )
    }
  }

  ## -------- standardize gene_results ----------
  needed_cols <- c(gene_col, p_col)
  if (!all(needed_cols %in% names(gene_results))) {
    stop("gene_results must contain columns '", gene_col,
         "' and '", p_col, "'.",
         call. = FALSE)
  }
  if (!is_onetail && !effect_col %in% names(gene_results)) {
    stop("effect_col = '", effect_col,
         "' not found in gene_results. Needed when is_onetail = FALSE.",
         call. = FALSE)
  }

  gr <- gene_results
  genes_all <- as.character(gr[[gene_col]])
  p_all     <- gr[[p_col]]

  genes_all_norm <- tolower(genes_all)
  gene_p_vec     <- stats::setNames(p_all, genes_all_norm)

  # effect signs if needed
  if (!is_onetail) {
    eff_all <- gr[[effect_col]]
    eff_vec <- stats::setNames(eff_all, genes_all_norm)
  } else {
    eff_vec <- NULL
  }

  # weights if requested
  if (!is.null(weight_col)) {
    if (!weight_col %in% names(gene_results)) {
      stop("weight_col = '", weight_col,
           "' not found in gene_results.",
           call. = FALSE)
    }
    w_all <- gr[[weight_col]]
    w_vec <- stats::setNames(w_all, genes_all_norm)
  } else {
    w_vec <- NULL
  }

  gene_map <- tapply(genes_all, genes_all_norm, function(x) x[1])

  ## -------- pathways -> named list + names ----------
  if (is.data.frame(pathways)) {
    if (!"pathway_id" %in% names(pathways) ||
        !"gene_id"    %in% names(pathways)) {
      stop("If 'pathways' is a data.frame it must have columns ",
           "'pathway_id' and 'gene_id'.",
           call. = FALSE)
    }
    if (!"pathway_name" %in% names(pathways)) {
      pathways$pathway_name <- pathways$pathway_id
    }

    p_list  <- split(pathways$gene_id, pathways$pathway_id)
    p_names <- tapply(
      pathways$pathway_name,
      pathways$pathway_id,
      FUN = function(x) x[1]
    )
  } else if (is.list(pathways)) {
    p_list  <- pathways
    p_names <- names(pathways)
    if (is.null(p_names)) {
      p_names <- paste0("PWY_", seq_along(pathways))
      names(p_list) <- p_names
    }
  } else {
    stop("'pathways' must be either a list or a data.frame.",
         call. = FALSE)
  }

  p_list <- lapply(p_list, function(g) tolower(as.character(g)))

  ## -------- prepare result container ----------
  n_pw <- length(p_list)
  res  <- data.frame(
    pathway_id           = names(p_list),
    pathway_name         = if (is.null(p_names)) names(p_list) else unname(p_names),
    n_genes              = NA_integer_,
    gene_names           = NA_character_,
    wfisher_p            = NA_real_,
    overall_eff_direction= NA_character_,
    stringsAsFactors     = FALSE
  )

  ## -------- loop over pathways ----------
  for (i in seq_len(n_pw)) {
    genes_i_norm <- p_list[[i]]
    p_i          <- gene_p_vec[genes_i_norm]

    # drop NAs
    keep <- !is.na(p_i)
    p_i  <- p_i[keep]
    genes_used_norm <- genes_i_norm[keep]

    d <- length(p_i)
    res$n_genes[i] <- d
    if (d == 0L) next

    # canonical gene IDs
    canon_ids <- unname(gene_map[genes_used_norm])
    canon_ids <- unique(canon_ids[!is.na(canon_ids)])
    res$gene_names[i] <- paste(canon_ids, collapse = ";")

    if (do_fix) {
      p_i <- fix_p_for_acat(p_i, min_p = min_p)
    }

    # weights
    if (!is.null(w_vec)) {
      w_i <- w_vec[genes_used_norm]
      # basic cleanup: replace non-positive/NA with median of positives or 1
      if (any(is.na(w_i) | w_i <= 0)) {
        w_pos <- w_i[!is.na(w_i) & w_i > 0]
        repl  <- if (length(w_pos)) stats::median(w_pos) else 1
        w_i[is.na(w_i) | w_i <= 0] <- repl
      }
    } else {
      w_i <- rep(1, d)
    }

    if (!is_onetail) {
      eff_i <- eff_vec[genes_used_norm]
      sgn   <- sign(eff_i)
      sgn[is.na(sgn) | sgn == 0] <- 1
      out <- wFisher(
        p         = p_i,
        weight    = w_i,
        is.onetail= FALSE,
        eff.sign  = sgn
      )
    } else {
      out <- wFisher(
        p         = p_i,
        weight    = w_i,
        is.onetail= TRUE
      )
    }

    res$wfisher_p[i]            <- out$p
    res$overall_eff_direction[i]<- out$overall.eff.direction
  }

  ## -------- sort by wfisher_p (ascending) ----------
  ord <- order(res$wfisher_p, decreasing = FALSE, na.last = TRUE)
  res <- res[ord, , drop = FALSE]

  ## -------- optional CSV output ----------
  if (output) {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    species_tag <- if (is.null(species)) "custom" else species
    out_path <- file.path(
      out_dir,
      paste0("magcat_wfisher_pathways_", species_tag, ".csv")
    )
    utils::write.csv(res, out_path, row.names = FALSE)
    attr(res, "file") <- out_path
  }

  res
}
