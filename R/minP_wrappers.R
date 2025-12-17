#' Pathway-level minP on MAGMA gene p-values
#'
#' For each pathway:
#'   * take gene-level p-values (e.g. MAGMA .genes.out),
#'   * compute the minimum p-value statistic: min(p_i),
#'   * optionally get analytic combined p via metap::minimump()
#'     (Wilkinson / Tippett type minP),
#'   * optionally calibrate with gene-set permutations for empirical p.
#'
#' All non-NA p-values in the pathway are used (no truncation / filtering),
#' except for cleaning to keep them strictly in (0,1) via fix_p_for_acat().
#'
#' @param gene_results data.frame with at least gene + p-value columns.
#'   Typically a MAGMA `.genes.out` file, with "GENE" and "P".
#' @param pathways either:
#'   * a named list: each element is a character vector of gene IDs
#'     in that pathway, OR
#'   * a data.frame with columns `pathway_id`, `gene_id`
#'     (and optional `pathway_name`).
#' @param species optional; one of "maize", "sorghum", "arabidopsis", "plant".
#'   If provided, built-in PMN pathways are loaded via `magcat_load_pathways()`.
#'   You must provide either `pathways` OR `species`, but not both.
#' @param pmn_gene_col optional; passed to `magcat_load_pathways(gene_col=...)`
#'   when `species` is used.
#' @param gene_col column name in `gene_results` containing gene IDs
#'   (default "GENE").
#' @param p_col column name in `gene_results` containing gene-level p-values
#'   (default "P").
#' @param min_p lower cap for very small p-values (default 1e-15).
#' @param do_fix logical; if TRUE, clean p-values with `fix_p_for_acat()`
#'   (to keep them strictly in (0,1)).
#' @param B_perm integer; number of permutations for empirical p
#'   (default 1000L = no permutations).
#' @param seed optional RNG seed for permutations (default NULL).
#' @param analytic_logical logical; if TRUE, compute analytic p via
#'   `metap::minimump()`. If `metap` is not available, falls back to
#'   Sidák: 1 - (1 - min_p)^K. (default TRUE)
#' @param output logical; if TRUE, write a CSV of results to `out_dir`.
#' @param out_dir directory to write CSV when `output = TRUE`
#'   (default "magcat_minp").
#'
#' @return data.frame with columns:
#'   * pathway_id
#'   * pathway_name
#'   * n_genes
#'   * gene_names
#'   * gene_pvals        (semicolon-separated gene-level p's used in the test)
#'   * minp_stat         (min(p_i) for that pathway)
#'   * minp_p_perm       (permutation p; NA if B_perm = 0)
#'   * minp_p_analytic   (analytic p from minimump / Sidák; NA if
#'                        analytic_logical = FALSE)
#'   Sorted by `minp_p_perm` ascending if B_perm > 0,
#'   otherwise by `minp_p_analytic` (ascending) if available,
#'   else by `minp_stat` (ascending).
#'   If `output = TRUE`, an attribute `"file"` is attached with the CSV path.
#' @export
magcat_minp_pathways <- function(gene_results,
                                 pathways         = NULL,
                                 species          = NULL,
                                 pmn_gene_col     = NULL,
                                 gene_col         = "GENE",
                                 p_col            = "P",
                                 min_p            = 1e-15,
                                 do_fix           = TRUE,
                                 B_perm           = 0L,
                                 seed             = NULL,
                                 analytic_logical = TRUE,
                                 output           = FALSE,
                                 out_dir          = "magcat_minp") {

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

  gr         <- gene_results
  genes_all  <- as.character(gr[[gene_col]])
  p_all      <- gr[[p_col]]
  genes_norm <- tolower(genes_all)

  # map from normalized -> canonical gene ID
  gene_map <- tapply(genes_all, genes_norm, function(x) x[1])

  # named vector of p's by normalized ID
  gene_p_vec <- stats::setNames(p_all, genes_norm)

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

  # normalize IDs inside pathways
  p_list <- lapply(p_list, function(g) tolower(as.character(g)))

  ## -------- result container ----------
  n_pw <- length(p_list)
  res  <- data.frame(
    pathway_id       = names(p_list),
    pathway_name     = if (is.null(p_names)) names(p_list) else unname(p_names),
    n_genes          = NA_integer_,
    gene_names       = NA_character_,
    gene_pvals       = NA_character_,
    minp_stat        = NA_real_,
    minp_p_perm      = NA_real_,
    minp_p_analytic  = NA_real_,
    stringsAsFactors = FALSE
  )

  ## -------- pool for permutations ----------
  idx_pool <- which(!is.na(p_all))
  if (length(idx_pool) == 0L) {
    stop("No non-NA gene-level p-values found in 'gene_results'.",
         call. = FALSE)
  }

  if (B_perm > 0L && !is.null(seed)) {
    set.seed(seed)
  }

  ## -------- helper: analytic minP combiner ----------
  .analytic_minp <- function(p_vec) {
    p_vec <- as.numeric(p_vec)
    p_vec <- p_vec[!is.na(p_vec) & p_vec > 0 & p_vec < 1]
    if (!length(p_vec)) return(NA_real_)

    if (requireNamespace("metap", quietly = TRUE)) {
      out <- tryCatch(
        metap::minimump(p_vec)$p,
        error = function(e) NA_real_
      )
      return(out)
    } else {
      # Sidák under independence
      k     <- length(p_vec)
      p_min <- min(p_vec)
      return(1 - (1 - p_min)^k)
    }
  }

  ## -------- main loop over pathways ----------
  for (i in seq_len(n_pw)) {
    genes_i_norm <- p_list[[i]]

    # pull gene-level p's
    p_i <- gene_p_vec[genes_i_norm]

    # drop NAs
    keep <- !is.na(p_i)
    p_i  <- as.numeric(p_i[keep])
    genes_used_norm <- genes_i_norm[keep]

    d <- length(p_i)
    res$n_genes[i] <- d
    if (d == 0L) next

    # canonical gene IDs for reporting
    canon_ids <- unname(gene_map[genes_used_norm])
    canon_ids <- unique(canon_ids[!is.na(canon_ids)])
    res$gene_names[i] <- paste(canon_ids, collapse = ";")

    # keep p's strictly inside (0,1) if requested
    if (do_fix) {
      p_i <- fix_p_for_acat(p_i, min_p = min_p)
    }

    # store the gene-level p's we actually used
    res$gene_pvals[i] <- paste(p_i, collapse = ";")

    ## ---- minP statistic ----
    min_stat <- min(p_i)
    res$minp_stat[i] <- min_stat

    ## ---- analytic p (optional) ----
    if (analytic_logical) {
      res$minp_p_analytic[i] <- .analytic_minp(p_i)
    }

    ## ---- permutations (optional, left-tail: smaller = more extreme) ----
    if (B_perm > 0L && is.finite(min_stat)) {
      min_perm <- numeric(B_perm)

      for (b in seq_len(B_perm)) {
        idx_perm <- sample(idx_pool, size = d, replace = FALSE)
        p_perm   <- as.numeric(p_all[idx_perm])

        # drop NAs in permuted set
        keep_perm <- !is.na(p_perm)
        p_perm    <- p_perm[keep_perm]

        if (!length(p_perm)) {
          min_perm[b] <- NA_real_
          next
        }

        if (do_fix) {
          p_perm <- fix_p_for_acat(p_perm, min_p = min_p)
        }

        min_perm[b] <- min(p_perm)
      }

      # left-tail: smaller min(p) = stronger signal
      res$minp_p_perm[i] <-
        (1 + sum(min_perm <= min_stat, na.rm = TRUE)) / (B_perm + 1)
    }
  }

  ## -------- sort ----------
  if (B_perm > 0L) {
    ord <- order(res$minp_p_perm, decreasing = FALSE, na.last = TRUE)
  } else if (analytic_logical) {
    ord <- order(res$minp_p_analytic, decreasing = FALSE, na.last = TRUE)
  } else {
    ord <- order(res$minp_stat, decreasing = FALSE, na.last = TRUE)
  }
  res <- res[ord, , drop = FALSE]

  ## -------- optional CSV output ----------
  if (output) {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    species_tag <- if (is.null(species)) "custom" else species
    out_path <- file.path(
      out_dir,
      paste0("magcat_minp_pathways_", species_tag, ".csv")
    )
    utils::write.csv(res, out_path, row.names = FALSE)
    attr(res, "file") <- out_path
  }

  res
}
