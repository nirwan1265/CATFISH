#' Pathway-level soft TFisher (stat.soft / p.soft) on MAGMA gene p-values
#'
#' For each pathway:
#'   * take gene-level p-values (e.g. MAGMA .genes.out),
#'   * compute soft-threshold TFisher statistic via TFisher::stat.soft(),
#'   * optionally get analytic p via TFisher::p.soft() under independence
#'     (converted to a right-tail test p = 1 - F(q)),
#'   * optionally calibrate with gene-set permutations for empirical p.
#'
#' All non-NA p-values in the pathway are used (no truncation / filtering).
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
#' @param tau1 soft TFisher truncation/normalization parameter (default 0.05);
#'   passed to stat.soft/p.soft.
#' @param min_p lower cap for very small p-values (default 1e-15).
#' @param do_fix logical; if TRUE, clean p-values with `fix_p_for_acat()`
#'   (to keep them strictly in (0,1)).
#' @param B_perm integer; number of permutations for empirical p
#'   (default 0L = no permutations).
#' @param seed optional RNG seed for permutations (default NULL).
#' @param analytic_logical logical; if TRUE, compute analytic p via p.soft()
#'   and convert to right-tail p = 1 - F(q) (default TRUE).
#' @param output logical; if TRUE, write a CSV of results to `out_dir`.
#' @param out_dir directory to write CSV when `output = TRUE`
#'   (default "magcat_tfisher_soft").
#'
#' @return data.frame with columns:
#'   * pathway_id
#'   * pathway_name
#'   * n_genes
#'   * gene_names
#'   * gene_pvals          (semicolon-separated gene-level p's used in the test)
#'   * tfisher_stat        (soft TFisher statistic)
#'   * tfisher_p_perm      (permutation p; NA if B_perm = 0)
#'   * tfisher_p_analytic  (right-tail analytic p = 1 - p.soft; NA if analytic_logical = FALSE)
#'   Sorted by `tfisher_p_perm` ascending if B_perm > 0,
#'   otherwise by `tfisher_stat` (descending).
#'   If `output = TRUE`, an attribute `"file"` is attached with the CSV path.
#' @export
magcat_soft_tfisher_pathways <- function(gene_results,
                                         pathways         = NULL,
                                         species          = NULL,
                                         pmn_gene_col     = NULL,
                                         gene_col         = "GENE",
                                         p_col            = "P",
                                         tau1             = 0.05,
                                         min_p            = 1e-15,
                                         do_fix           = TRUE,
                                         B_perm           = 0L,
                                         seed             = NULL,
                                         analytic_logical = TRUE,
                                         output           = FALSE,
                                         out_dir          = "magcat_tfisher_soft") {

  if (!requireNamespace("TFisher", quietly = TRUE)) {
    stop("magcat_soft_tfisher_pathways(): package 'TFisher' is required.",
         call. = FALSE)
  }

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
    pathway_id         = names(p_list),
    pathway_name       = if (is.null(p_names)) names(p_list) else unname(p_names),
    n_genes            = NA_integer_,
    gene_names         = NA_character_,
    gene_pvals         = NA_character_,
    tfisher_stat       = NA_real_,
    tfisher_p_perm     = NA_real_,
    tfisher_p_analytic = NA_real_,
    stringsAsFactors   = FALSE
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

    ## ---- soft TFisher statistic ----
    stat <- TFisher::stat.soft(p = p_i, tau1 = tau1)
    res$tfisher_stat[i] <- stat

    ## ---- analytic p (optional): right-tail p = 1 - F(q) ----
    if (analytic_logical) {
      Fq <- TFisher::p.soft(
        q    = stat,
        n    = length(p_i),
        tau1 = tau1,
        M    = NULL
      )
      res$tfisher_p_analytic[i] <- 1 - Fq
    }

    ## ---- permutations (optional, right-tail) ----
    if (B_perm > 0L && !is.na(stat)) {
      tf_perm <- numeric(B_perm)

      for (b in seq_len(B_perm)) {
        idx_perm <- sample(idx_pool, size = d, replace = FALSE)
        p_perm   <- as.numeric(p_all[idx_perm])

        # drop NAs in permuted set
        keep_perm <- !is.na(p_perm)
        p_perm    <- p_perm[keep_perm]

        if (!length(p_perm)) {
          tf_perm[b] <- NA_real_
          next
        }

        if (do_fix) {
          p_perm <- fix_p_for_acat(p_perm, min_p = min_p)
        }

        tf_perm[b] <- TFisher::stat.soft(p = p_perm, tau1 = tau1)
      }

      # right-tail: large statistic = more small p's = more signal
      res$tfisher_p_perm[i] <-
        (1 + sum(tf_perm >= stat, na.rm = TRUE)) / (B_perm + 1)
    }
  }

  ## -------- sort ----------
  if (B_perm > 0L) {
    ord <- order(res$tfisher_p_perm, decreasing = FALSE, na.last = TRUE)
  } else {
    ord <- order(res$tfisher_stat, decreasing = TRUE, na.last = TRUE)
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
      paste0("magcat_tfisher_soft_pathways_", species_tag, ".csv")
    )
    utils::write.csv(res, out_path, row.names = FALSE)
    attr(res, "file") <- out_path
  }

  res
}
