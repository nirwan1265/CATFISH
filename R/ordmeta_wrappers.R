## ---------- helper: Beta CDF for ith order statistic ----------
#' @title Beta probability for ith order statistic
#' @param p p-value (scalar)
#' @param i rank (1..n)
#' @param n number of inputs
#' @importFrom stats pbeta
F_i <- function(p, i, n) {
  a <- i
  b <- n - i + 1
  pbeta(q = p, shape1 = a, shape2 = b, lower.tail = TRUE)
}

## ---------- internal: compute MMP statistic --------------------
ordmeta_stat <- function(p2) {
  ord  <- order(p2, decreasing = FALSE)
  pord <- sort(p2, decreasing = FALSE)
  N    <- length(p2)

  Fi_vals <- numeric(N)
  for (i in seq_len(N)) {
    Fi_vals[i] <- F_i(pord[i], i, N)
  }

  alpha       <- min(Fi_vals)
  idx_minimum <- which.min(Fi_vals)
  eff.p.idx   <- ord[1:idx_minimum]

  list(
    alpha        = alpha,
    optimal_rank = idx_minimum,
    eff.p.idx    = eff.p.idx
  )
}

## ---------- internal: Monte Carlo null for MMP -----------------
ordmeta_core <- function(p2, B = 10000L) {
  N <- length(p2)
  if (N < 1L) {
    stop("Need at least one p-value for ordmeta.")
  }
  if (B < 1L) {
    stop("B must be >= 1.")
  }

  stat      <- ordmeta_stat(p2)
  alpha_obs <- stat$alpha

  # Monte Carlo null distribution of MMP
  alpha_sim <- numeric(B)
  for (b in seq_len(B)) {
    u  <- sort(stats::runif(N))
    Fi <- numeric(N)
    for (i in seq_len(N)) {
      Fi[i] <- F_i(u[i], i, N)
    }
    alpha_sim[b] <- min(Fi)
  }

  # p-value = P(alpha_null <= alpha_obs); add 1/(B+1) for stability
  pval <- (sum(alpha_sim <= alpha_obs) + 1) / (B + 1)

  list(
    p            = min(pval, 1),
    optimal_rank = stat$optimal_rank,
    eff.p.idx    = stat$eff.p.idx,
    MMP          = alpha_obs
  )
}

## ---------- exported ordmeta() without rSymPy ------------------
#' @title ordmeta
#' @description
#' Combine p-values using the minimum marginal p-value (MMP)
#' in the joint order distribution. This implementation uses
#' Monte Carlo simulation instead of symbolic integration, so
#' it does not depend on \pkg{rSymPy}.
#'
#' @param p A numeric vector of p-values.
#' @param is.onetail Logical. If TRUE, directions are ignored.
#'   If FALSE, `eff.sign` is used to consider effect directions.
#' @param eff.sign A numeric vector of effect signs (e.g. 1 or -1),
#'   required when `is.onetail = FALSE`.
#' @param B Integer; number of Monte Carlo simulations for the
#'   null distribution (default 10000).
#'
#' @return A list with elements:
#' \describe{
#'   \item{p}{Combined p-value.}
#'   \item{optimal_rank}{The optimal rank where the minimum marginal p-value occurs.}
#'   \item{eff.p.idx}{Indices of effective p-values (1..optimal_rank in the ordered set).}
#'   \item{MMP}{The minimum marginal p-value (alpha).}
#'   \item{overall.eff.direction}{If `is.onetail = FALSE`, the overall direction ("+" or "-").}
#' }
#'
#' @examples
#' \donttest{
#' ordmeta(p = c(0.01, 0.02, 0.8, 0.25),
#'         is.onetail = FALSE,
#'         eff.sign   = c(1, 1, 1, -1),
#'         B          = 2000)
#' }
#'
#' @seealso \code{\link{magcat_ordmeta_pathways}}
#' @export
ordmeta <- function(p,
                    is.onetail = TRUE,
                    eff.sign   = NULL,
                    B          = 10000L) {
  if (is.null(p)) {
    stop("Input p-values are required.")
  }
  if (!is.onetail && is.null(eff.sign)) {
    stop("Input the direction of effects via eff.sign when is.onetail = FALSE.")
  }

  # Drop NAs, like original
  idx_na <- which(is.na(p))
  if (length(idx_na) > 0L) {
    p <- p[-idx_na]
    if (!is.onetail) {
      eff.sign <- eff.sign[-idx_na]
    }
  }

  if (length(p) == 0L) {
    stop("No non-NA p-values supplied.")
  }

  if (is.onetail) {
    RES <- ordmeta_core(p2 = p, B = B)
    return(RES)
  } else {
    p1 <- p2 <- p
    idx_pos <- which(eff.sign >= 0)
    idx_neg <- which(eff.sign < 0)

    # map to one-sided p's (same as original code)
    p1[idx_pos] <- p[idx_pos] / 2
    p1[idx_neg] <- 1 - p[idx_neg] / 2

    p2[idx_pos] <- 1 - p[idx_pos] / 2
    p2[idx_neg] <- p[idx_neg] / 2

    RES1 <- ordmeta_core(p2 = p1, B = B)
    RES2 <- ordmeta_core(p2 = p2, B = B)

    if (RES1$p <= RES2$p) {
      RES <- RES1
      RES$overall.eff.direction <- "+"
    } else {
      RES <- RES2
      RES$overall.eff.direction <- "-"
    }

    # two-sided adjustment
    RES$p <- RES$p * 2
    if (RES$p > 1) RES$p <- 1

    return(RES)
  }
}


#' Pathway-level ordmeta test on MAGMA gene p-values
#'
#' @param gene_results data.frame with at least gene + p-value columns.
#'   Typically a MAGMA `.genes.out` file, with "GENE" and "P" (and "ZSTAT").
#' @param pathways either:
#'   * a named list: each element is a character vector of gene IDs
#'     in that pathway, OR
#'   * a data.frame with columns `pathway_id`, `gene_id`
#'     (and optional `pathway_name`).
#' @param species optional; one of "maize", "sorghum", "arabidopsis", "plant".
#'   If provided, built-in PMN pathways are loaded via `magcat_load_pathways()`.
#'   You must provide either `pathways` OR `species`, but not both.
#' @param pmn_gene_col optional; passed to `magcat_load_pathways(gene_col=...)`
#'   when `species` is used. If NULL, PMN loader prefers "Gene-name"
#'   and falls back to "Gene-id".
#' @param gene_col column name in `gene_results` containing gene IDs
#'   (default "GENE").
#' @param p_col column name in `gene_results` containing gene-level p-values
#'   (default "P").
#' @param effect_col character. Column name containing effect sizes (beta/log-odds
#'   or Z-statistic) to derive effect directions (+/-). Default "ZSTAT". Required
#'   when `is_onetail = FALSE` for directional ordmeta combination.
#' @param is_onetail logical; passed to `ordmeta(is.onetail=...)`.
#'   Default FALSE (two-sided p-values + directions).
#' @param min_p lower cap for very small p-values (default 1e-15).
#' @param do_fix logical; if TRUE, clean p-values with `fix_p_for_acat()`
#'   (0 -> min_p, 1 -> 1 - 1/d).
#' @param B Integer; number of Monte Carlo simulations used inside `ordmeta`
#'   (default 10000).
#' @param output logical; if TRUE, write a CSV of results to `out_dir`.
#' @param out_dir directory to write CSV when `output = TRUE`
#'   (default "magcat_ordmeta").
#'
#' @return data.frame with columns:
#'   * pathway_id
#'   * pathway_name
#'   * n_genes
#'   * gene_names
#'   * ordmeta_p
#'   * optimal_rank
#'   * eff_p_idx  (semicolon-separated indices)
#'   * MMP
#'   * overall_eff_direction
#'   Sorted by `ordmeta_p` ascending (most significant first).
#'   If `output = TRUE`, an attribute `"file"` is attached with the CSV path.
#'
#' @examples
#' \dontrun{
#' # Load gene results from MAGMA output
#' gene_results <- read.delim("magma_output.genes.out")
#'
#' # Run ordmeta on maize PMN pathways
#' ordmeta_res <- magcat_ordmeta_pathways(
#'   gene_results = gene_results,
#'   species = "maize",
#'   gene_col = "GENE",
#'   p_col = "P",
#'   effect_col = "ZSTAT",
#'   is_onetail = FALSE,
#'   B = 5000
#' )
#' head(ordmeta_res)
#' }
#'
#' @seealso \code{\link{ordmeta}}
#' @export
magcat_ordmeta_pathways <- function(gene_results,
                                    pathways     = NULL,
                                    species      = NULL,
                                    pmn_gene_col = NULL,
                                    gene_col     = "GENE",
                                    p_col        = "P",
                                    effect_col   = "ZSTAT",
                                    is_onetail   = FALSE,
                                    min_p        = 1e-15,
                                    do_fix       = TRUE,
                                    B            = 10000L,
                                    output       = FALSE,
                                    out_dir      = "magcat_ordmeta") {

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

  # map from normalized -> canonical gene ID (first occurrence)
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

  # normalize gene IDs in pathways
  p_list <- lapply(p_list, function(g) tolower(as.character(g)))

  ## -------- prepare result container ----------
  n_pw <- length(p_list)
  res  <- data.frame(
    pathway_id            = names(p_list),
    pathway_name          = if (is.null(p_names)) names(p_list) else unname(p_names),
    n_genes               = NA_integer_,
    gene_names            = NA_character_,
    ordmeta_p             = NA_real_,
    optimal_rank          = NA_integer_,
    eff_p_idx             = NA_character_,
    MMP                   = NA_real_,
    overall_eff_direction = NA_character_,
    stringsAsFactors      = FALSE
  )

  ## -------- loop over pathways ----------
  for (i in seq_len(n_pw)) {
    genes_i_norm <- p_list[[i]]
    p_i          <- gene_p_vec[genes_i_norm]

    # drop genes with NA p-values
    keep <- !is.na(p_i)
    p_i  <- p_i[keep]
    genes_used_norm <- genes_i_norm[keep]

    d <- length(p_i)
    res$n_genes[i] <- d

    if (d == 0L) {
      next
    }

    # canonical gene IDs (original case) actually used
    canon_ids <- unname(gene_map[genes_used_norm])
    canon_ids <- unique(canon_ids[!is.na(canon_ids)])
    res$gene_names[i] <- paste(canon_ids, collapse = ";")

    # clean p-values (same logic as ACAT wrapper)
    if (do_fix) {
      p_i <- fix_p_for_acat(p_i, min_p = min_p)
    }

    if (!is_onetail) {
      eff_i <- eff_vec[genes_used_norm]
      sgn   <- sign(eff_i)
      # treat 0 or NA as +1 (neutral-ish, like "no direction")
      sgn[is.na(sgn) | sgn == 0] <- 1
      out <- ordmeta(p = p_i,
                     is.onetail = FALSE,
                     eff.sign   = sgn,
                     B          = B)
    } else {
      out <- ordmeta(p = p_i,
                     is.onetail = TRUE,
                     eff.sign   = NULL,
                     B          = B)
    }

    res$ordmeta_p[i]             <- out$p
    res$optimal_rank[i]          <- out$optimal_rank
    res$MMP[i]                   <- out$MMP
    if (!is_onetail) {
      res$overall_eff_direction[i] <- out$overall.eff.direction
    }

    # eff.p.idx is a vector of indices: save as "1;2;5"
    res$eff_p_idx[i] <- paste(out$eff.p.idx, collapse = ";")
  }

  ## -------- sort by ordmeta_p (ascending) ----------
  ord <- order(res$ordmeta_p, decreasing = FALSE, na.last = TRUE)
  res <- res[ord, , drop = FALSE]

  ## -------- optional CSV output ----------
  if (output) {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    species_tag <- if (is.null(species)) "custom" else species
    out_path <- file.path(
      out_dir,
      paste0("magcat_ordmeta_pathways_", species_tag, ".csv")
    )
    utils::write.csv(res, out_path, row.names = FALSE)
    attr(res, "file") <- out_path
  }

  res
}
