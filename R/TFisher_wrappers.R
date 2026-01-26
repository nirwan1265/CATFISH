#' Internal: truncated Fisher (TPM) statistic + approximate p-value
#'
#' Uses:
#'   stat = -2 * sum(log p_i), over p_i < ptrunc
#'   df   = 2 * nincl
#' If no p-values pass truncation (nincl = 0), returns stat = 0, pval = 1.
#'
#' @param p numeric vector of p-values.
#' @param ptrunc scalar in (0, 1]; p-values < ptrunc are included.
#'
#' @return list(stat, nincl, nexcl, validp, ptrunc, pval)
magcat_tpm_stat <- function(p, ptrunc = 0.5) {
  p <- as.numeric(p)

  # keep valid
  keep   <- (p > 0) & (p <= 1)
  validp <- p[keep]

  if (length(validp) < 1L) {
    warning("magcat_tpm_stat(): no valid p values")
    return(list(
      stat   = NA_real_,
      nincl  = 0L,
      nexcl  = 0L,
      validp = validp,
      ptrunc = ptrunc,
      pval   = NA_real_
    ))
  }

  # sanity for ptrunc
  if ((ptrunc <= 0) || (ptrunc > 1)) {
    warning("The value of ptrunc must be between 0 and 1; setting to 0.5")
    ptrunc <- 0.5
  }

  nincl <- sum(validp < ptrunc)
  nexcl <- sum(validp >= ptrunc)

  # no p-values below ptrunc → define “no evidence”: p = 1
  if (nincl < 1L) {
    warning("No p values left after truncation")
    return(list(
      stat   = 0,
      nincl  = 0L,
      nexcl  = length(validp),
      validp = validp,
      ptrunc = ptrunc,
      pval   = 1
    ))
  }

  p_use <- validp[validp < ptrunc]
  stat  <- -2 * sum(log(p_use))

  df   <- 2L * nincl
  pval <- stats::pchisq(stat, df = df, lower.tail = FALSE)

  list(
    stat   = stat,
    nincl  = nincl,
    nexcl  = nexcl,
    validp = validp,
    ptrunc = ptrunc,
    pval   = pval
  )
}




#' Truncated Fisher (TPM) pathway test with gene-set permutations
#'
#' For each pathway, this function:
#' \itemize{
#'   \item takes gene-level p-values (e.g. from MAGMA .genes.out),
#'   \item computes the truncated Fisher (TPM) statistic with threshold \code{ptrunc}:
#'         \deqn{T = -2 \sum_{p_i < p_{trunc}} \log(p_i)}
#'   \item (optionally) reports an approximate analytic p-value via a chi-square
#'         distribution with \code{df = 2 * nincl},
#'   \item calibrates the test using gene-set permutations to obtain an
#'         empirical TPM p-value.
#' }
#'
#' This implementation does NOT depend on TFisher or metap; everything is done
#' in base R. The permutation p-value is the one you should trust for inference.
#'
#' @param gene_results data.frame with at least gene + p-value columns.
#'   Typically a MAGMA \code{.genes.out} file with "GENE" and "P".
#' @param pathways either:
#'   \itemize{
#'     \item a named list: each element is a character vector of gene IDs
#'           in that pathway, OR
#'     \item a data.frame with columns \code{pathway_id}, \code{gene_id}
#'           (and optional \code{pathway_name}).
#'   }
#' @param species optional; one of "maize", "sorghum", "arabidopsis", "plant".
#'   If provided, built-in PMN pathways are loaded via \code{magcat_load_pathways()}.
#'   You must provide either \code{pathways} OR \code{species}, but not both.
#' @param pmn_gene_col optional; passed to \code{magcat_load_pathways(gene_col=...)}
#'   when \code{species} is used.
#' @param gene_col column name in \code{gene_results} containing gene IDs
#'   (default "GENE").
#' @param p_col column name in \code{gene_results} containing gene-level p-values
#'   (default "P").
#' @param ptrunc truncation threshold; p-values \eqn{p_i < p_{trunc}} are included
#'   in the product (default 0.05).
#' @param min_p lower cap for very small p-values (default 1e-15).
#'   Values \code{<= 0} are set to \code{min_p}; values \code{>= 1} are set to
#'   \code{1 - 1/d}, where \code{d} is the number of p-values in that pathway.
#' @param do_fix logical; if TRUE, clean p-values with \code{fix_p_for_acat()}
#'   before computing the stat.
#' @param B_perm Integer; number of gene-set permutations used to obtain
#'   empirical TPM p-values (default 1000L). If 0, only the statistic
#'   (and optional approximate analytic p) are returned.
#' @param seed optional integer random seed for permutations (ignored if
#'   \code{B_perm == 0}).
#' @param analytic_logical logical; if TRUE, also compute \code{tpm_p_analytic}
#'   using a chi-square approximation under independence
#'   (\code{df = 2 * nincl}). This is an approximation, mainly for reference.
#' @param output logical; if TRUE, write a CSV of results to \code{out_dir}.
#' @param out_dir directory to write CSV when \code{output = TRUE}
#'   (default "magcat_tpm").
#'
#' @return data.frame with columns:
#'   \describe{
#'     \item{pathway_id}{pathway identifier}
#'     \item{pathway_name}{pathway name}
#'     \item{n_genes}{number of genes used (with non-NA p-values)}
#'     \item{gene_names}{semicolon-separated gene IDs}
#'     \item{tpm_stat}{truncated Fisher statistic for that pathway}
#'     \item{nincl}{number of p-values included after truncation}
#'     \item{tpm_p_perm}{empirical p-value from gene-set permutations
#'                       (NA if \code{B_perm == 0})}
#'     \item{tpm_p_analytic}{approximate analytic p-value via
#'                           chi-square(df = 2*nincl);
#'                           equals 1 when \code{nincl == 0} and
#'                           NA only if there are no valid p's.}
#'   }
#'   If \code{B_perm > 0}, rows are sorted by \code{tpm_p_perm} ascending;
#'   otherwise by \code{tpm_stat} descending (larger stat = stronger signal).
#'   If \code{output = TRUE}, an attribute \code{"file"} is attached with the CSV path.
#'
#' @examples
#' \dontrun{
#' # Load gene results from MAGMA output
#' gene_results <- read.delim("magma_output.genes.out")
#'
#' # Run truncated Fisher (TPM) on maize PMN pathways
#' tfisher_res <- magcat_tfisher_pathways(
#'   gene_results = gene_results,
#'   species = "maize",
#'   gene_col = "GENE",
#'   p_col = "P",
#'   ptrunc = 0.05,
#'   B_perm = 1000
#' )
#' head(tfisher_res)
#' }
#'
#' @seealso \code{\link{magcat_soft_tfisher_pathways}},
#'   \code{\link{magcat_soft_tfisher_adaptive_pathways}}
#' @export
magcat_tfisher_pathways <- function(gene_results,
                                    pathways         = NULL,
                                    species          = NULL,
                                    pmn_gene_col     = NULL,
                                    gene_col         = "GENE",
                                    p_col            = "P",
                                    ptrunc           = 0.05,
                                    min_p            = 1e-15,
                                    do_fix           = TRUE,
                                    B_perm           = 1000L,
                                    seed             = NULL,
                                    analytic_logical = TRUE,
                                    output           = FALSE,
                                    out_dir          = "magcat_tpm") {

  ## --- pathway source: pathways vs species ---
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

  ## --- standardize gene_results ---
  needed_cols <- c(gene_col, p_col)
  if (!all(needed_cols %in% names(gene_results))) {
    stop("gene_results must contain columns '", gene_col,
         "' and '", p_col, "'.",
         call. = FALSE)
  }

  gr        <- gene_results
  genes_all <- as.character(gr[[gene_col]])
  p_all     <- gr[[p_col]]

  genes_all_norm <- tolower(genes_all)

  ## --- pathways -> named list + names ---
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

  ## --- result container ---
  n_pw <- length(p_list)
  res  <- data.frame(
    pathway_id      = names(p_list),
    pathway_name    = if (is.null(p_names)) names(p_list) else unname(p_names),
    n_genes         = NA_integer_,
    gene_names      = NA_character_,
    tpm_stat        = NA_real_,
    nincl           = NA_integer_,
    tpm_p_perm      = NA_real_,
    tpm_p_analytic  = NA_real_,
    stringsAsFactors = FALSE
  )

  ## --- index pool for permutations ---
  idx_pool <- which(!is.na(p_all))
  if (length(idx_pool) == 0L) {
    stop("No non-NA gene-level p-values found in 'gene_results'.",
         call. = FALSE)
  }

  if (B_perm > 0L && !is.null(seed)) {
    set.seed(seed)
  }

  ## --- main loop over pathways ---
  for (i in seq_len(n_pw)) {
    genes_i_norm <- p_list[[i]]

    # match pathway genes to global gene indices
    idx_i <- match(genes_i_norm, genes_all_norm)
    idx_i <- idx_i[!is.na(idx_i)]

    d <- length(idx_i)
    res$n_genes[i] <- d
    if (d == 0L) {
      next
    }

    # canonical gene IDs for reporting
    canon_ids <- unique(genes_all[idx_i])
    res$gene_names[i] <- paste(canon_ids, collapse = ";")

    # observed p-values
    p_i <- p_all[idx_i]

    # drop NA p's within pathway
    keep  <- !is.na(p_i)
    p_i   <- p_i[keep]
    idx_i <- idx_i[keep]
    d     <- length(p_i)
    res$n_genes[i] <- d

    if (d == 0L) {
      next
    }

    # clean p-values if requested
    if (do_fix) {
      p_i <- fix_p_for_acat(p_i, min_p = min_p)
    }

    # TPM stat (truncated Fisher) + analytic p
    tpm_out <- magcat_tpm_stat(p_i, ptrunc = ptrunc)
    tpm_stat <- tpm_out$stat
    nincl    <- tpm_out$nincl

    res$tpm_stat[i] <- tpm_stat
    res$nincl[i]    <- nincl

    # analytic p: use magcat_tpm_stat pval directly
    if (analytic_logical) {
      res$tpm_p_analytic[i] <- tpm_out$pval
    }

    # permutation p-value (gene-set null)
    if (B_perm > 0L) {
      if (is.na(tpm_stat) || nincl == 0L) {
        # no signal / no truncated p's → define perm p = 1
        res$tpm_p_perm[i] <- 1
      } else {
        tpm_perm <- numeric(B_perm)

        for (b in seq_len(B_perm)) {
          idx_perm <- sample(idx_pool, size = d, replace = FALSE)
          p_perm   <- p_all[idx_perm]
          p_perm   <- p_perm[!is.na(p_perm)]

          if (!length(p_perm)) {
            tpm_perm[b] <- NA_real_
            next
          }

          if (do_fix) {
            p_perm <- fix_p_for_acat(p_perm, min_p = min_p)
          }

          tpm_perm[b] <- magcat_tpm_stat(p_perm, ptrunc = ptrunc)$stat
        }

        res$tpm_p_perm[i] <-
          (1 + sum(tpm_perm >= tpm_stat, na.rm = TRUE)) / (B_perm + 1)
      }
    }
  }

  ## --- sorting ---
  if (B_perm > 0L) {
    ord <- order(res$tpm_p_perm, decreasing = FALSE, na.last = TRUE)
  } else {
    ord <- order(res$tpm_stat, decreasing = TRUE, na.last = TRUE)
  }
  res <- res[ord, , drop = FALSE]

  ## --- optional CSV output ---
  if (output) {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    species_tag <- if (is.null(species)) "custom" else species
    out_path <- file.path(
      out_dir,
      paste0("magcat_tpm_pathways_", species_tag, ".csv")
    )
    utils::write.csv(res, out_path, row.names = FALSE)
    attr(res, "file") <- out_path
  }

  res
}
