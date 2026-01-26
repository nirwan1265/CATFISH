## ---------- PMN pathway helpers ------------------------------------

#' Path to built-in PMN pathway file
#'
#' @param species Which built-in DB to use. One of:
#' \itemize{
#'   \item{"maize"}{"corncyc_pathways.20230103"}
#'   \item{"sorghum"}{"sorghumbicolorcyc_pathways.20230103"}
#'   \item{"arabidopsis"}{"aracyc_pathways.20230103"}
#'   \item{"plant"}{"plantcyc_pathways.20230103"}
#' }
#'
#' @return A file path to the selected PMN pathway file in the MAGCAT package.
#' @keywords internal
magcat_pmn_file <- function(
  species = c("maize", "sorghum", "arabidopsis", "plant")
) {
  species <- match.arg(species)

  fname <- switch(
    species,
    maize       = "corncyc_pathways.20230103",
    sorghum     = "sorghumbicolorcyc_pathways.20230103",
    arabidopsis = "aracyc_pathways.20230103",
    plant       = "plantcyc_pathways.20230103"
  )

  path <- system.file("extdata/pathway", fname, package = "MAGCAT")
  if (identical(path, "")) {
    stop("Could not find ", fname, " in inst/extdata/pathway of MAGCAT.",
         call. = FALSE)
  }
  path
}

#' Load PMN pathways as a long data.frame
#'
#' Reads a Plant Metabolic Network (PMN) pathway file from
#' `inst/extdata/pathway` and returns a long table of (pathway, gene) membership.
#'
#' @param species One of "maize", "sorghum", "arabidopsis", "plant".
#' @param gene_col Which column from the PMN file to use for gene IDs.
#'   If NULL (default), the function will:
#'   \itemize{
#'     \item use "Gene-name" if present;
#'     \item otherwise fall back to "Gene-id".
#'   }
#' @param drop_unknown If TRUE, drop rows with missing/unknown gene IDs.
#'
#' @return A data.frame with columns:
#' \describe{
#'   \item{pathway_id}{PMN pathway ID (e.g., "PWY-7634")}
#'   \item{pathway_name}{Human-readable pathway name}
#'   \item{gene_id}{Gene identifier from `gene_col`}
#' }
#'
#' @examples
#' \dontrun{
#' # Load maize PMN pathways
#' maize_pathways <- magcat_load_pathways(species = "maize")
#' head(maize_pathways)
#'
#' # Load arabidopsis pathways with specific gene column
#' arab_pathways <- magcat_load_pathways(
#'   species = "arabidopsis",
#'   gene_col = "Gene-id"
#' )
#' }
#'
#' @seealso \code{\link{magcat_acat_pathways}}, \code{\link{magcat_fisher_pathways}}
#' @export
magcat_load_pathways <- function(
  species      = c("maize", "sorghum", "arabidopsis", "plant"),
  gene_col     = NULL,
  drop_unknown = TRUE
) {
  species <- match.arg(species)
  fpath <- magcat_pmn_file(species)

  x <- suppressWarnings(
    utils::read.delim(
      fpath,
      header           = TRUE,
      stringsAsFactors = FALSE,
      check.names      = FALSE
    )
  )

  # 1) pick gene_col
  if (is.null(gene_col)) {
    candidates <- c("Gene-name", "Gene-id")
    gene_col <- intersect(candidates, names(x))[1]
    if (is.na(gene_col)) {
      stop(
        "Pathway file ", basename(fpath),
        " does not contain any of: ", paste(candidates, collapse = ", "),
        call. = FALSE
      )
    }
  } else if (!gene_col %in% names(x)) {
    stop(
      "Requested gene_col = '", gene_col,
      "' is not a column in ", basename(fpath), ".",
      call. = FALSE
    )
  }

  needed <- c("Pathway-id", "Pathway-name", gene_col)
  if (!all(needed %in% names(x))) {
    stop("Pathway file ", basename(fpath), " does not contain columns: ",
         paste(needed, collapse = ", "),
         call. = FALSE)
  }

  df <- data.frame(
    pathway_id   = x[["Pathway-id"]],
    pathway_name = x[["Pathway-name"]],
    gene_id      = x[[gene_col]],
    stringsAsFactors = FALSE
  )

  if (isTRUE(drop_unknown)) {
    bad <- is.na(df$gene_id) | df$gene_id == "" | df$gene_id == "unknown"
    df <- df[!bad, , drop = FALSE]
  }

  df
}

## ---------- p-value cleaner -----------------------------------------

#' Clean p-values for ACAT
#'
#' ACAT uses Cauchy transforms that can behave poorly at exactly 0 or 1.
#' This helper removes NA values and caps extreme p-values.
#'
#' @param p Numeric vector of p-values (NA allowed). Names are preserved.
#' @param min_p Lower cap for very small p-values (default 1e-15).
#'
#' @return Cleaned numeric vector (NAs removed, extremes capped).
#' @keywords internal
fix_p_for_acat <- function(p, min_p = 1e-15) {
  p <- p[!is.na(p)]
  d <- length(p)
  if (d == 0L) return(p)

  # lower cap
  p[p <= 0] <- min_p

  # upper cap: exact 1's (or >1) become 1 - 1/d
  p[p >= 1] <- 1 - 1 / d

  p
}

## ---------- internal helpers ----------------------------------------

#' @keywords internal
.magcat_assert_cols <- function(df, cols, df_name = "data.frame") {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0L) {
    stop(df_name, " is missing required column(s): ",
         paste(missing, collapse = ", "),
         call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
.magcat_pathways_to_list <- function(pathways) {
  if (is.data.frame(pathways)) {
    .magcat_assert_cols(pathways, c("pathway_id", "gene_id"), "pathways")

    if (!"pathway_name" %in% names(pathways)) {
      pathways$pathway_name <- pathways$pathway_id
    }

    p_list  <- split(as.character(pathways$gene_id), pathways$pathway_id)

    # pathway_name: first name per id (named character vector)
    p_names <- tapply(pathways$pathway_name, pathways$pathway_id, function(x) x[1])

    list(p_list = p_list, p_names = p_names)
  } else if (is.list(pathways)) {
    p_list <- pathways
    nm <- names(p_list)
    if (is.null(nm)) {
      nm <- paste0("PWY_", seq_along(p_list))
      names(p_list) <- nm
    }
    # For list input, names are used as pathway_name (and pathway_id).
    p_names <- stats::setNames(nm, nm)

    list(p_list = p_list, p_names = p_names)
  } else {
    stop("'pathways' must be either a list or a data.frame.", call. = FALSE)
  }
}

## ---------- Pathway-level ACAT wrapper ------------------------------

#' ACAT pathway test on gene-level p-values
#'
#' Computes a self-contained pathway p-value per gene set by combining the
#' member gene p-values using ACAT (Cauchy combination).
#'
#' Gene IDs are matched case-insensitively: both `gene_results[[gene_col]]` and
#' `pathways` gene IDs are normalized to lower-case before matching. Genes that
#' do not appear in `gene_results` are dropped for that pathway.
#'
#' @param gene_results data.frame with at least a gene ID column and a p-value column.
#'   For MAGMA `.genes.out`, these are commonly `GENE` and `P`.
#' @param pathways Either:
#' \itemize{
#'   \item A named list: each element is a character vector of gene IDs in that pathway, OR
#'   \item A data.frame with columns `pathway_id`, `gene_id` (and optional `pathway_name`).
#' }
#' @param species Optional; one of "maize", "sorghum", "arabidopsis", "plant".
#'   If provided, built-in PMN pathways are loaded via `magcat_load_pathways()`.
#'   Provide either `pathways` OR `species`, but not both.
#' @param pmn_gene_col Optional; passed to `magcat_load_pathways(gene_col=...)`
#'   when `species` is used. If NULL, PMN loader prefers "Gene-name" then "Gene-id".
#' @param gene_col Column name in `gene_results` containing gene IDs (default "GENE").
#' @param p_col Column name in `gene_results` containing gene-level p-values (default "P").
#' @param min_p Lower cap for very small p-values before ACAT (default 1e-15).
#' @param do_fix If TRUE (default), clean/cap p-values with `fix_p_for_acat()`.
#' @param output If TRUE, write a CSV of results to `out_dir`.
#' @param out_dir Directory to write CSV when `output = TRUE` (default "magcat_acat").
#'
#' @return A data.frame with columns:
#' \describe{
#'   \item{pathway_id}{Pathway identifier}
#'   \item{pathway_name}{Pathway name (if available; otherwise equals `pathway_id`)}
#'   \item{n_genes}{Number of genes used in ACAT for this pathway}
#'   \item{gene_names}{Semicolon-separated gene IDs actually used (from `gene_results`)}
#'   \item{acat_p}{ACAT p-value}
#' }
#' Rows are sorted by `acat_p` in increasing order (smallest p first; NA at bottom).
#' If `output = TRUE`, an attribute `"file"` is attached with the CSV path.
#'
#' @examples
#' \dontrun{
#' # Load gene results from MAGMA output
#' gene_results <- read.delim("magma_output.genes.out")
#'
#' # Run ACAT on maize PMN pathways
#' acat_res <- magcat_acat_pathways(
#'   gene_results = gene_results,
#'   species = "maize",
#'   gene_col = "GENE",
#'   p_col = "P"
#' )
#' head(acat_res)
#'
#' # Run with custom pathways
#' my_pathways <- list(
#'   pathway1 = c("gene1", "gene2", "gene3"),
#'   pathway2 = c("gene4", "gene5")
#' )
#' acat_custom <- magcat_acat_pathways(
#'   gene_results = gene_results,
#'   pathways = my_pathways
#' )
#' }
#'
#' @seealso \code{\link{magcat_fisher_pathways}}, \code{\link{magcat_minp_pathways}},
#'   \code{\link{magcat_omni2_pathways}}, \code{\link{magcat_load_pathways}}
#' @export
magcat_acat_pathways <- function(gene_results,
                                 pathways     = NULL,
                                 species      = NULL,
                                 pmn_gene_col = NULL,
                                 gene_col     = "GENE",
                                 p_col        = "P",
                                 min_p        = 1e-15,
                                 do_fix       = TRUE,
                                 output       = FALSE,
                                 out_dir      = "magcat_acat") {
  if (!requireNamespace("ACAT", quietly = TRUE)) {
    stop("Package 'ACAT' is required for magcat_acat_pathways(). ",
         "Please install it with install.packages('ACAT').",
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

  # If user gave species, load PMN pathways
  if (is.null(pathways) && !is.null(species)) {
    pathways <- if (is.null(pmn_gene_col)) {
      magcat_load_pathways(species = species)
    } else {
      magcat_load_pathways(species = species, gene_col = pmn_gene_col)
    }
  }

  ## -------- standardize gene_results ----------
  .magcat_assert_cols(gene_results, c(gene_col, p_col), "gene_results")

  genes_all <- as.character(gene_results[[gene_col]])
  p_all     <- as.numeric(gene_results[[p_col]])

  if (anyNA(p_all) && all(is.na(p_all))) {
    stop("All values in gene_results[['", p_col, "']] are NA.", call. = FALSE)
  }

  # normalize (lowercase) gene IDs for matching
  genes_all_norm <- tolower(genes_all)

  # handle possible duplicate gene IDs robustly: keep the minimum p-value per gene
  ok <- !is.na(genes_all_norm) & genes_all_norm != "" & !is.na(p_all)
  if (!any(ok)) {
    stop("No valid (non-NA) gene IDs and p-values found in gene_results.", call. = FALSE)
  }

  p_min <- tapply(p_all[ok], genes_all_norm[ok], min, na.rm = TRUE)
  gene_p_vec <- p_min

  # map from normalized -> canonical gene ID (first occurrence)
  gene_map <- tapply(genes_all[ok], genes_all_norm[ok], function(x) x[1])

  ## -------- pathways -> named list + names ----------
  pw <- .magcat_pathways_to_list(pathways)
  p_list  <- pw$p_list
  p_names <- pw$p_names

  # normalize pathway gene_ids to lower-case
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
    acat_p       = NA_real_,
    stringsAsFactors = FALSE
  )

  ## -------- loop over pathways ----------
  for (i in seq_along(p_list)) {
    genes_i_norm <- p_list[[i]]

    # vector of p-values for genes present in gene_results
    p_i <- gene_p_vec[genes_i_norm]
    p_i <- p_i[!is.na(p_i)]

    d <- length(p_i)
    res$n_genes[i] <- d

    if (d == 0L) {
      res$acat_p[i]     <- NA_real_
      res$gene_names[i] <- NA_character_
      next
    }

    # canonical gene IDs (original case) for genes actually used
    used_norm <- names(p_i)
    canon_ids <- unname(gene_map[used_norm])
    canon_ids <- unique(canon_ids[!is.na(canon_ids) & canon_ids != ""])
    res$gene_names[i] <- paste(canon_ids, collapse = ";")

    if (isTRUE(do_fix)) {
      p_i <- fix_p_for_acat(p_i, min_p = min_p)
    }

    # ACAT p-value
    res$acat_p[i] <- ACAT::ACAT(Pvals = p_i)
  }

  ## -------- sort by acat_p (smallest first) ----------
  ord <- order(res$acat_p, decreasing = FALSE, na.last = TRUE)
  res <- res[ord, , drop = FALSE]

  ## -------- optional CSV output ----------
  if (isTRUE(output)) {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    species_tag <- if (is.null(species)) "custom" else species
    out_path <- file.path(out_dir, paste0("magcat_acat_pathways_", species_tag, ".csv"))
    utils::write.csv(res, out_path, row.names = FALSE)
    attr(res, "file") <- out_path
  }

  res
}
