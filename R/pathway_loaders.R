## ==========================================================================
## Pathway Loading Utilities
## ==========================================================================
##
## Functions for loading built-in pathway databases (PMN for plants, FlyCyc
## for Drosophila). These are used by all pathway testing wrappers.
##
## Data files are stored in: inst/extdata/pathway/
## ==========================================================================

#' Path to built-in pathway file
#'
#' Returns the file path to a built-in pathway database bundled with MAGCAT.
#'
#' @param species Which built-in database to use:
#' \describe{
#'   \item{"maize"}{CornCyc pathways (corncyc_pathways.20230103)}
#'   \item{"sorghum"}{SorghumBicolorCyc pathways (sorghumbicolorcyc_pathways.20230103)}
#'   \item{"arabidopsis"}{AraCyc pathways (aracyc_pathways.20230103)}
#'   \item{"plant"}{PlantCyc general plant pathways (plantcyc_pathways.20230103)}
#'   \item{"fly"}{FlyCyc Drosophila pathways (Fly_Cyc.tsv)}
#' }
#'
#' @return Character. File path to the selected pathway file in the MAGCAT
#'   package installation.
#'
#' @examples
#' \dontrun{
#' # Get path to maize pathway file
#' maize_file <- magcat_pmn_file("maize")
#'
#' # Get path to fly pathway file
#' fly_file <- magcat_pmn_file("fly")
#' }
#'
#' @seealso \code{\link{magcat_load_pathways}} for loading pathways as a
#'   data.frame
#' @export
magcat_pmn_file <- function(
    species = c("maize", "sorghum", "arabidopsis", "plant", "fly")
) {

  species <- match.arg(species)

  fname <- switch(
    species,
    maize       = "corncyc_pathways.20230103",
    sorghum     = "sorghumbicolorcyc_pathways.20230103",
    arabidopsis = "aracyc_pathways.20230103",
    plant       = "plantcyc_pathways.20230103",
    fly         = "Fly_Cyc.tsv"
  )

  path <- system.file("extdata/pathway", fname, package = "MAGCAT")
  if (identical(path, "")) {
    stop(
      "Could not find ", fname, " in inst/extdata/pathway of MAGCAT.",
      call. = FALSE
    )
  }
  path
}


#' Load pathways as a data.frame
#'
#' Reads a pathway file from the MAGCAT package and returns a long-format
#' data.frame of pathway-gene memberships. Supports Plant Metabolic Network
#' (PMN) pathway files for plants and FlyCyc for Drosophila.
#'
#' @param species Character. One of "maize", "sorghum", "arabidopsis",
#'   "plant", or "fly".
#' @param gene_col Character or NULL. Column name to use for gene IDs.
#'   If NULL (default), the function will:
#'   \itemize{
#'     \item For PMN files: use "Gene-name" if present, else "Gene-id"
#'     \item For fly: use "gene"
#'   }
#' @param drop_unknown Logical. If TRUE (default), drop rows where gene_id
#'   is NA, empty, or "unknown".
#'
#' @return A data.frame with columns:
#' \describe{
#'   \item{pathway_id}{Pathway identifier (e.g., "PWY-7634")}
#'   \item{pathway_name}{Human-readable pathway name}
#'   \item{gene_id}{Gene identifier from the selected gene column}
#' }
#'
#' @details
#' Pathway data files are stored in \code{inst/extdata/pathway/} within
#' the MAGCAT package. The supported species and their source databases are:
#' \itemize{
#'   \item \strong{maize}: CornCyc from Plant Metabolic Network
#'   \item \strong{sorghum}: SorghumBicolorCyc from PMN
#'   \item \strong{arabidopsis}: AraCyc from PMN
#'   \item \strong{plant}: PlantCyc (general plant pathways) from PMN
#'   \item \strong{fly}: FlyCyc (Drosophila melanogaster pathways)
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
#'
#' # Load fly pathways (FlyCyc)
#' fly_pathways <- magcat_load_pathways(species = "fly")
#' head(fly_pathways)
#'
#' # Check how many pathways and genes
#' length(unique(maize_pathways$pathway_id))
#' length(unique(maize_pathways$gene_id))
#' }
#'
#' @seealso
#' \code{\link{magcat_pmn_file}} for getting the raw file path
#'
#' \code{\link{magcat_acat_pathways}}, \code{\link{magcat_fisher_pathways}},
#' \code{\link{magcat_minp_pathways}} for pathway enrichment tests
#' @export
magcat_load_pathways <- function(
    species      = c("maize", "sorghum", "arabidopsis", "plant", "fly"),
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

  # Handle fly format (different column names)
  if (species == "fly") {
    # Fly_Cyc.tsv has: pathway_name, pathway_id, gene
    if (is.null(gene_col)) gene_col <- "gene"
    if (!gene_col %in% names(x)) {
      stop("Fly pathway file does not contain column: ", gene_col,
           call. = FALSE)
    }
    df <- data.frame(
      pathway_id   = x[["pathway_id"]],
      pathway_name = x[["pathway_name"]],
      gene_id      = x[[gene_col]],
      stringsAsFactors = FALSE
    )
  } else {
    # PMN format (maize, sorghum, arabidopsis, plant)
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
      stop(
        "Pathway file ", basename(fpath), " does not contain columns: ",
        paste(needed, collapse = ", "),
        call. = FALSE
      )
    }

    df <- data.frame(
      pathway_id   = x[["Pathway-id"]],
      pathway_name = x[["Pathway-name"]],
      gene_id      = x[[gene_col]],
      stringsAsFactors = FALSE
    )
  }

  if (isTRUE(drop_unknown)) {
    bad <- is.na(df$gene_id) | df$gene_id == "" | df$gene_id == "unknown"
    df <- df[!bad, , drop = FALSE]
  }

  df
}
