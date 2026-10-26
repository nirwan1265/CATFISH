## ==========================================================================
## Generic pathway loaders
## ==========================================================================
##
## The built-in loaders in pathway_loaders.R are restricted to the five
## bundled databases (maize, sorghum, arabidopsis, plant, fly). These
## functions let CATFISH read pathway definitions for any organism from
## either a BioCyc/Pathway Tools export or a GMT file, returning the same
## long data.frame the component tests expect:
##
##     pathway_id | pathway_name | gene_id
## ==========================================================================

#' Load pathways from a GMT file
#'
#' Reads a Broad-style GMT file (one pathway per line: set name, description,
#' then gene identifiers) and returns it in CATFISH's long format. Works with
#' MSigDB, Reactome, WikiPathways or any other GMT source.
#'
#' @param file Path to a .gmt file.
#' @param min_size,max_size Integer. Keep only pathways whose gene count falls
#'   within these bounds after de-duplication. Defaults keep everything;
#'   10 and 200 are common choices for GWAS pathway analysis.
#' @param drop_unknown Logical. Drop blank and "unknown" gene identifiers.
#'
#' @return data.frame with columns \code{pathway_id}, \code{pathway_name},
#'   \code{gene_id}, one row per pathway-gene pair.
#'
#' @examples
#' \dontrun{
#' pw <- catfish_load_gmt("c2.cp.v2026.2.Hs.entrez.gmt",
#'                        min_size = 10, max_size = 200)
#' }
#' @export
catfish_load_gmt <- function(
    file,
    min_size     = 1L,
    max_size     = Inf,
    drop_unknown = TRUE
) {
  if (!file.exists(file)) {
    stop("GMT file not found: ", file, call. = FALSE)
  }

  lines <- readLines(file, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]

  parts <- strsplit(lines, "\t", fixed = TRUE)

  out <- lapply(parts, function(p) {
    if (length(p) < 3L) return(NULL)
    # GMT column 2 is a description/URL, not a gene
    genes <- unique(trimws(p[-(1:2)]))
    genes <- genes[nzchar(genes)]
    if (isTRUE(drop_unknown)) genes <- genes[genes != "unknown"]
    if (!length(genes)) return(NULL)
    data.frame(
      pathway_id       = p[1],
      pathway_name     = p[1],
      gene_id          = genes,
      stringsAsFactors = FALSE
    )
  })

  df <- do.call(rbind, out)
  if (is.null(df)) {
    stop("No usable pathways parsed from ", basename(file), call. = FALSE)
  }

  .catfish_filter_sizes(df, min_size, max_size)
}


#' Load pathways from a BioCyc / Pathway Tools SmartTable export
#'
#' BioCyc SmartTable CSV exports (HumanCyc, MetaCyc, PlantCyc, ...) are
#' "sparse": the pathway identifier and name appear only on the first row of
#' each block, and subsequent rows carry just one more gene with the leading
#' columns left blank. This reader forward-fills those columns.
#'
#' @param file Path to the exported .csv.
#' @param id_col,name_col Column names holding the pathway object ID and the
#'   pathway name.
#' @param gene_col Column holding the gene identifier to use. For human
#'   analyses with MAGMA this should be an NCBI/Entrez gene ID column, since
#'   MAGMA's gene location files are Entrez-keyed.
#' @param min_size,max_size Pathway size bounds, applied after de-duplication.
#' @param drop_unknown Logical. Drop blank and "unknown" gene identifiers.
#'
#' @return data.frame with columns \code{pathway_id}, \code{pathway_name},
#'   \code{gene_id}.
#'
#' @examples
#' \dontrun{
#' pw <- catfish_load_biocyc("humancyc_pathways.csv",
#'                           gene_col = "NCBI Gene ID")
#' }
#' @export
catfish_load_biocyc <- function(
    file,
    id_col       = "Object ID",
    name_col     = "Pathways",
    gene_col     = "Genes of Pathways",
    min_size     = 1L,
    max_size     = Inf,
    drop_unknown = TRUE
) {
  if (!file.exists(file)) {
    stop("BioCyc export not found: ", file, call. = FALSE)
  }

  x <- utils::read.csv(
    file,
    stringsAsFactors = FALSE,
    check.names      = FALSE,
    colClasses       = "character",
    fileEncoding     = "UTF-8-BOM"
  )

  missing_cols <- setdiff(c(id_col, name_col, gene_col), names(x))
  if (length(missing_cols)) {
    stop(
      "Export is missing column(s): ", paste(missing_cols, collapse = ", "),
      "\nColumns found: ", paste(names(x), collapse = ", "),
      call. = FALSE
    )
  }

  # Forward-fill the sparse pathway columns down each block.
  ff <- function(v) {
    v <- trimws(v)
    v[!nzchar(v)] <- NA_character_
    idx <- cumsum(!is.na(v))
    idx[idx == 0L] <- NA_integer_
    v[!is.na(v)][idx]
  }

  df <- data.frame(
    pathway_id       = ff(x[[id_col]]),
    pathway_name     = ff(x[[name_col]]),
    gene_id          = trimws(x[[gene_col]]),
    stringsAsFactors = FALSE
  )

  df <- df[!is.na(df$pathway_id), , drop = FALSE]

  if (isTRUE(drop_unknown)) {
    bad <- is.na(df$gene_id) | df$gene_id == "" | df$gene_id == "unknown"
    df  <- df[!bad, , drop = FALSE]
  }

  .catfish_filter_sizes(df, min_size, max_size)
}


#' Write a CATFISH pathway data.frame as a MAGMA --set-annot file
#'
#' @param pathways data.frame from \code{catfish_load_gmt} or
#'   \code{catfish_load_biocyc}.
#' @param file Output path. One row per pathway: set name followed by its
#'   gene identifiers, tab separated.
#'
#' @return \code{file}, invisibly.
#' @export
catfish_write_set_annot <- function(pathways, file) {
  stopifnot(all(c("pathway_id", "gene_id") %in% names(pathways)))
  sp <- split(pathways$gene_id, pathways$pathway_id)
  lines <- vapply(
    names(sp),
    function(p) paste(c(p, unique(sp[[p]])), collapse = "\t"),
    character(1)
  )
  writeLines(lines, file)
  invisible(file)
}


## internal -----------------------------------------------------------------
.catfish_filter_sizes <- function(df, min_size, max_size) {
  df <- unique(df)
  n  <- table(df$pathway_id)
  keep <- names(n)[n >= min_size & n <= max_size]
  df <- df[df$pathway_id %in% keep, , drop = FALSE]
  rownames(df) <- NULL
  df
}
