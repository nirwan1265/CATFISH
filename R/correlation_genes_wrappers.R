#' Extract gene-gene correlations from MAGMA .genes.raw files
#'
#' Parses MAGMA's banded correlation output (.genes.raw) and extracts
#' pairwise gene correlations into a simple 3-column format (gene1, gene2, r).
#' This file is required for MVN-based resampling in \code{\link{magcat_omni2_pathways}}.
#'
#' @param genes_raw_file Path to MAGMA .genes.raw file (one chromosome).
#' @param out_pairs_file Output path for the correlation pairs file.
#' @param fixed_fields Integer. Number of fixed fields before correlation values
#'   in MAGMA output (default 9 for multi=snp-wise model).
#' @param gene_regex Regular expression to identify gene IDs. Use organism-specific
#'   patterns:
#'   \itemize{
#'     \item Maize: \code{"^Zm"}
#'     \item Fly: \code{"^FBgn"}
#'     \item Arabidopsis: \code{"^AT"}
#'     \item Sorghum: \code{"^SORBI"}
#'   }
#' @param keep_abs_r_ge Numeric. Only keep correlations with |r| >= this value.
#'   Default 0 keeps all correlations.
#' @param overwrite Logical. Overwrite output file if it exists.
#' @param verbose Logical. Print progress messages.
#'
#' @return Invisibly returns the path to the output file.
#'
#' @details
#' MAGMA stores gene-gene correlations in a banded format within .genes.raw files.
#' This function parses that format and outputs a simple tab-separated file with
#' columns: gene1, gene2, r.
#'
#' For multi-chromosome analyses, run this function on each chromosome's .raw file
#' and concatenate the results (dropping headers after the first file).
#'
#' @examples
#' \dontrun{
#' # Process a single chromosome
#' magma_genesraw_to_cor_pairs_banded(
#'   genes_raw_file = "magma_output_chr1.genes.raw",
#'   out_pairs_file = "gene_cor_chr1.txt",
#'   gene_regex = "^Zm"  # Maize genes
#' )
#'
#' # Combine multiple chromosomes - see usage2.R for full example
#' }
#'
#' @seealso \code{\link{magcat_omni2_pathways}} for using the correlation file
#' @export
magma_genesraw_to_cor_pairs_banded <- function(genes_raw_file,
                                               out_pairs_file,
                                               fixed_fields = 9L,
                                               gene_regex = "^FBgn",
                                               keep_abs_r_ge = 0,
                                               overwrite = TRUE,
                                               verbose = TRUE) {
  stopifnot(file.exists(genes_raw_file))
  if (file.exists(out_pairs_file) && !overwrite) {
    stop("out_pairs_file exists and overwrite=FALSE: ", out_pairs_file, call. = FALSE)
  }

  dir.create(dirname(out_pairs_file), recursive = TRUE, showWarnings = FALSE)
  writeLines("gene1\tgene2\tr", out_pairs_file, useBytes = TRUE)

  con <- file(genes_raw_file, open = "r")
  on.exit(close(con), add = TRUE)

  genes_seen <- character()
  last_chr   <- NA_character_

  n_phys  <- 0L
  n_rec   <- 0L
  n_pairs <- 0L

  append_pairs <- function(g1, g2, r) {
    out <- paste(g1, g2, format(r, digits = 10, scientific = TRUE), sep = "\t")
    cat(out, sep = "\n", file = out_pairs_file, append = TRUE)
    cat("\n", file = out_pairs_file, append = TRUE)
    invisible(length(r))
  }

  is_int_chr <- function(x) grepl("^[0-9]+$", x)
  is_new_record <- function(tok) {
    if (length(tok) < fixed_fields) return(FALSE)
    if (!grepl(gene_regex, tok[1])) return(FALSE)
    if (!is_int_chr(tok[2])) return(FALSE)
    TRUE
  }

  cur_gene <- NULL
  cur_chr  <- NULL
  cur_corr <- character()

  flush_cur <- function() {
    if (is.null(cur_gene)) return(invisible(NULL))

    K <- length(cur_corr)

    # Map to *last K* previously-seen genes on this chromosome block
    if (K > 0L) {
      if (K > length(genes_seen)) {
        stop("Banded parse error: K correlations but only ", length(genes_seen),
             " prior genes seen. gene=", cur_gene, " chr=", cur_chr,
             " K=", K, call. = FALSE)
      }

      g2 <- tail(genes_seen, K)
      r  <- suppressWarnings(as.numeric(cur_corr))

      ok <- is.finite(r) & !is.na(r)
      if (keep_abs_r_ge > 0) ok <- ok & (abs(r) >= keep_abs_r_ge)

      if (any(ok)) {
        n_pairs <<- n_pairs + append_pairs(cur_gene, g2[ok], r[ok])
      }
    }

    genes_seen <<- c(genes_seen, cur_gene)
    n_rec <<- n_rec + 1L

    cur_gene <<- NULL
    cur_chr  <<- NULL
    cur_corr <<- character()
    invisible(NULL)
  }

  repeat {
    lines <- readLines(con, n = 5000L, warn = FALSE)
    if (!length(lines)) break

    for (ln in lines) {
      if (!nzchar(trimws(ln))) next
      if (grepl("^\\s*#", ln)) next

      n_phys <- n_phys + 1L
      tok <- strsplit(trimws(ln), "\\s+")[[1]]
      if (!length(tok)) next

      if (is_new_record(tok)) {
        # new gene record begins => flush previous record
        flush_cur()

        cur_gene <- tok[1]
        cur_chr  <- tok[2]

        # reset band history when chr changes
        if (is.na(last_chr) || cur_chr != last_chr) {
          genes_seen <- character()
          last_chr <- cur_chr
        }

        # trailing correlations after fixed fields
        cur_corr <- if (length(tok) > fixed_fields) tok[(fixed_fields + 1L):length(tok)] else character()

      } else {
        # continuation line: just more correlations
        if (is.null(cur_gene)) next
        cur_corr <- c(cur_corr, tok)
      }
    }

    if (verbose) {
      message(sprintf("...%d physical lines; %d gene records; %d pairs", n_phys, n_rec, n_pairs))
    }
  }

  # flush last record
  flush_cur()

  if (verbose) {
    message(sprintf("DONE: %d physical lines, %d gene records, %d pairs -> %s",
                    n_phys, n_rec, n_pairs, out_pairs_file))
  }

  invisible(out_pairs_file)
}
