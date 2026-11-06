#' Extract gene-gene correlations from MAGMA .genes.raw files
#'
#' Parses MAGMA's banded correlation output (.genes.raw) and extracts
#' pairwise gene correlations into a simple 3-column format (gene1, gene2, r).
#' This file supplies the gene-gene correlations for Westfall--Young calibration in
#' \code{\link{catfish_omni2_pathways}} and \code{\link{catfish_wy_pathways}}.
#'
#' @param genes_raw_file Path(s) to MAGMA \code{.genes.raw} files. Give one file
#'   per chromosome when MAGMA was run chromosome by chromosome; all are written
#'   to one output file. MAGMA computes gene--gene correlations only within a
#'   chromosome, so the files together hold the full correlation structure.
#' @param out_pairs_file Output path for the correlation pairs file.
#' @param fixed_fields Integer, or \code{NULL} (default) to read it from each
#'   file's \code{# COVAR} header line: the seven fixed MAGMA fields (GENE, CHR,
#'   START, STOP, NSNPS, NPARAM, N) plus one per covariate. For MAGMA 1.10 output
#'   (\code{# COVAR = NSAMP MAC}) this is 9 for every gene model. Used when a
#'   file has no \code{# COVAR} line.
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
#' For per-chromosome MAGMA runs, pass all the chromosome files at once.
#' The function stops if no gene record matches \code{gene_regex}, which
#' usually means the pattern does not fit the gene IDs.
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
#' # All chromosomes of a per-chromosome run, into one file
#' raw <- sprintf("trait_chr%02d.genes.raw", 1:10)
#' magma_genesraw_to_cor_pairs_banded(raw, "trait_gene_cor_pairs.txt",
#'                                    gene_regex = "^SORBI")
#' }
#'
#' @seealso \code{\link{catfish_omni2_pathways}} for using the correlation file
#' @export
magma_genesraw_to_cor_pairs_banded <- function(genes_raw_file,
                                               out_pairs_file,
                                               fixed_fields = NULL,
                                               gene_regex = "^FBgn",
                                               keep_abs_r_ge = 0,
                                               overwrite = TRUE,
                                               verbose = TRUE) {
  genes_raw_file <- as.character(genes_raw_file)
  if (!length(genes_raw_file)) stop("genes_raw_file is empty.", call. = FALSE)
  miss <- genes_raw_file[!file.exists(genes_raw_file)]
  if (length(miss)) stop("genes_raw_file not found: ", paste(miss, collapse = ", "), call. = FALSE)
  if (file.exists(out_pairs_file) && !overwrite) {
    stop("out_pairs_file exists and overwrite=FALSE: ", out_pairs_file, call. = FALSE)
  }

  dir.create(dirname(out_pairs_file), recursive = TRUE, showWarnings = FALSE)
  writeLines("gene1\tgene2\tr", out_pairs_file, useBytes = TRUE)

  fixed_arg  <- fixed_fields
  fixed_from_header <- function(f) {
    hdr <- readLines(f, n = 20L, warn = FALSE)
    cv  <- grep("^#\\s*COVAR\\s*=", hdr, value = TRUE)
    if (!length(cv)) return(NA_integer_)
    v <- trimws(sub("^#\\s*COVAR\\s*=", "", cv[1]))
    7L + if (nzchar(v)) length(strsplit(v, "\\s+")[[1]]) else 0L
  }

  genes_seen <- character()
  last_chr   <- NA_character_

  n_phys  <- 0L
  n_rec   <- 0L
  n_pairs <- 0L

  append_pairs <- function(g1, g2, r) {
    out <- paste(g1, g2, format(r, digits = 10, scientific = TRUE), sep = "\t")
    cat(out, file = out_pairs_file, sep = "\n", append = TRUE)
    invisible(length(r))
  }

  # Accept any non-empty chromosome (Fly has 2L, 2R, 3L, 3R, X, Y, etc.)
  is_new_record <- function(tok) {
    if (length(tok) < fixed_fields) return(FALSE)
    if (!grepl(gene_regex, tok[1])) return(FALSE)
    if (!nzchar(tok[2])) return(FALSE)  # just check non-empty, not numeric
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

  for (raw_f in genes_raw_file) {
  fixed_fields <- if (is.null(fixed_arg)) fixed_from_header(raw_f) else as.integer(fixed_arg)
  if (is.na(fixed_fields)) fixed_fields <- 9L
  n_rec_before <- n_rec
  genes_seen <- character()
  last_chr   <- NA_character_
  con <- file(raw_f, open = "r")
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
      message(sprintf("...%s: %d physical lines; %d gene records; %d pairs",
                      basename(raw_f), n_phys, n_rec, n_pairs))
    }
  }

  close(con)
  # flush last record of this file
  flush_cur()
  if (n_rec == n_rec_before)
    stop("No gene record in ", raw_f, " matches gene_regex = '", gene_regex,
         "'. Set gene_regex to fit the gene IDs (e.g. '^SORBI', '^Zm', '^FBgn', '^[0-9]+$').",
         call. = FALSE)
  }

  if (verbose) {
    message(sprintf("DONE: %d physical lines, %d gene records, %d pairs -> %s",
                    n_phys, n_rec, n_pairs, out_pairs_file))
  }

  invisible(out_pairs_file)
}
