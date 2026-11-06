#' Read and combine MAGMA gene results
#'
#' Reads one or more MAGMA \code{.genes.out} files and returns them as one
#' data.frame. Use it after a per-chromosome MAGMA run (one file per
#' chromosome) to get the gene table that \code{\link{catfish_wy_pathways}},
#' \code{\link{catfish_omni2_pathways}} and \code{\link{gene_prioritize}} take.
#'
#' Every file must have the same columns. A gene that appears in more than one
#' file is kept once (its first occurrence) and reported.
#'
#' With \code{verbose = TRUE} the median of every gene \eqn{p}-value column is
#' printed. Where most genes carry no signal it should be near 0.5. The
#' Westfall--Young null assumes null gene \eqn{p}-values are uniform, so a
#' column with a clearly smaller median (MAGMA's SNP-wise top model, and the
#' multi model that includes it, are often inflated) will be read as pathway
#' signal; choose the gene model, and \code{p_col} downstream, accordingly.
#'
#' @param files Character vector of \code{.genes.out} paths. Provide this or
#'   \code{dir}.
#' @param dir Directory to search for files matching \code{pattern}.
#' @param pattern Regular expression for the file names when \code{dir} is
#'   used.
#' @param verbose Print the number of genes and chromosomes and the median of
#'   each gene \eqn{p}-value column.
#'
#' @return data.frame of gene results (MAGMA columns, e.g. \code{GENE},
#'   \code{CHR}, \code{START}, \code{STOP}, \code{NSNPS}, \code{N},
#'   \code{ZSTAT}, \code{P} or \code{P_MULTI}, \code{P_SNPWISE_MEAN}, ...),
#'   ordered by chromosome and start position. The files read are attached as
#'   attribute \code{files}.
#'
#' @examples
#' f1 <- tempfile(fileext = ".genes.out"); f2 <- tempfile(fileext = ".genes.out")
#' write.table(data.frame(GENE = c("g1", "g2"), CHR = 1, START = c(1, 50),
#'                        STOP = c(40, 90), NSNPS = 5, NPARAM = 2, N = 300,
#'                        ZSTAT = c(2.1, -0.3), P = c(0.018, 0.62)),
#'             f1, sep = "\t", quote = FALSE, row.names = FALSE)
#' write.table(data.frame(GENE = "g3", CHR = 2, START = 10, STOP = 70,
#'                        NSNPS = 8, NPARAM = 3, N = 300, ZSTAT = 0.4, P = 0.34),
#'             f2, sep = "\t", quote = FALSE, row.names = FALSE)
#' genes <- magma_read_genes(c(f1, f2))
#' genes
#'
#' \dontrun{
#' # all chromosomes of a per-chromosome MAGMA run
#' genes <- magma_read_genes(dir = "magma", pattern = "^trait_chr.*\\.genes\\.out$")
#' }
#' @seealso \code{\link{magma_gene}}, \code{\link{magma_genesraw_to_cor_pairs_banded}},
#'   \code{\link{catfish_wy_pathways}}
#' @export
magma_read_genes <- function(files   = NULL,
                             dir     = NULL,
                             pattern = "\\.genes\\.out$",
                             verbose = TRUE) {
  if (is.null(files) == is.null(dir))
    stop("Provide either 'files' or 'dir'.", call. = FALSE)
  if (!is.null(dir)) {
    if (!dir.exists(dir)) stop("dir does not exist: ", dir, call. = FALSE)
    files <- sort(list.files(dir, pattern = pattern, full.names = TRUE))
    if (!length(files))
      stop("No file in ", dir, " matches pattern '", pattern, "'.", call. = FALSE)
  }
  files <- as.character(files)
  miss <- files[!file.exists(files)]
  if (length(miss)) stop("File(s) not found: ", paste(miss, collapse = ", "), call. = FALSE)

  tabs <- lapply(files, function(f)
    utils::read.table(f, header = TRUE, stringsAsFactors = FALSE,
                      comment.char = "#", check.names = FALSE))
  cols <- names(tabs[[1]])
  for (i in seq_along(tabs)[-1]) {
    if (!identical(names(tabs[[i]]), cols))
      stop("Columns differ between files:\n  ", basename(files[1]), ": ",
           paste(cols, collapse = " "), "\n  ", basename(files[i]), ": ",
           paste(names(tabs[[i]]), collapse = " "), call. = FALSE)
  }
  g <- do.call(rbind, tabs)
  if (!"GENE" %in% names(g)) stop("No GENE column in the MAGMA output.", call. = FALSE)

  dup <- duplicated(g$GENE)
  if (any(dup)) {
    message(sprintf("magma_read_genes: %d gene(s) appear in more than one file; first occurrence kept (e.g. %s).",
                    sum(dup), paste(utils::head(unique(g$GENE[dup]), 3), collapse = ", ")))
    g <- g[!dup, , drop = FALSE]
  }
  if (all(c("CHR", "START") %in% names(g))) {
    chr_num <- suppressWarnings(as.numeric(as.character(g$CHR)))
    o <- if (all(is.finite(chr_num))) order(chr_num, g$START) else order(as.character(g$CHR), g$START)
    g <- g[o, , drop = FALSE]
  }
  rownames(g) <- NULL

  if (isTRUE(verbose)) {
    n_chr <- if ("CHR" %in% names(g)) length(unique(g$CHR)) else NA
    message(sprintf("magma_read_genes: %d genes on %s chromosome(s) from %d file(s)",
                    nrow(g), format(n_chr), length(files)))
    pc <- intersect(c("P", "P_MULTI", "P_SNPWISE_MEAN", "P_SNPWISE_TOP1"), names(g))
    for (cc in pc)
      message(sprintf("  median %-15s %.3f", cc, stats::median(as.numeric(g[[cc]]), na.rm = TRUE)))
    if (length(pc))
      message("  (about 0.5 is expected when most genes are null)")
  }
  attr(g, "files") <- files
  g
}
