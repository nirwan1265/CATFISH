#' Extract gene lengths from a GFF3 file
#'
#' Reads a (possibly huge) GFF3 and returns a table of gene coordinates + length.
#' Uses a proper GFF3 parser (rtracklayer) if available; otherwise uses fread(cmd=...)
#' to strip comment lines for speed and compatibility with older data.table.
#'
#' @param gff3_file Path to a GFF3 file (.gff3 or .gff3.gz).
#' @param feature Which feature type to treat as "gene" (default "gene").
#' @param id_key Attribute key used for gene ID (default "ID").
#' @param keep_chr_prefix If non-NULL, enforce a prefix in CHR (e.g. "chr").
#'   If NULL, leaves seqid as-is.
#' @param output Logical; if TRUE write TSV to output_dir/file_name.
#' @param output_dir Directory for output if output=TRUE.
#' @param file_name Output TSV file name if output=TRUE.
#'
#' @return data.frame with columns: gene_id, chr, start, end, strand, length
#'
#' @examples
#' \dontrun{
#' # Extract gene lengths from a GFF3 file
#' gene_lengths <- get_gene_lengths(
#'   gff3_file = "reference.gff3",
#'   feature = "gene",
#'   id_key = "ID"
#' )
#' head(gene_lengths)
#'
#' # Write output to file
#' gene_lengths <- get_gene_lengths(
#'   gff3_file = "reference.gff3",
#'   output = TRUE,
#'   output_dir = "output",
#'   file_name = "my_gene_lengths.tsv"
#' )
#' }
#'
#' @seealso \code{\link{gff3_to_geneloc}}, \code{\link{magcat_adjust_gene_p}}
#' @export
get_gene_lengths <- function(gff3_file,
                             feature         = "gene",
                             id_key          = "ID",
                             keep_chr_prefix = NULL,
                             output          = FALSE,
                             output_dir      = ".",
                             file_name       = "gene_lengths.tsv") {

  if (!is.character(gff3_file) || length(gff3_file) != 1L) {
    stop("get_gene_lengths(): 'gff3_file' must be a single file path.", call. = FALSE)
  }
  if (!file.exists(gff3_file)) {
    stop("get_gene_lengths(): file not found: ", gff3_file, call. = FALSE)
  }

  # -----------------------------
  # 1) Preferred: true GFF3 reader
  # -----------------------------
  if (requireNamespace("rtracklayer", quietly = TRUE) &&
      requireNamespace("GenomicRanges", quietly = TRUE)) {

    gr <- rtracklayer::import(gff3_file)

    # rtracklayer stores feature type in 'type'
    if (!("type" %in% names(S4Vectors::mcols(gr)))) {
      stop("get_gene_lengths(): rtracklayer import did not yield 'type' metadata.", call. = FALSE)
    }

    gr <- gr[S4Vectors::mcols(gr)$type == feature]
    if (length(gr) == 0L) {
      stop("get_gene_lengths(): no rows with type == '", feature, "'.", call. = FALSE)
    }

    m <- S4Vectors::mcols(gr)

    # gene id key: ID is usually present as m$ID; else try Name; else parse attributes if present
    gene_id <- NULL
    if (id_key %in% names(m)) {
      gene_id <- as.character(m[[id_key]])
    } else if ("ID" %in% names(m)) {
      gene_id <- as.character(m[["ID"]])
    } else if ("Name" %in% names(m)) {
      gene_id <- as.character(m[["Name"]])
    } else if ("attributes" %in% names(m)) {
      # last-ditch: parse attributes string if present
      pat <- paste0("(^|;)", id_key, "=([^;]+)")
      gene_id <- sub(pat, "\\2", as.character(m[["attributes"]]))
      gene_id[gene_id == as.character(m[["attributes"]])] <- NA_character_
      if (all(is.na(gene_id)) && "attributes" %in% names(m)) {
        gene_id2 <- sub("(^|;)Name=([^;]+)", "\\2", as.character(m[["attributes"]]))
        gene_id2[gene_id2 == as.character(m[["attributes"]])] <- NA_character_
        gene_id <- gene_id2
      }
    } else {
      stop(
        "get_gene_lengths(): couldn't find gene IDs in rtracklayer metadata (tried ID/Name/attributes).",
        call. = FALSE
      )
    }

    chr <- as.character(GenomicRanges::seqnames(gr))
    if (!is.null(keep_chr_prefix)) {
      chr <- ifelse(grepl(paste0("^", keep_chr_prefix), chr), chr, paste0(keep_chr_prefix, chr))
    }

    start <- as.integer(GenomicRanges::start(gr))
    end   <- as.integer(GenomicRanges::end(gr))
    strand <- as.character(GenomicRanges::strand(gr))
    strand[strand == "*"] <- "."

    out <- data.frame(
      gene_id = gene_id,
      chr     = chr,
      start   = start,
      end     = end,
      strand  = strand,
      length  = end - start + 1L,
      stringsAsFactors = FALSE
    )

  } else {
    # ------------------------------------------
    # 2) Fast fallback: fread(cmd=...) compatible
    # ------------------------------------------
    use_dt <- requireNamespace("data.table", quietly = TRUE)

    # Build a command that strips comment lines starting with '#'
    # Works for both .gff3 and .gff3.gz on mac/linux.
    strip_cmd <- if (grepl("\\.gz$", gff3_file, ignore.case = TRUE)) {
      paste0("gzip -dc ", shQuote(gff3_file), " | grep -v '^#'")
    } else {
      paste0("grep -v '^#' ", shQuote(gff3_file))
    }

    if (use_dt) {
      gff <- data.table::fread(
        cmd = strip_cmd,
        sep = "\t",
        header = FALSE,
        quote = "",
        data.table = FALSE,
        fill = TRUE,
        showProgress = FALSE
      )
    } else {
      # Base R fallback: read stripped text via pipe
      con <- pipe(strip_cmd, open = "r")
      on.exit(close(con), add = TRUE)
      gff <- utils::read.delim(
        con,
        header = FALSE,
        sep = "\t",
        quote = "",
        stringsAsFactors = FALSE
      )
    }

    if (ncol(gff) < 9L) {
      stop("get_gene_lengths(): expected >=9 columns in GFF3 after parsing.", call. = FALSE)
    }

    gff <- gff[, 1:9]
    colnames(gff) <- c("seqid","source","type","start","end","score","strand","phase","attributes")

    gff <- gff[gff$type == feature, , drop = FALSE]
    if (!nrow(gff)) {
      stop("get_gene_lengths(): no rows with type == '", feature, "'.", call. = FALSE)
    }

    pat <- paste0("(^|;)", id_key, "=([^;]+)")
    gene_id <- sub(pat, "\\2", gff$attributes)
    gene_id[gene_id == gff$attributes] <- NA_character_

    if (all(is.na(gene_id))) {
      gene_id2 <- sub("(^|;)Name=([^;]+)", "\\2", gff$attributes)
      gene_id2[gene_id2 == gff$attributes] <- NA_character_
      gene_id <- gene_id2
    }

    start <- suppressWarnings(as.integer(gff$start))
    end   <- suppressWarnings(as.integer(gff$end))

    swap <- which(!is.na(start) & !is.na(end) & start > end)
    if (length(swap)) {
      tmp <- start[swap]; start[swap] <- end[swap]; end[swap] <- tmp
    }

    chr <- as.character(gff$seqid)
    if (!is.null(keep_chr_prefix)) {
      chr <- ifelse(grepl(paste0("^", keep_chr_prefix), chr), chr, paste0(keep_chr_prefix, chr))
    }

    out <- data.frame(
      gene_id = as.character(gene_id),
      chr     = chr,
      start   = start,
      end     = end,
      strand  = as.character(gff$strand),
      length  = ifelse(is.na(start) | is.na(end), NA_integer_, end - start + 1L),
      stringsAsFactors = FALSE
    )
  }

  # Common cleanup (both paths)
  out <- out[!is.na(out$gene_id) & nzchar(out$gene_id) & !is.na(out$start) & !is.na(out$end), , drop = FALSE]

  if (anyDuplicated(out$gene_id)) {
    o <- order(out$gene_id, -out$length)
    out <- out[o, , drop = FALSE]
    out <- out[!duplicated(out$gene_id), , drop = FALSE]
  }

  out <- out[order(out$chr, out$start), , drop = FALSE]

  if (output) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    out_path <- file.path(output_dir, file_name)
    utils::write.table(
      out,
      file = out_path,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE
    )
    attr(out, "file") <- out_path
  }

  out
}
