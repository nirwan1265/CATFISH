## =========================================
## Helpers: chr mapping (package-internal)
## =========================================

#' Make a chromosome recode map (original -> "1","2",...)
#'
#' @param chr_vec Character vector of chromosome/contig labels observed in data.
#' @param chr_order Optional character vector specifying desired order of labels;
#'   will be mapped to 1..length(chr_order).
#' @param chr_map Optional named character vector mapping original->new
#'   (e.g. c("2L"="1","2R"="2",...)). Names are original labels.
#' @param strict If TRUE, error on labels not in chr_order/chr_map.
#'   If FALSE, append unseen labels after known ones.
#' @param auto_sort_extra If TRUE, extra labels are appended in sorted order.
#' @return A named character vector: names = original labels, values = new labels.
magcat_make_chr_map <- function(chr_vec,
                                chr_order = NULL,
                                chr_map = NULL,
                                strict = FALSE,
                                auto_sort_extra = TRUE) {
  chr_vec <- as.character(chr_vec)
  chr_vec <- chr_vec[!is.na(chr_vec) & nzchar(chr_vec)]
  all_chr <- unique(chr_vec)

  # If user gives an explicit map, trust it (and optionally append extras)
  if (!is.null(chr_map)) {
    if (is.null(names(chr_map)) || any(!nzchar(names(chr_map)))) {
      stop("chr_map must be a *named* character vector: names=original, values=new.",
           call. = FALSE)
    }
    chr_map <- as.character(chr_map)
    names(chr_map) <- as.character(names(chr_map))

    extra <- setdiff(all_chr, names(chr_map))
    if (length(extra) && strict) {
      stop("Found CHR labels not in chr_map: ", paste(extra, collapse = ", "),
           call. = FALSE)
    }
    if (length(extra)) {
      if (auto_sort_extra) extra <- sort(extra)
      start <- length(chr_map) + 1L
      chr_map <- c(chr_map, setNames(as.character(seq.int(start, length.out = length(extra))), extra))
    }
    return(chr_map)
  }

  # If user gives an order, map it to 1..K (and optionally append extras)
  if (!is.null(chr_order)) {
    chr_order <- as.character(chr_order)
    chr_order <- chr_order[!is.na(chr_order) & nzchar(chr_order)]
    if (!length(chr_order)) stop("chr_order provided but empty.", call. = FALSE)

    m <- setNames(as.character(seq_along(chr_order)), chr_order)

    extra <- setdiff(all_chr, names(m))
    if (length(extra) && strict) {
      stop("Found CHR labels not in chr_order: ", paste(extra, collapse = ", "),
           call. = FALSE)
    }
    if (length(extra)) {
      if (auto_sort_extra) extra <- sort(extra)
      start <- length(m) + 1L
      m <- c(m, setNames(as.character(seq.int(start, length.out = length(extra))), extra))
    }
    return(m)
  }

  # AUTO mode: build a stable mapping from observed labels:
  # - numeric-looking first (1,2,3,...)
  # - then common sex/mt chromosomes (X,Y,M/MT/mito)
  # - then everything else alphabetically
  is_int <- function(x) !is.na(suppressWarnings(as.integer(x)))
  num_chr <- all_chr[is_int(all_chr)]
  oth_chr <- setdiff(all_chr, num_chr)

  num_chr <- num_chr[order(as.integer(num_chr))]

  # Put common special chroms in a predictable order if present
  special_order <- c("X","Y","W","Z","MT","M","mitochondrion","mitochondrion_genome","chrM","chrMT")
  specials <- oth_chr[match(toupper(oth_chr), toupper(special_order), nomatch = 0) > 0]
  rest     <- setdiff(oth_chr, specials)

  specials <- specials[order(match(toupper(specials), toupper(special_order)))]
  rest     <- sort(rest)

  ord <- c(num_chr, specials, rest)
  setNames(as.character(seq_along(ord)), ord)
}

#' Apply a chromosome recode map
#'
#' @param chr_vec vector of chromosome labels
#' @param chr_map named character vector original->new
#' @param strict if TRUE error on unknown labels; if FALSE keep unknown as-is
#' @return character vector
magcat_apply_chr_map <- function(chr_vec, chr_map, strict = TRUE) {
  x <- as.character(chr_vec)
  y <- unname(chr_map[x])

  if (anyNA(y)) {
    bad <- unique(x[is.na(y)])
    if (strict) {
      stop("CHR recode map missing entries for: ", paste(bad, collapse = ", "),
           call. = FALSE)
    } else {
      # keep originals for unknowns
      y[is.na(y)] <- x[is.na(y)]
    }
  }
  y
}


## =========================================
## Main: GFF3 -> MAGMA gene-loc (numeric CHR)
## =========================================

#' Convert a GFF3 file to a MAGMA gene location file (optionally recoding CHR to 1..K)
#'
#' @param gff Path to a GFF3 file.
#' @param out Path to write the MAGMA gene location file (tab-delimited).
#' @param feature_type GFF3 feature type(s) to use (default "gene").
#' @param id_fields Candidate attribute tags to use as gene IDs
#'   (searched in order), e.g. c("gene_id","ID","Name").
#' @param drop_scaffolds Logical; if TRUE, drop seqids that look like scaffolds.
#' @param scaffold_regex Regex for scaffold/unplaced contigs to drop when drop_scaffolds=TRUE.
#' @param chr_prefix Optional string to strip from start of seqid when reading (e.g. "chr").
#' @param recode_chr One of "none","auto","order","map".
#' @param chr_order If recode_chr="order", the desired order of original chromosome labels.
#' @param chr_map If recode_chr="map", a named vector mapping original->new numeric labels.
#' @param strict_chr If TRUE, error when encountering CHR labels not in chr_order/chr_map.
#' @param write_chr_map If TRUE, write a mapping file.
#' @param chr_map_out Optional explicit path for the mapping file. Default: paste0(out, ".chr_map.tsv")
#'
#' @return Invisibly returns a list with:
#'   - gene_loc_df: data.frame with GENE, CHR_ORIG, CHR, START, STOP
#'   - chr_map: named character vector original->new (or NULL if no recode)
#'   - chr_map_out: path (or NULL)
#'
#' @examples
#' \dontrun{
#' # Convert GFF3 to MAGMA gene location file
#' result <- gff3_to_geneloc(
#'   gff = "reference.gff3",
#'   out = "output/genes.loc",
#'   feature_type = "gene",
#'   id_fields = c("gene_id", "ID", "Name")
#' )
#'
#' # With automatic chromosome recoding for non-numeric chromosomes
#' result <- gff3_to_geneloc(
#'   gff = "drosophila.gff3",
#'   out = "output/fly_genes.loc",
#'   recode_chr = "auto",
#'   write_chr_map = TRUE
#' )
#' }
#'
#' @seealso \code{\link{magma_gene}}, \code{\link{get_gene_lengths}}
#' @export
gff3_to_geneloc <- function(gff,
                            out,
                            feature_type   = "gene",
                            id_fields      = c("gene_id", "ID", "Name"),
                            drop_scaffolds = TRUE,
                            scaffold_regex = "^(scaf|scaffold|un|unplaced|chrUn|KI|GL|NW_|NT_)",
                            chr_prefix     = "chr",
                            recode_chr     = c("auto","none","order","map"),
                            chr_order      = NULL,
                            chr_map        = NULL,
                            strict_chr     = FALSE,
                            write_chr_map  = TRUE,
                            chr_map_out    = NULL) {

  recode_chr <- match.arg(recode_chr)

  if (!file.exists(gff)) stop("GFF3 not found: ", gff, call. = FALSE)

  x <- utils::read.delim(
    gff,
    header           = FALSE,
    comment.char     = "#",
    stringsAsFactors = FALSE,
    sep              = "\t"
  )
  if (ncol(x) < 9L) stop("GFF3 must have >= 9 columns.", call. = FALSE)

  colnames(x)[1:9] <- c("seqid","source","type","start","end","score","strand","phase","attr")

  x <- x[x$type %in% feature_type, , drop = FALSE]
  if (!nrow(x)) stop("No rows found for feature_type = ", paste(feature_type, collapse = ", "), call. = FALSE)

  # drop scaffold-like seqids (optional)
  if (drop_scaffolds) {
    keep <- !grepl(scaffold_regex, x$seqid, ignore.case = TRUE)
    x <- x[keep, , drop = FALSE]
    if (!nrow(x)) stop("All rows removed by drop_scaffolds/scaffold_regex.", call. = FALSE)
  }

  parse_id <- function(attr_string, id_fields) {
    if (is.na(attr_string) || attr_string == "") return(NA_character_)
    parts <- strsplit(attr_string, ";", fixed = TRUE)[[1]]
    kv <- strsplit(parts, "=", fixed = TRUE)
    keys <- vapply(kv, function(z) if (length(z) >= 1L) z[1] else "", character(1))
    vals <- vapply(kv, function(z) if (length(z) >= 2L) z[2] else NA_character_, character(1))
    for (tag in id_fields) {
      idx <- which(keys == tag)
      if (length(idx) > 0L && !is.na(vals[idx[1]]) && nzchar(vals[idx[1]])) return(vals[idx[1]])
    }
    NA_character_
  }

  gene_ids <- vapply(x$attr, parse_id, character(1), id_fields = id_fields)

  res <- data.frame(
    GENE     = gene_ids,
    CHR_ORIG = x$seqid,
    START    = as.integer(x$start),
    STOP     = as.integer(x$end),
    stringsAsFactors = FALSE
  )
  res <- res[!is.na(res$GENE) & nzchar(res$GENE), , drop = FALSE]
  if (!nrow(res)) stop("No gene IDs extracted using id_fields: ", paste(id_fields, collapse = ", "), call. = FALSE)

  # normalize chromosome labels (prefix stripping only; keeps 2L/2R/etc intact)
  chr_norm <- as.character(res$CHR_ORIG)
  if (!is.null(chr_prefix) && nzchar(chr_prefix)) {
    chr_norm <- sub(paste0("^", chr_prefix), "", chr_norm, ignore.case = TRUE)
  }
  res$CHR_ORIG <- chr_norm

  # decide whether to recode CHR to 1..K
  chr_map_used <- NULL
  if (recode_chr == "none") {
    res$CHR <- res$CHR_ORIG
  } else if (recode_chr == "order") {
    if (is.null(chr_order)) stop("recode_chr='order' requires chr_order.", call. = FALSE)
    chr_map_used <- magcat_make_chr_map(res$CHR_ORIG, chr_order = chr_order, strict = strict_chr)
    res$CHR <- magcat_apply_chr_map(res$CHR_ORIG, chr_map_used, strict = TRUE)
  } else if (recode_chr == "map") {
    if (is.null(chr_map)) stop("recode_chr='map' requires chr_map.", call. = FALSE)
    chr_map_used <- magcat_make_chr_map(res$CHR_ORIG, chr_map = chr_map, strict = strict_chr)
    res$CHR <- magcat_apply_chr_map(res$CHR_ORIG, chr_map_used, strict = TRUE)
  } else {
    # auto: if already all numeric, keep; else recode
    chr_num <- suppressWarnings(as.integer(res$CHR_ORIG))
    if (all(!is.na(chr_num))) {
      res$CHR <- res$CHR_ORIG
    } else {
      chr_map_used <- magcat_make_chr_map(res$CHR_ORIG, strict = strict_chr)
      res$CHR <- magcat_apply_chr_map(res$CHR_ORIG, chr_map_used, strict = TRUE)
    }
  }

  # order rows by CHR (numeric if possible), then START
  chr_num2 <- suppressWarnings(as.integer(res$CHR))
  if (all(!is.na(chr_num2))) {
    o <- order(chr_num2, res$START, res$STOP)
  } else {
    o <- order(res$CHR, res$START, res$STOP)
  }
  res <- res[o, , drop = FALSE]

  # write gene-loc (MAGMA expects: GENE CHR START STOP; keep only those 4 columns)
  out_dir <- dirname(out)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  gene_loc_out <- res[, c("GENE","CHR","START","STOP"), drop = FALSE]

  utils::write.table(
    gene_loc_out,
    file      = out,
    quote     = FALSE,
    sep       = "\t",
    row.names = FALSE,
    col.names = TRUE
  )

  # write mapping file (original -> magma_chr)
  map_path <- NULL
  if (isTRUE(write_chr_map) && !is.null(chr_map_used)) {
    if (is.null(chr_map_out)) chr_map_out <- paste0(out, ".chr_map.tsv")
    map_path <- chr_map_out
    map_df <- data.frame(
      chr_original = names(chr_map_used),
      chr_magma    = unname(chr_map_used),
      stringsAsFactors = FALSE
    )
    utils::write.table(map_df, map_path, quote = FALSE, sep = "\t",
                       row.names = FALSE, col.names = TRUE)
  }

  invisible(list(
    gene_loc_df  = res,
    chr_map      = chr_map_used,
    chr_map_out  = map_path
  ))
}










# #' Convert a GFF3 file to a MAGMA gene location file
# #'
# #' @param gff Path to a GFF3 file.
# #' @param out Path to write the MAGMA gene location file (tab-delimited).
# #' @param feature_type GFF3 feature type to use (default "gene").
# #' @param id_fields Candidate attribute tags to use as gene IDs
# #'   (searched in order), e.g. "gene_id", "ID", "Name".
# #' @param drop_scaffolds Logical; if TRUE (default), drop seqids that look
# #'   like scaffolds (seqid starting with "scaf", case-insensitive).
# #' @param chr_prefix Optional string to strip from the start of chromosome
# #'   names (seqid) when generating the CHR column, e.g. "chr" will turn
# #'   "chr1" into "1", "chr10" into "10". Use NULL to keep seqids as-is.
# #'
# #' @return Invisibly returns a data.frame with columns GENE, CHR, START, STOP.
# #' @export
# gff3_to_geneloc <- function(gff,
#                             out,
#                             feature_type   = "gene",
#                             id_fields      = c("gene_id", "ID", "Name"),
#                             drop_scaffolds = TRUE,
#                             chr_prefix     = "chr") {

#   # read GFF3, skipping comment lines starting with '#'
#   x <- utils::read.delim(
#     gff,
#     header           = FALSE,
#     comment.char     = "#",
#     stringsAsFactors = FALSE,
#     sep              = "\t"
#   )

#   if (ncol(x) < 9L) {
#     stop("GFF3 file must have at least 9 columns.")
#   }

#   colnames(x)[1:9] <- c(
#     "seqid", "source", "type", "start",
#     "end", "score", "strand", "phase", "attr"
#   )

#   # keep only the desired feature type(s)
#   keep <- x$type %in% feature_type
#   x <- x[keep, , drop = FALSE]

#   if (!nrow(x)) {
#     stop(sprintf("No rows with type %s found in GFF3.", feature_type))
#   }

#   # optionally drop scaffold seqids
#   if (drop_scaffolds) {
#     not_scaf <- !grepl("^scaf", x$seqid, ignore.case = TRUE)
#     x <- x[not_scaf, , drop = FALSE]

#     if (!nrow(x)) {
#       stop("All entries were removed by drop_scaffolds=TRUE; ",
#            "no non-scaffold seqids left.")
#     }
#   }

#   # helper to parse the attributes string and extract ID
#   parse_id <- function(attr_string, id_fields) {
#     if (is.na(attr_string) || attr_string == "") return(NA_character_)

#     parts <- strsplit(attr_string, ";", fixed = TRUE)[[1]]
#     kv <- strsplit(parts, "=", fixed = TRUE)

#     keys <- vapply(kv, function(z) if (length(z) >= 1L) z[1] else "", character(1))
#     vals <- vapply(kv, function(z) if (length(z) >= 2L) z[2] else NA_character_, character(1))

#     for (tag in id_fields) {
#       idx <- which(keys == tag)
#       if (length(idx) > 0L && !is.na(vals[idx[1]])) {
#         return(vals[idx[1]])
#       }
#     }

#     NA_character_
#   }

#   gene_ids <- vapply(x$attr, parse_id, character(1), id_fields = id_fields)

#   res <- data.frame(
#     GENE  = gene_ids,
#     CHR   = x$seqid,
#     START = x$start,
#     STOP  = x$end,
#     stringsAsFactors = FALSE
#   )

#   # drop rows without a gene ID
#   res <- res[!is.na(res$GENE) & res$GENE != "", , drop = FALSE]

#   if (!nrow(res)) {
#     stop("Could not extract any gene IDs from attributes using tags: ",
#          paste(id_fields, collapse = ", "))
#   }

#   # strip chromosome prefix if asked (e.g. "chr1" -> "1")
#   if (!is.null(chr_prefix)) {
#     pattern <- paste0("^", chr_prefix)
#     res$CHR <- sub(pattern, "", res$CHR, ignore.case = TRUE)
#   }

#   ## ----- order by CHR (numeric if possible), then START -----
#   chr_num <- suppressWarnings(as.integer(res$CHR))
#   if (all(!is.na(chr_num))) {
#     o <- order(chr_num, res$START)
#   } else {
#     # fallback: alphabetical CHR, then START
#     o <- order(res$CHR, res$START)
#   }

#   res <- res[o, , drop = FALSE]
#   ## ---------------------------------------------------------

#   # write MAGMA gene-loc file
#   utils::write.table(
#     res,
#     file      = out,
#     quote     = FALSE,
#     sep       = "\t",
#     row.names = FALSE,
#     col.names = TRUE
#   )

#   invisible(res)
# }
