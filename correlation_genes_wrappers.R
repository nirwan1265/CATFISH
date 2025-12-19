magma_genesraw_to_cor_pairs_banded <- function(genes_raw_file,
                                               out_pairs_file,
                                               fixed_fields = 9L,
                                               gene_regex = "^Zm",   # adjust if needed
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


# Run correlation calculations
raw_dir <- "/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_multi_snp_wise_genes_by_chr_N_maize/"
chr_files <- Sys.glob(file.path(raw_dir, "N_maize_MLM_chr*.multi_snp_wise.genes.raw"))

chr_num <- as.integer(sub(".*_chr([0-9]+)\\..*$", "\\1", basename(chr_files)))
chr_files <- chr_files[order(chr_num)]

out_pairs <- file.path(raw_dir, "magma_gene_cor_pairs_MLM.txt")
if (file.exists(out_pairs)) file.remove(out_pairs)

first <- TRUE
for (f in chr_files) {
  tmp <- paste0(tempfile(), ".txt")

  magma_genesraw_to_cor_pairs_banded(
    genes_raw_file = f,
    out_pairs_file = tmp,
    keep_abs_r_ge  = 0,
    overwrite      = TRUE,
    verbose        = FALSE
  )

  x <- readLines(tmp, warn = FALSE)
  if (!length(x)) next
  if (!first) x <- x[-1]  # drop header
  writeLines(x, out_pairs, useBytes = TRUE, sep = "\n", append = !first)
  first <- FALSE
}

head(x)
out_pairs

x <- x[nzchar(x)]          # <-- drop blank lines

if (!first) x <- x[-1]     # drop header after first file

cat(paste(x, collapse = "\n"), "\n",
    file = out_pairs, append = !first)
first <- FALSE
head(x)

writeLines(x, out_pairs, useBytes = TRUE, sep = "\n")
