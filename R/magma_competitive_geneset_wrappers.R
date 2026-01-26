# #' Run MAGMA competitive gene-set (pathway) analysis
# #'
# #' MAGMA is picky about the --gene-results file: it expects the standard gene output
# #' columns and will error if extra columns exist (like P_SNPWISE_MEAN/TOP1).
# #' This wrapper can auto-create a MAGMA-clean gene results file with the required
# #' columns (GENE, CHR, START, STOP, NSNPS, NPARAM, N, ZSTAT, P), and optionally
# #' rename user-provided columns to those MAGMA names.
# #'
# #' @param gene_results_raw Path to gene results table (TSV/space-delimited) OR a data.frame.
# #' @param set_annot Path to MAGMA set annotation file (*.genes.annot).
# #' @param out_prefix Output prefix (no extension).
# #' @param out_dir Optional directory for outputs.
# #' @param clean_gene_results Logical; if TRUE, write a cleaned file for MAGMA and use it.
# #' @param required_cols The exact columns MAGMA expects for --gene-results.
# #' @param col_map Optional named character vector mapping YOUR columns -> MAGMA columns.
# #'   Example: c(gene="GENE", chromosome="CHR", bp_start="START", bp_end="STOP", nsnps="NSNPS",
# #'              nparam="NPARAM", n="N", z="ZSTAT", p="P")
# #'   Names are "your column names", values are MAGMA standard names.
# #' @param sep Output separator for cleaned file (default tab).
# #'
# #' @return Invisible output prefix path used by MAGMA.
# #' @export
# #' Run MAGMA competitive gene-set (pathway) analysis
# #' @export
# magma_geneset_competitive <- function(gene_results_raw,
#                                       set_annot,
#                                       out_prefix,
#                                       out_dir = NULL,
#                                       clean_gene_results = FALSE) {

#   mp <- magma_path()

#   if (!file.exists(gene_results_raw)) stop("gene_results_raw not found: ", gene_results_raw, call.=FALSE)
#   if (!file.exists(set_annot))        stop("set_annot not found: ", set_annot, call.=FALSE)

#   prefix_full <- if (is.null(out_dir)) out_prefix else {
#     if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#     file.path(out_dir, out_prefix)
#   }

#   # Auto-detect MAGMA .genes.raw style (starts with "# VERSION")
#   first_line <- readLines(gene_results_raw, n = 1L, warn = FALSE)
#   is_genes_raw <- length(first_line) && grepl("^#\\s*VERSION\\s*=", first_line)

#   if (is_genes_raw) {
#     clean_gene_results <- FALSE
#   }

#   if (isTRUE(clean_gene_results)) {
#     stop("This MAGMA build expects .genes.raw style input (with COVAR=NSAMP MAC). ",
#          "Pass the MAGMA-produced *.genes.raw file, or set clean_gene_results=FALSE.",
#          call.=FALSE)
#   }

#   args <- c(
#     "--gene-results", gene_results_raw,
#     "--set-annot",    set_annot,
#     "--out",          prefix_full
#   )

#   system2(mp, args)
#   invisible(prefix_full)
# }




# #' Run MAGMA competitive gene-set (pathway) analysis + return tidy table
# #'
# #' This keeps your wrapper simple (runs MAGMA), then immediately reads the
# #' resulting *.gsa.out and returns a CATFISH-ready table:
# #' pathway_id, pathway_name, n_genes, gene_names, gene_pvals, magma_pvalue
# #'
# #' IMPORTANT:
# #' - This custom MAGMA build expects *.genes.raw for --gene-results.
# #' - Your *.gsa.out has 4 metadata lines starting with '#', so we drop the first 4 lines.
# #'
# #' @param gene_results_raw Path to merged MAGMA *.genes.raw (all chr).
# #' @param set_annot Path to MAGMA set annotation (WIDE format).
# #' @param out_prefix Output prefix (no extension).
# #' @param out_dir Optional directory for outputs.
# #' @param pathways Pathway mapping table (e.g. maize_pw) with columns:
# #'   pathway_id, pathway_name, gene_id
# #' @param genes_all Gene results table with columns: GENE, P
# #' @param write_tidy Logical; write tidy TSV to disk (default TRUE)
# #' @param tidy_suffix Suffix for tidy file (default ".tidy.tsv")
# #'
# #' @return data.table with columns:
# #'   pathway_id, pathway_name, n_genes, gene_names, gene_pvals, magma_pvalue
# #' @export
# magma_geneset_competitive <- function(gene_results_raw,
#                                       set_annot,
#                                       out_prefix,
#                                       out_dir = NULL,
#                                       pathways,
#                                       genes_all,
#                                       write_tidy = TRUE,
#                                       tidy_suffix = ".tidy.tsv") {

#   # ---------- run MAGMA ----------
#   mp <- magma_path()

#   if (!file.exists(gene_results_raw)) stop("gene_results_raw not found: ", gene_results_raw, call.=FALSE)
#   if (!file.exists(set_annot))        stop("set_annot not found: ", set_annot, call.=FALSE)

#   prefix_full <- if (is.null(out_dir)) out_prefix else {
#     if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#     file.path(out_dir, out_prefix)
#   }

#   first_line <- readLines(gene_results_raw, n = 1L, warn = FALSE)
#   is_genes_raw <- length(first_line) && grepl("^#\\s*VERSION\\s*=", first_line)
#   if (!is_genes_raw) {
#     stop("This MAGMA build expects *.genes.raw (starts with '# VERSION ='). Provided: ", gene_results_raw, call.=FALSE)
#   }

#   args <- c(
#     "--gene-results", gene_results_raw,
#     "--set-annot",    set_annot,
#     "--out",          prefix_full
#   )

#   status <- system2(mp, args)
#   if (!is.null(status) && status != 0) {
#     stop("MAGMA failed (exit status ", status, "). Check log: ", paste0(prefix_full, ".log"), call.=FALSE)
#   }

#   # ---------- read gsa.out ----------
#   if (!requireNamespace("data.table", quietly = TRUE)) {
#     stop("magma_geneset_competitive(): requires the 'data.table' package.", call.=FALSE)
#   }
#   dt <- data.table::as.data.table

#   gsa_file <- paste0(prefix_full, ".gsa.out")
#   if (!file.exists(gsa_file)) stop("MAGMA .gsa.out not found: ", gsa_file, call.=FALSE)

#   lns  <- readLines(gsa_file, warn = FALSE)
#   if (length(lns) < 6) stop("gsa.out too short: ", gsa_file, call.=FALSE)

#   lns2 <- lns[-(1:4)]  # drop 4 metadata lines

#   gsa <- data.table::fread(
#     text = paste(lns2, collapse = "\n"),
#     header = TRUE,
#     data.table = TRUE
#   )

#   if (!all(c("VARIABLE","NGENES","P") %in% names(gsa))) {
#     stop("Unexpected gsa.out columns. Expected VARIABLE, NGENES, P. File: ", gsa_file, call.=FALSE)
#   }

#   # ---------- build tidy output ----------
#   # 1) MAGMA results
#   gsa_tidy <- gsa[, list(
#     pathway_id    = as.character(VARIABLE),
#     n_genes_magma = as.integer(NGENES),
#     magma_pvalue  = as.numeric(P)
#   )]

#   # checks for inputs
#   pathways <- dt(pathways)
#   genes_all <- dt(genes_all)

#   if (!all(c("pathway_id","pathway_name","gene_id") %in% names(pathways))) {
#     stop("pathways must have columns: pathway_id, pathway_name, gene_id", call.=FALSE)
#   }
#   if (!all(c("GENE","P") %in% names(genes_all))) {
#     stop("genes_all must have columns: GENE, P", call.=FALSE)
#   }

#   # 2) pathway_id → pathway_name
#   pw_map <- unique(pathways[, list(
#     pathway_id   = as.character(pathway_id),
#     pathway_name = as.character(pathway_name)
#   )])

#   # 3) pathway → genes
#   pw_genes <- unique(pathways[, list(
#     pathway_id = as.character(pathway_id),
#     gene_id    = sub("^gene:", "", as.character(gene_id))
#   )])

#   # 4) gene → pvalue
#   gene_p <- unique(genes_all[, list(
#     gene_id = as.character(GENE),
#     gene_p  = as.numeric(P)
#   )], by = "gene_id")

#   # 5) attach gene pvals to pathways
#   pw_gene_p <- merge(pw_genes, gene_p, by = "gene_id", all.x = TRUE)

#   pw_summ <- pw_gene_p[, list(
#     n_genes    = sum(!is.na(gene_p)),
#     gene_names = paste(gene_id[!is.na(gene_p)], collapse = ";"),
#     gene_pvals = paste(format(gene_p[!is.na(gene_p)], scientific = TRUE, digits = 6), collapse = ";")
#   ), by = pathway_id]

#   # 6) merge everything
#   out <- merge(gsa_tidy, pw_map, by = "pathway_id", all.x = TRUE)
#   out <- merge(out, pw_summ, by = "pathway_id", all.x = TRUE)

#   # fallbacks
#   out[is.na(pathway_name), pathway_name := pathway_id]
#   out[is.na(n_genes), n_genes := 0L]
#   out[is.na(gene_names), gene_names := ""]
#   out[is.na(gene_pvals), gene_pvals := ""]

#   # final shape
#   out <- out[, list(
#     pathway_id,
#     pathway_name,
#     n_genes,
#     gene_names,
#     gene_pvals,
#     magma_pvalue
#   )]

#   data.table::setorder(out, magma_pvalue)

#   # write tidy file next to MAGMA outputs (optional)
#   if (isTRUE(write_tidy)) {
#     tidy_file <- paste0(prefix_full, tidy_suffix)
#     data.table::fwrite(out, tidy_file, sep = "\t")
#     attr(out, "tidy_file") <- tidy_file
#   }

#   attr(out, "prefix_full") <- prefix_full
#   attr(out, "gsa_file") <- gsa_file
#   out
# }


# #' Run MAGMA competitive gene-set analysis and return a CATFISH-ready tidy table
# #'
# #' NEW: If set_annot is NULL, this function will *build* a MAGMA WIDE set-annotation
# #' file from `pathways` (data.frame/data.table OR a file path) using the same
# #' "post-editing" steps you were doing manually (sanitize gene ids, unique, and
# #' optionally keep only genes present in your gene results file).
# #'
# #' @param gene_results_raw Path to merged MAGMA gene results file used by --gene-results
# #'   (your build may accept *.genes.raw or your merged *.genes file with '#' metadata).
# #' @param set_annot Path to MAGMA set annotation (WIDE format). If NULL, it is created.
# #' @param out_prefix Output prefix (no extension).
# #' @param out_dir Optional directory for outputs.
# #' @param pathways Either:
# #'   (i) a data.frame/data.table with columns pathway_id, pathway_name, gene_id; OR
# #'   (ii) a file path to a TSV/CSV that contains those columns.
# #' @param genes_all Gene results table with columns: GENE, P (used to attach per-gene pvals).
# #' @param filter_to_gene_results Logical; keep only pathway genes that exist in gene_results_raw.
# #' @param write_tidy Logical; write tidy TSV next to MAGMA outputs.
# #' @param tidy_suffix Suffix for tidy file (default ".tidy.tsv")
# #' @param set_annot_suffix Suffix for auto-built set annotation (default ".sets.annot.WIDE")
# #'
# #' @return data.table with columns:
# #'   pathway_id, pathway_name, n_genes, gene_names, gene_pvals, magma_pvalue
# #' @export
# magma_geneset_competitive <- function(gene_results_raw,
#                                       set_annot = NULL,
#                                       out_prefix,
#                                       out_dir = NULL,
#                                       pathways,
#                                       genes_all,
#                                       filter_to_gene_results = TRUE,
#                                       write_tidy = TRUE,
#                                       tidy_suffix = ".tidy.tsv",
#                                       set_annot_suffix = ".sets.annot.WIDE") {

#   if (!requireNamespace("data.table", quietly = TRUE)) {
#     stop("magma_geneset_competitive(): requires the 'data.table' package.", call. = FALSE)
#   }

#   asDT <- data.table::as.data.table

#   # ---------- helpers ----------
#   read_gene_ids_from_gene_results <- function(path) {
#     if (!file.exists(path)) stop("gene_results_raw not found: ", path, call. = FALSE)
#     lns <- readLines(path, warn = FALSE)
#     lns <- lns[!is.na(lns) & nzchar(lns)]
#     # drop comment/metadata lines
#     dat <- lns[!grepl("^\\s*#", lns)]
#     if (!length(dat)) return(character(0))
#     # first token per line = gene id (your current approach)
#     unique(sub("\\s+.*$", "", dat))
#   }

#   load_pathways_any <- function(x) {
#     if (is.data.frame(x) || data.table::is.data.table(x)) return(asDT(x))
#     if (is.character(x) && length(x) == 1L && file.exists(x)) {
#       # fread will handle TSV/CSV reasonably; user can also give sep via file extension
#       return(asDT(data.table::fread(x, data.table = TRUE)))
#     }
#     stop("pathways must be a data.frame/data.table OR an existing file path.", call. = FALSE)
#   }

#   build_set_annot_wide <- function(pathways_dt, gene_ids_ok = NULL, out_path) {
#     if (!all(c("pathway_id", "gene_id") %in% names(pathways_dt))) {
#       stop("pathways must have columns: pathway_id, gene_id (and ideally pathway_name).", call. = FALSE)
#     }

#     dt <- asDT(pathways_dt)[, .(
#       set_id  = as.character(pathway_id),
#       gene_id = sub("^gene:", "", as.character(gene_id))
#     )]

#     dt <- dt[!is.na(set_id) & nzchar(set_id) & !is.na(gene_id) & nzchar(gene_id)]
#     dt <- unique(dt)

#     if (!is.null(gene_ids_ok) && length(gene_ids_ok)) {
#       dt <- dt[gene_id %in% gene_ids_ok]
#     }

#     setlist <- split(dt$gene_id, dt$set_id)
#     setlist <- lapply(setlist, unique)

#     dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

#     con <- file(out_path, open = "wt")
#     on.exit(close(con), add = TRUE)

#     for (sid in names(setlist)) {
#       genes <- setlist[[sid]]
#       if (!length(genes)) next
#       writeLines(paste(c(sid, genes), collapse = "\t"), con = con)
#     }

#     out_path
#   }

#   # ---------- validate gene results ----------
#   if (!file.exists(gene_results_raw)) stop("gene_results_raw not found: ", gene_results_raw, call. = FALSE)

#   first_line <- readLines(gene_results_raw, n = 1L, warn = FALSE)
#   if (!length(first_line)) stop("gene_results_raw is empty: ", gene_results_raw, call. = FALSE)
#   # accept either "# VERSION =" OR "# MEAN_SAMPLE_SIZE =" OR any '#' metadata style
#   if (!grepl("^\\s*#", first_line)) {
#     stop("gene_results_raw does not look like a MAGMA gene results file (expected '#' metadata header). Provided: ",
#          gene_results_raw, call. = FALSE)
#   }

#   # ---------- resolve output prefix ----------
#   prefix_full <- if (is.null(out_dir)) out_prefix else {
#     if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#     file.path(out_dir, out_prefix)
#   }

#   # ---------- build set annotation if needed ----------
#   pathways_dt <- load_pathways_any(pathways)

#   if (is.null(set_annot)) {
#     gene_ids_ok <- NULL
#     if (isTRUE(filter_to_gene_results)) {
#       gene_ids_ok <- read_gene_ids_from_gene_results(gene_results_raw)
#       if (!length(gene_ids_ok)) {
#         stop("Could not extract any gene IDs from gene_results_raw for filtering.", call. = FALSE)
#       }
#     }
#     set_annot <- paste0(prefix_full, set_annot_suffix)
#     build_set_annot_wide(pathways_dt, gene_ids_ok = gene_ids_ok, out_path = set_annot)
#   } else {
#     if (!file.exists(set_annot)) stop("set_annot not found: ", set_annot, call. = FALSE)
#   }

#   # ---------- run MAGMA ----------
#   mp <- magma_path()

#   args <- c(
#     "--gene-results", gene_results_raw,
#     "--set-annot",    set_annot,
#     "--out",          prefix_full
#   )

#   status <- system2(mp, args)
#   if (!is.null(status) && status != 0) {
#     stop("MAGMA failed (exit status ", status, "). Check log: ", paste0(prefix_full, ".log"), call. = FALSE)
#   }

#   # ---------- read gsa.out ----------
#   gsa_file <- paste0(prefix_full, ".gsa.out")
#   if (!file.exists(gsa_file)) stop("MAGMA .gsa.out not found: ", gsa_file, call. = FALSE)

#   lns  <- readLines(gsa_file, warn = FALSE)
#   if (length(lns) < 6) stop("gsa.out too short: ", gsa_file, call. = FALSE)

#   # your build writes 4 metadata lines starting with '#'
#   lns2 <- lns[-(1:4)]

#   gsa <- data.table::fread(
#     text = paste(lns2, collapse = "\n"),
#     header = TRUE,
#     data.table = TRUE
#   )

#   if (!all(c("VARIABLE", "NGENES", "P") %in% names(gsa))) {
#     stop("Unexpected gsa.out columns. Expected VARIABLE, NGENES, P. File: ", gsa_file, call. = FALSE)
#   }

#   # ---------- tidy output ----------
#   gsa_tidy <- gsa[, .(
#     pathway_id    = as.character(VARIABLE),
#     n_genes_magma = as.integer(NGENES),
#     magma_pvalue  = as.numeric(P)
#   )]

#   genes_all <- asDT(genes_all)
#   if (!all(c("GENE", "P") %in% names(genes_all))) {
#     stop("genes_all must have columns: GENE, P", call. = FALSE)
#   }

#   # pathway_id -> pathway_name
#   pw_map <- unique(pathways_dt[, .(
#     pathway_id   = as.character(pathway_id),
#     pathway_name = if ("pathway_name" %in% names(pathways_dt)) as.character(pathway_name) else as.character(pathway_id)
#   )])

#   # pathway -> genes
#   pw_genes <- unique(pathways_dt[, .(
#     pathway_id = as.character(pathway_id),
#     gene_id    = sub("^gene:", "", as.character(gene_id))
#   )])

#   # gene -> pvalue
#   gene_p <- unique(genes_all[, .(
#     gene_id = as.character(GENE),
#     gene_p  = as.numeric(P)
#   )], by = "gene_id")

#   pw_gene_p <- merge(pw_genes, gene_p, by = "gene_id", all.x = TRUE)

#   pw_summ <- pw_gene_p[, .(
#     n_genes    = sum(!is.na(gene_p)),
#     gene_names = paste(gene_id[!is.na(gene_p)], collapse = ";"),
#     gene_pvals = paste(format(gene_p[!is.na(gene_p)], scientific = TRUE, digits = 6), collapse = ";")
#   ), by = pathway_id]

#   out <- merge(gsa_tidy, pw_map, by = "pathway_id", all.x = TRUE)
#   out <- merge(out, pw_summ, by = "pathway_id", all.x = TRUE)

#   out[is.na(pathway_name), pathway_name := pathway_id]
#   out[is.na(n_genes), n_genes := 0L]
#   out[is.na(gene_names), gene_names := ""]
#   out[is.na(gene_pvals), gene_pvals := ""]

#   out <- out[, .(
#     pathway_id,
#     pathway_name,
#     n_genes,
#     gene_names,
#     gene_pvals,
#     magma_pvalue
#   )]

#   data.table::setorder(out, magma_pvalue)

#   if (isTRUE(write_tidy)) {
#     tidy_file <- paste0(prefix_full, tidy_suffix)
#     data.table::fwrite(out, tidy_file, sep = "\t")
#     attr(out, "tidy_file") <- tidy_file
#   }

#   attr(out, "prefix_full") <- prefix_full
#   attr(out, "gsa_file") <- gsa_file
#   attr(out, "set_annot") <- set_annot
#   out
# }

# -------------------------
# Example usage (your case):
# -------------------------
# maize_pw <- magcat_load_pathways("maize", gene_col = "Gene-name")
#
# out <- magma_geneset_competitive(
#   gene_results_raw = "/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_multi_snp_wise_genes_by_chr_N_maize/N_maize_MLM_ALLCHR.multi_snp_wise.genes",
#   set_annot        = "annot/maize_pathways.sets.annot.WIDE",  # OR set_annot = NULL to auto-build
#   out_prefix       = "N_maize_MLM_ALLCHR.PMN_COMP",
#   out_dir          = "magma_geneset",
#   pathways         = maize_pw,
#   genes_all        = genes_all,
#   filter_to_gene_results = TRUE
# )







# Yes — you can ship built-in MAGMA WIDE set-annotation files and use those.
# BUT: the WIDE file only contains:
#   set_id <tab> gene1 <tab> gene2 <tab> ...
# It does NOT contain pathway_name (human-readable name) and it does NOT contain
# a long membership table for building gene_names/gene_pvals lists.
#
# So:
# - If you only need MAGMA competitive p-values: WIDE is enough.
# - If you want your CATFISH-ready tidy table WITH pathway_name + gene_names + gene_pvals:
#   you still need a mapping table somewhere (at least pathway_id -> pathway_name, and ideally
#   pathway_id -> gene_id to list genes/pvals).
#
# Best practice for your package:
#   inst/extdata/MAGMA_pathway/
#     maize_pathways.sets.annot.WIDE        # for MAGMA
#     maize_pathways.long.tsv               # pathway_id, pathway_name, gene_id  (for tidy output)
#
# If you insist on ONLY ONE file: you can’t recover pathway_name from WIDE.

# ----------------------------
# 1) Built-in file path helper
# ----------------------------
pmn_gene_col     = "Gene-name"
builtin_magma_set_annot_path <- function(
  species = c("maize","sorghum","arabidopsis","plant","fly")
) {
  species <- match.arg(tolower(species),
                       c("maize","sorghum","arabidopsis","plant","fly"))

  fname <- switch(
    species,
    maize       = "maize_pathway.sets.annot.WIDE.txt",
    sorghum     = "sorghum_pathway.sets.annot.WIDE",
    arabidopsis = "arabidopsis_pathway.sets.annot.WIDE",
    plant       = "plant_pathway.sets.annot.WIDE"
  )

  # 1) DEV: robust path under load_all()
  root <- tryCatch(pkgload::pkg_path(), error = function(e) NULL)
  if (!is.null(root)) {
    dev_path <- file.path(root, "inst", "extdata", "MAGMA_pathway", fname)
    if (file.exists(dev_path)) {
      return(normalizePath(dev_path, winslash = "/", mustWork = TRUE))
    }
  }

  # 2) INSTALLED package
  path <- system.file("extdata", "MAGMA_pathway", fname, package = "MAGCAT")
  if (path == "") stop("Built-in set-annot not found: ", fname, call. = FALSE)
  path
}


builtin_pathway_long_path <- function(
  species = c("maize","sorghum","arabidopsis","plant","fly")
) {
  species <- match.arg(tolower(species),
                       c("maize","sorghum","arabidopsis","plant","fly"))

  fname <- switch(
    species,
    maize       = "maize_pathways.long.tsv",
    sorghum     = "sorghum_pathways.long.tsv",
    arabidopsis = "arabidopsis_pathways.long.tsv",
    plant       = "plant_pathways.long.tsv"
  )

  dev_path <- file.path("inst", "extdata", "MAGMA_pathway", fname)
  if (file.exists(dev_path)) return(normalizePath(dev_path, winslash = "/", mustWork = TRUE))

  path <- system.file("extdata", "MAGMA_pathway", fname, package = "MAGCAT")
  if (path == "") stop("Built-in pathway long table not found: ", fname, call. = FALSE)
  path
}


# ---------------------------------------------------------
# 2) Minimal change to YOUR working function (species mode)
# ---------------------------------------------------------
# Add two NEW args:
#   species = NULL
#   use_builtin_pathways = TRUE
#
# Logic:
#   - If set_annot is NULL and species provided: set_annot <- builtin_magma_set_annot_path(species)
#   - If pathways is missing and species provided: pathways <- fread(builtin_pathway_long_path(species))
#
# That way ONE call works and stays stable.
magma_geneset_competitive <- function(gene_results_raw,
                                      set_annot = NULL,
                                      out_prefix,
                                      out_dir = NULL,
                                      pathways = NULL,
                                      species = NULL,
                                      pmn_gene_col = "Gene-name",
                                      genes_all,
                                      write_tidy = TRUE,
                                      tidy_suffix = ".tidy.tsv") {

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("magma_geneset_competitive(): requires the 'data.table' package.", call. = FALSE)
  }

  if (!file.exists(gene_results_raw)) {
    stop("gene_results_raw not found: ", gene_results_raw, call. = FALSE)
  }

  prefix_full <- if (is.null(out_dir)) out_prefix else {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    file.path(out_dir, out_prefix)
  }

  # If user wants built-ins, fill them in here
  if (is.null(set_annot) && !is.null(species)) {
    set_annot <- builtin_magma_set_annot_path(species)
  }
  if (is.null(pathways) && !is.null(species)) {
    pathways <- magcat_load_pathways(species = species, gene_col = pmn_gene_col)
  }

  if (is.null(set_annot) || !file.exists(set_annot)) {
    stop("set_annot missing/not found. Provide set_annot or species.", call. = FALSE)
  }

  # ---------- run MAGMA ----------
  mp <- magma_path()
  args <- c(
    "--gene-results", gene_results_raw,
    "--set-annot",    set_annot,
    "--out",          prefix_full
  )

  status <- system2(mp, args)
  if (!is.null(status) && status != 0) {
    stop("MAGMA failed (exit status ", status, "). Check log: ", paste0(prefix_full, ".log"), call. = FALSE)
  }

  
  # ---------- read gsa.out ----------
  # ---------- read gsa.out (SIMPLE + STABLE) ----------
  gsa_file <- paste0(prefix_full, ".gsa.out")
  if (!file.exists(gsa_file)) stop("MAGMA .gsa.out not found: ", gsa_file, call. = FALSE)

  lns <- readLines(gsa_file, warn = FALSE)
  if (length(lns) < 6) stop("gsa.out too short: ", gsa_file, call. = FALSE)

  # your build writes exactly 4 metadata lines at top
  lns2 <- lns[-(1:4)]

  # parse as data.frame (avoid data.table NSE issues)
  gsa <- utils::read.table(
    text = paste(lns2, collapse = "\n"),
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  need <- c("VARIABLE", "NGENES", "P")
  if (!all(need %in% names(gsa))) {
    stop("Unexpected gsa.out columns. Found: ",
         paste(names(gsa), collapse = ", "),
         "\nNeed at least: ", paste(need, collapse = ", "),
         call. = FALSE)
  }

  # build the minimal MAGMA tidy table (base df)
  gsa_tidy <- data.frame(
    pathway_id    = as.character(gsa[["VARIABLE"]]),
    n_genes_magma = as.integer(gsa[["NGENES"]]),
    magma_pvalue  = as.numeric(gsa[["P"]]),
    stringsAsFactors = FALSE
  )

  # If no pathways mapping, return minimal tidy
  if (is.null(pathways)) {
    out <- data.frame(
      pathway_id   = gsa_tidy$pathway_id,
      pathway_name = gsa_tidy$pathway_id,
      n_genes      = gsa_tidy$n_genes_magma,
      gene_names   = "",
      gene_pvals   = "",
      magma_pvalue = gsa_tidy$magma_pvalue,
      stringsAsFactors = FALSE
    )

    out <- out[order(out$magma_pvalue), , drop = FALSE]
    if (isTRUE(write_tidy)) {
      utils::write.table(out, file = paste0(prefix_full, tidy_suffix),
                         sep = "\t", quote = FALSE, row.names = FALSE)
    }
    return(out)
  }

  # ---------- attach names + gene lists (BASE R) ----------
  pathways  <- as.data.frame(pathways, stringsAsFactors = FALSE)
  genes_all <- as.data.frame(genes_all, stringsAsFactors = FALSE)

  if (!all(c("pathway_id","pathway_name","gene_id") %in% names(pathways))) {
    stop("pathways must have columns: pathway_id, pathway_name, gene_id.\nFound: ",
         paste(names(pathways), collapse = ", "), call. = FALSE)
  }
  if (!all(c("GENE","P") %in% names(genes_all))) {
    stop("genes_all must have columns: GENE, P.\nFound: ",
         paste(names(genes_all), collapse = ", "), call. = FALSE)
  }

  # pathway map
  pw_map <- unique(pathways[, c("pathway_id","pathway_name"), drop = FALSE])

  # pathway -> genes (clean)
  pw_genes <- unique(data.frame(
    pathway_id = as.character(pathways[["pathway_id"]]),
    gene_id    = sub("^gene:", "", as.character(pathways[["gene_id"]])),
    stringsAsFactors = FALSE
  ))
  pw_genes <- pw_genes[nzchar(pw_genes$pathway_id) & nzchar(pw_genes$gene_id), , drop = FALSE]

  # gene -> p
  gene_p <- unique(data.frame(
    gene_id = as.character(genes_all[["GENE"]]),
    gene_p  = as.numeric(genes_all[["P"]]),
    stringsAsFactors = FALSE
  ))

  # join genes to p
  pw_gene_p <- merge(pw_genes, gene_p, by = "gene_id", all.x = TRUE)

  # summarize per pathway
  split_list <- split(pw_gene_p, pw_gene_p$pathway_id)
  pw_summ <- do.call(rbind, lapply(names(split_list), function(pid) {
    dd <- split_list[[pid]]
    ok <- !is.na(dd$gene_p)
    data.frame(
      pathway_id = pid,
      n_genes    = sum(ok),
      gene_names = paste(dd$gene_id[ok], collapse = ";"),
      gene_pvals = paste(format(dd$gene_p[ok], scientific = TRUE, digits = 6), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }))

  # merge everything onto MAGMA results
  out <- merge(gsa_tidy, pw_map,  by = "pathway_id", all.x = TRUE)
  out <- merge(out,     pw_summ, by = "pathway_id", all.x = TRUE)

  # fill
  out$pathway_name <- ifelse(
    is.na(out$pathway_name) | out$pathway_name == "",
    out$pathway_id,
    out$pathway_name
  )

  out$n_genes     <- ifelse(is.na(out$n_genes), 0L, out$n_genes)
  out$gene_names  <- ifelse(is.na(out$gene_names), "", out$gene_names)
  out$gene_pvals  <- ifelse(is.na(out$gene_pvals), "", out$gene_pvals)

  # final columns (NO fancy subset)
  out <- out[, c("pathway_id","pathway_name","n_genes","gene_names","gene_pvals","magma_pvalue"), drop = FALSE]
  out <- out[order(out$magma_pvalue), , drop = FALSE]

  if (isTRUE(write_tidy)) {
    utils::write.table(out, file = paste0(prefix_full, tidy_suffix),
                       sep = "\t", quote = FALSE, row.names = FALSE)
  }

  return(out)


}

