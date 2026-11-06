# ----------------------------
# 1) Built-in file path helper
# ----------------------------
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
  path <- system.file("extdata", "MAGMA_pathway", fname, package = "CATFISH")
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

  path <- system.file("extdata", "MAGMA_pathway", fname, package = "CATFISH")
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
    pathways <- catfish_load_pathways(species = species, gene_col = pmn_gene_col)
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

