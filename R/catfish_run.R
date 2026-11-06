#' Run the whole CATFISH analysis in one call
#'
#' Runs every step from GWAS summary statistics, or from MAGMA output you
#' already have, to pathway results and candidate genes:
#' \enumerate{
#'   \item MAGMA SNP-to-gene annotation (\code{\link{magma_annotate}}) and gene
#'     analysis (\code{\link{magma_gene}}), when \code{gwas_file} is given;
#'   \item reading the gene results of every chromosome
#'     (\code{\link{magma_read_genes}});
#'   \item gene--gene correlations from the \code{.genes.raw} files
#'     (\code{\link{magma_genesraw_to_cor_pairs_banded}});
#'   \item Westfall--Young pathway testing (\code{\link{catfish_wy_pathways}});
#'   \item optionally, the five component tests (ACAT, Fisher, soft TFisher,
#'     minP, Stouffer) for every pathway, added to the same table
#'     (\code{\link{catfish_omni2_pathways}} with \code{perm_mode = "none"});
#'   \item candidate gene ranking (\code{\link{gene_prioritize}}).
#' }
#' Each step is also an exported function, so any part can be run on its own
#' (for example MAGMA one chromosome per cluster job, then \code{catfish_run()}
#' with \code{magma_dir}).
#'
#' @section Where the analysis starts:
#' Give exactly one of
#' \describe{
#'   \item{\code{gwas_file} (without \code{magma_dir})}{GWAS summary statistics;
#'     MAGMA is run (needs \code{bfile}, \code{rename_columns} and
#'     \code{gene_loc} or a \code{species} with built-in gene locations).}
#'   \item{\code{magma_dir}}{a folder of MAGMA \code{.genes.out} files and their
#'     matching \code{.genes.raw} files; MAGMA is not run. \code{gwas_file} may
#'     still be given, then it is used only for the SNP layer of the candidate
#'     gene ranking.}
#' }
#'
#' @param gwas_file GWAS summary statistics file (one SNP per row).
#' @param rename_columns Named character vector mapping the standard names to
#'   the columns of \code{gwas_file}: \code{CHR}, \code{SNP}, \code{POS},
#'   \code{PVALUE}, and optionally a per-SNP sample size (\code{N},
#'   \code{NOBS}, or \code{NMISS} with \code{n_total}). Example:
#'   \code{c(CHR = "chr", SNP = "rs", POS = "ps", PVALUE = "p_wald", NMISS = "n_miss")}.
#' @param bfile PLINK prefix of the LD reference genotypes (\code{.bed/.bim/.fam}).
#' @param gene_loc MAGMA gene-location file. If \code{NULL}, the built-in file of
#'   \code{species} is used (maize, sorghum, arabidopsis).
#' @param n_total Total sample size, when \code{rename_columns} has no
#'   per-SNP \code{N}/\code{NOBS}, or with \code{NMISS}.
#' @param window SNP-to-gene window in kb, \code{c(upstream, downstream)}.
#' @param gene_model MAGMA gene model (one). \code{"snp-wise=mean"} is
#'   recommended: see Details.
#' @param chroms Chromosomes to analyse one MAGMA run at a time (e.g.
#'   \code{1:10}). \code{NULL} runs MAGMA once on the whole genome.
#' @param n_threads Parallel MAGMA runs when \code{chroms} is given.
#' @param nonhuman Passed to MAGMA annotation (\code{TRUE} for non-human data).
#' @param sep Field separator of \code{gwas_file}; inferred from the extension
#'   if \code{NULL}.
#' @param chr_map_path Optional chromosome-renaming map (see
#'   \code{\link{gff3_to_geneloc}}).
#' @param magma_dir Folder with existing MAGMA output (\code{.genes.out} and
#'   \code{.genes.raw}).
#' @param magma_pattern Regular expression selecting the \code{.genes.out}
#'   files in \code{magma_dir}; default all of them.
#' @param gene_annot Optional MAGMA \code{.genes.annot} for the SNP layer of the
#'   candidate gene ranking when starting from \code{magma_dir} with a
#'   \code{gwas_file}. Without it, SNPs are assigned to genes by position.
#' @param pathways Pathways: named list of gene vectors, or data.frame with
#'   \code{pathway_id}, \code{gene_id} (optional \code{pathway_name}). Give this
#'   or \code{species}.
#' @param species Built-in pathways (and gene locations), e.g. \code{"sorghum"};
#'   see \code{\link{catfish_load_pathways}}.
#' @param pmn_gene_col Gene column of a built-in pathway file.
#' @param p_col Gene \eqn{p}-value column. Default: \code{P_SNPWISE_MEAN} when
#'   present (a \code{multi=snp-wise} run), otherwise \code{P}.
#' @param gene_regex Regular expression identifying gene IDs in the
#'   \code{.genes.raw} files; worked out from the gene IDs if \code{NULL}.
#' @param B Westfall--Young null replicates.
#' @param use_tfisher Include soft TFisher in the omnibus.
#' @param tau_grid Fixed soft TFisher grid (also used for the component table).
#' @param min_genes,max_genes Pathway size limits (genes present).
#' @param seed Seed for the null draws.
#' @param components Add the five analytic component \eqn{p}-values to the
#'   pathway table.
#' @param prioritize Rank candidate genes.
#' @param targets Optional known genes; their ranks are reported.
#' @param w_gwas,w_magma,w_pathway,w_hub Layer weights of
#'   \code{\link{gene_prioritize}}.
#' @param out_dir Output folder.
#' @param out_prefix Prefix of the output files.
#' @param verbose Print progress.
#'
#' @details
#' \strong{Gene model.} Westfall--Young calibration assumes gene \eqn{p}-values
#' are uniform where there is no association. MAGMA's SNP-wise mean model meets
#' this when the GWAS is well calibrated; the SNP-wise top model, and
#' \code{multi=snp-wise}, are often inflated in panels with long-range LD.
#' The median of every gene \eqn{p}-value column is printed so this can be
#' checked (about 0.5 is expected).
#'
#' \strong{Files written to \code{out_dir}}: \code{magma/} (MAGMA output, when
#' MAGMA is run), \code{<out_prefix>_gene_cor_pairs.txt} (reused on a rerun if
#' newer than the \code{.genes.raw} files), \code{<out_prefix>_pathways.csv},
#' \code{<out_prefix>_candidate_genes.csv}.
#'
#' @return A list:
#' \describe{
#'   \item{pathways}{Westfall--Young results, one row per pathway (columns of
#'     \code{\link{catfish_wy_pathways}}, plus \code{acat_p}, \code{fisher_p},
#'     \code{tfisher_p}, \code{minp_p}, \code{stouffer_p} when
#'     \code{components = TRUE})}
#'   \item{candidates}{\code{\link{gene_prioritize}} table, or \code{NULL}}
#'   \item{genes}{gene-level results used}
#'   \item{pathway_sets}{pathway memberships used}
#'   \item{settings}{\code{p_col}, \code{gene_regex}, \code{B}, \code{tau_grid}, ...}
#'   \item{files}{paths of everything written}
#' }
#'
#' @examples
#' \dontrun{
#' options(magma.path = "/path/to/magma")
#'
#' # from GWAS summary statistics
#' res <- catfish_run(
#'   gwas_file      = "gwas.assoc.txt",
#'   rename_columns = c(CHR = "chr", SNP = "rs", POS = "ps", PVALUE = "p_wald",
#'                      NMISS = "n_miss"),
#'   n_total        = 327,
#'   bfile          = "genotypes",
#'   species        = "sorghum",
#'   chroms         = 1:10, n_threads = 4,
#'   out_dir        = "dry_tons")
#' head(res$pathways)
#' head(res$candidates)
#'
#' # from MAGMA output computed elsewhere (e.g. on a cluster)
#' res <- catfish_run(magma_dir = "magma", magma_pattern = "^dry_tons_chr",
#'                    species = "sorghum", out_dir = "dry_tons")
#' }
#' @seealso \code{\link{catfish_wy_pathways}}, \code{\link{gene_prioritize}},
#'   \code{\link{magma_gene}}
#' @export
catfish_run <- function(gwas_file      = NULL,
                        rename_columns = NULL,
                        bfile          = NULL,
                        gene_loc       = NULL,
                        n_total        = NULL,
                        window         = c(25, 25),
                        gene_model     = "snp-wise=mean",
                        chroms         = NULL,
                        n_threads      = 1L,
                        nonhuman       = TRUE,
                        sep            = NULL,
                        chr_map_path   = NULL,
                        magma_dir      = NULL,
                        magma_pattern  = NULL,
                        gene_annot     = NULL,
                        pathways       = NULL,
                        species        = NULL,
                        pmn_gene_col   = NULL,
                        p_col          = NULL,
                        gene_regex     = NULL,
                        B              = 10000L,
                        use_tfisher    = TRUE,
                        tau_grid       = c(0.01, 0.05, 0.5),
                        min_genes      = 10L,
                        max_genes      = 200L,
                        seed           = 20260706L,
                        components     = TRUE,
                        prioritize     = TRUE,
                        targets        = NULL,
                        w_gwas         = 1,
                        w_magma        = 1,
                        w_pathway      = 0.5,
                        w_hub          = 0.3,
                        out_dir        = "catfish_results",
                        out_prefix     = "catfish",
                        verbose        = TRUE) {
  say <- function(...) if (isTRUE(verbose)) message(...)
  run_magma <- is.null(magma_dir)
  if (run_magma && is.null(gwas_file))
    stop("Give 'gwas_file' (MAGMA is run) or 'magma_dir' (existing MAGMA output).", call. = FALSE)
  if (is.null(pathways) && is.null(species))
    stop("Give 'pathways' or 'species'.", call. = FALSE)
  if (!is.null(gwas_file) && !file.exists(gwas_file))
    stop("gwas_file not found: ", gwas_file, call. = FALSE)
  if (!is.null(gwas_file)) {
    need <- c("CHR", "SNP", "POS", "PVALUE")
    if (is.null(rename_columns) || !all(need %in% names(rename_columns)))
      stop("rename_columns must map CHR, SNP, POS and PVALUE to columns of gwas_file.", call. = FALSE)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  files <- list()
  ann <- NULL

  ## ---- 1. MAGMA ------------------------------------------------------------
  if (run_magma) {
    if (is.null(bfile)) stop("bfile (PLINK LD reference) is needed to run MAGMA.", call. = FALSE)
    if (length(gene_model) != 1L) stop("Give one gene_model.", call. = FALSE)
    if (is.null(gene_loc)) {
      if (is.null(species)) stop("Give gene_loc, or a species with built-in gene locations.", call. = FALSE)
      gene_loc <- builtin_geneloc_path(species)
    }
    magma_dir <- file.path(out_dir, "magma")
    say("[CATFISH 1/6] MAGMA annotation")
    ann <- magma_annotate(stats_file = gwas_file, rename_columns = rename_columns,
                          gene_loc = gene_loc, out_prefix = out_prefix,
                          out_dir = magma_dir, window = window, sep = sep,
                          nonhuman = nonhuman, chr_map_path = chr_map_path,
                          verbose = verbose)
    if (!file.exists(ann$gene_annot))
      stop("MAGMA annotation wrote no ", ann$gene_annot, "; see MAGMA's .log in ", magma_dir, call. = FALSE)
    gene_annot <- ann$gene_annot
    files$gene_annot <- gene_annot
    say("[CATFISH 1/6] MAGMA gene analysis (", gene_model, ")")
    args <- list(bfile = bfile, gene_annot = gene_annot, stats_file = gwas_file,
                 rename_columns = rename_columns, out_prefix = out_prefix,
                 out_dir = magma_dir, gene_model = gene_model, chroms = chroms,
                 n_threads = n_threads, sep = sep, chr_map_path = chr_map_path,
                 verbose = verbose)
    if (!is.null(n_total)) args$n_total <- n_total
    do.call(magma_gene, args)
    tag <- gsub("[^A-Za-z0-9]+", "_", gene_model)
    magma_pattern <- paste0("^", gsub("([.|()^${}+*?\\[\\]\\\\])", "\\\\\\1", out_prefix),
                            "(_chr[^.]+)?\\.", tag, "\\.genes\\.out$")
    files$magma_dir <- magma_dir
  } else {
    say("[CATFISH 1/6] MAGMA not run: reading ", magma_dir)
  }

  ## ---- 2. gene results ------------------------------------------------------
  say("[CATFISH 2/6] gene results")
  genes <- magma_read_genes(dir = magma_dir,
                            pattern = if (is.null(magma_pattern)) "\\.genes\\.out$" else magma_pattern,
                            verbose = verbose)
  for (cl in c("GENE", "CHR"))
    if (!cl %in% names(genes)) stop("MAGMA output has no ", cl, " column.", call. = FALSE)
  if (is.null(p_col)) p_col <- if ("P_SNPWISE_MEAN" %in% names(genes)) "P_SNPWISE_MEAN" else "P"
  if (!p_col %in% names(genes))
    stop("p_col '", p_col, "' not in the MAGMA output; columns: ",
         paste(names(genes), collapse = ", "), call. = FALSE)
  say(sprintf("  gene p-value column: %s (median %.3f)", p_col,
              stats::median(as.numeric(genes[[p_col]]), na.rm = TRUE)))

  ## ---- 3. gene-gene correlations ------------------------------------------
  say("[CATFISH 3/6] gene-gene correlations")
  raw <- sub("\\.genes\\.out$", ".genes.raw", attr(genes, "files"))
  miss <- raw[!file.exists(raw)]
  if (length(miss))
    stop("Missing MAGMA .genes.raw file(s), which hold the gene correlations:\n  ",
         paste(miss, collapse = "\n  "), call. = FALSE)
  if (is.null(gene_regex)) gene_regex <- .catfish_gene_regex(genes$GENE)
  cor_file <- file.path(out_dir, paste0(out_prefix, "_gene_cor_pairs.txt"))
  if (file.exists(cor_file) && file.mtime(cor_file) >= max(file.mtime(raw))) {
    say("  reusing ", cor_file, " (newer than the .genes.raw files)")
  } else {
    magma_genesraw_to_cor_pairs_banded(raw, cor_file, gene_regex = gene_regex,
                                       verbose = FALSE)
    say("  ", cor_file)
  }
  files$cor_pairs <- cor_file

  ## ---- 4. Westfall-Young ---------------------------------------------------
  say("[CATFISH 4/6] Westfall-Young pathway test (B = ", B, ")")
  pw <- if (!is.null(pathways)) pathways else if (is.null(pmn_gene_col))
    catfish_load_pathways(species = species) else
    catfish_load_pathways(species = species, gene_col = pmn_gene_col)
  wy <- catfish_wy_pathways(gene_results = genes, pathways = pw, gene_col = "GENE",
                            p_col = p_col, chr_col = "CHR", cor_file = cor_file,
                            B = B, use_tfisher = use_tfisher, tau_grid = tau_grid,
                            min_genes = min_genes, max_genes = max_genes,
                            seed = seed, verbose = verbose)

  ## ---- 5. component tests -------------------------------------------------
  if (isTRUE(components)) {
    say("[CATFISH 5/6] component tests (analytic)")
    g2 <- genes
    pp <- pmin(pmax(as.numeric(g2[[p_col]]), 1e-300), 1)
    g2$CATFISH_Z <- stats::qnorm(pp, lower.tail = FALSE)   # z of the same p column
    om <- suppressWarnings(catfish_omni2_pathways(
      gene_results = g2, pathways = pw, gene_col = "GENE", p_raw_col = p_col,
      z_col = "CATFISH_Z", tau_grid = tau_grid, use_tfisher = use_tfisher,
      min_genes = min_genes, max_genes = max_genes, perm_mode = "none"))
    take <- c(acat_p = "acat_p", fisher_p = "fisher_p", tfisher_p = "tfisher_p_analytic",
              minp_p = "minp_p_analytic", stouffer_p = "stouffer_p_analytic")
    m <- match(wy$pathway_id, om$pathway_id)
    for (k in names(take))
      wy[[k]] <- if (take[[k]] %in% names(om)) om[[take[[k]]]][m] else NA_real_
  } else {
    say("[CATFISH 5/6] component tests skipped")
  }
  files$pathways <- file.path(out_dir, paste0(out_prefix, "_pathways.csv"))
  utils::write.csv(wy, files$pathways, row.names = FALSE)

  ## ---- 6. candidate genes -------------------------------------------------
  cand <- NULL
  if (isTRUE(prioritize)) {
    say("[CATFISH 6/6] candidate genes")
    gp <- list(magma = genes, magma_p_col = p_col, pathway_results = wy,
               pathways = pw, targets = targets, w_gwas = w_gwas,
               w_magma = w_magma, w_pathway = w_pathway, w_hub = w_hub,
               verbose = verbose)
    if (!is.null(gwas_file)) {
      gp$gwas <- gwas_file
      gp$gwas_snp_col <- unname(rename_columns[["SNP"]])
      gp$gwas_p_col   <- unname(rename_columns[["PVALUE"]])
      gp$gwas_chr_col <- unname(rename_columns[["CHR"]])
      gp$gwas_pos_col <- unname(rename_columns[["POS"]])
      if (!is.null(gene_annot)) {
        gp$gene_annot <- gene_annot
      } else {
        ## by position, with the gene coordinates MAGMA reports (window included)
        gp$gene_loc <- data.frame(gene = genes$GENE, chr = genes$CHR,
                                  start = genes$START, end = genes$STOP)
        gp$window <- 0L
      }
    }
    cand <- do.call(gene_prioritize, gp)
    files$candidates <- file.path(out_dir, paste0(out_prefix, "_candidate_genes.csv"))
    utils::write.csv(cand, files$candidates, row.names = FALSE)
  } else {
    say("[CATFISH 6/6] candidate genes skipped")
  }

  if (isTRUE(verbose)) {
    message(sprintf("\nCATFISH: %d pathways tested; %d with Westfall-Young FWER <= 0.05",
                    nrow(wy), sum(wy$p_fwer_std <= 0.05, na.rm = TRUE)))
    message("  ", paste(unlist(files), collapse = "\n  "))
  }
  list(pathways = wy, candidates = cand, genes = genes, pathway_sets = pw,
       settings = list(p_col = p_col, gene_regex = gene_regex, gene_model =
                         if (run_magma) gene_model else NA_character_,
                       B = B, tau_grid = tau_grid, use_tfisher = use_tfisher,
                       min_genes = min_genes, max_genes = max_genes, seed = seed),
       files = files)
}

## Regular expression matching every gene ID, for the .genes.raw parser:
## the shared leading letters (SORBI_, Zm, FBgn, AT, ...) or digits only (Entrez).
.catfish_gene_regex <- function(ids) {
  ids <- unique(as.character(ids))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (!length(ids)) stop("No gene IDs.", call. = FALSE)
  if (all(grepl("^[0-9]+$", ids))) return("^[0-9]+$")
  pre <- unique(sub("^([A-Za-z_]*).*$", "\\1", ids))
  if (length(pre) == 1L && nzchar(pre)) return(paste0("^", pre))
  pre <- pre[nzchar(pre)]
  if (length(pre) && length(pre) <= 20L && all(grepl("^[A-Za-z_]", ids)))
    return(paste0("^(", paste(sort(pre), collapse = "|"), ")"))
  stop("Could not work out a pattern for the gene IDs (e.g. ", paste(utils::head(ids, 3), collapse = ", "),
       "); give gene_regex.", call. = FALSE)
}
