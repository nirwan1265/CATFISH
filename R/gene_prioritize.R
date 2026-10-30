# =============================================================================
# gene_prioritize(): three-layer candidate-gene prioritisation
#
#   score(g) = w_gwas * S_gwas(g) + w_magma * S_magma(g)
#            + w_pathway * S_pathway(g) + w_hub * B_hub(g)
#
#   S_layer(g) = -log10(rank / N), rank 1 = strongest, N = genes with evidence
#                of that kind; genes without that kind of evidence get 0
#   B_hub(g)   = log10(1 + number of tested pathways containing g)
#
# Same formula and defaults as catfish2/scripts/06_candidate_genes.R.
# =============================================================================

## ---------------------------------------------------------------- helpers ---

## data.table syntax is used without importing data.table (it is in Suggests)
.datatable.aware <- TRUE
if (getRversion() >= "2.15.1") utils::globalVariables("NL")

## rank 1 = smallest value = strongest evidence; NA stays NA
.gp_rank <- function(x) {
  r <- rep(NA_real_, length(x))
  ok <- !is.na(x)
  r[ok] <- rank(x[ok], ties.method = "average")
  r
}

## -log10(rank / N) in [0, log10(N)]; no evidence -> 0
.gp_rank_score <- function(x) {
  out <- -log10(.gp_rank(x) / sum(!is.na(x)))
  out[is.na(out)] <- 0
  out
}

.gp_key <- function(x) tolower(trimws(as.character(x)))

## a data.frame passes through; a single existing file path is read
.gp_table <- function(x, what, header = TRUE) {
  if (is.data.frame(x)) return(as.data.frame(x, stringsAsFactors = FALSE))
  if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) stop(what, " file not found: ", x, call. = FALSE)
    if (requireNamespace("data.table", quietly = TRUE))
      return(data.table::fread(x, header = header, data.table = FALSE,
                               showProgress = FALSE))
    return(utils::read.table(x, header = header, stringsAsFactors = FALSE,
                             comment.char = "", quote = "", fill = TRUE))
  }
  stop(what, " must be a data.frame or the path to a file.", call. = FALSE)
}

.gp_need_cols <- function(df, cols, what, arg) {
  miss <- setdiff(cols, names(df))
  if (length(miss))
    stop(sprintf("%s is missing column(s) %s. Columns present: %s. Set %s to the right name.",
                 what, paste(miss, collapse = ", "),
                 paste(names(df), collapse = ", "), arg), call. = FALSE)
}

## first candidate column name present in df, or NA
.gp_pick <- function(df, cand) {
  hit <- intersect(cand, names(df))
  if (length(hit)) hit[1] else NA_character_
}

## pathway membership from a list, a data.frame (pathway_id, gene_id) or a set
## file (one pathway per line: id then genes; .gmt files carry a description
## in field 2, which is dropped)
.gp_pathways <- function(pathways) {
  if (is.character(pathways) && length(pathways) == 1L) {
    if (!file.exists(pathways))
      stop("pathways file not found: ", pathways, call. = FALSE)
    ln <- readLines(pathways, warn = FALSE)
    ln <- ln[nzchar(trimws(ln)) & !startsWith(trimws(ln), "#")]
    gmt <- grepl("\\.gmt$", pathways, ignore.case = TRUE)
    f <- strsplit(trimws(ln), if (gmt) "\t" else "[ \t]+")
    id <- vapply(f, `[`, character(1), 1L)
    mem <- lapply(f, function(x) x[-seq_len(if (gmt) 2L else 1L)])
    names(mem) <- id
    return(mem)
  }
  .catfish_pathways_to_list(pathways)$p_list
}

## ------------------------------------------------------------------ export ---

#' Prioritise candidate genes from GWAS, MAGMA and pathway evidence
#'
#' Ranks every gene in the MAGMA gene-level results by combining three layers
#' of evidence: the strongest SNP association near the gene (GWAS layer), the
#' MAGMA gene-level \eqn{p}-value (MAGMA layer), and the strongest tested
#' pathway that contains the gene (pathway layer), plus a small bonus for genes
#' that belong to many tested pathways.
#'
#' Each layer is turned into a rank score \eqn{S = -\log_{10}(\mathrm{rank}/N)},
#' where rank 1 is the strongest gene and \eqn{N} is the number of genes with
#' evidence of that kind; genes without that kind of evidence score 0. The
#' combined score is
#' \deqn{\mathrm{score} = w_{gwas} S_{gwas} + w_{magma} S_{magma} +
#'       w_{pathway} S_{pathway} + w_{hub} \log_{10}(1 + n_{pathways}).}
#' Working on ranks makes the layers comparable whatever their scale, and keeps
#' the order of genes whose \eqn{p}-values are far below double precision when
#' the GWAS layer is given as \eqn{-\log_{10} p}.
#'
#' @section Inputs and required columns:
#' Every input can be an R object or the path to a file. Gene IDs must be in
#' the same ID space in all inputs (matching ignores case and surrounding
#' spaces).
#' \describe{
#'   \item{\code{magma}}{MAGMA gene results (e.g. \code{.genes.out}). Required
#'     columns: gene ID (\code{magma_gene_col}, default \code{"GENE"}) and gene
#'     \eqn{p}-value (\code{magma_p_col}; by default \code{"P_MULTI"} if present,
#'     otherwise \code{"P"}).}
#'   \item{\code{pathway_results}}{One row per tested pathway, e.g. the output
#'     of \code{\link{catfish_wy_pathways}} or \code{\link{catfish_omni2_pathways}}.
#'     Required columns: \code{pathway_id} (\code{pathway_id_col}) and a
#'     statistic for which \strong{smaller means stronger}
#'     (\code{pathway_stat_col}). By default the first of \code{z_std},
#'     \code{omni_z_std}, \code{log10_p_omnibus}, \code{omni_log10_p_wy_stat},
#'     \code{p_omnibus}, \code{omni_p_final} that is present. The default is
#'     the Westfall--Young standardised statistic (each pathway's observed
#'     omnibus standardised against its own null), so the pathway layer uses
#'     the same ordering as the family-wise inference. It orders pathways
#'     exactly as the standardised family-wise \eqn{p}-value does, without that
#'     \eqn{p}-value's floor at \eqn{1/(B+1)}, which ties the strongest
#'     pathways. Standardisation does not remove the dependence of pathway
#'     evidence on pathway size.}
#'   \item{\code{pathways} or \code{species}}{Pathway membership, which must
#'     be the same sets that produced \code{pathway_results}. Either
#'     \code{pathways}: a named list of gene-ID vectors, a data.frame with
#'     \code{pathway_id} and \code{gene_id}, or a set file with one pathway
#'     per line (\code{id gene gene ...}, as used by MAGMA; \code{.gmt} files
#'     are also read); or \code{species}: a built-in set loaded with
#'     \code{\link{catfish_load_pathways}} (\code{"maize"}, \code{"sorghum"},
#'     \code{"arabidopsis"}, \code{"plant"}, \code{"fly"}, \code{"human"}
#'     for HumanCyc, \code{"human_reactome"} for Reactome). Only pathways
#'     present in \code{pathway_results} are used.}
#'   \item{\code{gwas} (optional)}{SNP-level results. Required columns: SNP ID
#'     (\code{gwas_snp_col}, default \code{"SNP"}) and either
#'     \eqn{-\log_{10} p} (\code{gwas_neglog10p_col}; default the first of
#'     \code{NEGLOG10P}, \code{neglog10p}) or \eqn{p} (\code{gwas_p_col};
#'     default the first of \code{P}, \code{p}, \code{PVAL}). \eqn{-\log_{10} p}
#'     is used when both are present, because a \eqn{p} column is usually
#'     floored (e.g. at \eqn{10^{-300}}) and then ties the strongest genes. For
#'     the position route, also chromosome and base-pair position columns
#'     (\code{gwas_chr_col}, \code{gwas_pos_col}).}
#'   \item{\code{gene_annot} (optional)}{How SNPs are assigned to genes for the
#'     GWAS layer: MAGMA's \code{.genes.annot} file, or a data.frame with
#'     columns \code{gene} and \code{SNP}. Using the same annotation as the
#'     MAGMA run keeps the GWAS and MAGMA layers consistent.}
#'   \item{\code{gene_loc} (optional)}{Gene locations, as a MAGMA
#'     \code{.gene.loc} file (columns: gene ID, chromosome, start, end, and
#'     optionally strand and symbol) or a data.frame with columns \code{gene},
#'     \code{chr}, \code{start}, \code{end} and optionally \code{symbol}. Used
#'     to assign SNPs to genes by position (\eqn{\pm} \code{window} bp) when no
#'     \code{gene_annot} is given, and to add gene symbols.}
#' }
#' If \code{gwas} is \code{NULL}, or neither \code{gene_annot} nor
#' \code{gene_loc} can place its SNPs, the GWAS layer is 0 for every gene and
#' the ranking uses the MAGMA and pathway layers only.
#'
#' @param magma MAGMA gene results: data.frame or file path.
#' @param pathway_results Pathway results: data.frame or file path (CSV/TSV).
#' @param pathways Pathway membership: list, data.frame or set-file path.
#'   Give this or \code{species}.
#' @param species Built-in pathway set (see \code{\link{catfish_load_pathways}});
#'   give this or \code{pathways}.
#' @param pmn_gene_col Gene-ID column of the built-in set, passed to
#'   \code{catfish_load_pathways(gene_col = ...)} when \code{species} is used.
#' @param gwas Optional SNP-level results: data.frame or file path.
#' @param gene_annot Optional SNP-to-gene annotation: MAGMA \code{.genes.annot}
#'   path or data.frame with \code{gene} and \code{SNP}.
#' @param gene_loc Optional gene locations: \code{.gene.loc} path or data.frame.
#' @param magma_gene_col,magma_p_col Gene ID and \eqn{p}-value columns of
#'   \code{magma}.
#' @param pathway_id_col,pathway_stat_col Pathway ID and statistic columns of
#'   \code{pathway_results}.
#' @param gwas_snp_col,gwas_p_col,gwas_neglog10p_col,gwas_chr_col,gwas_pos_col
#'   Columns of \code{gwas}.
#' @param symbol_col Column of \code{gene_loc} holding gene symbols: a name, or
#'   a number for a header-less \code{.gene.loc} file (MAGMA's NCBI files keep
#'   the symbol in column 6). \code{NULL} uses a \code{symbol} column if
#'   present, otherwise column 6 of a 6-column file.
#' @param window Base pairs added on each side of a gene for the position
#'   route. Default 10000.
#' @param w_gwas,w_magma,w_pathway,w_hub Layer weights. Defaults 1, 1, 0.5, 0.3.
#' @param targets Optional character vector of known genes (IDs or symbols) to
#'   report as positive controls.
#' @param chunk Genes read at a time from a \code{.genes.annot} file.
#' @param verbose Print progress.
#' @param output If \code{TRUE}, write the table to \code{out_file}.
#' @param out_file CSV path used when \code{output = TRUE}.
#'
#' @return A data.frame with one row per MAGMA gene, sorted by \code{score}:
#' \describe{
#'   \item{rank}{Overall rank (1 = top candidate)}
#'   \item{gene}{Gene ID as given in \code{magma}}
#'   \item{symbol}{Gene symbol, if \code{gene_loc} supplies one}
#'   \item{score}{Combined score}
#'   \item{n_layers}{Number of the three layers with evidence for the gene}
#'   \item{gwas_score, magma_score, pathway_score, hub_bonus}{The four terms
#'     before weighting}
#'   \item{gwas_neglog10p, gwas_rank}{Strongest SNP association near the gene,
#'     as \eqn{-\log_{10} p}, and its rank among genes with a SNP}
#'   \item{magma_p, magma_rank}{MAGMA gene \eqn{p}-value and its rank}
#'   \item{best_pathway, pathway_stat, pathway_rank}{Strongest tested pathway
#'     containing the gene, its statistic (\code{pathway_stat_col}) and the
#'     gene's rank among pathway-annotated genes}
#'   \item{n_pathways}{Number of tested pathways containing the gene}
#' }
#' Attributes: \code{layer_sizes} (genes with evidence in each layer),
#' \code{pathway_stat_col}, \code{gwas_route}, \code{weights}, and
#' \code{targets} (a table of the positive controls) when \code{targets} is
#' given.
#'
#' @examples
#' ## Toy data: 200 genes, 2000 SNPs, 15 pathways; gene g1 is the true signal
#' set.seed(1)
#' genes <- paste0("g", 1:200)
#' magma <- data.frame(GENE = genes, P = runif(200))
#' magma$P[1:3] <- c(1e-12, 1e-6, 1e-4)
#' gwas <- data.frame(SNP = paste0("rs", 1:2000), P = runif(2000))
#' gwas$P[1:3] <- c(1e-20, 1e-8, 1e-5)          # the SNPs annotated to g1-g3
#' annot <- data.frame(gene = rep(genes, each = 10), SNP = gwas$SNP)
#' pathways <- c(list(PWY_signal = genes[1:12]),
#'               setNames(lapply(1:14, function(i) sample(genes[-(1:12)], 12)),
#'                        paste0("PWY_", 1:14)))
#' pw_res <- data.frame(pathway_id = names(pathways),
#'                      log10_p_omnibus = c(-15, log10(runif(14))))
#'
#' cand <- gene_prioritize(magma, pw_res, pathways, gwas = gwas,
#'                         gene_annot = annot, verbose = FALSE)
#' head(cand[, c("rank", "gene", "score", "gwas_rank", "magma_rank",
#'               "pathway_rank", "best_pathway")])
#'
#' \dontrun{
#' ## From files: MAGMA output, a Westfall-Young result and the set file used
#' cand <- gene_prioritize(
#'   magma           = "ldl_eur.genes.out",
#'   pathway_results = "ldl_eur_reactome_5comp_westfall_young.csv",
#'   species         = "human_reactome",             # built-in Reactome sets
#'   gwas            = "ldl_eur.neglog10.txt",       # SNP, NEGLOG10P
#'   gene_annot      = "g1000_eur.genes.annot",      # from magma --annotate
#'   gene_loc        = "NCBI37.3.gene.loc",          # for gene symbols
#'   targets         = c("LDLR", "PCSK9", "APOB", "APOE"))
#' attr(cand, "targets")
#' }
#' @export
gene_prioritize <- function(magma,
                            pathway_results,
                            pathways           = NULL,
                            species            = NULL,
                            pmn_gene_col       = NULL,
                            gwas               = NULL,
                            gene_annot         = NULL,
                            gene_loc           = NULL,
                            magma_gene_col     = "GENE",
                            magma_p_col        = NULL,
                            pathway_id_col     = "pathway_id",
                            pathway_stat_col   = NULL,
                            gwas_snp_col       = "SNP",
                            gwas_p_col         = NULL,
                            gwas_neglog10p_col = NULL,
                            gwas_chr_col       = "CHR",
                            gwas_pos_col       = "BP",
                            symbol_col         = NULL,
                            window             = 10000L,
                            w_gwas             = 1,
                            w_magma            = 1,
                            w_pathway          = 0.5,
                            w_hub              = 0.3,
                            targets            = NULL,
                            chunk              = 2000L,
                            verbose            = TRUE,
                            output             = FALSE,
                            out_file           = "gene_prioritize.csv") {
  say <- function(...) if (isTRUE(verbose)) message(sprintf(...))

  ## ---------------------------------------------------------- MAGMA layer --
  mg <- .gp_table(magma, "magma")
  if (is.null(magma_p_col)) magma_p_col <- if ("P_MULTI" %in% names(mg)) "P_MULTI" else "P"
  .gp_need_cols(mg, c(magma_gene_col, magma_p_col), "magma",
                "magma_gene_col / magma_p_col")
  res <- data.frame(gene    = as.character(mg[[magma_gene_col]]),
                    magma_p = suppressWarnings(as.numeric(mg[[magma_p_col]])),
                    stringsAsFactors = FALSE)
  res <- res[!is.na(res$magma_p) & nzchar(trimws(res$gene)), , drop = FALSE]
  res <- res[!duplicated(.gp_key(res$gene)), , drop = FALSE]
  if (!nrow(res)) stop("no genes with a MAGMA p-value.", call. = FALSE)
  key <- .gp_key(res$gene)
  res$magma_score <- .gp_rank_score(res$magma_p)
  res$magma_rank  <- .gp_rank(res$magma_p)
  say("MAGMA layer: %d genes (p column '%s')", nrow(res), magma_p_col)

  ## -------------------------------------------------------- pathway layer --
  pr <- .gp_table(pathway_results, "pathway_results")
  if (is.null(pathway_stat_col)) {
    stat_cand <- c("z_std", "omni_z_std", "log10_p_omnibus",
                   "omni_log10_p_wy_stat", "p_omnibus", "omni_p_final")
    pathway_stat_col <- NA_character_
    for (cc in stat_cand)
      if (cc %in% names(pr) && any(is.finite(suppressWarnings(as.numeric(pr[[cc]]))))) {
        pathway_stat_col <- cc; break
      }
    if (is.na(pathway_stat_col))
      stop("pathway_results has none of ", paste(stat_cand, collapse = ", "),
           " with finite values; set pathway_stat_col (smaller = stronger).",
           call. = FALSE)
  }
  .gp_need_cols(pr, c(pathway_id_col, pathway_stat_col), "pathway_results",
                "pathway_id_col / pathway_stat_col")
  pstat <- suppressWarnings(as.numeric(pr[[pathway_stat_col]]))
  pid   <- as.character(pr[[pathway_id_col]])
  ok    <- !is.na(pstat) & !is.na(pid)
  pw_stat <- stats::setNames(pstat[ok], pid[ok])
  pw_stat <- pw_stat[!duplicated(names(pw_stat))]

  if (!is.null(pathways) && !is.null(species))
    stop("Provide either 'pathways' or 'species', not both.", call. = FALSE)
  if (is.null(pathways) && is.null(species))
    stop("Provide 'pathways' (list, data.frame or set file) or 'species'.",
         call. = FALSE)
  if (is.null(pathways))
    pathways <- if (is.null(pmn_gene_col)) catfish_load_pathways(species = species)
                else catfish_load_pathways(species = species, gene_col = pmn_gene_col)
  mem <- .gp_pathways(pathways)
  mem <- mem[names(mem) %in% names(pw_stat)]
  if (!length(mem))
    stop("no pathway in 'pathways' matches a pathway_id in pathway_results.",
         call. = FALSE)
  long <- data.frame(key     = .gp_key(unlist(mem, use.names = FALSE)),
                     pathway = rep(names(mem), lengths(mem)),
                     stringsAsFactors = FALSE)
  long <- long[nzchar(long$key) & long$key %in% key, , drop = FALSE]
  long <- long[!duplicated(long), , drop = FALSE]
  if (!nrow(long))
    warning("no pathway member matches a MAGMA gene ID; check that both use ",
            "the same ID space. The pathway layer is 0 for every gene.",
            call. = FALSE)
  long$stat <- unname(pw_stat[long$pathway])
  long <- long[order(long$key, long$stat), , drop = FALSE]
  first <- !duplicated(long$key)
  best_stat <- stats::setNames(long$stat[first], long$key[first])
  best_pw   <- stats::setNames(long$pathway[first], long$key[first])
  n_pw      <- table(long$key)

  res$best_pathway <- unname(best_pw[key])
  res$pathway_stat <- unname(best_stat[key])
  res$n_pathways   <- as.integer(ifelse(is.na(n_pw[key]), 0L, n_pw[key]))
  ann <- !is.na(res$pathway_stat)
  res$pathway_score <- 0
  res$pathway_score[ann] <- .gp_rank_score(res$pathway_stat[ann])
  res$pathway_rank <- NA_real_
  res$pathway_rank[ann] <- .gp_rank(res$pathway_stat[ann])
  res$hub_bonus <- log10(res$n_pathways + 1)
  say("pathway layer: %d tested pathways, %d genes annotated (statistic '%s')",
      length(mem), sum(ann), pathway_stat_col)

  ## ----------------------------------------------------------- GWAS layer --
  res$gwas_neglog10p <- NA_real_
  route <- "none"
  loc <- NULL
  if (!is.null(gene_loc)) {
    loc <- .gp_table(gene_loc, "gene_loc",
                     header = !(is.character(gene_loc) && length(gene_loc) == 1L))
    if (!all(c("gene", "chr", "start", "end") %in% names(loc))) {
      if (ncol(loc) < 4L) stop("gene_loc needs gene, chr, start and end columns.",
                               call. = FALSE)
      names(loc)[1:4] <- c("gene", "chr", "start", "end")
    }
  }

  if (!is.null(gwas)) {
    gw <- .gp_table(gwas, "gwas")
    if (is.null(gwas_neglog10p_col))
      gwas_neglog10p_col <- .gp_pick(gw, c("NEGLOG10P", "neglog10p", "NEG_LOG10_P",
                                           "pvalue_neg_log10"))
    if (is.null(gwas_p_col))
      gwas_p_col <- .gp_pick(gw, c("P", "p", "PVAL", "pval", "P_VALUE"))
    if (!is.na(gwas_neglog10p_col) && gwas_neglog10p_col %in% names(gw)) {
      nl <- suppressWarnings(as.numeric(gw[[gwas_neglog10p_col]]))
      say("GWAS layer: ranking on -log10(p) column '%s'", gwas_neglog10p_col)
    } else if (!is.na(gwas_p_col) && gwas_p_col %in% names(gw)) {
      pp <- suppressWarnings(as.numeric(gw[[gwas_p_col]]))
      nl <- -log10(pp)
      say("GWAS layer: ranking on p column '%s' (check for floor ties)", gwas_p_col)
    } else {
      stop("gwas needs a -log10(p) or a p column; set gwas_neglog10p_col or ",
           "gwas_p_col. Columns present: ", paste(names(gw), collapse = ", "),
           call. = FALSE)
    }
    gw_ok <- is.finite(nl)

    if (!is.null(gene_annot)) {
      ## ---- route 1: SNP-to-gene annotation (MAGMA .genes.annot)
      .gp_need_cols(gw, gwas_snp_col, "gwas", "gwas_snp_col")
      if (!requireNamespace("data.table", quietly = TRUE))
        stop("Package 'data.table' is required for the gene_annot route.",
             call. = FALSE)
      route <- "gene_annot"
      pv <- data.table::data.table(SNP = as.character(gw[[gwas_snp_col]])[gw_ok],
                                   NL  = nl[gw_ok])
      pv <- pv[, list(NL = max(NL)), by = "SNP"]
      data.table::setkeyv(pv, "SNP")
      best <- NULL
      if (is.data.frame(gene_annot)) {
        .gp_need_cols(gene_annot, c("gene", "SNP"), "gene_annot", "gene_annot")
        dt <- data.table::data.table(gene = .gp_key(gene_annot$gene),
                                     SNP  = as.character(gene_annot$SNP))
        dt <- pv[dt, on = "SNP", nomatch = 0L]
        best <- dt[, list(NL = max(NL)), by = "gene"]
      } else {
        if (!(is.character(gene_annot) && length(gene_annot) == 1L &&
              file.exists(gene_annot)))
          stop("gene_annot must be a data.frame (gene, SNP) or an existing ",
               ".genes.annot file.", call. = FALSE)
        con <- file(gene_annot, "r"); on.exit(close(con), add = TRUE)
        acc <- list(); n_read <- 0L
        repeat {
          ln <- readLines(con, n = chunk, warn = FALSE)
          if (!length(ln)) break
          ln <- ln[!startsWith(ln, "#")]
          if (!length(ln)) next
          n_read <- n_read + length(ln)
          parts <- strsplit(ln, "[ \t]+")
          gid  <- vapply(parts, `[`, character(1), 1L)
          snps <- lapply(parts, function(x) if (length(x) > 2L) x[-(1:2)] else character(0))
          k <- lengths(snps) > 0L
          if (!any(k)) next
          dt <- data.table::data.table(gene = .gp_key(rep(gid[k], lengths(snps[k]))),
                                       SNP  = unlist(snps[k], use.names = FALSE))
          dt <- pv[dt, on = "SNP", nomatch = 0L]
          if (nrow(dt)) acc[[length(acc) + 1L]] <- dt[, list(NL = max(NL)), by = "gene"]
        }
        say("  annotation genes read: %d", n_read)
        if (length(acc))
          best <- data.table::rbindlist(acc)[, list(NL = max(NL)), by = "gene"]
      }
      if (!is.null(best))
        res$gwas_neglog10p <- best$NL[match(key, best$gene)]
    } else if (!is.null(loc)) {
      ## ---- route 2: SNP positions within +/- window of the gene
      .gp_need_cols(gw, c(gwas_chr_col, gwas_pos_col), "gwas",
                    "gwas_chr_col / gwas_pos_col (or give gene_annot)")
      route <- "position"
      g_chr <- as.character(gw[[gwas_chr_col]])[gw_ok]
      g_pos <- suppressWarnings(as.numeric(gw[[gwas_pos_col]]))[gw_ok]
      g_nl  <- nl[gw_ok]
      L <- data.frame(key = .gp_key(loc$gene), chr = as.character(loc$chr),
                      start = as.numeric(loc$start), end = as.numeric(loc$end),
                      stringsAsFactors = FALSE)
      best <- rep(NA_real_, nrow(L))
      for (cc in unique(L$chr)) {
        li <- which(L$chr == cc)
        si <- which(g_chr == cc & is.finite(g_pos))
        if (!length(si)) next
        o <- si[order(g_pos[si])]
        pos <- g_pos[o]; v <- g_nl[o]
        lo <- findInterval(L$start[li] - window - 0.5, pos) + 1L
        hi <- findInterval(L$end[li] + window, pos)
        for (j in seq_along(li))
          if (hi[j] >= lo[j]) best[li[j]] <- max(v[lo[j]:hi[j]])
      }
      res$gwas_neglog10p <- best[match(key, L$key)]
    } else {
      warning("gwas was given but neither gene_annot nor gene_loc, so its SNPs ",
              "cannot be assigned to genes; the GWAS layer is 0 for every gene.",
              call. = FALSE)
    }
  }
  hit <- !is.na(res$gwas_neglog10p)
  res$gwas_score <- 0
  res$gwas_rank  <- NA_real_
  if (any(hit)) {
    res$gwas_score[hit] <- .gp_rank_score(-res$gwas_neglog10p[hit])
    res$gwas_rank[hit]  <- .gp_rank(-res$gwas_neglog10p[hit])
    top <- max(res$gwas_neglog10p[hit])
    nt  <- sum(res$gwas_neglog10p[hit] >= top * (1 - 1e-9))
    if (nt > 1L)
      warning(sprintf(paste0("%d genes tie at the strongest GWAS value (-log10 p = %.4g); ",
                             "if the p-values were floored, supply -log10(p) instead."),
                      nt, top), call. = FALSE)
  }
  say("GWAS layer: %d genes with a SNP (route: %s)", sum(hit), route)

  ## ------------------------------------------------------------ composite --
  res$score <- w_gwas * res$gwas_score + w_magma * res$magma_score +
               w_pathway * res$pathway_score + w_hub * res$hub_bonus
  res$n_layers <- as.integer((res$gwas_score > 0) + (res$magma_score > 0) +
                             (res$pathway_score > 0))

  res$symbol <- NA_character_
  if (!is.null(loc)) {
    sc <- symbol_col
    if (is.null(sc)) sc <- if ("symbol" %in% names(loc)) "symbol" else if (ncol(loc) >= 6L) 6L else NA
    if (!is.na(sc[1])) {
      sv <- if (is.numeric(sc) && sc > ncol(loc)) NULL else loc[[sc]]
      if (!is.null(sv))
        res$symbol <- unname(stats::setNames(as.character(sv), .gp_key(loc$gene))[key])
    }
  }

  res <- res[order(-res$score, res$magma_p), , drop = FALSE]
  res$rank <- seq_len(nrow(res))
  res <- res[, c("rank", "gene", "symbol", "score", "n_layers",
                 "gwas_score", "magma_score", "pathway_score", "hub_bonus",
                 "gwas_neglog10p", "gwas_rank", "magma_p", "magma_rank",
                 "best_pathway", "pathway_stat", "pathway_rank", "n_pathways")]
  rownames(res) <- NULL

  attr(res, "layer_sizes") <- c(gwas = sum(hit), magma = nrow(res),
                                pathway = sum(ann))
  attr(res, "pathway_stat_col") <- pathway_stat_col
  attr(res, "gwas_route") <- route
  attr(res, "weights") <- c(gwas = w_gwas, magma = w_magma,
                            pathway = w_pathway, hub = w_hub)

  if (length(targets)) {
    tk <- .gp_key(targets)
    i  <- match(tk, .gp_key(res$symbol))
    j  <- match(tk, .gp_key(res$gene))
    i[is.na(i)] <- j[is.na(i)]
    attr(res, "targets") <- data.frame(
      target       = targets,
      rank         = res$rank[i],
      gwas_rank    = res$gwas_rank[i],
      magma_rank   = res$magma_rank[i],
      pathway_rank = res$pathway_rank[i],
      best_pathway = res$best_pathway[i],
      stringsAsFactors = FALSE)
    if (isTRUE(verbose)) {
      tt <- attr(res, "targets")
      say("positive controls found: %d of %d; median overall rank %s of %d",
          sum(!is.na(tt$rank)), nrow(tt),
          format(stats::median(tt$rank, na.rm = TRUE)), nrow(res))
    }
  }

  if (isTRUE(output)) {
    d <- dirname(out_file)
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(res, out_file, row.names = FALSE)
    attr(res, "file") <- out_file
  }
  res
}
