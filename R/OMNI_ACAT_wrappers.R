# #' Omnibus (ACAT-O) pathway test
# #'
# #' Two usage modes:
# #'
# #' \strong{1) From MAGMA gene-level results (full pipeline)}\cr
# #'   - Provide \code{gene_results} + \code{pathways} or \code{species}. \cr
# #'   - For each pathway, compute:
# #'       * ACAT on gene p-values,
# #'       * wFisher (optionally weighted by \code{weight_col}),
# #'       * truncated Fisher TPM (ptrunc),
# #'     then combine them with ACAT-O (via \code{sumFREGAT::ACATO}).
# #'   - Optionally perform gene-set permutations for empirical omnibus p-values.
# #'
# #' \strong{2) From an existing table of pathway-level p-values}\cr
# #'   - Provide \code{tab} (data.frame or CSV path), \code{id_col}, and \code{p_cols}. \cr
# #'   - For each row (pathway), apply ACAT-O across the p-columns in \code{p_cols}.
# #'
# #' Multiple-testing correction:
# #'   - BH FDR (Benjamini–Hochberg) is computed for omni and component p-values.
# #'   - If the \code{qvalue} package is installed, q-values are also reported.
# #'
# #' @param gene_results data.frame with at least gene + p-value columns.
# #'   Typically a MAGMA \code{.genes.out} file, with "GENE" and "P" (and "ZSTAT").
# #'   Only used when \code{tab = NULL}.
# #' @param pathways either:
# #'   \itemize{
# #'     \item a named list: each element is a character vector of gene IDs
# #'           in that pathway, OR
# #'     \item a data.frame with columns \code{pathway_id}, \code{gene_id}
# #'           (and optional \code{pathway_name}).
# #'   }
# #'   Only used when \code{tab = NULL}.
# #' @param species optional; one of "maize", "sorghum", "arabidopsis", "plant".
# #'   If provided, built-in PMN pathways are loaded via \code{magcat_load_pathways()}.
# #'   You must provide either \code{pathways} OR \code{species}, but not both.
# #'   Only used when \code{tab = NULL}.
# #' @param pmn_gene_col optional; passed to \code{magcat_load_pathways(gene_col=...)}
# #'   when \code{species} is used.
# #' @param gene_col column name in \code{gene_results} containing gene IDs
# #'   (default "GENE").
# #' @param p_col column name in \code{gene_results} containing gene-level p-values
# #'   (default "P").
# #' @param effect_col column with effect size / Z-statistic to derive direction
# #'   (default "ZSTAT"). Only used when \code{is_onetail = FALSE}.
# #' @param weight_col optional column in \code{gene_results} to use as weight
# #'   for wFisher (e.g. "NSNPS"). If NULL, all weights = 1.
# #' @param is_onetail logical; passed to wFisher (\code{is.onetail=...}).
# #'   Default FALSE (two-sided p-values + directions).
# #' @param min_p lower cap for very small p-values (default 1e-15).
# #' @param do_fix logical; if TRUE, clean p-values with \code{fix_p_for_acat()}
# #'   (0 -> \code{min_p}, 1 -> \code{1 - 1/d}).
# #' @param ptrunc truncation threshold for truncated Fisher TPM (default 0.05).
# #' @param B_perm Integer; number of gene-set permutations for empirical
# #'   omnibus p-values. Default 0 (no permutations). Only used when \code{tab = NULL}.
# #' @param seed optional integer random seed for permutations (ignored if
# #'   \code{B_perm == 0}).
# #' @param output logical; if TRUE, write a CSV of results to \code{out_dir}.
# #' @param out_dir directory to write CSV when \code{output = TRUE}
# #'   (default "magcat_omni").
# #'
# #' @param tab Optional: if non-NULL, use table mode. Either a data.frame
# #'   or a single CSV file path with columns including \code{id_col} and \code{p_cols}.
# #' @param id_col Name of pathway ID column in \code{tab} (default "pathway_id").
# #' @param p_cols Character vector of column names in \code{tab} that contain
# #'   p-values to be combined by ACAT-O (length >= 1).
# #'
# #' @return
# #' \strong{Table mode (\code{tab} not NULL)}: \cr
# #'   data.frame with \code{id_col}, \code{p_cols}, \code{omni_p},
# #'   and BH/qvalue-adjusted omni p-values. \cr
# #'
# #' \strong{Pipeline mode (\code{tab} = NULL)}: \cr
# #'   data.frame with columns:
# #'   \describe{
# #'     \item{pathway_id}{pathway identifier}
# #'     \item{pathway_name}{pathway name}
# #'     \item{n_genes}{number of genes used}
# #'     \item{gene_names}{semicolon-separated gene IDs}
# #'     \item{acat_p}{ACAT pathway p-value}
# #'     \item{wfisher_p}{wFisher pathway p-value}
# #'     \item{tpm_p}{truncated Fisher (TPM) pathway p-value (chi-square approx)}
# #'     \item{omni_p}{ACAT-O combined p-value of the three methods}
# #'     \item{omni_perm_p}{empirical omnibus p-value from permutations
# #'                        (NA if \code{B_perm == 0})}
# #'     \item{acat_p_BH}{BH FDR for \code{acat_p}}
# #'     \item{wfisher_p_BH}{BH FDR for \code{wfisher_p}}
# #'     \item{tpm_p_BH}{BH FDR for \code{tpm_p}}
# #'     \item{omni_p_BH}{BH FDR for \code{omni_p}}
# #'     \item{omni_perm_p_BH}{BH FDR for \code{omni_perm_p} (if \code{B_perm > 0})}
# #'     \item{acat_p_q}{qvalue for \code{acat_p} (if \code{qvalue} installed)}
# #'     \item{wfisher_p_q}{qvalue for \code{wfisher_p}}
# #'     \item{tpm_p_q}{qvalue for \code{tpm_p}}
# #'     \item{omni_p_q}{qvalue for \code{omni_p}}
# #'     \item{omni_perm_p_q}{qvalue for \code{omni_perm_p} (if \code{B_perm > 0})}
# #'   }
# #'
# #' @export
# magcat_omni_pathways <- function(gene_results = NULL,
#                                  pathways     = NULL,
#                                  species      = NULL,
#                                  pmn_gene_col = NULL,
#                                  gene_col     = "GENE",
#                                  p_col        = "P",
#                                  effect_col   = "ZSTAT",
#                                  weight_col   = NULL,
#                                  is_onetail   = FALSE,
#                                  min_p        = 1e-15,
#                                  do_fix       = TRUE,
#                                  ptrunc       = 0.05,
#                                  B_perm       = 0L,
#                                  seed         = NULL,
#                                  output       = FALSE,
#                                  out_dir      = "magcat_omni",
#                                  tab          = NULL,
#                                  id_col       = "pathway_id",
#                                  p_cols       = NULL) {

#   ## small helpers for BH and qvalue that respect NAs
#   .p_adjust_BH <- function(p) {
#     p <- as.numeric(p)
#     out <- rep(NA_real_, length(p))
#     idx <- which(!is.na(p))
#     if (length(idx) > 0L) {
#       out[idx] <- stats::p.adjust(p[idx], method = "BH")
#     }
#     out
#   }

#   .qvalue_vec <- function(p) {
#     p <- as.numeric(p)
#     out <- rep(NA_real_, length(p))
#     idx <- which(!is.na(p) & p > 0 & p <= 1)
#     if (length(idx) > 1L && requireNamespace("qvalue", quietly = TRUE)) {
#       qq <- tryCatch(
#         qvalue::qvalue(p[idx])$qvalues,
#         error = function(e) rep(NA_real_, length(idx))
#       )
#       out[idx] <- qq
#     }
#     out
#   }

#   ## ==================== MODE 2: TABLE OF P-VALUES =====================
#   if (!is.null(tab)) {
#     if (!requireNamespace("sumFREGAT", quietly = TRUE)) {
#       stop("magcat_omni_pathways(): table mode requires sumFREGAT (ACATO).",
#            call. = FALSE)
#     }
#     if (missing(p_cols) || length(p_cols) == 0L) {
#       stop("magcat_omni_pathways(): in table mode you must supply 'p_cols'.",
#            call. = FALSE)
#     }

#     if (is.character(tab) && length(tab) == 1L) {
#       if (!file.exists(tab)) {
#         stop("magcat_omni_pathways(): file not found: ", tab, call. = FALSE)
#       }
#       tab <- utils::read.csv(tab, stringsAsFactors = FALSE, check.names = FALSE)
#     } else if (!is.data.frame(tab)) {
#       stop("magcat_omni_pathways(): 'tab' must be a data.frame or a single CSV path.",
#            call. = FALSE)
#     }

#     if (!id_col %in% names(tab)) {
#       stop("magcat_omni_pathways(): id_col = '", id_col, "' not found in 'tab'.",
#            call. = FALSE)
#     }
#     missing_p <- setdiff(p_cols, names(tab))
#     if (length(missing_p)) {
#       stop("magcat_omni_pathways(): p_cols not found in 'tab': ",
#            paste(missing_p, collapse = ", "),
#            call. = FALSE)
#     }

#     res <- tab[, c(id_col, p_cols), drop = FALSE]

#     P_mat <- as.matrix(res[, p_cols, drop = FALSE])
#     n_pw  <- nrow(P_mat)
#     omni  <- rep(NA_real_, n_pw)

#     for (i in seq_len(n_pw)) {
#       p_vec <- as.numeric(P_mat[i, ])
#       p_vec <- p_vec[!is.na(p_vec)]
#       if (!length(p_vec)) {
#         omni[i] <- NA_real_
#       } else {
#         p_vec   <- fix_p_for_acat(p_vec, min_p = min_p)
#         omni[i] <- as.numeric(sumFREGAT::ACATO(p_vec))
#       }
#     }

#     res$omni_p    <- omni
#     res$omni_p_BH <- .p_adjust_BH(res$omni_p)
#     res$omni_p_q  <- .qvalue_vec(res$omni_p)

#     ord <- order(res$omni_p, decreasing = FALSE, na.last = TRUE)
#     res <- res[ord, , drop = FALSE]

#     if (output) {
#       if (!dir.exists(out_dir)) {
#         dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#       }
#       out_path <- file.path(out_dir, "magcat_omni_from_pvals.csv")
#       utils::write.csv(res, out_path, row.names = FALSE)
#       attr(res, "file") <- out_path
#     }

#     return(res)
#   }

#   ## ==================== MODE 1: FULL PIPELINE =========================

#   if (is.null(gene_results)) {
#     stop("magcat_omni_pathways(): in pipeline mode you must provide 'gene_results' ",
#          "and either 'pathways' or 'species'.",
#          call. = FALSE)
#   }

#   if (!requireNamespace("ACAT", quietly = TRUE)) {
#     stop("Package 'ACAT' is required for magcat_omni_pathways() (gene-level ACAT).",
#          call. = FALSE)
#   }
#   if (!requireNamespace("sumFREGAT", quietly = TRUE)) {
#     stop("Package 'sumFREGAT' is required for ACAT-O combination (ACATO).",
#          call. = FALSE)
#   }
#   if (!exists("wFisher")) {
#     stop("Function 'wFisher' not found. Make sure it is defined/loaded.",
#          call. = FALSE)
#   }

#   ## -------- decide pathway source: pathways vs species ----------
#   if (!is.null(pathways) && !is.null(species)) {
#     stop("Provide either 'pathways' OR 'species', not both.", call. = FALSE)
#   }
#   if (is.null(pathways) && is.null(species)) {
#     stop("You must provide either:\n",
#          "  * 'pathways' (list/data.frame), OR\n",
#          "  * 'species' = 'maize' | 'sorghum' | 'arabidopsis' | 'plant'.",
#          call. = FALSE)
#   }

#   if (is.null(pathways) && !is.null(species)) {
#     if (is.null(pmn_gene_col)) {
#       pathways <- magcat_load_pathways(species = species)
#     } else {
#       pathways <- magcat_load_pathways(
#         species  = species,
#         gene_col = pmn_gene_col
#       )
#     }
#   }

#   ## -------- standardize gene_results ----------
#   needed_cols <- c(gene_col, p_col)
#   if (!all(needed_cols %in% names(gene_results))) {
#     stop("gene_results must contain columns '", gene_col,
#          "' and '", p_col, "'.",
#          call. = FALSE)
#   }
#   if (!is_onetail && !effect_col %in% names(gene_results)) {
#     stop("effect_col = '", effect_col,
#          "' not found in gene_results. Needed when is_onetail = FALSE.",
#          call. = FALSE)
#   }

#   gr        <- gene_results
#   genes_all <- as.character(gr[[gene_col]])
#   p_all     <- gr[[p_col]]

#   genes_all_norm <- tolower(genes_all)

#   # effect signs if needed
#   if (!is_onetail) {
#     eff_all <- gr[[effect_col]]
#     eff_vec <- stats::setNames(eff_all, genes_all_norm)
#   } else {
#     eff_vec <- NULL
#   }

#   # weights if requested
#   if (!is.null(weight_col)) {
#     if (!weight_col %in% names(gene_results)) {
#       stop("weight_col = '", weight_col,
#            "' not found in gene_results.",
#            call. = FALSE)
#     }
#     w_all <- gr[[weight_col]]
#     w_vec <- stats::setNames(w_all, genes_all_norm)
#   } else {
#     w_vec <- NULL
#   }

#   ## -------- pathways -> named list + names ----------
#   if (is.data.frame(pathways)) {
#     if (!"pathway_id" %in% names(pathways) ||
#         !"gene_id"    %in% names(pathways)) {
#       stop("If 'pathways' is a data.frame it must have columns ",
#            "'pathway_id' and 'gene_id'.",
#            call. = FALSE)
#     }
#     if (!"pathway_name" %in% names(pathways)) {
#       pathways$pathway_name <- pathways$pathway_id
#     }

#     p_list  <- split(pathways$gene_id, pathways$pathway_id)
#     p_names <- tapply(
#       pathways$pathway_name,
#       pathways$pathway_id,
#       FUN = function(x) x[1]
#     )
#   } else if (is.list(pathways)) {
#     p_list  <- pathways
#     p_names <- names(pathways)
#     if (is.null(p_names)) {
#       p_names <- paste0("PWY_", seq_along(pathways))
#       names(p_list) <- p_names
#     }
#   } else {
#     stop("'pathways' must be either a list or a data.frame.",
#          call. = FALSE)
#   }

#   # normalize gene IDs in pathways
#   p_list <- lapply(p_list, function(g) tolower(as.character(g)))

#   ## -------- result container ----------
#   n_pw <- length(p_list)
#   res  <- data.frame(
#     pathway_id      = names(p_list),
#     pathway_name    = if (is.null(p_names)) names(p_list) else unname(p_names),
#     n_genes         = NA_integer_,
#     gene_names      = NA_character_,
#     acat_p          = NA_real_,
#     wfisher_p       = NA_real_,
#     tpm_p           = NA_real_,
#     omni_p          = NA_real_,
#     omni_perm_p     = NA_real_,
#     acat_p_BH       = NA_real_,
#     wfisher_p_BH    = NA_real_,
#     tpm_p_BH        = NA_real_,
#     omni_p_BH       = NA_real_,
#     omni_perm_p_BH  = NA_real_,
#     acat_p_q        = NA_real_,
#     wfisher_p_q     = NA_real_,
#     tpm_p_q         = NA_real_,
#     omni_p_q        = NA_real_,
#     omni_perm_p_q   = NA_real_,
#     stringsAsFactors = FALSE
#   )

#   ## -------- helper: per-pathway stats (obs or perm) --------------
#   compute_pathway_stats <- function(idx) {
#     p_i <- p_all[idx]

#     # drop NA p's
#     keep <- !is.na(p_i)
#     p_i  <- p_i[keep]
#     idx  <- idx[keep]

#     d <- length(p_i)
#     if (d == 0L) {
#       return(list(
#         acat_p    = NA_real_,
#         wfisher_p = NA_real_,
#         tpm_p     = NA_real_,
#         omni_p    = NA_real_
#       ))
#     }

#     # clean gene-level p-values
#     if (do_fix) {
#       p_i <- fix_p_for_acat(p_i, min_p = min_p)
#     }

#     # directions
#     if (!is_onetail) {
#       eff_i <- eff_vec[genes_all_norm[idx]]
#       sgn   <- sign(eff_i)
#       sgn[is.na(sgn) | sgn == 0] <- 1
#     } else {
#       sgn <- NULL
#     }

#     # weights
#     if (!is.null(w_vec)) {
#       w_i <- w_vec[genes_all_norm[idx]]
#       if (any(is.na(w_i) | w_i <= 0)) {
#         w_pos <- w_i[!is.na(w_i) & w_i > 0]
#         repl  <- if (length(w_pos)) stats::median(w_pos) else 1
#         w_i[is.na(w_i) | w_i <= 0] <- repl
#       }
#     } else {
#       w_i <- rep(1, d)
#     }

#     ## --- ACAT (gene-level) ---
#     p_acat <- ACAT::ACAT(Pvals = p_i)

#     ## --- wFisher ---
#     if (!is_onetail) {
#       out_wf <- wFisher(
#         p          = p_i,
#         weight     = w_i,
#         is.onetail = FALSE,
#         eff.sign   = sgn
#       )
#     } else {
#       out_wf <- wFisher(
#         p          = p_i,
#         weight     = w_i,
#         is.onetail = TRUE
#       )
#     }
#     p_wf <- out_wf$p

#     ## --- truncated Fisher (TPM) ---
#     tpm <- magcat_tpm_stat(p_i, ptrunc = ptrunc)
#     p_tpm <- tpm$pval

#     ## --- ACAT-O on the three pathway p's ---
#     p_vec <- c(p_acat, p_wf, p_tpm)
#     p_vec <- p_vec[!is.na(p_vec)]
#     if (!length(p_vec)) {
#       p_omni <- NA_real_
#     } else {
#       p_vec  <- fix_p_for_acat(p_vec, min_p = min_p)
#       p_omni <- as.numeric(sumFREGAT::ACATO(p_vec))
#     }

#     list(
#       acat_p    = p_acat,
#       wfisher_p = p_wf,
#       tpm_p     = p_tpm,
#       omni_p    = p_omni
#     )
#   }

#   ## -------- observed stats for each pathway ----------------------
#   idx_pool <- which(!is.na(p_all))

#   for (i in seq_len(n_pw)) {
#     genes_i_norm <- p_list[[i]]

#     idx_i <- match(genes_i_norm, genes_all_norm)
#     idx_i <- idx_i[!is.na(idx_i)]

#     d <- length(idx_i)
#     res$n_genes[i] <- d
#     if (d == 0L) {
#       next
#     }

#     canon_ids <- unique(genes_all[idx_i])
#     res$gene_names[i] <- paste(canon_ids, collapse = ";")

#     stats_i <- compute_pathway_stats(idx_i)

#     res$acat_p[i]    <- stats_i$acat_p
#     res$wfisher_p[i] <- stats_i$wfisher_p
#     res$tpm_p[i]     <- stats_i$tpm_p
#     res$omni_p[i]    <- stats_i$omni_p
#   }

#   ## -------- permutations for omnibus p (optional) ----------------
#   if (B_perm > 0L) {
#     if (!is.null(seed)) {
#       set.seed(seed)
#     }

#     n_pool <- length(idx_pool)
#     if (n_pool == 0L) {
#       stop("No non-NA gene-level p-values available for permutations.",
#            call. = FALSE)
#     }

#     omni_perm_p <- rep(NA_real_, n_pw)

#     for (i in seq_len(n_pw)) {
#       d <- res$n_genes[i]
#       if (is.na(d) || d == 0L || is.na(res$omni_p[i])) {
#         next
#       }

#       omni_perm <- numeric(B_perm)

#       for (b in seq_len(B_perm)) {
#         idx_perm <- sample(idx_pool, size = d, replace = FALSE)
#         stats_b  <- compute_pathway_stats(idx_perm)
#         omni_perm[b] <- stats_b$omni_p
#       }

#       omni_obs <- res$omni_p[i]
#       omni_perm_p[i] <- (1 + sum(omni_perm <= omni_obs, na.rm = TRUE)) /
#                         (B_perm + 1)
#     }

#     res$omni_perm_p <- omni_perm_p
#   }

#   ## -------- multiple-testing correction (BH + qvalue) -------------
#   res$acat_p_BH      <- .p_adjust_BH(res$acat_p)
#   res$wfisher_p_BH   <- .p_adjust_BH(res$wfisher_p)
#   res$tpm_p_BH       <- .p_adjust_BH(res$tpm_p)
#   res$omni_p_BH      <- .p_adjust_BH(res$omni_p)
#   if (B_perm > 0L) {
#     res$omni_perm_p_BH <- .p_adjust_BH(res$omni_perm_p)
#   }

#   res$acat_p_q       <- .qvalue_vec(res$acat_p)
#   res$wfisher_p_q    <- .qvalue_vec(res$wfisher_p)
#   res$tpm_p_q        <- .qvalue_vec(res$tpm_p)
#   res$omni_p_q       <- .qvalue_vec(res$omni_p)
#   if (B_perm > 0L) {
#     res$omni_perm_p_q <- .qvalue_vec(res$omni_perm_p)
#   }

#   ## -------- sort by omni_p (ascending) ---------------------------
#   ord <- order(res$omni_p, decreasing = FALSE, na.last = TRUE)
#   res <- res[ord, , drop = FALSE]

#   ## -------- optional CSV output ----------------------------------
#   if (output) {
#     if (!dir.exists(out_dir)) {
#       dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#     }
#     species_tag <- if (is.null(species)) "custom" else species
#     out_path <- file.path(
#       out_dir,
#       paste0("magcat_omni_pathways_", species_tag, ".csv")
#     )
#     utils::write.csv(res, out_path, row.names = FALSE)
#     attr(res, "file") <- out_path
#   }

#   res
# }
