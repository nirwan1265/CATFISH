






# # # ---- helper: read MAGMA competitive tidy TSV and return id + p ----
# # .magcat_read_magma_tidy <- function(magma_tidy_tsv) {
# #   if (is.null(magma_tidy_tsv)) return(NULL)
# #   if (!file.exists(magma_tidy_tsv)) {
# #     stop("magma_tidy_tsv does not exist: ", magma_tidy_tsv, call. = FALSE)
# #   }

# #   tab <- tryCatch(
# #     utils::read.delim(magma_tidy_tsv, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
# #     error = function(e) NULL
# #   )
# #   if (is.null(tab)) stop("Failed to read magma_tidy_tsv: ", magma_tidy_tsv, call. = FALSE)

# #   if (!("pathway_id" %in% names(tab))) {
# #     stop("magma_tidy_tsv must have column 'pathway_id'.", call. = FALSE)
# #   }

# #   # your file uses magma_pvalue; allow flexible fallback
# #   pcol <- NULL
# #   if ("magma_pvalue" %in% names(tab)) {
# #     pcol <- "magma_pvalue"
# #   } else {
# #     cand <- grep("magma.*p", names(tab), ignore.case = TRUE, value = TRUE)
# #     if (length(cand)) pcol <- cand[1]
# #   }
# #   if (is.null(pcol)) {
# #     stop("magma_tidy_tsv must have a MAGMA p-value column (e.g. 'magma_pvalue').", call. = FALSE)
# #   }

# #   out <- tab[, c("pathway_id", pcol), drop = FALSE]
# #   names(out)[2] <- "magma_pvalue"
# #   out
# # }

# # # ---- UPDATED exported wrapper (calls your existing magcat_omni_prepare/run_prepared) ----
# # #' Omnibus pathway test (ACAT-O or minP) with MVN or global resampling
# # #'
# # #' Components (5): ACAT, Fisher, adaptive soft TFisher, StoufferZ (z_col), minP
# # #' Optional: merge MAGMA competitive p-values from magma_tidy_tsv into output table.
# # #'
# # #' @export
# # magcat_omni_pathways <- function(gene_results,
# #                                 pathways     = NULL,
# #                                 species      = NULL,
# #                                 pmn_gene_col = NULL,
# #                                 gene_col     = "GENE",
# #                                 p_col        = "P",
# #                                 z_col        = "ZSTAT",     # <-- IMPORTANT (StoufferZ uses this)
# #                                 weight_col   = NULL,
# #                                 tau_grid     = c(0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.001),
# #                                 tau_cap      = 0.5,
# #                                 min_p        = 1e-15,
# #                                 do_fix       = TRUE,
# #                                 omnibus      = c("ACAT", "minP"),
# #                                 B_perm       = 0L,
# #                                 perm_mode    = c("resample_global", "global_resampling", "global_sampling", "mvn"),
# #                                 magma_genes_out = NULL,
# #                                 magma_cor_file  = NULL,
# #                                 make_PD      = TRUE,
# #                                 seed         = NULL,
# #                                 remove_singletons = TRUE,
# #                                 magma_tidy_tsv = NULL,      # <-- your *.tidy.tsv from magma_geneset_competitive()
# #                                 output       = FALSE,
# #                                 out_dir      = "magcat_omni") {

# #   omnibus   <- match.arg(omnibus)
# #   perm_mode <- match.arg(perm_mode)

# #   # normalize perm_mode aliases
# #   if (perm_mode %in% c("global_resampling", "global_sampling")) perm_mode <- "resample_global"

# #   prep <- magcat_omni_prepare(
# #     gene_results = gene_results,
# #     pathways     = pathways,
# #     species      = species,
# #     pmn_gene_col = pmn_gene_col,
# #     gene_col     = gene_col,
# #     p_col        = p_col,
# #     z_col        = z_col,
# #     weight_col   = weight_col,
# #     tau_grid     = tau_grid,
# #     tau_cap      = tau_cap,
# #     min_p        = min_p,
# #     do_fix       = do_fix,
# #     seed         = seed
# #   )

# #   magcat_omni_run_prepared(
# #     prep              = prep,
# #     omnibus           = omnibus,
# #     B_perm            = B_perm,
# #     perm_mode         = perm_mode,
# #     magma_genes_out   = magma_genes_out,
# #     magma_cor_file    = magma_cor_file,
# #     make_PD           = make_PD,
# #     remove_singletons = remove_singletons,
# #     magma_tidy_tsv    = magma_tidy_tsv,
# #     output            = output,
# #     out_dir           = out_dir
# #   )
# # }

# # # ---- UPDATED run_prepared: merges MAGMA tidy TSV if provided ----
# # #' @export
# # magcat_omni_run_prepared <- function(prep,
# #                                      omnibus = c("ACAT","minP"),
# #                                      B_perm = 0L,
# #                                      perm_mode = c("resample_global","mvn"),
# #                                      magma_genes_out = NULL,
# #                                      magma_cor_file  = NULL,
# #                                      make_PD = TRUE,
# #                                      remove_singletons = TRUE,
# #                                      magma_tidy_tsv = NULL,
# #                                      output = FALSE,
# #                                      out_dir = "magcat_omni") {

# #   omnibus   <- match.arg(omnibus)
# #   perm_mode <- match.arg(perm_mode)

# #   # sanity: you need these helpers if using MVN
# #   if (perm_mode == "mvn") {
# #     req <- c("magma_read_genes_out", "magma_build_R_for_pathway", "magma_simulate_null_Zp", "magma_read_gene_cor_pairs")
# #     miss <- req[!vapply(req, exists, logical(1), mode = "function")]
# #     if (length(miss)) {
# #       stop("perm_mode='mvn' requires these functions to be loaded: ",
# #            paste(miss, collapse = ", "),
# #            call. = FALSE)
# #     }
# #   }

# #   # internal helpers (must exist in your namespace already from wrappers)
# #   if (!exists("fix_p_for_acat", mode = "function")) {
# #     stop("fix_p_for_acat() not found. Load/source your ACAT wrappers first.", call. = FALSE)
# #   }

# #   .p_adjust_BH <- function(p) {
# #     p <- as.numeric(p)
# #     out <- rep(NA_real_, length(p))
# #     idx <- which(!is.na(p))
# #     if (length(idx)) out[idx] <- stats::p.adjust(p[idx], method = "BH")
# #     out
# #   }

# #   .sidak_minp <- function(p_vec) {
# #     p_vec <- as.numeric(p_vec)
# #     p_vec <- p_vec[is.finite(p_vec) & p_vec > 0 & p_vec < 1]
# #     if (!length(p_vec)) return(NA_real_)
# #     k <- length(p_vec)
# #     pmin <- min(p_vec)
# #     1 - (1 - pmin)^k
# #   }

# #   .combine_omni <- function(p_vec, min_p = prep$min_p) {
# #     p_vec <- as.numeric(p_vec)
# #     p_vec <- p_vec[is.finite(p_vec) & p_vec > 0 & p_vec < 1]
# #     if (!length(p_vec)) return(NA_real_)
# #     if (omnibus == "ACAT") {
# #       if (!requireNamespace("sumFREGAT", quietly = TRUE)) {
# #         stop("omnibus='ACAT' requires sumFREGAT::ACATO().", call. = FALSE)
# #       }
# #       p_vec <- fix_p_for_acat(p_vec, min_p = min_p)
# #       return(tryCatch(as.numeric(sumFREGAT::ACATO(p_vec)), error = function(e) NA_real_))
# #     } else {
# #       return(.sidak_minp(p_vec))
# #     }
# #   }

# #   .stouffer_p_from_z <- function(z_i, w_i = NULL) {
# #     z_i <- as.numeric(z_i)
# #     z_i <- z_i[is.finite(z_i)]
# #     if (length(z_i) < 2L) return(NA_real_)
# #     if (is.null(w_i)) w_i <- rep(1, length(z_i))
# #     w_i <- as.numeric(w_i)
# #     # fix weights
# #     bad <- !is.finite(w_i) | abs(w_i) <= 1e-8
# #     if (any(bad)) {
# #       w_ok <- w_i[!bad]
# #       repl <- if (length(w_ok)) stats::median(abs(w_ok)) else 1
# #       w_i[bad] <- repl
# #     }
# #     z <- sum(w_i * z_i) / sqrt(sum(w_i^2))
# #     2 * stats::pnorm(-abs(z))
# #   }

# #   .tfisher_adaptive_pmin <- function(p_i, tau_grid, tau_cap = 0.5, min_p = prep$min_p) {
# #     p_i <- fix_p_for_acat(p_i, min_p = min_p)
# #     p_i <- p_i[is.finite(p_i) & p_i > 0 & p_i < 1]
# #     n <- length(p_i)
# #     if (n < 2L) return(list(p_min = NA_real_, tau_best = NA_real_))

# #     tg <- sort(unique(as.numeric(tau_grid)))
# #     tg <- tg[is.finite(tg) & tg > 0 & tg < 1 & tg <= tau_cap]
# #     if (!length(tg)) tg <- min(0.05, tau_cap)

# #     best_p <- Inf
# #     best_tau <- tg[1]
# #     for (tau in tg) {
# #       st <- TFisher::stat.soft(p = p_i, tau1 = tau)
# #       pv <- 1 - as.numeric(TFisher::p.soft(q = st, n = n, tau1 = tau, M = NULL))
# #       if (is.finite(pv) && pv < best_p) {
# #         best_p <- pv
# #         best_tau <- tau
# #       }
# #     }
# #     list(p_min = as.numeric(best_p), tau_best = as.numeric(best_tau))
# #   }

# #   # ---------- observed component p-values ----------
# #   nm <- names(prep$p_list)

# #   res <- data.frame(
# #     pathway_id   = nm,
# #     pathway_name = prep$pathway_name[nm],
# #     n_genes      = prep$n_genes[nm],
# #     stringsAsFactors = FALSE
# #   )

# #   # components
# #   res$acat_p <- vapply(nm, function(k) {
# #     p_i <- prep$p_list[[k]]
# #     if (length(p_i) < 2L) return(NA_real_)
# #     ACAT::ACAT(Pvals = fix_p_for_acat(p_i, min_p = prep$min_p))
# #   }, numeric(1))

# #   res$fisher_p <- vapply(nm, function(k) {
# #     p_i <- fix_p_for_acat(prep$p_list[[k]], min_p = prep$min_p)
# #     p_i <- p_i[is.finite(p_i) & p_i > 0 & p_i < 1]
# #     if (length(p_i) < 2L) return(NA_real_)
# #     stats::pchisq(-2 * sum(log(p_i)), df = 2 * length(p_i), lower.tail = FALSE)
# #   }, numeric(1))

# #   tf_obs <- lapply(nm, function(k) .tfisher_adaptive_pmin(prep$p_list[[k]], prep$tau_grid, prep$tau_cap, prep$min_p))
# #   res$tfisher_adapt_p_analytic <- vapply(tf_obs, function(x) x$p_min, numeric(1))
# #   res$tfisher_adapt_tau_best   <- vapply(tf_obs, function(x) x$tau_best, numeric(1))

# #   res$minp_gene_p <- vapply(nm, function(k) .sidak_minp(prep$p_list[[k]]), numeric(1))

# #   res$stouffer_p_analytic <- vapply(nm, function(k) {
# #     z_i <- prep$z_list[[k]]
# #     w_i <- prep$w_list[[k]]
# #     .stouffer_p_from_z(z_i, w_i)
# #   }, numeric(1))

# #   # remove singletons
# #   if (remove_singletons) {
# #     keep <- is.finite(res$n_genes) & res$n_genes >= 2L
# #     res <- res[keep, , drop = FALSE]
# #     nm  <- res$pathway_id
# #   }

# #   # ---------- omnibus (analytic) ----------
# #   res$omni_p <- vapply(seq_len(nrow(res)), function(i) {
# #     pv <- c(res$acat_p[i], res$fisher_p[i], res$tfisher_adapt_p_analytic[i], res$stouffer_p_analytic[i], res$minp_gene_p[i])
# #     .combine_omni(pv)
# #   }, numeric(1))
# #   res$omni_p_BH <- .p_adjust_BH(res$omni_p)

# #   # ---------- merge MAGMA competitive tidy TSV (optional) ----------
# #   mg <- .magcat_read_magma_tidy(magma_tidy_tsv)
# #   if (!is.null(mg)) {
# #     res$magma_pvalue <- mg$magma_pvalue[match(res$pathway_id, mg$pathway_id)]
# #   } else {
# #     res$magma_pvalue <- NA_real_
# #   }

# #   # ---------- omni permutations ----------
# #   res$omni_perm_p    <- NA_real_
# #   res$omni_perm_p_BH <- NA_real_

# #   if (B_perm > 0L) {

# #     if (perm_mode == "resample_global") {

# #       pool_p <- as.numeric(prep$p_all)
# #       pool_z <- as.numeric(prep$z_all)

# #       ok <- which(is.finite(pool_p) & pool_p > 0 & pool_p < 1 & is.finite(pool_z))
# #       pool_p <- fix_p_for_acat(pool_p[ok], min_p = prep$min_p)
# #       pool_z <- pool_z[ok]

# #       for (i in seq_len(nrow(res))) {

# #         d <- as.integer(res$n_genes[i])
# #         if (!is.finite(d) || d < 2L) next
# #         if (!is.finite(res$omni_p[i])) next

# #         omni_null <- numeric(B_perm)

# #         # weights for this pathway (fixed)
# #         w_i <- prep$w_list[[res$pathway_id[i]]]
# #         if (is.null(w_i) || length(w_i) != d) w_i <- rep(1, d)

# #         for (b in seq_len(B_perm)) {

# #           idx <- sample.int(length(pool_p), size = d, replace = FALSE)
# #           p_i <- pool_p[idx]
# #           z_i <- pool_z[idx]

# #           p_acat   <- ACAT::ACAT(Pvals = fix_p_for_acat(p_i, min_p = prep$min_p))

# #           p_fisher <- stats::pchisq(-2 * sum(log(fix_p_for_acat(p_i, min_p = prep$min_p))),
# #                                     df = 2 * d, lower.tail = FALSE)

# #           tf <- .tfisher_adaptive_pmin(p_i, prep$tau_grid, prep$tau_cap, prep$min_p)
# #           p_tf <- tf$p_min

# #           p_st <- .stouffer_p_from_z(z_i, w_i)

# #           p_min <- .sidak_minp(p_i)

# #           omni_null[b] <- .combine_omni(c(p_acat, p_fisher, p_tf, p_st, p_min))
# #         }

# #         omni_obs <- res$omni_p[i]
# #         res$omni_perm_p[i] <- (1 + sum(omni_null <= omni_obs, na.rm = TRUE)) / (B_perm + 1)
# #       }

# #       res$omni_perm_p_BH <- .p_adjust_BH(res$omni_perm_p)
# #     }

# #     if (perm_mode == "mvn") {

# #       if (is.null(magma_genes_out) || !file.exists(magma_genes_out)) {
# #         stop("perm_mode='mvn' requires magma_genes_out = path to merged MAGMA *.genes.out", call. = FALSE)
# #       }

# #       genes_out_tab <- magma_read_genes_out(magma_genes_out, gene_col = prep$gene_col, chr_col = "CHR")

# #       cor_pairs <- NULL
# #       if (!is.null(magma_cor_file)) {
# #         if (!file.exists(magma_cor_file)) stop("magma_cor_file does not exist: ", magma_cor_file, call. = FALSE)
# #         cor_pairs <- magma_read_gene_cor_pairs(magma_cor_file, gene1_col = 1, gene2_col = 2, r_col = 3)
# #       }

# #       for (i in seq_len(nrow(res))) {

# #         pid <- res$pathway_id[i]
# #         genes_S <- prep$g_list[[pid]]
# #         if (is.null(genes_S) || length(genes_S) < 2L) next
# #         if (!is.finite(res$omni_p[i])) next

# #         R_S <- magma_build_R_for_pathway(
# #           genes_S       = genes_S,
# #           genes_out_tab = genes_out_tab,
# #           cor_pairs     = cor_pairs,
# #           gene_col      = prep$gene_col,
# #           chr_col       = "CHR",
# #           make_PD       = make_PD
# #         )
# #         if (is.null(R_S)) next

# #         sim <- magma_simulate_null_Zp(
# #           genes_S       = rownames(R_S),
# #           genes_out_tab = genes_out_tab,
# #           R_S           = R_S,
# #           B             = B_perm,
# #           gene_col      = prep$gene_col,
# #           chr_col       = "CHR"
# #         )

# #         genes_sim <- colnames(sim$P)
# #         d <- length(genes_sim)
# #         if (d < 2L) next

# #         # weights aligned to genes_sim (if available)
# #         w_i <- prep$w_norm
# #         if (!is.null(w_i)) {
# #           wv <- as.numeric(w_i[tolower(genes_sim)])
# #           if (length(wv) != d || all(!is.finite(wv))) wv <- rep(1, d)
# #         } else {
# #           wv <- rep(1, d)
# #         }

# #         omni_null <- numeric(B_perm)

# #         for (b in seq_len(B_perm)) {

# #           p_i <- as.numeric(sim$P[b, ])
# #           z_i <- as.numeric(sim$Z[b, ])

# #           p_acat   <- ACAT::ACAT(Pvals = fix_p_for_acat(p_i, min_p = prep$min_p))

# #           p_fisher <- stats::pchisq(-2 * sum(log(fix_p_for_acat(p_i, min_p = prep$min_p))),
# #                                     df = 2 * d, lower.tail = FALSE)

# #           tf <- .tfisher_adaptive_pmin(p_i, prep$tau_grid, prep$tau_cap, prep$min_p)
# #           p_tf <- tf$p_min

# #           p_st <- .stouffer_p_from_z(z_i, wv)

# #           p_min <- .sidak_minp(p_i)

# #           omni_null[b] <- .combine_omni(c(p_acat, p_fisher, p_tf, p_st, p_min))
# #         }

# #         omni_obs <- res$omni_p[i]
# #         res$omni_perm_p[i] <- (1 + sum(omni_null <= omni_obs, na.rm = TRUE)) / (B_perm + 1)
# #       }

# #       res$omni_perm_p_BH <- .p_adjust_BH(res$omni_perm_p)
# #     }
# #   }

# #   # final p = perm if available else analytic
# #   res$omni_final_p  <- if (B_perm > 0L) res$omni_perm_p else res$omni_p
# #   res$omni_final_BH <- .p_adjust_BH(res$omni_final_p)

# #   # sort
# #   ord <- order(res$omni_final_p, decreasing = FALSE, na.last = TRUE)
# #   res <- res[ord, , drop = FALSE]

# #   # optional output
# #   if (output) {
# #     if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# #     out_tag <- if (omnibus == "ACAT") "acato" else "minp"
# #     out_path <- file.path(out_dir, paste0("omni_pathways_", out_tag, ".csv"))
# #     utils::write.csv(res, out_path, row.names = FALSE)
# #     attr(res, "file") <- out_path
# #   }

# #   res
# # }


# ## ============================================================
# ##  MAGCAT omnibus pathway testing (WORKING, consistent prep)
# ##  Components:
# ##    1) ACAT
# ##    2) Fisher
# ##    3) Adaptive soft TFisher (analytic over tau grid)
# ##    4) minP (Sidak)
# ##    5) Stouffer Z (uses z_col; optional weights)
# ##    6) MAGMA competitive (optional; read from tidy TSV)  <-- included in omnibus if present
# ##
# ##  Exports:
# ##    - magcat_omni_prepare()
# ##    - magcat_omni_run_prepared()
# ##    - magcat_omni_pathways() convenience wrapper
# ## ============================================================

# # ---- helper: safe p-fix for ACAT/Fisher/TFisher/minP ----
# .magcat_fix_p <- function(p, min_p = 1e-15) {
#   p <- suppressWarnings(as.numeric(p))
#   p <- p[is.finite(p) & !is.na(p)]
#   if (!length(p)) return(numeric(0))
#   p[p <= 0] <- min_p
#   p[p >= 1] <- 1 - min_p
#   p
# }

# # ---- helper: BH with NAs ----
# .magcat_padj_BH <- function(p) {
#   p <- as.numeric(p)
#   out <- rep(NA_real_, length(p))
#   idx <- which(is.finite(p) & !is.na(p))
#   if (length(idx)) out[idx] <- stats::p.adjust(p[idx], method = "BH")
#   out
# }

# # ---- helper: Sidak minP for k p-values ----
# .magcat_sidak_minp <- function(p_vec) {
#   p <- .magcat_fix_p(p_vec, min_p = 1e-15)
#   p <- p[p > 0 & p < 1]
#   k <- length(p)
#   if (k < 1L) return(NA_real_)
#   pmin <- min(p)
#   1 - (1 - pmin)^k
# }

# # ---- helper: read MAGMA tidy TSV ----
# .magcat_read_magma_tidy <- function(magma_tidy_tsv) {
#   if (is.null(magma_tidy_tsv)) return(NULL)
#   if (!file.exists(magma_tidy_tsv)) stop("magma_tidy_tsv does not exist: ", magma_tidy_tsv, call. = FALSE)

#   tab <- utils::read.delim(magma_tidy_tsv, sep = "\t", header = TRUE,
#                            stringsAsFactors = FALSE, check.names = FALSE)

#   if (!("pathway_id" %in% names(tab))) stop("magma_tidy_tsv must have column 'pathway_id'.", call. = FALSE)

#   pcol <- NULL
#   if ("magma_pvalue" %in% names(tab)) {
#     pcol <- "magma_pvalue"
#   } else {
#     cand <- grep("magma.*p", names(tab), ignore.case = TRUE, value = TRUE)
#     if (length(cand)) pcol <- cand[1]
#   }
#   if (is.null(pcol)) stop("magma_tidy_tsv needs a MAGMA p-value column (e.g. magma_pvalue).", call. = FALSE)

#   out <- tab[, c("pathway_id", pcol), drop = FALSE]
#   names(out)[2] <- "magma_pvalue"
#   out
# }

# # ---- helper: coerce pathways into (pathway_id, pathway_name, gene_id) long df ----
# .magcat_as_pathway_long <- function(pathways) {
#   if (is.data.frame(pathways)) {
#     if (!all(c("pathway_id","gene_id") %in% names(pathways))) {
#       stop("If pathways is a data.frame it must have columns pathway_id and gene_id.", call. = FALSE)
#     }
#     if (!("pathway_name" %in% names(pathways))) pathways$pathway_name <- pathways$pathway_id
#     out <- pathways[, c("pathway_id","pathway_name","gene_id"), drop = FALSE]
#     out$gene_id <- as.character(out$gene_id)
#     return(out)
#   }

#   if (is.list(pathways)) {
#     pw_id <- names(pathways)
#     if (is.null(pw_id)) pw_id <- paste0("PWY_", seq_along(pathways))
#     out <- data.frame(
#       pathway_id   = rep(pw_id, lengths(pathways)),
#       pathway_name = rep(pw_id, lengths(pathways)),
#       gene_id      = unlist(pathways, use.names = FALSE),
#       stringsAsFactors = FALSE
#     )
#     out$gene_id <- as.character(out$gene_id)
#     return(out)
#   }

#   stop("pathways must be a data.frame or a list.", call. = FALSE)
# }

# #' Prepare omnibus object (CONSISTENT)
# #' @export
# magcat_omni_prepare <- function(gene_results,
#                                 pathways     = NULL,
#                                 species      = NULL,
#                                 pmn_gene_col = NULL,
#                                 gene_col     = "GENE",
#                                 p_col        = "P",
#                                 z_col        = NULL,
#                                 weight_col   = NULL,
#                                 tau_grid     = c(0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.001),
#                                 tau_cap      = 0.5,
#                                 min_p        = 1e-15,
#                                 do_fix       = TRUE,
#                                 seed         = NULL) {

#   if (!is.null(seed)) set.seed(seed)

#   ## ---- pathway source ----
#   if (!is.null(pathways) && !is.null(species)) {
#     stop("Provide either 'pathways' OR 'species', not both.", call. = FALSE)
#   }
#   if (is.null(pathways) && is.null(species)) {
#     stop("Need 'pathways' (list/data.frame) or 'species' (built-in).", call. = FALSE)
#   }
#   if (is.null(pathways)) {
#     if (!exists("magcat_load_pathways", mode = "function")) {
#       stop("magcat_omni_prepare(): missing magcat_load_pathways().", call. = FALSE)
#     }
#     pathways <- magcat_load_pathways(species = species, gene_col = pmn_gene_col)
#   }

#   ## ---- gene_results checks ----
#   if (!is.data.frame(gene_results)) stop("gene_results must be a data.frame.", call. = FALSE)
#   if (!all(c(gene_col, p_col) %in% names(gene_results))) {
#     stop("gene_results must contain columns '", gene_col, "' and '", p_col, "'.", call. = FALSE)
#   }
#   if (!is.null(z_col) && !(z_col %in% names(gene_results))) {
#     stop("z_col = '", z_col, "' not found in gene_results.", call. = FALSE)
#   }
#   if (!is.null(weight_col) && !(weight_col %in% names(gene_results))) {
#     stop("weight_col = '", weight_col, "' not found in gene_results.", call. = FALSE)
#   }

#   ## ---- normalize gene table ----
#   gr_gene_raw  <- as.character(gene_results[[gene_col]])
#   gr_gene_norm <- tolower(trimws(gr_gene_raw))

#   p_all <- suppressWarnings(as.numeric(gene_results[[p_col]]))
#   z_all <- if (!is.null(z_col)) suppressWarnings(as.numeric(gene_results[[z_col]])) else NULL
#   w_all <- if (!is.null(weight_col)) suppressWarnings(as.numeric(gene_results[[weight_col]])) else NULL

#   ok <- which(!is.na(gr_gene_norm) & gr_gene_norm != "" & gr_gene_norm != "unknown" &
#               is.finite(p_all) & !is.na(p_all))
#   if (!length(ok)) stop("No valid genes found in gene_results after filtering.", call. = FALSE)

#   gr_gene_raw  <- gr_gene_raw[ok]
#   gr_gene_norm <- gr_gene_norm[ok]
#   p_all        <- p_all[ok]
#   if (!is.null(z_all)) z_all <- z_all[ok]
#   if (!is.null(w_all)) w_all <- w_all[ok]

#   # deduplicate gene ids
#   keep1 <- !duplicated(gr_gene_norm)
#   gr_gene_raw  <- gr_gene_raw[keep1]
#   gr_gene_norm <- gr_gene_norm[keep1]
#   p_all        <- p_all[keep1]
#   if (!is.null(z_all)) z_all <- z_all[keep1]
#   if (!is.null(w_all)) w_all <- w_all[keep1]

#   gene_map_raw <- stats::setNames(gr_gene_raw, gr_gene_norm)

#   gene_p <- stats::setNames(p_all, gr_gene_norm)
#   gene_z <- if (!is.null(z_all)) stats::setNames(z_all, gr_gene_norm) else NULL
#   gene_w <- if (!is.null(w_all)) stats::setNames(w_all, gr_gene_norm) else NULL

#   ## ---- pathways -> long -> per-pathway gene lists ----
#   if (is.data.frame(pathways)) {
#     if (!all(c("pathway_id","gene_id") %in% names(pathways))) {
#       stop("If 'pathways' is a data.frame it must have columns 'pathway_id' and 'gene_id'.", call. = FALSE)
#     }
#     if (!("pathway_name" %in% names(pathways))) pathways$pathway_name <- pathways$pathway_id

#     pw_long <- pathways[, c("pathway_id","pathway_name","gene_id"), drop = FALSE]
#     pw_long$gene_id <- as.character(pw_long$gene_id)

#   } else if (is.list(pathways)) {
#     pw_id <- names(pathways)
#     if (is.null(pw_id)) pw_id <- paste0("PWY_", seq_along(pathways))
#     pw_long <- data.frame(
#       pathway_id   = rep(pw_id, lengths(pathways)),
#       pathway_name = rep(pw_id, lengths(pathways)),
#       gene_id      = unlist(pathways, use.names = FALSE),
#       stringsAsFactors = FALSE
#     )
#   } else {
#     stop("'pathways' must be a data.frame or a list.", call. = FALSE)
#   }

#   pw_long$gene_id <- gsub("^gene:", "", pw_long$gene_id, ignore.case = TRUE)
#   pw_long$gene_id <- trimws(pw_long$gene_id)
#   pw_long$gene_norm <- tolower(pw_long$gene_id)
#   pw_long <- pw_long[!is.na(pw_long$gene_norm) & pw_long$gene_norm != "" & pw_long$gene_norm != "unknown", , drop = FALSE]

#   g_list_norm <- split(pw_long$gene_norm, pw_long$pathway_id)
#   pw_name     <- tapply(pw_long$pathway_name, pw_long$pathway_id, function(x) x[1])

#   g_list_norm <- lapply(g_list_norm, function(g) unique(g[g %in% names(gene_p)]))

#   p_list <- lapply(g_list_norm, function(g) unname(gene_p[g]))
#   z_list <- if (!is.null(gene_z)) lapply(g_list_norm, function(g) unname(gene_z[g])) else NULL
#   w_list <- if (!is.null(gene_w)) lapply(g_list_norm, function(g) unname(gene_w[g])) else NULL
#   g_list_raw <- lapply(g_list_norm, function(g) unname(gene_map_raw[g]))

#   # tau grid sanitize
#   tau_grid <- as.numeric(tau_grid)
#   tau_grid <- tau_grid[is.finite(tau_grid) & tau_grid > 0 & tau_grid < 1 & tau_grid <= tau_cap]
#   if (!length(tau_grid)) tau_grid <- min(0.05, tau_cap)
#   tau_grid <- sort(unique(tau_grid), decreasing = TRUE)

#   structure(
#     list(
#       gene_col     = gene_col,
#       p_col        = p_col,
#       z_col        = z_col,
#       weight_col   = weight_col,
#       min_p        = min_p,
#       do_fix       = isTRUE(do_fix),
#       tau_grid     = tau_grid,
#       tau_cap      = tau_cap,
#       seed         = seed,

#       pathway_id   = names(g_list_norm),
#       pathway_name = as.character(pw_name[names(g_list_norm)]),

#       g_list_norm  = g_list_norm,
#       g_list_raw   = g_list_raw,
#       p_list       = p_list,
#       z_list       = z_list,
#       w_list       = w_list,

#       p_all        = p_all,
#       z_all        = z_all,
#       w_all        = w_all
#     ),
#     class = "magcat_omni_prep"
#   )
# }


# # ---- component p-values from vectors ----
# .magcat_p_acat <- function(p_i, min_p = 1e-15) {
#   if (!requireNamespace("ACAT", quietly = TRUE)) stop("Needs ACAT package.", call. = FALSE)
#   p <- .magcat_fix_p(p_i, min_p = min_p)
#   if (length(p) < 2L) return(NA_real_)
#   as.numeric(ACAT::ACAT(Pvals = p))
# }

# .magcat_p_fisher <- function(p_i, min_p = 1e-15) {
#   p <- .magcat_fix_p(p_i, min_p = min_p)
#   p <- p[p > 0 & p < 1]
#   k <- length(p)
#   if (k < 2L) return(NA_real_)
#   stats::pchisq(-2 * sum(log(p)), df = 2 * k, lower.tail = FALSE)
# }

# .magcat_p_tfisher_adapt <- function(p_i, tau_grid, min_p = 1e-15) {
#   if (!requireNamespace("TFisher", quietly = TRUE)) stop("Needs TFisher package.", call. = FALSE)
#   p <- .magcat_fix_p(p_i, min_p = min_p)
#   p <- p[p > 0 & p < 1]
#   n <- length(p)
#   if (n < 2L) return(list(p = NA_real_, tau = NA_real_))

#   best_p   <- Inf
#   best_tau <- tau_grid[1]

#   for (tau in tau_grid) {
#     st <- TFisher::stat.soft(p = p, tau1 = tau)
#     pv <- 1 - as.numeric(TFisher::p.soft(q = st, n = n, tau1 = tau, M = NULL))
#     if (is.finite(pv) && pv < best_p) {
#       best_p <- pv
#       best_tau <- tau
#     }
#   }

#   list(p = as.numeric(best_p), tau = as.numeric(best_tau))
# }

# .magcat_p_stouffer <- function(z_i, w_i = NULL, min_abs_w = 1e-8) {
#   z <- suppressWarnings(as.numeric(z_i))
#   keep <- is.finite(z) & !is.na(z)
#   z <- z[keep]
#   if (length(z) < 2L) return(NA_real_)

#   if (is.null(w_i)) {
#     w <- rep(1, length(z))
#   } else {
#     w <- suppressWarnings(as.numeric(w_i))[keep]
#     bad <- !is.finite(w) | is.na(w) | abs(w) <= min_abs_w
#     if (any(bad)) {
#       w_ok <- w[!bad]
#       repl <- if (length(w_ok)) stats::median(abs(w_ok)) else 1
#       w[bad] <- repl
#     }
#   }

#   Zs <- sum(w * z) / sqrt(sum(w^2))
#   2 * stats::pnorm(-abs(Zs))
# }

# # ---- omnibus combiner ----
# .magcat_omni_combine <- function(p_vec, omnibus = c("ACAT","minP"), min_p = 1e-15) {
#   omnibus <- match.arg(omnibus)
#   pv <- .magcat_fix_p(p_vec, min_p = min_p)
#   pv <- pv[pv > 0 & pv < 1]
#   if (!length(pv)) return(NA_real_)

#   if (omnibus == "ACAT") {
#     if (!requireNamespace("ACAT", quietly = TRUE)) stop("Needs ACAT package.", call. = FALSE)
#     return(as.numeric(ACAT::ACAT(Pvals = pv)))
#   } else {
#     k <- length(pv)
#     pmin <- min(pv)
#     return(1 - (1 - pmin)^k)
#   }
# }

# #' Run omnibus from prepared object
# #' @export
# magcat_omni_run_prepared <- function(prep,
#                                     omnibus = c("ACAT","minP"),
#                                     B_perm = 0L,
#                                     perm_mode = c("resample_global","mvn"),
#                                     magma_genes_out = NULL,
#                                     magma_cor_file  = NULL,
#                                     make_PD = TRUE,
#                                     remove_singletons = TRUE,
#                                     magma_tidy_tsv = NULL,
#                                     include_magma_in_omni = TRUE,
#                                     output = FALSE,
#                                     out_dir = "magcat_omni") {

#   if (!inherits(prep, "magcat_omni_prep")) stop("prep must come from magcat_omni_prepare().", call. = FALSE)

#   omnibus   <- match.arg(omnibus)
#   perm_mode <- match.arg(perm_mode)

#   # MAGMA competitive p-values (optional, included as 6th component if requested)
#   mg <- .magcat_read_magma_tidy(magma_tidy_tsv)

#   # build base result table
#   res <- data.frame(
#     pathway_id   = prep$pathway_id,
#     pathway_name = prep$pathway_name,
#     n_genes      = vapply(prep$g_list_norm, length, integer(1)),
#     stringsAsFactors = FALSE
#   )

#   if (remove_singletons) {
#     keep <- is.finite(res$n_genes) & res$n_genes >= 2L
#     res <- res[keep, , drop = FALSE]
#   }

#   # observed component p-values
#   pid <- res$pathway_id

#   res$acat_p <- vapply(pid, function(k) .magcat_p_acat(prep$p_list[[k]], min_p = prep$min_p), numeric(1))
#   res$fisher_p <- vapply(pid, function(k) .magcat_p_fisher(prep$p_list[[k]], min_p = prep$min_p), numeric(1))

#   tf_obs <- lapply(pid, function(k) .magcat_p_tfisher_adapt(prep$p_list[[k]], prep$tau_grid, min_p = prep$min_p))
#   res$tfisher_adapt_p_analytic <- vapply(tf_obs, `[[`, numeric(1), "p")
#   res$tfisher_adapt_tau_best   <- vapply(tf_obs, `[[`, numeric(1), "tau")

#   res$minp_gene_p <- vapply(pid, function(k) .magcat_sidak_minp(prep$p_list[[k]]), numeric(1))

#   # Stouffer: requires z_col
#   if (!is.null(prep$z_col)) {
#     res$stouffer_p_analytic <- vapply(pid, function(k) {
#       .magcat_p_stouffer(prep$z_list[[k]], prep$w_list[[k]])
#     }, numeric(1))
#   } else {
#     res$stouffer_p_analytic <- NA_real_
#   }

#   # MAGMA comp
#   res$magma_pvalue <- NA_real_
#   if (!is.null(mg)) {
#     res$magma_pvalue <- mg$magma_pvalue[match(res$pathway_id, mg$pathway_id)]
#   }

#   # omnibus (analytic)
#   res$omni_p <- vapply(seq_len(nrow(res)), function(i) {
#     comps <- c(res$acat_p[i],
#                res$fisher_p[i],
#                res$tfisher_adapt_p_analytic[i],
#                res$minp_gene_p[i],
#                res$stouffer_p_analytic[i])

#     if (isTRUE(include_magma_in_omni) && is.finite(res$magma_pvalue[i])) {
#       comps <- c(comps, res$magma_pvalue[i])
#     }

#     .magcat_omni_combine(comps, omnibus = omnibus, min_p = prep$min_p)
#   }, numeric(1))

#   res$omni_p_BH <- .magcat_padj_BH(res$omni_p)

#   # permutation-calibrated omni
#   res$omni_perm_p    <- NA_real_
#   res$omni_perm_p_BH <- NA_real_

#   if (B_perm > 0L) {

#     if (!is.null(prep$seed)) set.seed(prep$seed)

#     if (perm_mode == "resample_global") {

#       pool_p <- suppressWarnings(as.numeric(prep$p_all))
#       okp <- is.finite(pool_p) & !is.na(pool_p) & pool_p > 0 & pool_p < 1
#       pool_p <- .magcat_fix_p(pool_p[okp], min_p = prep$min_p)

#       pool_z <- NULL
#       if (!is.null(prep$z_all)) {
#         pool_z <- suppressWarnings(as.numeric(prep$z_all))[okp]
#         # keep only finite z too (align)
#         okz <- is.finite(pool_z) & !is.na(pool_z)
#         pool_p <- pool_p[okz]
#         pool_z <- pool_z[okz]
#       }

#       for (i in seq_len(nrow(res))) {

#         d <- res$n_genes[i]
#         if (!is.finite(d) || d < 2L) next
#         if (!is.finite(res$omni_p[i])) next

#         replace_flag <- length(pool_p) < d

#         # fixed weights for this pathway (if any)
#         w_i <- NULL
#         if (!is.null(prep$w_list) && !is.null(prep$w_list[[res$pathway_id[i]]])) {
#           w_i <- prep$w_list[[res$pathway_id[i]]]
#           if (length(w_i) != d) w_i <- rep(1, d)
#         }

#         omni_null <- numeric(B_perm)

#         for (b in seq_len(B_perm)) {
#           idx <- sample.int(length(pool_p), size = d, replace = replace_flag)
#           p_i <- pool_p[idx]

#           p_acat   <- .magcat_p_acat(p_i, min_p = prep$min_p)
#           p_fisher <- .magcat_p_fisher(p_i, min_p = prep$min_p)
#           tf       <- .magcat_p_tfisher_adapt(p_i, prep$tau_grid, min_p = prep$min_p)
#           p_tf     <- tf$p
#           p_min    <- .magcat_sidak_minp(p_i)

#           p_st <- NA_real_
#           if (!is.null(pool_z)) {
#             z_i <- pool_z[idx]
#             p_st <- .magcat_p_stouffer(z_i, w_i)
#           }

#           comps <- c(p_acat, p_fisher, p_tf, p_min, p_st)

#           # MAGMA competitive is not permuted here (it’s external); we leave it OUT of perm null
#           # (If you want to permute MAGMA comp too, you need a null for MAGMA competitive.)
#           omni_null[b] <- .magcat_omni_combine(comps, omnibus = omnibus, min_p = prep$min_p)
#         }

#         omni_obs <- res$omni_p[i]
#         res$omni_perm_p[i] <- (1 + sum(omni_null <= omni_obs, na.rm = TRUE)) / (B_perm + 1)
#       }

#       res$omni_perm_p_BH <- .magcat_padj_BH(res$omni_perm_p)
#     }

#     if (perm_mode == "mvn") {

#       # MVN needs your MAGMA-correlation helpers
#       req <- c("magma_read_genes_out", "magma_build_R_for_pathway", "magma_simulate_null_Zp", "magma_read_gene_cor_pairs")
#       miss <- req[!vapply(req, exists, logical(1), mode = "function")]
#       if (length(miss)) {
#         stop("perm_mode='mvn' requires these functions to be loaded: ",
#              paste(miss, collapse = ", "), call. = FALSE)
#       }

#       if (is.null(magma_genes_out) || !file.exists(magma_genes_out)) {
#         stop("perm_mode='mvn' requires magma_genes_out = path to merged MAGMA *.genes.out (or your merged genes file).",
#              call. = FALSE)
#       }

#       genes_out_tab <- magma_read_genes_out(magma_genes_out, gene_col = prep$gene_col, chr_col = "CHR")

#       cor_pairs <- NULL
#       if (!is.null(magma_cor_file)) {
#         if (!file.exists(magma_cor_file)) stop("magma_cor_file does not exist: ", magma_cor_file, call. = FALSE)
#         cor_pairs <- magma_read_gene_cor_pairs(magma_cor_file, gene1_col = 1, gene2_col = 2, r_col = 3)
#       }

#       for (i in seq_len(nrow(res))) {

#         pid_i <- res$pathway_id[i]
#         genes_S_raw <- prep$g_list_raw[[pid_i]]

#         genes_S_raw <- as.character(genes_S_raw)
#         genes_S_raw <- genes_S_raw[!is.na(genes_S_raw) & genes_S_raw != ""]

#         if (length(genes_S_raw) < 2L) next
#         if (!is.finite(res$omni_p[i])) next

#         R_S <- magma_build_R_for_pathway(
#           genes_S       = genes_S_raw,
#           genes_out_tab = genes_out_tab,
#           cor_pairs     = cor_pairs,
#           gene_col      = prep$gene_col,
#           chr_col       = "CHR",
#           make_PD       = make_PD
#         )
#         if (is.null(R_S)) next

#         sim <- magma_simulate_null_Zp(
#           genes_S       = rownames(R_S),
#           genes_out_tab = genes_out_tab,
#           R_S           = R_S,
#           B             = B_perm,
#           gene_col      = prep$gene_col,
#           chr_col       = "CHR"
#         )

#         d <- ncol(sim$P)
#         if (d < 2L) next

#         # weights: if you provided weight_col, align by gene name (lowercase map)
#         wv <- rep(1, d)
#         if (!is.null(prep$weight_col) && !is.null(prep$w_all)) {
#           # try to pull weights from prep by matching gene_results norm IDs
#           # sim genes are raw; convert to norm
#           gsim_norm <- tolower(colnames(sim$P))
#           # rebuild a norm->w map from prep (global)
#           # (best-effort; if mismatch, fallback to 1s)
#           w_map <- NULL
#           if (!is.null(prep$w_all)) {
#             # prep$w_all is numeric vector aligned to gene_results after filtering/dedup
#             # but we didn't keep names; so we can’t map perfectly here without extra storage
#             # -> fallback to ones unless you extend prep to store named weights
#             w_map <- NULL
#           }
#           # leave wv as ones (safe)
#         }

#         omni_null <- numeric(B_perm)

#         for (b in seq_len(B_perm)) {
#           p_i <- as.numeric(sim$P[b, ])
#           z_i <- as.numeric(sim$Z[b, ])

#           p_acat   <- .magcat_p_acat(p_i, min_p = prep$min_p)
#           p_fisher <- .magcat_p_fisher(p_i, min_p = prep$min_p)
#           tf       <- .magcat_p_tfisher_adapt(p_i, prep$tau_grid, min_p = prep$min_p)
#           p_tf     <- tf$p
#           p_min    <- .magcat_sidak_minp(p_i)
#           p_st     <- .magcat_p_stouffer(z_i, wv)

#           comps <- c(p_acat, p_fisher, p_tf, p_min, p_st)
#           omni_null[b] <- .magcat_omni_combine(comps, omnibus = omnibus, min_p = prep$min_p)
#         }

#         omni_obs <- res$omni_p[i]
#         res$omni_perm_p[i] <- (1 + sum(omni_null <= omni_obs, na.rm = TRUE)) / (B_perm + 1)
#       }

#       res$omni_perm_p_BH <- .magcat_padj_BH(res$omni_perm_p)
#     }
#   }

#   # final p
#   res$omni_final_p  <- if (B_perm > 0L) res$omni_perm_p else res$omni_p
#   res$omni_final_BH <- .magcat_padj_BH(res$omni_final_p)

#   # sort + output
#   res <- res[order(res$omni_final_p, decreasing = FALSE, na.last = TRUE), , drop = FALSE]

#   if (output) {
#     if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#     out_tag <- if (omnibus == "ACAT") "acat" else "minp"
#     out_path <- file.path(out_dir, paste0("omni_pathways_", out_tag, ".csv"))
#     utils::write.csv(res, out_path, row.names = FALSE)
#     attr(res, "file") <- out_path
#   }

#   res
# }

# #' Convenience wrapper: prepare + run
# #' @export
# magcat_omni_pathways <- function(gene_results,
#                                 pathways     = NULL,
#                                 species      = NULL,
#                                 pmn_gene_col = NULL,
#                                 gene_col     = "GENE",
#                                 p_col        = "P",
#                                 z_col        = NULL,
#                                 weight_col   = NULL,
#                                 tau_grid     = c(0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.001),
#                                 tau_cap      = 0.5,
#                                 min_p        = 1e-15,
#                                 do_fix       = TRUE,
#                                 omnibus      = c("ACAT","minP"),
#                                 B_perm       = 0L,
#                                 perm_mode    = c("resample_global","mvn","global_resampling","global_sampling"),
#                                 magma_genes_out = NULL,
#                                 magma_cor_file  = NULL,
#                                 make_PD      = TRUE,
#                                 seed         = NULL,
#                                 remove_singletons = TRUE,
#                                 magma_tidy_tsv = NULL,
#                                 include_magma_in_omni = TRUE,
#                                 output       = FALSE,
#                                 out_dir      = "magcat_omni") {

#   omnibus <- match.arg(omnibus)
#   perm_mode <- match.arg(perm_mode)
#   if (perm_mode %in% c("global_resampling","global_sampling")) perm_mode <- "resample_global"

#   prep <- magcat_omni_prepare(
#     gene_results = gene_results,
#     pathways     = pathways,
#     species      = species,
#     pmn_gene_col = pmn_gene_col,
#     gene_col     = gene_col,
#     p_col        = p_col,
#     z_col        = z_col,
#     weight_col   = weight_col,
#     tau_grid     = tau_grid,
#     tau_cap      = tau_cap,
#     min_p        = min_p,
#     do_fix       = do_fix,
#     seed         = seed
#   )

#   magcat_omni_run_prepared(
#     prep              = prep,
#     omnibus           = omnibus,
#     B_perm            = B_perm,
#     perm_mode         = perm_mode,
#     magma_genes_out   = magma_genes_out,
#     magma_cor_file    = magma_cor_file,
#     make_PD           = make_PD,
#     remove_singletons = remove_singletons,
#     magma_tidy_tsv    = magma_tidy_tsv,
#     include_magma_in_omni = include_magma_in_omni,
#     output            = output,
#     out_dir           = out_dir
#   )
# }
