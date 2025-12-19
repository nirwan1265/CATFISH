# ## ============================================================
# ##  MAGCAT: prepared (cached) pathway-testing wrappers
# ##  Methods supported here (6):
# ##    1) ACAT
# ##    2) Fisher
# ##    3) Adaptive soft TFisher
# ##    4) minP
# ##    5) Stouffer Z (optionally weighted)
# ##    6) MAGMA competitive gene-set test (requires MAGMA + set.annot)
# ##
# ##  This file defines:
# ##    - magcat_prepare_core()
# ##    - *_prepare() / *_run_prepared() pairs for each method
# ##
# ##  Notes:
# ##    * These prepared runners compute ANALYTIC p-values only.
# ##      (If you want permutation calibration, use your existing wrapper
# ##       functions that implement perm_mode = 'resample_global'/'mvn'.)
# ##    * MAGMA competitive is inherently "external" (calls MAGMA), so the
# ##      prepared runner just freezes paths/inputs and then runs MAGMA.
# ## ============================================================

# #' Prepare common pathway + gene maps for fast reuse
# #'
# #' Freezes repeated setup steps: pathway list building, gene normalization,
# #' canonical gene-id mapping, and named vectors for p / Z / weights.
# #'
# #' Requires:
# #'   - magcat_load_pathways()  (from ACAT_wrappers.R)
# #'   - fix_p_for_acat()        (from ACAT_wrappers.R)
# #'
# #' @keywords internal
# magcat_prepare_core <- function(gene_results,
#                                 pathways     = NULL,
#                                 species      = NULL,
#                                 pmn_gene_col = NULL,
#                                 gene_col     = "GENE",
#                                 p_col        = "P",
#                                 z_col        = NULL,
#                                 weight_col   = NULL) {

#   ## ---- pathway source ----
#   if (!is.null(pathways) && !is.null(species)) {
#     stop("Provide either 'pathways' OR 'species', not both.", call. = FALSE)
#   }
#   if (is.null(pathways) && is.null(species)) {
#     stop("Need 'pathways' (list/data.frame) or 'species' (built-in).", call. = FALSE)
#   }
#   if (is.null(pathways)) {
#     if (!exists("magcat_load_pathways", mode = "function")) {
#       stop("magcat_prepare_core(): missing magcat_load_pathways().", call. = FALSE)
#     }
#     pathways <- magcat_load_pathways(
#       species  = species,
#       gene_col = pmn_gene_col
#     )
#   }

#   ## ---- standardize gene_results ----
#   if (!is.data.frame(gene_results)) {
#     stop("gene_results must be a data.frame.", call. = FALSE)
#   }
#   if (!all(c(gene_col, p_col) %in% names(gene_results))) {
#     stop("gene_results must contain columns '", gene_col, "' and '", p_col, "'.", call. = FALSE)
#   }
#   if (!is.null(z_col) && !(z_col %in% names(gene_results))) {
#     stop("z_col = '", z_col, "' not found in gene_results.", call. = FALSE)
#   }
#   if (!is.null(weight_col) && !(weight_col %in% names(gene_results))) {
#     stop("weight_col = '", weight_col, "' not found in gene_results.", call. = FALSE)
#   }

#   gr <- gene_results

#   genes_all_u <- as.character(gr[[gene_col]])
#   genes_norm  <- tolower(genes_all_u)

#   p_all <- suppressWarnings(as.numeric(gr[[p_col]]))

#   z_all <- NULL
#   if (!is.null(z_col)) {
#     z_all <- suppressWarnings(as.numeric(gr[[z_col]]))
#   }

#   w_all <- NULL
#   if (!is.null(weight_col)) {
#     w_all <- suppressWarnings(as.numeric(gr[[weight_col]]))
#   }

#   ## ---- drop obviously bad gene IDs ----
#   ok_gene <- which(!is.na(genes_norm) & genes_norm != "" & genes_norm != "unknown")
#   if (!length(ok_gene)) stop("No valid genes found in gene_results.", call. = FALSE)

#   genes_all_u <- genes_all_u[ok_gene]
#   genes_norm  <- genes_norm[ok_gene]
#   p_all       <- p_all[ok_gene]
#   if (!is.null(z_all)) z_all <- z_all[ok_gene]
#   if (!is.null(w_all)) w_all <- w_all[ok_gene]

#   ## ---- deduplicate by normalized gene id (keep first) ----
#   keep1 <- !duplicated(genes_norm)
#   genes_all_u <- genes_all_u[keep1]
#   genes_norm  <- genes_norm[keep1]
#   p_all       <- p_all[keep1]
#   if (!is.null(z_all)) z_all <- z_all[keep1]
#   if (!is.null(w_all)) w_all <- w_all[keep1]

#   ## ---- canonical map: norm -> original case (first) ----
#   gene_map <- stats::setNames(genes_all_u, genes_norm)

#   ## ---- named p-vector ----
#   gene_p_vec <- stats::setNames(p_all, genes_norm)

#   ## ---- named z-vector (optional) ----
#   gene_z_vec <- NULL
#   if (!is.null(z_all)) {
#     gene_z_vec <- stats::setNames(z_all, genes_norm)
#   }

#   ## ---- named weights (optional) ----
#   w_norm <- NULL
#   if (!is.null(w_all)) {
#     w_norm <- stats::setNames(w_all, genes_norm)
#   }

#   ## ---- pathways -> p_list + p_names ----
#   if (is.data.frame(pathways)) {
#     if (!all(c("pathway_id", "gene_id") %in% names(pathways))) {
#       stop("If 'pathways' is a data.frame it must have columns 'pathway_id' and 'gene_id'.", call. = FALSE)
#     }
#     if (!("pathway_name" %in% names(pathways))) pathways$pathway_name <- pathways$pathway_id

#     gene_id <- as.character(pathways$gene_id)
#     gene_id <- gsub("^gene:", "", gene_id, ignore.case = TRUE)
#     gene_id <- trimws(gene_id)

#     p_list <- split(gene_id, pathways$pathway_id)
#     p_names <- tapply(pathways$pathway_name, pathways$pathway_id, FUN = function(x) x[1])

#   } else if (is.list(pathways)) {
#     p_list <- pathways
#     if (is.null(names(p_list))) {
#       names(p_list) <- paste0("PWY_", seq_along(p_list))
#     }
#     p_names <- names(p_list)

#   } else {
#     stop("'pathways' must be a data.frame or a list.", call. = FALSE)
#   }

#   ## normalize pathway gene ids
#   p_list <- lapply(p_list, function(g) {
#     g <- as.character(g)
#     g <- gsub("^gene:", "", g, ignore.case = TRUE)
#     g <- trimws(g)
#     g <- tolower(g)
#     g[!is.na(g) & g != "" & g != "unknown"]
#   })

#   structure(
#     list(
#       gene_col     = gene_col,
#       p_col        = p_col,
#       z_col        = z_col,
#       weight_col   = weight_col,
#       genes_all_u  = genes_all_u,
#       genes_norm   = genes_norm,
#       p_all        = p_all,
#       z_all        = z_all,
#       w_all        = w_all,
#       gene_p_vec   = gene_p_vec,
#       gene_z_vec   = gene_z_vec,
#       w_norm       = w_norm,
#       gene_map     = gene_map,
#       p_list       = p_list,
#       p_names      = p_names
#     ),
#     class = "magcat_core_prep"
#   )
# }

# ## ============================================================
# ##  1) ACAT (prepared)
# ## ============================================================

# #' @export
# magcat_acat_prepare <- function(gene_results,
#                                 pathways     = NULL,
#                                 species      = NULL,
#                                 pmn_gene_col = NULL,
#                                 gene_col     = "GENE",
#                                 p_col        = "P",
#                                 min_p        = 1e-15,
#                                 do_fix       = TRUE) {

#   if (!requireNamespace("ACAT", quietly = TRUE)) {
#     stop("magcat_acat_prepare(): requires ACAT package.", call. = FALSE)
#   }
#   if (!exists("fix_p_for_acat", mode = "function")) {
#     stop("magcat_acat_prepare(): missing fix_p_for_acat().", call. = FALSE)
#   }

#   core <- magcat_prepare_core(
#     gene_results  = gene_results,
#     pathways      = pathways,
#     species       = species,
#     pmn_gene_col  = pmn_gene_col,
#     gene_col      = gene_col,
#     p_col         = p_col,
#     z_col         = NULL,
#     weight_col    = NULL
#   )

#   structure(
#     c(core, list(min_p = min_p, do_fix = do_fix)),
#     class = c("magcat_acat_prep", class(core))
#   )
# }

# #' @export
# magcat_acat_run_prepared <- function(prep,
#                                      output  = FALSE,
#                                      out_dir = "magcat_acat") {

#   if (!inherits(prep, "magcat_acat_prep")) {
#     stop("magcat_acat_run_prepared(): 'prep' must come from magcat_acat_prepare().", call. = FALSE)
#   }

#   p_list  <- prep$p_list
#   p_names <- prep$p_names
#   pw_ids <- names(p_list)
#   pw_names <- pw_ids
#   if (!is.null(p_names)) {
#     tmp <- p_names[pw_ids]
#     tmp <- as.character(tmp)
#     pw_names <- ifelse(is.na(tmp) | tmp == "", pw_ids, unname(tmp))
#   }
#   gene_p  <- prep$gene_p_vec
#   gmap    <- prep$gene_map

#   n_pw <- length(p_list)
#   res  <- data.frame(
#     pathway_id   = pw_ids,
#     pathway_name = pw_names,
#     n_genes      = NA_integer_,
#     gene_names   = NA_character_,
#     gene_pvals   = NA_character_,
#     acat_p       = NA_real_,
#     stringsAsFactors = FALSE
#   )

#   for (i in seq_len(n_pw)) {

#     g_i <- p_list[[i]]
#     p_i <- gene_p[g_i]
#     keep <- !is.na(p_i) & is.finite(p_i) & p_i > 0 & p_i <= 1
#     p_i  <- as.numeric(p_i[keep])
#     g_use <- g_i[keep]

#     d <- length(p_i)
#     res$n_genes[i] <- d

#     if (d == 0L) next

#     canon <- unname(gmap[g_use])
#     canon <- canon[!is.na(canon) & canon != ""]
#     if (length(canon)) {
#       res$gene_names[i] <- paste(canon, collapse = ";")
#       res$gene_pvals[i] <- paste(format(p_i, digits = 6, scientific = TRUE), collapse = ";")
#     }

#     if (prep$do_fix) {
#       p_i <- fix_p_for_acat(p_i, min_p = prep$min_p)
#     }

#     res$acat_p[i] <- tryCatch(
#       as.numeric(ACAT::ACAT(Pvals = p_i)),
#       error = function(e) NA_real_
#     )
#   }

#   ord <- order(res$acat_p, decreasing = FALSE, na.last = TRUE)
#   res <- res[ord, , drop = FALSE]

#   if (output) {
#     if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#     out_path <- file.path(out_dir, "magcat_acat_prepared.csv")
#     utils::write.csv(res, out_path, row.names = FALSE)
#     attr(res, "file") <- out_path
#   }

#   res
# }

# ## ============================================================
# ##  2) Fisher (prepared)
# ## ============================================================

# #' @export
# magcat_fisher_prepare <- function(gene_results,
#                                   pathways     = NULL,
#                                   species      = NULL,
#                                   pmn_gene_col = NULL,
#                                   gene_col     = "GENE",
#                                   p_col        = "P",
#                                   min_p        = 1e-15,
#                                   do_fix       = TRUE) {

#   if (!exists("fix_p_for_acat", mode = "function")) {
#     stop("magcat_fisher_prepare(): missing fix_p_for_acat().", call. = FALSE)
#   }

#   core <- magcat_prepare_core(
#     gene_results  = gene_results,
#     pathways      = pathways,
#     species       = species,
#     pmn_gene_col  = pmn_gene_col,
#     gene_col      = gene_col,
#     p_col         = p_col,
#     z_col         = NULL,
#     weight_col    = NULL
#   )

#   structure(
#     c(core, list(min_p = min_p, do_fix = do_fix)),
#     class = c("magcat_fisher_prep", class(core))
#   )
# }

# #' @export
# magcat_fisher_run_prepared <- function(prep,
#                                        output  = FALSE,
#                                        out_dir = "magcat_fisher") {

#   if (!inherits(prep, "magcat_fisher_prep")) {
#     stop("magcat_fisher_run_prepared(): 'prep' must come from magcat_fisher_prepare().", call. = FALSE)
#   }

#   p_list  <- prep$p_list
#   p_names <- prep$p_names
#   pw_ids <- names(p_list)
#   pw_names <- pw_ids
#   if (!is.null(p_names)) {
#     tmp <- p_names[pw_ids]
#     tmp <- as.character(tmp)
#     pw_names <- ifelse(is.na(tmp) | tmp == "", pw_ids, unname(tmp))
#   }
#   gene_p  <- prep$gene_p_vec
#   gmap    <- prep$gene_map

#   n_pw <- length(p_list)
#   res  <- data.frame(
#     pathway_id   = pw_ids,
#     pathway_name = pw_names,
#     n_genes      = NA_integer_,
#     gene_names   = NA_character_,
#     gene_pvals   = NA_character_,
#     fisher_p     = NA_real_,
#     stringsAsFactors = FALSE
#   )

#   for (i in seq_len(n_pw)) {
#     g_i <- p_list[[i]]
#     p_i <- gene_p[g_i]
#     keep <- !is.na(p_i) & is.finite(p_i) & p_i > 0 & p_i <= 1
#     p_i  <- as.numeric(p_i[keep])
#     g_use <- g_i[keep]

#     d <- length(p_i)
#     res$n_genes[i] <- d
#     if (d < 2L) next

#     canon <- unname(gmap[g_use])
#     canon <- canon[!is.na(canon) & canon != ""]
#     res$gene_names[i] <- paste(canon, collapse = ";")
#     res$gene_pvals[i] <- paste(format(p_i, digits = 6, scientific = TRUE), collapse = ";")

#     if (prep$do_fix) p_i <- fix_p_for_acat(p_i, min_p = prep$min_p)

#     stat <- -2 * sum(log(p_i))
#     res$fisher_p[i] <- stats::pchisq(stat, df = 2 * d, lower.tail = FALSE)
#   }

#   ord <- order(res$fisher_p, decreasing = FALSE, na.last = TRUE)
#   res <- res[ord, , drop = FALSE]

#   if (output) {
#     if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#     out_path <- file.path(out_dir, "magcat_fisher_prepared.csv")
#     utils::write.csv(res, out_path, row.names = FALSE)
#     attr(res, "file") <- out_path
#   }

#   res
# }

# ## ============================================================
# ##  3) Adaptive soft TFisher (prepared)
# ## ============================================================

# #' @export
# magcat_soft_tfisher_adaptive_prepare <- function(gene_results,
#                                                  pathways     = NULL,
#                                                  species      = NULL,
#                                                  pmn_gene_col = NULL,
#                                                  gene_col     = "GENE",
#                                                  p_col        = "P",
#                                                  tau_grid     = c(0.20, 0.10, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001),
#                                                  min_p        = 1e-15,
#                                                  do_fix       = TRUE) {

#   if (!requireNamespace("TFisher", quietly = TRUE)) {
#     stop("magcat_soft_tfisher_adaptive_prepare(): requires TFisher package.", call. = FALSE)
#   }
#   if (!exists("fix_p_for_acat", mode = "function")) {
#     stop("magcat_soft_tfisher_adaptive_prepare(): missing fix_p_for_acat().", call. = FALSE)
#   }

#   tau_grid <- as.numeric(tau_grid)
#   tau_grid <- tau_grid[is.finite(tau_grid) & !is.na(tau_grid) & tau_grid > 0 & tau_grid < 1]
#   if (!length(tau_grid)) stop("tau_grid must contain values in (0,1).", call. = FALSE)
#   tau_grid <- sort(unique(tau_grid), decreasing = TRUE)

#   core <- magcat_prepare_core(
#     gene_results  = gene_results,
#     pathways      = pathways,
#     species       = species,
#     pmn_gene_col  = pmn_gene_col,
#     gene_col      = gene_col,
#     p_col         = p_col,
#     z_col         = NULL,
#     weight_col    = NULL
#   )

#   structure(
#     c(core, list(tau_grid = tau_grid, min_p = min_p, do_fix = do_fix)),
#     class = c("magcat_soft_tfisher_adaptive_prep", class(core))
#   )
# }

# #' @export
# magcat_soft_tfisher_adaptive_run_prepared <- function(prep,
#                                                       output  = FALSE,
#                                                       out_dir = "magcat_soft_tfisher_adaptive") {

#   if (!inherits(prep, "magcat_soft_tfisher_adaptive_prep")) {
#     stop("..._run_prepared(): 'prep' must come from magcat_soft_tfisher_adaptive_prepare().", call. = FALSE)
#   }

#   p_list  <- prep$p_list
#   p_names <- prep$p_names
#   pw_ids <- names(p_list)
#   pw_names <- pw_ids
#   if (!is.null(p_names)) {
#     tmp <- p_names[pw_ids]
#     tmp <- as.character(tmp)
#     pw_names <- ifelse(is.na(tmp) | tmp == "", pw_ids, unname(tmp))
#   }
#   gene_p  <- prep$gene_p_vec
#   gmap    <- prep$gene_map
#   tau_grid <- prep$tau_grid

#   n_pw <- length(p_list)
#   res  <- data.frame(
#     pathway_id        = pw_ids,
#     pathway_name      = pw_names,
#     n_genes           = NA_integer_,
#     gene_names        = NA_character_,
#     gene_pvals        = NA_character_,
#     tfisher_stat_hat  = NA_real_,
#     tfisher_p_analytic= NA_real_,
#     tau_hat           = NA_real_,
#     stringsAsFactors  = FALSE
#   )

#   for (i in seq_len(n_pw)) {

#     g_i <- p_list[[i]]
#     p_i <- gene_p[g_i]
#     keep <- !is.na(p_i) & is.finite(p_i) & p_i > 0 & p_i <= 1
#     p_i  <- as.numeric(p_i[keep])
#     g_use <- g_i[keep]

#     d <- length(p_i)
#     res$n_genes[i] <- d
#     if (d < 2L) next

#     canon <- unname(gmap[g_use])
#     canon <- canon[!is.na(canon) & canon != ""]
#     res$gene_names[i] <- paste(canon, collapse = ";")
#     res$gene_pvals[i] <- paste(format(p_i, digits = 6, scientific = TRUE), collapse = ";")

#     if (prep$do_fix) p_i <- fix_p_for_acat(p_i, min_p = prep$min_p)

#     # evaluate tau grid
#     best_p <- NA_real_
#     best_tau <- NA_real_
#     best_stat <- NA_real_

#     for (tau in tau_grid) {
#       st <- TFisher::stat.soft(p = p_i, tau1 = tau)

#       # p.soft returns CDF; right-tail p = 1 - CDF
#       p_tau <- 1 - as.numeric(TFisher::p.soft(q = st, n = d, tau1 = tau, M = NULL))

#       if (!is.finite(p_tau) || is.na(p_tau)) next
#       if (is.na(best_p) || p_tau < best_p) {
#         best_p <- p_tau
#         best_tau <- tau
#         best_stat <- st
#       }
#     }

#     res$tfisher_stat_hat[i]   <- best_stat
#     res$tfisher_p_analytic[i] <- best_p
#     res$tau_hat[i]            <- best_tau
#   }

#   ord <- order(res$tfisher_p_analytic, decreasing = FALSE, na.last = TRUE)
#   res <- res[ord, , drop = FALSE]

#   if (output) {
#     if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#     out_path <- file.path(out_dir, "magcat_soft_tfisher_adaptive_prepared.csv")
#     utils::write.csv(res, out_path, row.names = FALSE)
#     attr(res, "file") <- out_path
#   }

#   res
# }

# ## ============================================================
# ##  4) minP (prepared)
# ## ============================================================

# #' @export
# magcat_minp_prepare <- function(gene_results,
#                                 pathways     = NULL,
#                                 species      = NULL,
#                                 pmn_gene_col = NULL,
#                                 gene_col     = "GENE",
#                                 p_col        = "P",
#                                 min_p        = 1e-15,
#                                 do_fix       = TRUE) {

#   if (!exists("fix_p_for_acat", mode = "function")) {
#     stop("magcat_minp_prepare(): missing fix_p_for_acat().", call. = FALSE)
#   }

#   core <- magcat_prepare_core(
#     gene_results  = gene_results,
#     pathways      = pathways,
#     species       = species,
#     pmn_gene_col  = pmn_gene_col,
#     gene_col      = gene_col,
#     p_col         = p_col,
#     z_col         = NULL,
#     weight_col    = NULL
#   )

#   structure(
#     c(core, list(min_p = min_p, do_fix = do_fix)),
#     class = c("magcat_minp_prep", class(core))
#   )
# }

# #' @export
# magcat_minp_run_prepared <- function(prep,
#                                      output  = FALSE,
#                                      out_dir = "magcat_minp") {

#   if (!inherits(prep, "magcat_minp_prep")) {
#     stop("magcat_minp_run_prepared(): 'prep' must come from magcat_minp_prepare().", call. = FALSE)
#   }

#   .analytic_minp <- function(p_vec) {
#     p_vec <- as.numeric(p_vec)
#     p_vec <- p_vec[!is.na(p_vec) & p_vec > 0 & p_vec <= 1]
#     if (!length(p_vec)) return(NA_real_)
#     if (requireNamespace("metap", quietly = TRUE)) {
#       out <- tryCatch(metap::minimump(p_vec)$p, error = function(e) NA_real_)
#       return(out)
#     }
#     k <- length(p_vec)
#     pmin <- min(p_vec)
#     1 - (1 - pmin)^k
#   }

#   p_list  <- prep$p_list
#   p_names <- prep$p_names
#   pw_ids <- names(p_list)
#   pw_names <- pw_ids
#   if (!is.null(p_names)) {
#     tmp <- p_names[pw_ids]
#     tmp <- as.character(tmp)
#     pw_names <- ifelse(is.na(tmp) | tmp == "", pw_ids, unname(tmp))
#   }
#   gene_p  <- prep$gene_p_vec
#   gmap    <- prep$gene_map

#   n_pw <- length(p_list)
#   res  <- data.frame(
#     pathway_id      = pw_ids,
#     pathway_name    = pw_names,
#     n_genes         = NA_integer_,
#     gene_names      = NA_character_,
#     gene_pvals      = NA_character_,
#     minp_stat       = NA_real_,
#     minp_p_analytic = NA_real_,
#     stringsAsFactors = FALSE
#   )

#   for (i in seq_len(n_pw)) {
#     g_i <- p_list[[i]]
#     p_i <- gene_p[g_i]
#     keep <- !is.na(p_i) & is.finite(p_i) & p_i > 0 & p_i <= 1
#     p_i  <- as.numeric(p_i[keep])
#     g_use <- g_i[keep]

#     d <- length(p_i)
#     res$n_genes[i] <- d
#     if (d < 2L) next

#     canon <- unname(gmap[g_use])
#     canon <- canon[!is.na(canon) & canon != ""]
#     res$gene_names[i] <- paste(canon, collapse = ";")
#     res$gene_pvals[i] <- paste(format(p_i, digits = 6, scientific = TRUE), collapse = ";")

#     if (prep$do_fix) p_i <- fix_p_for_acat(p_i, min_p = prep$min_p)

#     res$minp_stat[i]       <- min(p_i, na.rm = TRUE)
#     res$minp_p_analytic[i] <- .analytic_minp(p_i)
#   }

#   ord <- order(res$minp_p_analytic, decreasing = FALSE, na.last = TRUE)
#   res <- res[ord, , drop = FALSE]

#   if (output) {
#     if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#     out_path <- file.path(out_dir, "magcat_minp_prepared.csv")
#     utils::write.csv(res, out_path, row.names = FALSE)
#     attr(res, "file") <- out_path
#   }

#   res
# }

# ## ============================================================
# ##  5) Stouffer Z (prepared)
# ## ============================================================

# #' @export
# magcat_stoufferZ_prepare <- function(gene_results,
#                                      pathways     = NULL,
#                                      species      = NULL,
#                                      pmn_gene_col = NULL,
#                                      gene_col     = "GENE",
#                                      p_col        = "P",
#                                      z_col        = "ZSTAT",
#                                      weight_col   = NULL,
#                                      min_abs_w    = 1e-8) {

#   core <- magcat_prepare_core(
#     gene_results  = gene_results,
#     pathways      = pathways,
#     species       = species,
#     pmn_gene_col  = pmn_gene_col,
#     gene_col      = gene_col,
#     p_col         = p_col,
#     z_col         = z_col,
#     weight_col    = weight_col
#   )

#   structure(
#     c(core, list(min_abs_w = min_abs_w)),
#     class = c("magcat_stoufferZ_prep", class(core))
#   )
# }

# #' @export
# magcat_stoufferZ_run_prepared <- function(prep,
#                                           output  = FALSE,
#                                           out_dir = "magcat_stouffer_z") {

#   if (!inherits(prep, "magcat_stoufferZ_prep")) {
#     stop("magcat_stoufferZ_run_prepared(): 'prep' must come from magcat_stoufferZ_prepare().", call. = FALSE)
#   }

#   z_list  <- prep$p_list
#   p_names <- prep$p_names
#   pw_ids <- names(z_list)
#   pw_names <- pw_ids
#   if (!is.null(p_names)) {
#     tmp <- p_names[pw_ids]
#     tmp <- as.character(tmp)
#     pw_names <- ifelse(is.na(tmp) | tmp == "", pw_ids, unname(tmp))
#   }
#   z_vec   <- prep$gene_z_vec
#   w_vec   <- prep$w_norm
#   gmap    <- prep$gene_map
#   min_abs_w <- prep$min_abs_w

#   if (is.null(z_vec)) stop("Prepared object has no gene_z_vec; check z_col.", call. = FALSE)

#   .fix_w <- function(w) {
#     w <- as.numeric(w)
#     bad <- !is.finite(w) | is.na(w) | w <= 0
#     if (any(bad)) w[bad] <- min_abs_w
#     w
#   }

#   .stouffer_from_z <- function(z_i, w_i = NULL) {
#     z_i <- as.numeric(z_i)
#     keep <- is.finite(z_i) & !is.na(z_i)
#     z_i <- z_i[keep]
#     if (!length(z_i)) return(list(Z = NA_real_, p = NA_real_))

#     if (is.null(w_i)) {
#       w_i <- rep(1, length(z_i))
#     } else {
#       w_i <- .fix_w(w_i)
#       if (length(w_i) != length(z_i)) w_i <- rep(1, length(z_i))
#     }

#     Zs <- sum(w_i * z_i) / sqrt(sum(w_i^2))
#     ps <- 2 * stats::pnorm(-abs(Zs))
#     list(Z = Zs, p = ps)
#   }

#   n_pw <- length(z_list)
#   res  <- data.frame(
#     pathway_id          = pw_ids,
#     pathway_name        = pw_names,
#     n_genes             = NA_integer_,
#     gene_names          = NA_character_,
#     gene_zvals          = NA_character_,
#     stouffer_Z          = NA_real_,
#     stouffer_p_analytic = NA_real_,
#     stringsAsFactors    = FALSE
#   )

#   for (i in seq_len(n_pw)) {

#     g_i <- z_list[[i]]
#     z_i <- z_vec[g_i]
#     keep <- !is.na(z_i) & is.finite(z_i)
#     z_i  <- as.numeric(z_i[keep])
#     g_use <- g_i[keep]

#     d <- length(z_i)
#     res$n_genes[i] <- d
#     if (d < 2L) next

#     canon <- unname(gmap[g_use])
#     canon <- canon[!is.na(canon) & canon != ""]
#     res$gene_names[i] <- paste(canon, collapse = ";")
#     res$gene_zvals[i] <- paste(format(z_i, digits = 5, scientific = TRUE), collapse = ";")

#     w_i <- NULL
#     if (!is.null(w_vec)) {
#       w_i <- as.numeric(w_vec[g_use])
#     }

#     st <- .stouffer_from_z(z_i, w_i)
#     res$stouffer_Z[i]          <- st$Z
#     res$stouffer_p_analytic[i] <- st$p
#   }

#   ord <- order(res$stouffer_p_analytic, decreasing = FALSE, na.last = TRUE)
#   res <- res[ord, , drop = FALSE]

#   if (output) {
#     if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#     out_path <- file.path(out_dir, "magcat_stoufferZ_prepared.csv")
#     utils::write.csv(res, out_path, row.names = FALSE)
#     attr(res, "file") <- out_path
#   }

#   res
# }

# ## ============================================================
# ##  6) MAGMA competitive gene-set (prepared)
# ## ============================================================

# #' @export
# magcat_magma_competitive_prepare <- function(gene_results,
#                                              gene_results_raw,
#                                              set_annot        = NULL,
#                                              pathways         = NULL,
#                                              species          = NULL,
#                                              pmn_gene_col     = NULL,
#                                              gene_col         = "GENE",
#                                              p_col            = "P",
#                                              out_prefix       = "MAGMA_GENESET",
#                                              out_dir          = "magma_geneset",
#                                              write_tidy       = TRUE,
#                                              tidy_suffix      = ".tidy.tsv") {

#   if (!exists("magma_geneset_competitive", mode = "function")) {
#     stop("magcat_magma_competitive_prepare(): missing magma_geneset_competitive().", call. = FALSE)
#   }
#   if (!file.exists(gene_results_raw)) {
#     stop("gene_results_raw file does not exist: ", gene_results_raw, call. = FALSE)
#   }

#   # pathways long table
#   if (!is.null(pathways) && !is.null(species)) {
#     stop("Provide either 'pathways' OR 'species', not both.", call. = FALSE)
#   }
#   if (is.null(pathways) && is.null(species)) {
#     stop("Need pathways (data.frame/list) or species for MAGMA competitive.", call. = FALSE)
#   }
#   if (is.null(pathways)) {
#     if (!exists("magcat_load_pathways", mode = "function")) {
#       stop("magcat_magma_competitive_prepare(): missing magcat_load_pathways().", call. = FALSE)
#     }
#     pathways <- magcat_load_pathways(species = species, gene_col = pmn_gene_col)
#   }

#   # ensure long data.frame for magma wrapper
#   if (is.list(pathways) && !is.data.frame(pathways)) {
#     pw_id <- names(pathways)
#     if (is.null(pw_id)) pw_id <- paste0("PWY_", seq_along(pathways))
#     pathways <- data.frame(
#       pathway_id   = rep(pw_id, lengths(pathways)),
#       pathway_name = rep(pw_id, lengths(pathways)),
#       gene_id      = unlist(pathways, use.names = FALSE),
#       stringsAsFactors = FALSE
#     )
#   }

#   if (!is.data.frame(pathways) || !all(c("pathway_id","gene_id") %in% names(pathways))) {
#     stop("pathways must be a data.frame with pathway_id and gene_id (and optional pathway_name).", call. = FALSE)
#   }
#   if (!("pathway_name" %in% names(pathways))) pathways$pathway_name <- pathways$pathway_id

#   # genes_all table used for nice tidy output (gene_names + gene_pvals)
#   if (!all(c(gene_col, p_col) %in% names(gene_results))) {
#     stop("gene_results must have columns '", gene_col, "' and '", p_col, "' for genes_all.", call. = FALSE)
#   }
#   genes_all <- data.frame(
#     GENE = as.character(gene_results[[gene_col]]),
#     P    = suppressWarnings(as.numeric(gene_results[[p_col]])),
#     stringsAsFactors = FALSE
#   )

#   structure(
#     list(
#       gene_results_raw = gene_results_raw,
#       set_annot        = set_annot,
#       pathways         = pathways,
#       species          = species,
#       pmn_gene_col     = pmn_gene_col,
#       genes_all        = genes_all,
#       out_prefix       = out_prefix,
#       out_dir          = out_dir,
#       write_tidy       = write_tidy,
#       tidy_suffix      = tidy_suffix
#     ),
#     class = "magcat_magma_competitive_prep"
#   )
# }

# #' @export
# magcat_magma_competitive_run_prepared <- function(prep,
#                                                   output_dir    = NULL,
#                                                   output_prefix = NULL,
#                                                   write_tidy    = NULL) {

#   if (!inherits(prep, "magcat_magma_competitive_prep")) {
#     stop("magcat_magma_competitive_run_prepared(): 'prep' must come from magcat_magma_competitive_prepare().",
#          call. = FALSE)
#   }
#   if (!exists("magma_geneset_competitive", mode = "function")) {
#     stop("magcat_magma_competitive_run_prepared(): missing magma_geneset_competitive().", call. = FALSE)
#   }

#   # allow overrides
#   out_dir  <- if (is.null(output_dir))    prep$out_dir    else output_dir
#   out_pref <- if (is.null(output_prefix)) prep$out_prefix else output_prefix
#   wt       <- if (is.null(write_tidy))    prep$write_tidy else write_tidy

#   out <- magma_geneset_competitive(
#     gene_results_raw = prep$gene_results_raw,
#     set_annot        = prep$set_annot,
#     out_prefix       = out_pref,
#     out_dir          = out_dir,
#     species          = prep$species,
#     pmn_gene_col     = prep$pmn_gene_col,
#     pathways         = prep$pathways,
#     genes_all        = prep$genes_all,
#     write_tidy       = isTRUE(wt),
#     tidy_suffix      = prep$tidy_suffix
#   )

#   out
# }







# # ## ============================================================
# # ##  Shared "freeze" core used by Fisher / minP / softTFisher / Stouffer
# # ## ============================================================

# # #' Prepare common pathway + gene maps for fast reuse
# # #'
# # #' This freezes the expensive, repeated setup steps you currently redo
# # #' inside every wrapper call: pathway list building, gene normalization,
# # #' named p-vector creation, canonical gene-id map, and permutation pools.
# # #'
# # #' Requires your existing:
# # #'   - magcat_load_pathways()
# # #'   - fix_p_for_acat()
# # #'
# # #' @keywords internal
# # magcat_prepare_core <- function(gene_results,
# #                                 pathways     = NULL,
# #                                 species      = NULL,
# #                                 pmn_gene_col = NULL,
# #                                 gene_col     = "GENE",
# #                                 p_col        = "P",
# #                                 weight_col   = NULL) {

# #   ## pathway source
# #   if (!is.null(pathways) && !is.null(species)) {
# #     stop("Provide either 'pathways' OR 'species', not both.", call. = FALSE)
# #   }
# #   if (is.null(pathways) && is.null(species)) {
# #     stop("You must provide either 'pathways' or 'species'.", call. = FALSE)
# #   }
# #   if (is.null(pathways) && !is.null(species)) {
# #     if (!exists("magcat_load_pathways", mode = "function")) {
# #       stop("magcat_prepare_core(): missing magcat_load_pathways().", call. = FALSE)
# #     }
# #     pathways <- magcat_load_pathways(
# #       species  = species,
# #       gene_col = pmn_gene_col
# #     )
# #   }

# #   ## gene_results checks
# #   if (!all(c(gene_col, p_col) %in% names(gene_results))) {
# #     stop("gene_results must contain columns '", gene_col, "' and '", p_col, "'.",
# #          call. = FALSE)
# #   }

# #   gr         <- gene_results
# #   genes_all  <- as.character(gr[[gene_col]])
# #   p_all      <- as.numeric(gr[[p_col]])
# #   genes_norm <- tolower(genes_all)

# #   ## named vectors
# #   gene_p_vec <- stats::setNames(p_all, genes_norm)
# #   gene_map   <- tapply(genes_all, genes_norm, function(x) x[1])

# #   ## optional weights (for Stouffer)
# #   w_norm <- NULL
# #   if (!is.null(weight_col)) {
# #     if (!weight_col %in% names(gr)) {
# #       stop("weight_col = '", weight_col, "' not found in gene_results.", call. = FALSE)
# #     }
# #     w_all  <- as.numeric(gr[[weight_col]])
# #     w_norm <- stats::setNames(w_all, genes_norm)
# #   }

# #   ## pathways -> named list + names
# #   if (is.data.frame(pathways)) {
# #     if (!all(c("pathway_id", "gene_id") %in% names(pathways))) {
# #       stop("If 'pathways' is a data.frame it must have columns pathway_id and gene_id.",
# #            call. = FALSE)
# #     }
# #     if (!"pathway_name" %in% names(pathways)) {
# #       pathways$pathway_name <- pathways$pathway_id
# #     }
# #     p_list  <- split(pathways$gene_id, pathways$pathway_id)
# #     p_names <- tapply(pathways$pathway_name, pathways$pathway_id, function(x) x[1])
# #   } else if (is.list(pathways)) {
# #     p_list  <- pathways
# #     p_names <- names(pathways)
# #     if (is.null(p_names)) {
# #       p_names <- paste0("PWY_", seq_along(pathways))
# #       names(p_list) <- p_names
# #     }
# #   } else {
# #     stop("'pathways' must be either a list or a data.frame.", call. = FALSE)
# #   }

# #   ## normalize IDs inside pathways
# #   p_list <- lapply(p_list, function(g) tolower(as.character(g)))

# #   ## permutation pools (store both styles)
# #   idx_pool_na     <- which(!is.na(p_all))
# #   idx_pool_valid  <- which(is.finite(p_all) & p_all > 0 & p_all < 1)

# #   if (!length(idx_pool_na)) {
# #     stop("No non-NA gene-level p-values found in gene_results.", call. = FALSE)
# #   }
# #   if (!length(idx_pool_valid)) {
# #     ## some wrappers use valid01 pool (e.g. Stouffer); keep NA pool anyway
# #     idx_pool_valid <- idx_pool_na
# #   }

# #   structure(
# #     list(
# #       gene_col     = gene_col,
# #       p_col        = p_col,
# #       weight_col   = weight_col,
# #       genes_all    = genes_all,
# #       genes_norm   = genes_norm,
# #       p_all        = p_all,
# #       gene_p_vec   = gene_p_vec,
# #       gene_map     = gene_map,
# #       w_norm       = w_norm,
# #       p_list       = p_list,
# #       p_names      = p_names,
# #       idx_pool_na  = idx_pool_na,
# #       idx_pool_01  = idx_pool_valid
# #     ),
# #     class = "magcat_core_prep"
# #   )
# # }


# # ## ============================================================
# # ##  Fisher: prepare + run_prepared
# # ## ============================================================

# # #' Prepare Fisher pathway test objects
# # #' @export
# # magcat_fisher_prepare <- function(gene_results,
# #                                   pathways     = NULL,
# #                                   species      = NULL,
# #                                   pmn_gene_col = NULL,
# #                                   gene_col     = "GENE",
# #                                   p_col        = "P",
# #                                   min_p        = 1e-15,
# #                                   do_fix       = TRUE) {

# #   core <- magcat_prepare_core(
# #     gene_results  = gene_results,
# #     pathways      = pathways,
# #     species       = species,
# #     pmn_gene_col  = pmn_gene_col,
# #     gene_col      = gene_col,
# #     p_col         = p_col,
# #     weight_col    = NULL
# #   )

# #   structure(
# #     c(core, list(min_p = min_p, do_fix = do_fix)),
# #     class = c("magcat_fisher_prep", class(core))
# #   )
# # }

# # #' Run Fisher pathway test from a prepared object
# # #' @export
# # magcat_fisher_run_prepared <- function(prep,
# #                                        output  = FALSE,
# #                                        out_dir = "magcat_fisher") {

# #   if (!inherits(prep, "magcat_fisher_prep")) {
# #     stop("magcat_fisher_run_prepared(): 'prep' must come from magcat_fisher_prepare().",
# #          call. = FALSE)
# #   }
# #   if (!exists("fix_p_for_acat", mode = "function")) {
# #     stop("magcat_fisher_run_prepared(): missing fix_p_for_acat().", call. = FALSE)
# #   }

# #   p_list  <- prep$p_list
# #   p_names <- prep$p_names
# #   gene_p  <- prep$gene_p_vec
# #   gmap    <- prep$gene_map

# #   n_pw <- length(p_list)
# #   res  <- data.frame(
# #     pathway_id   = names(p_list),
# #     pathway_name = if (is.null(p_names)) names(p_list) else unname(p_names),
# #     n_genes      = NA_integer_,
# #     gene_names   = NA_character_,
# #     fisher_p     = NA_real_,
# #     stringsAsFactors = FALSE
# #   )

# #   for (i in seq_len(n_pw)) {
# #     g_i <- p_list[[i]]
# #     p_i <- gene_p[g_i]
# #     keep <- !is.na(p_i)
# #     p_i  <- as.numeric(p_i[keep])
# #     g_use <- g_i[keep]

# #     d <- length(p_i)
# #     res$n_genes[i] <- d
# #     if (d == 0L) next

# #     canon <- unname(gmap[g_use])
# #     canon <- unique(canon[!is.na(canon)])
# #     res$gene_names[i] <- paste(canon, collapse = ";")

# #     if (prep$do_fix) {
# #       p_i <- fix_p_for_acat(p_i, min_p = prep$min_p)
# #     }

# #     stat <- -2 * sum(log(p_i))
# #     res$fisher_p[i] <- stats::pchisq(stat, df = 2 * length(p_i), lower.tail = FALSE)
# #   }

# #   ord <- order(res$fisher_p, decreasing = FALSE, na.last = TRUE)
# #   res <- res[ord, , drop = FALSE]

# #   if (output) {
# #     if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# #     utils::write.csv(res, file.path(out_dir, "magcat_fisher_prepared.csv"), row.names = FALSE)
# #   }
# #   res
# # }


# # ## ============================================================
# # ##  minP: prepare + run_prepared (analytic only; perm handled in omni)
# # ## ============================================================

# # #' Prepare minP pathway test objects
# # #' @export
# # magcat_minp_prepare <- function(gene_results,
# #                                 pathways     = NULL,
# #                                 species      = NULL,
# #                                 pmn_gene_col = NULL,
# #                                 gene_col     = "GENE",
# #                                 p_col        = "P",
# #                                 min_p        = 1e-15,
# #                                 do_fix       = TRUE) {

# #   core <- magcat_prepare_core(
# #     gene_results  = gene_results,
# #     pathways      = pathways,
# #     species       = species,
# #     pmn_gene_col  = pmn_gene_col,
# #     gene_col      = gene_col,
# #     p_col         = p_col,
# #     weight_col    = NULL
# #   )

# #   structure(
# #     c(core, list(min_p = min_p, do_fix = do_fix)),
# #     class = c("magcat_minp_prep", class(core))
# #   )
# # }

# # #' Run minP pathway test from a prepared object (analytic)
# # #' @export
# # magcat_minp_run_prepared <- function(prep,
# #                                      output  = FALSE,
# #                                      out_dir = "magcat_minp") {

# #   if (!inherits(prep, "magcat_minp_prep")) {
# #     stop("magcat_minp_run_prepared(): 'prep' must come from magcat_minp_prepare().",
# #          call. = FALSE)
# #   }
# #   if (!exists("fix_p_for_acat", mode = "function")) {
# #     stop("magcat_minp_run_prepared(): missing fix_p_for_acat().", call. = FALSE)
# #   }

# #   .analytic_minp <- function(p_vec) {
# #     p_vec <- as.numeric(p_vec)
# #     p_vec <- p_vec[!is.na(p_vec) & p_vec > 0 & p_vec < 1]
# #     if (!length(p_vec)) return(NA_real_)
# #     if (requireNamespace("metap", quietly = TRUE)) {
# #       return(tryCatch(metap::minimump(p_vec)$p, error = function(e) NA_real_))
# #     }
# #     k <- length(p_vec)
# #     p_min <- min(p_vec)
# #     1 - (1 - p_min)^k
# #   }

# #   p_list  <- prep$p_list
# #   p_names <- prep$p_names
# #   gene_p  <- prep$gene_p_vec
# #   gmap    <- prep$gene_map

# #   n_pw <- length(p_list)
# #   res  <- data.frame(
# #     pathway_id      = names(p_list),
# #     pathway_name    = if (is.null(p_names)) names(p_list) else unname(p_names),
# #     n_genes         = NA_integer_,
# #     gene_names      = NA_character_,
# #     minp_stat       = NA_real_,
# #     minp_p_analytic = NA_real_,
# #     stringsAsFactors = FALSE
# #   )

# #   for (i in seq_len(n_pw)) {
# #     g_i <- p_list[[i]]
# #     p_i <- gene_p[g_i]
# #     keep <- !is.na(p_i)
# #     p_i  <- as.numeric(p_i[keep])
# #     g_use <- g_i[keep]

# #     d <- length(p_i)
# #     res$n_genes[i] <- d
# #     if (d == 0L) next

# #     canon <- unname(gmap[g_use])
# #     canon <- unique(canon[!is.na(canon)])
# #     res$gene_names[i] <- paste(canon, collapse = ";")

# #     if (prep$do_fix) {
# #       p_i <- fix_p_for_acat(p_i, min_p = prep$min_p)
# #     }

# #     res$minp_stat[i]       <- min(p_i)
# #     res$minp_p_analytic[i] <- .analytic_minp(p_i)
# #   }

# #   ord <- order(res$minp_p_analytic, decreasing = FALSE, na.last = TRUE)
# #   res <- res[ord, , drop = FALSE]

# #   if (output) {
# #     if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# #     utils::write.csv(res, file.path(out_dir, "magcat_minp_prepared.csv"), row.names = FALSE)
# #   }
# #   res
# # }


# # ## ============================================================
# # ##  soft TFisher: prepare + run_prepared (analytic only; perm handled in omni)
# # ## ============================================================

# # #' Prepare soft-TFisher pathway test objects
# # #' @export
# # magcat_soft_tfisher_prepare <- function(gene_results,
# #                                         pathways     = NULL,
# #                                         species      = NULL,
# #                                         pmn_gene_col = NULL,
# #                                         gene_col     = "GENE",
# #                                         p_col        = "P",
# #                                         tau1         = 0.05,
# #                                         min_p        = 1e-15,
# #                                         do_fix       = TRUE) {

# #   if (!requireNamespace("TFisher", quietly = TRUE)) {
# #     stop("magcat_soft_tfisher_prepare(): requires TFisher package.", call. = FALSE)
# #   }

# #   core <- magcat_prepare_core(
# #     gene_results  = gene_results,
# #     pathways      = pathways,
# #     species       = species,
# #     pmn_gene_col  = pmn_gene_col,
# #     gene_col      = gene_col,
# #     p_col         = p_col,
# #     weight_col    = NULL
# #   )

# #   structure(
# #     c(core, list(tau1 = tau1, min_p = min_p, do_fix = do_fix)),
# #     class = c("magcat_soft_tfisher_prep", class(core))
# #   )
# # }

# # #' Run soft-TFisher pathway test from a prepared object (analytic)
# # #' @export
# # magcat_soft_tfisher_run_prepared <- function(prep,
# #                                              output  = FALSE,
# #                                              out_dir = "magcat_tfisher_soft") {

# #   if (!inherits(prep, "magcat_soft_tfisher_prep")) {
# #     stop("magcat_soft_tfisher_run_prepared(): 'prep' must come from magcat_soft_tfisher_prepare().",
# #          call. = FALSE)
# #   }
# #   if (!exists("fix_p_for_acat", mode = "function")) {
# #     stop("magcat_soft_tfisher_run_prepared(): missing fix_p_for_acat().", call. = FALSE)
# #   }

# #   p_list  <- prep$p_list
# #   p_names <- prep$p_names
# #   gene_p  <- prep$gene_p_vec
# #   gmap    <- prep$gene_map

# #   n_pw <- length(p_list)
# #   res  <- data.frame(
# #     pathway_id         = names(p_list),
# #     pathway_name       = if (is.null(p_names)) names(p_list) else unname(p_names),
# #     n_genes            = NA_integer_,
# #     gene_names         = NA_character_,
# #     tfisher_stat       = NA_real_,
# #     tfisher_p_analytic = NA_real_,
# #     stringsAsFactors = FALSE
# #   )

# #   for (i in seq_len(n_pw)) {
# #     g_i <- p_list[[i]]
# #     p_i <- gene_p[g_i]
# #     keep <- !is.na(p_i)
# #     p_i  <- as.numeric(p_i[keep])
# #     g_use <- g_i[keep]

# #     d <- length(p_i)
# #     res$n_genes[i] <- d
# #     if (d == 0L) next

# #     canon <- unname(gmap[g_use])
# #     canon <- unique(canon[!is.na(canon)])
# #     res$gene_names[i] <- paste(canon, collapse = ";")

# #     if (prep$do_fix) {
# #       p_i <- fix_p_for_acat(p_i, min_p = prep$min_p)
# #     }

# #     stat <- TFisher::stat.soft(p = p_i, tau1 = prep$tau1)
# #     res$tfisher_stat[i] <- stat

# #     Fq <- TFisher::p.soft(q = stat, n = length(p_i), tau1 = prep$tau1, M = NULL)
# #     res$tfisher_p_analytic[i] <- 1 - as.numeric(Fq)
# #   }

# #   ord <- order(res$tfisher_p_analytic, decreasing = FALSE, na.last = TRUE)
# #   res <- res[ord, , drop = FALSE]

# #   if (output) {
# #     if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# #     utils::write.csv(res, file.path(out_dir, "magcat_tfisher_soft_prepared.csv"), row.names = FALSE)
# #   }
# #   res
# # }


# # ## ============================================================
# # ##  Stouffer: prepare + run_prepared (analytic only; perm handled in omni)
# # ## ============================================================

# # #' Prepare Stouffer pathway test objects
# # #' @export
# # magcat_stouffer_prepare <- function(gene_results,
# #                                     pathways     = NULL,
# #                                     species      = NULL,
# #                                     pmn_gene_col = NULL,
# #                                     gene_col     = "GENE",
# #                                     p_col        = "P",
# #                                     weight_col   = NULL,
# #                                     min_p        = 1e-15,
# #                                     do_fix       = TRUE) {

# #   if (!requireNamespace("metap", quietly = TRUE)) {
# #     stop("magcat_stouffer_prepare(): requires metap package.", call. = FALSE)
# #   }

# #   core <- magcat_prepare_core(
# #     gene_results  = gene_results,
# #     pathways      = pathways,
# #     species       = species,
# #     pmn_gene_col  = pmn_gene_col,
# #     gene_col      = gene_col,
# #     p_col         = p_col,
# #     weight_col    = weight_col
# #   )

# #   structure(
# #     c(core, list(min_p = min_p, do_fix = do_fix)),
# #     class = c("magcat_stouffer_prep", class(core))
# #   )
# # }

# # #' Run Stouffer pathway test from a prepared object (analytic)
# # #' @export
# # magcat_stouffer_run_prepared <- function(prep,
# #                                          output  = FALSE,
# #                                          out_dir = "magcat_stouffer") {

# #   if (!inherits(prep, "magcat_stouffer_prep")) {
# #     stop("magcat_stouffer_run_prepared(): 'prep' must come from magcat_stouffer_prepare().",
# #          call. = FALSE)
# #   }
# #   if (!exists("fix_p_for_acat", mode = "function")) {
# #     stop("magcat_stouffer_run_prepared(): missing fix_p_for_acat().", call. = FALSE)
# #   }

# #   p_list  <- prep$p_list
# #   p_names <- prep$p_names
# #   gene_p  <- prep$gene_p_vec
# #   gmap    <- prep$gene_map
# #   w_norm  <- prep$w_norm

# #   n_pw <- length(p_list)
# #   res  <- data.frame(
# #     pathway_id          = names(p_list),
# #     pathway_name        = if (is.null(p_names)) names(p_list) else unname(p_names),
# #     n_genes             = NA_integer_,
# #     gene_names          = NA_character_,
# #     stouffer_z          = NA_real_,
# #     stouffer_p_analytic = NA_real_,
# #     stringsAsFactors = FALSE
# #   )

# #   for (i in seq_len(n_pw)) {
# #     g_i <- p_list[[i]]
# #     p_i <- gene_p[g_i]
# #     keep <- !is.na(p_i)
# #     p_i  <- as.numeric(p_i[keep])
# #     g_use <- g_i[keep]

# #     d <- length(p_i)
# #     res$n_genes[i] <- d
# #     if (d < 2L) next

# #     canon <- unname(gmap[g_use])
# #     canon <- unique(canon[!is.na(canon)])
# #     res$gene_names[i] <- paste(canon, collapse = ";")

# #     if (prep$do_fix) {
# #       p_i <- fix_p_for_acat(p_i, min_p = prep$min_p)
# #     }

# #     ## weights
# #     if (!is.null(w_norm)) {
# #       w_i <- as.numeric(w_norm[g_use])
# #       if (any(is.na(w_i) | w_i <= 0)) {
# #         w_pos <- w_i[!is.na(w_i) & w_i > 0]
# #         repl  <- if (length(w_pos)) stats::median(w_pos) else 1
# #         w_i[is.na(w_i) | w_i <= 0] <- repl
# #       }
# #     } else {
# #       w_i <- rep(1, length(p_i))
# #     }

# #     sz <- tryCatch(metap::sumz(p = p_i, weights = w_i), error = function(e) NULL)
# #     if (!is.null(sz)) {
# #       res$stouffer_z[i]          <- as.numeric(sz$z)
# #       res$stouffer_p_analytic[i] <- as.numeric(sz$p)
# #     }
# #   }

# #   ord <- order(res$stouffer_p_analytic, decreasing = FALSE, na.last = TRUE)
# #   res <- res[ord, , drop = FALSE]

# #   if (output) {
# #     if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# #     utils::write.csv(res, file.path(out_dir, "magcat_stouffer_prepared.csv"), row.names = FALSE)
# #   }
# #   res
# # }


# # ## ============================================================
# # ##  ACAT
# # ## ============================================================

# # #' Run ACAT pathway test from a prepared object
# # #' @export
# # magcat_acat_run_prepared <- function(prep,
# #                                      output  = FALSE,
# #                                      out_dir = "magcat_acat") {

# #   if (is.null(prep$idx_list) || is.null(prep$p_all_u) || is.null(prep$pathway_id)) {
# #     stop("magcat_acat_run_prepared(): 'prep' does not look like output of magcat_acat_prepare().",
# #          call. = FALSE)
# #   }
# #   if (!requireNamespace("ACAT", quietly = TRUE)) {
# #     stop("magcat_acat_run_prepared(): requires ACAT package.", call. = FALSE)
# #   }

# #   idx_list <- prep$idx_list
# #   p_all_u  <- as.numeric(prep$p_all_u)

# #   # optional: original-case gene IDs (if you stored them)
# #   genes_all_u <- prep$genes_all_u
# #   if (is.null(genes_all_u)) {
# #     # fall back: use normalized ids
# #     genes_all_u <- prep$genes_norm_u
# #   }

# #   n_pw <- length(idx_list)

# #   res <- data.frame(
# #     pathway_id   = prep$pathway_id,
# #     pathway_name = prep$pathway_name,
# #     n_genes      = as.integer(prep$n_genes),
# #     gene_names   = NA_character_,
# #     acat_p       = NA_real_,
# #     stringsAsFactors = FALSE
# #   )

# #   for (i in seq_len(n_pw)) {
# #     idx <- idx_list[[i]]
# #     d   <- length(idx)
# #     if (d < 1L) next

# #     p_i <- p_all_u[idx]
# #     # (p_all_u is already fixed in prepare if do_fix=TRUE, but safe anyway)
# #     p_i <- fix_p_for_acat(p_i, min_p = prep$min_p)

# #     res$acat_p[i] <- ACAT::ACAT(Pvals = p_i)
# #     res$gene_names[i] <- paste(unique(genes_all_u[idx]), collapse = ";")
# #   }

# #   ord <- order(res$acat_p, decreasing = FALSE, na.last = TRUE)
# #   res <- res[ord, , drop = FALSE]

# #   if (output) {
# #     if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# #     utils::write.csv(res, file.path(out_dir, "magcat_acat_prepared.csv"), row.names = FALSE)
# #   }

# #   res
# # }


# # #' Prepare ACAT pathway test objects (fast reusable setup)
# # #' @export
# # magcat_acat_prepare <- function(gene_results,
# #                                 pathways     = NULL,
# #                                 species      = NULL,
# #                                 pmn_gene_col = NULL,
# #                                 gene_col     = "GENE",
# #                                 p_col        = "P",
# #                                 min_p        = 1e-15,
# #                                 do_fix       = TRUE) {

# #   # pathway source
# #   if (!is.null(pathways) && !is.null(species)) {
# #     stop("Provide either 'pathways' OR 'species', not both.", call. = FALSE)
# #   }
# #   if (is.null(pathways) && is.null(species)) {
# #     stop("You must provide either 'pathways' OR 'species'.", call. = FALSE)
# #   }
# #   if (is.null(pathways) && !is.null(species)) {
# #     if (!exists("magcat_load_pathways", mode = "function")) {
# #       stop("magcat_acat_prepare(): missing magcat_load_pathways().", call. = FALSE)
# #     }
# #     pathways <- magcat_load_pathways(
# #       species  = species,
# #       gene_col = pmn_gene_col
# #     )
# #   }

# #   # gene_results checks
# #   if (!all(c(gene_col, p_col) %in% names(gene_results))) {
# #     stop("gene_results must contain columns '", gene_col, "' and '", p_col, "'.", call. = FALSE)
# #   }
# #   if (!exists("fix_p_for_acat", mode = "function")) {
# #     stop("magcat_acat_prepare(): missing fix_p_for_acat().", call. = FALSE)
# #   }

# #   genes_all  <- as.character(gene_results[[gene_col]])
# #   p_all      <- as.numeric(gene_results[[p_col]])
# #   genes_norm <- tolower(genes_all)

# #   ok <- which(!is.na(genes_norm) & genes_norm != "" &
# #               !is.na(p_all) & is.finite(p_all) & p_all > 0 & p_all <= 1)

# #   genes_all  <- genes_all[ok]
# #   genes_norm <- genes_norm[ok]
# #   p_all      <- p_all[ok]

# #   # dedupe by normalized gene id, keep first (matches your earlier logic)
# #   keep_first <- !duplicated(genes_norm)
# #   genes_all_u  <- genes_all[keep_first]
# #   genes_norm_u <- genes_norm[keep_first]
# #   p_all_u      <- p_all[keep_first]

# #   if (do_fix) {
# #     p_all_u <- fix_p_for_acat(p_all_u, min_p = min_p)
# #   }

# #   # pathways -> list
# #   if (is.data.frame(pathways)) {
# #     if (!all(c("pathway_id","gene_id") %in% names(pathways))) {
# #       stop("If 'pathways' is a data.frame it must have columns pathway_id and gene_id.", call. = FALSE)
# #     }
# #     p_list <- split(pathways$gene_id, pathways$pathway_id)
# #     p_names <- if ("pathway_name" %in% names(pathways)) {
# #       tapply(pathways$pathway_name, pathways$pathway_id, function(x) x[1])
# #     } else {
# #       setNames(names(p_list), names(p_list))
# #     }
# #   } else if (is.list(pathways)) {
# #     p_list <- pathways
# #     if (is.null(names(p_list))) names(p_list) <- paste0("PWY_", seq_along(p_list))
# #     p_names <- setNames(names(p_list), names(p_list))
# #   } else {
# #     stop("'pathways' must be list or data.frame.", call. = FALSE)
# #   }

# #   # normalize pathway genes + precompute indices into genes_norm_u
# #   p_list_norm <- lapply(p_list, function(g) tolower(as.character(g)))
# #   idx_list <- lapply(p_list_norm, function(g) {
# #     idx <- match(g, genes_norm_u)
# #     idx <- idx[!is.na(idx)]
# #     unique(idx)
# #   })

# #   n_genes <- vapply(idx_list, length, integer(1))

# #   structure(
# #     list(
# #       genes_all_u  = genes_all_u,
# #       genes_norm_u = genes_norm_u,
# #       p_all_u      = p_all_u,
# #       idx_list     = idx_list,
# #       pathway_id   = names(idx_list),
# #       pathway_name = unname(p_names[names(idx_list)]),
# #       n_genes      = n_genes,
# #       min_p        = min_p,
# #       do_fix       = do_fix
# #     ),
# #     class = "magcat_acat_prep"
# #   )
# # }
