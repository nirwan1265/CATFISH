#' Pathway-level (weighted) Stouffer Z test on MAGMA gene Z statistics
#'
#' For each pathway:
#'   * take gene-level Z statistics (e.g. MAGMA `.genes.out` has a Z column),
#'   * optionally apply weights (e.g. sqrt(NSNPS), 1/SE, etc.),
#'   * compute the weighted Stouffer Z:
#'       Z_S = sum_i w_i Z_i / sqrt(sum_i w_i^2),
#'     and return a two-sided p-value p = 2 * Phi(-|Z_S|).
#'
#' Optional permutation calibration:
#'   - perm_mode="resample_global": resample Z's from the global pool (exchangeable null)
#'   - perm_mode="mvn": simulate correlated null Z's using MAGMA gene correlation
#'
#' IMPORTANT:
#'   - If your Z's have direction, this keeps direction (unlike p-value Stouffer-from-p).
#'   - Use `B_perm > 0` if you want empirical p-values that account for correlation.
#'
#' @param gene_results data.frame with at least gene + Z columns.
#'   Typically MAGMA `.genes.out` (or your merged file) with `gene_col` and `z_col`.
#' @param pathways either:
#'   * a named list: each element is a character vector of gene IDs in that pathway, OR
#'   * a data.frame with columns `pathway_id`, `gene_id` (and optional `pathway_name`).
#' @param species optional; one of "maize", "sorghum", "arabidopsis", "plant".
#'   If provided, built-in PMN pathways are loaded via `magcat_load_pathways()`.
#'   Provide either `pathways` OR `species`, but not both.
#' @param pmn_gene_col optional; passed to `magcat_load_pathways(gene_col=...)` when `species` is used.
#' @param gene_col column name in `gene_results` containing gene IDs (default "GENE").
#' @param z_col column name in `gene_results` containing gene-level Z statistics (default "ZSTAT").
#' @param weight_col optional column name in `gene_results` for per-gene weights (default NULL = equal weights).
#' @param min_abs_w replace non-finite / <=0 weights by this small positive value (default 1e-8).
#' @param B_perm integer; number of permutations for empirical p (default 0L = no permutations).
#' @param perm_mode either "resample_global" or "mvn" (default "resample_global").
#' @param magma_genes_out required if perm_mode="mvn"; path to merged MAGMA *.genes.out (must include CHR).
#' @param magma_cor_file optional if perm_mode="mvn"; 3-col file gene1 gene2 r of gene-gene correlations.
#' @param make_PD logical; if TRUE, force correlation matrix PSD (default TRUE).
#' @param alternative either "greater", "two.sided" or "less", greater, less = one-sided in the positive and negative direction respectively (default "greater")
#' @param seed optional RNG seed (default NULL).
#' @param output logical; if TRUE, write a CSV of results to `out_dir`.
#' @param out_dir directory to write CSV when `output = TRUE` (default "magcat_stouffer_z").
#'
#' @return data.frame with columns:
#'   * pathway_id
#'   * pathway_name
#'   * n_genes
#'   * gene_names
#'   * gene_zvals           (semicolon-separated Z's used)
#'   * stouffer_Z
#'   * stouffer_p_analytic  (two-sided analytic p)
#'   * stouffer_p_perm      (perm p; NA if B_perm = 0)
#' Sorted by `stouffer_p_perm` ascending if B_perm>0 else by `stouffer_p_analytic`.
#' If `output = TRUE`, attribute `"file"` contains the CSV path.
#'
#' @export
magcat_stoufferZ_pathways <- function(gene_results,
                                      pathways        = NULL,
                                      species         = NULL,
                                      pmn_gene_col    = NULL,
                                      gene_col        = "GENE",
                                      z_col           = "ZSTAT",
                                      weight_col      = NULL,
                                      min_abs_w       = 1e-8,
                                      B_perm          = 0L,
                                      perm_mode       = c("resample_global", "mvn"),
                                      magma_genes_out = NULL,
                                      magma_cor_file  = NULL,
                                      make_PD         = TRUE,
                                      seed            = NULL,
                                      alternative     = c("greater", "two.sided", "less"), # greater, less = one-sided in the positive and negative direction
                                      output          = FALSE,
                                      out_dir         = "magcat_stouffer_z") {

  perm_mode   <- match.arg(perm_mode)
  alternative <- match.arg(alternative)

  ## -------- decide pathway source: pathways vs species ----------
  if (!is.null(pathways) && !is.null(species)) {
    stop("Provide either 'pathways' OR 'species', not both.", call. = FALSE)
  }
  if (is.null(pathways) && is.null(species)) {
    stop("You must provide either 'pathways' (list/data.frame) OR 'species'.", call. = FALSE)
  }
  if (is.null(pathways) && !is.null(species)) {
    if (is.null(pmn_gene_col)) {
      pathways <- magcat_load_pathways(species = species)
    } else {
      pathways <- magcat_load_pathways(species = species, gene_col = pmn_gene_col)
    }
  }

  ## -------- standardize gene_results ----------
  needed_cols <- c(gene_col, z_col)
  if (!all(needed_cols %in% names(gene_results))) {
    stop("gene_results must contain columns '", gene_col, "' and '", z_col, "'.", call. = FALSE)
  }
  if (!is.null(weight_col) && !(weight_col %in% names(gene_results))) {
    stop("weight_col='", weight_col, "' not found in gene_results.", call. = FALSE)
  }

  gr         <- gene_results
  genes_all  <- as.character(gr[[gene_col]])
  z_all      <- as.numeric(gr[[z_col]])
  genes_norm <- tolower(genes_all)

  # map normalized -> canonical
  gene_map <- tapply(genes_all, genes_norm, function(x) x[1])

  # named vectors for fast lookup
  z_vec <- stats::setNames(z_all, genes_norm)

  w_vec <- NULL
  if (!is.null(weight_col)) {
    w_all <- as.numeric(gr[[weight_col]])
    w_vec <- stats::setNames(w_all, genes_norm)
  }

  ## -------- pathways -> named list + names ----------
  if (is.data.frame(pathways)) {
    if (!("pathway_id" %in% names(pathways)) || !("gene_id" %in% names(pathways))) {
      stop("If 'pathways' is a data.frame it must have columns 'pathway_id' and 'gene_id'.", call. = FALSE)
    }
    if (!("pathway_name" %in% names(pathways))) pathways$pathway_name <- pathways$pathway_id
    p_list  <- split(pathways$gene_id, pathways$pathway_id)
    p_names <- tapply(pathways$pathway_name, pathways$pathway_id, FUN = function(x) x[1])
  } else if (is.list(pathways)) {
    p_list  <- pathways
    p_names <- names(pathways)
    if (is.null(p_names)) {
      p_names <- paste0("PWY_", seq_along(pathways))
      names(p_list) <- p_names
    }
  } else {
    stop("'pathways' must be either a list or a data.frame.", call. = FALSE)
  }
  p_list <- lapply(p_list, function(g) tolower(as.character(g)))

  ## -------- helpers ----------
  .fix_w <- function(w) {
    w <- as.numeric(w)
    bad <- !is.finite(w) | is.na(w) | w <= 0
    if (any(bad)) w[bad] <- min_abs_w
    w
  }

  .p_from_Z <- function(Zs) {
    if (!is.finite(Zs) || is.na(Zs)) return(NA_real_)
    switch(alternative,
      "two.sided" = 2 * stats::pnorm(-abs(Zs)),
      "greater"   = stats::pnorm(Zs, lower.tail = FALSE),
      "less"      = stats::pnorm(Zs, lower.tail = TRUE)
    )
  }

  .stouffer_from_z <- function(z_i, w_i = NULL) {
    z_i <- as.numeric(z_i)
    z_i <- z_i[is.finite(z_i) & !is.na(z_i)]
    if (!length(z_i)) return(list(Z = NA_real_, p = NA_real_))

    if (is.null(w_i)) {
      w_i <- rep(1, length(z_i))
    } else {
      w_i <- .fix_w(w_i)
      if (length(w_i) != length(z_i)) w_i <- rep(1, length(z_i))
    }

    Zs <- sum(w_i * z_i) / sqrt(sum(w_i^2))
    ps <- .p_from_Z(Zs)
    list(Z = Zs, p = ps)
  }

  ## -------- global pools for permutations ----------
  ok_pool <- which(is.finite(z_all) & !is.na(z_all) & !is.na(genes_norm) & genes_norm != "")
  if (!length(ok_pool)) stop("No valid gene-level Z statistics found.", call. = FALSE)

  z_pool <- as.numeric(z_all[ok_pool])

  if (!is.null(seed)) set.seed(seed)

  ## -------- MVN prep (once) ----------
  genes_out_tab <- NULL
  cor_pairs     <- NULL
  if (B_perm > 0L && perm_mode == "mvn") {
    if (is.null(magma_genes_out) || !file.exists(magma_genes_out)) {
      stop("perm_mode='mvn' requires magma_genes_out = path to merged MAGMA *.genes.out", call. = FALSE)
    }
    genes_out_tab <- magma_read_genes_out(magma_genes_out, gene_col = gene_col, chr_col = "CHR")

    if (!is.null(magma_cor_file)) {
      if (!file.exists(magma_cor_file)) stop("magma_cor_file does not exist: ", magma_cor_file, call. = FALSE)
      cor_pairs <- magma_read_gene_cor_pairs(magma_cor_file, gene1_col = 1, gene2_col = 2, r_col = 3)
    }
  }

  ## -------- result container ----------
  n_pw <- length(p_list)
  res  <- data.frame(
    pathway_id          = names(p_list),
    pathway_name        = if (is.null(p_names)) names(p_list) else unname(p_names),
    n_genes             = NA_integer_,
    gene_names          = NA_character_,
    gene_zvals          = NA_character_,
    stouffer_Z          = NA_real_,
    stouffer_p_analytic = NA_real_,
    stouffer_p_perm     = NA_real_,
    stringsAsFactors    = FALSE
  )

  ## -------- main loop ----------
  for (i in seq_len(n_pw)) {
    genes_i <- p_list[[i]]

    z_i <- z_vec[genes_i]
    keep <- is.finite(z_i) & !is.na(z_i)
    z_i  <- as.numeric(z_i[keep])
    genes_used <- genes_i[keep]

    d <- length(z_i)
    res$n_genes[i] <- d
    if (!d) next

    canon_ids <- unname(gene_map[genes_used])
    canon_ids <- unique(canon_ids[!is.na(canon_ids)])
    res$gene_names[i] <- paste(canon_ids, collapse = ";")
    res$gene_zvals[i] <- paste(z_i, collapse = ";")

    w_i <- NULL
    if (!is.null(w_vec)) w_i <- as.numeric(w_vec[genes_used])

    st <- .stouffer_from_z(z_i, w_i)
    res$stouffer_Z[i]          <- st$Z
    res$stouffer_p_analytic[i] <- st$p

    ## ---- permutations ----
    if (B_perm > 0L && is.finite(st$Z)) {

      if (alternative == "two.sided") {
        T_obs <- abs(st$Z)
        T_of  <- function(Z) abs(Z)
        cmp   <- function(Tnull) Tnull >= T_obs
      } else if (alternative == "greater") {
        T_obs <- st$Z
        T_of  <- function(Z) Z
        cmp   <- function(Tnull) Tnull >= T_obs
      } else { # "less"
        T_obs <- st$Z
        T_of  <- function(Z) Z
        cmp   <- function(Tnull) Tnull <= T_obs
      }

      T_null <- rep(NA_real_, B_perm)

      if (perm_mode == "resample_global") {
        for (b in seq_len(B_perm)) {
          z_perm <- sample(z_pool, size = d, replace = FALSE)
          st_b <- .stouffer_from_z(z_perm, if (is.null(w_i)) NULL else w_i)
          T_null[b] <- T_of(st_b$Z)
        }
        res$stouffer_p_perm[i] <- (1 + sum(cmp(T_null), na.rm = TRUE)) / (B_perm + 1)

      } else if (perm_mode == "mvn") {

        R_S <- magma_build_R_for_pathway(
          genes_S       = genes_used,
          genes_out_tab = genes_out_tab,
          cor_pairs     = cor_pairs,
          gene_col      = gene_col,
          chr_col       = "CHR",
          make_PD       = make_PD
        )
        if (is.null(R_S)) next

        sim <- magma_simulate_null_Zp(
          genes_S       = rownames(R_S),
          genes_out_tab = genes_out_tab,
          R_S           = R_S,
          B             = B_perm,
          gene_col      = gene_col,
          chr_col       = "CHR"
        )

        Zmat <- sim$Z
        if (is.null(Zmat)) {
          # If only P is returned, convert P -> Z using MAGMA-style association-strength:
          #   z = qnorm(1 - p)
          Pmat <- sim$P
          if (is.null(Pmat)) next

          for (b in seq_len(B_perm)) {
            p_perm <- as.numeric(Pmat[b, ])
            p_perm <- p_perm[is.finite(p_perm) & !is.na(p_perm)]
            if (!length(p_perm)) next
            p_perm <- pmin(pmax(p_perm, 1e-15), 1 - 1e-15)

            z_perm <- stats::qnorm(1 - p_perm)

            st_b <- .stouffer_from_z(z_perm, if (is.null(w_i)) NULL else w_i)
            T_null[b] <- T_of(st_b$Z)
          }

        } else {
          for (b in seq_len(B_perm)) {
            z_perm <- as.numeric(Zmat[b, ])
            z_perm <- z_perm[is.finite(z_perm) & !is.na(z_perm)]
            if (!length(z_perm)) next
            st_b <- .stouffer_from_z(z_perm, if (is.null(w_i)) NULL else w_i)
            T_null[b] <- T_of(st_b$Z)
          }
        }

        res$stouffer_p_perm[i] <- (1 + sum(cmp(T_null), na.rm = TRUE)) / (B_perm + 1)
      }
    }
  }

  ## -------- sort ----------
  if (B_perm > 0L) {
    ord <- order(res$stouffer_p_perm, decreasing = FALSE, na.last = TRUE)
  } else {
    ord <- order(res$stouffer_p_analytic, decreasing = FALSE, na.last = TRUE)
  }
  res <- res[ord, , drop = FALSE]

  ## -------- optional CSV output ----------
  if (output) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    species_tag <- if (is.null(species)) "custom" else species
    out_path <- file.path(out_dir, paste0("magcat_stoufferZ_pathways_", species_tag, ".csv"))
    utils::write.csv(res, out_path, row.names = FALSE)
    attr(res, "file") <- out_path
  }

  res
}











# magcat_stoufferZ_pathways <- function(gene_results,
#                                       pathways        = NULL,
#                                       species         = NULL,
#                                       pmn_gene_col    = NULL,
#                                       gene_col        = "GENE",
#                                       z_col           = "ZSTAT",
#                                       weight_col      = NULL,
#                                       min_abs_w       = 1e-8,
#                                       B_perm          = 0L,
#                                       perm_mode       = c("resample_global", "mvn"),
#                                       magma_genes_out = NULL,
#                                       magma_cor_file  = NULL,
#                                       make_PD         = TRUE,
#                                       seed            = NULL,
#                                       output          = FALSE,
#                                       out_dir         = "magcat_stouffer_z") {

#   perm_mode <- match.arg(perm_mode)

#   ## -------- decide pathway source: pathways vs species ----------
#   if (!is.null(pathways) && !is.null(species)) {
#     stop("Provide either 'pathways' OR 'species', not both.", call. = FALSE)
#   }
#   if (is.null(pathways) && is.null(species)) {
#     stop("You must provide either 'pathways' (list/data.frame) OR 'species'.", call. = FALSE)
#   }
#   if (is.null(pathways) && !is.null(species)) {
#     if (is.null(pmn_gene_col)) {
#       pathways <- magcat_load_pathways(species = species)
#     } else {
#       pathways <- magcat_load_pathways(species = species, gene_col = pmn_gene_col)
#     }
#   }

#   ## -------- standardize gene_results ----------
#   needed_cols <- c(gene_col, z_col)
#   if (!all(needed_cols %in% names(gene_results))) {
#     stop("gene_results must contain columns '", gene_col, "' and '", z_col, "'.", call. = FALSE)
#   }
#   if (!is.null(weight_col) && !(weight_col %in% names(gene_results))) {
#     stop("weight_col='", weight_col, "' not found in gene_results.", call. = FALSE)
#   }

#   gr         <- gene_results
#   genes_all  <- as.character(gr[[gene_col]])
#   z_all      <- as.numeric(gr[[z_col]])
#   genes_norm <- tolower(genes_all)

#   # map normalized -> canonical
#   gene_map <- tapply(genes_all, genes_norm, function(x) x[1])

#   # named vectors for fast lookup
#   z_vec <- stats::setNames(z_all, genes_norm)

#   w_vec <- NULL
#   if (!is.null(weight_col)) {
#     w_all <- as.numeric(gr[[weight_col]])
#     w_vec <- stats::setNames(w_all, genes_norm)
#   }

#   ## -------- pathways -> named list + names ----------
#   if (is.data.frame(pathways)) {
#     if (!("pathway_id" %in% names(pathways)) || !("gene_id" %in% names(pathways))) {
#       stop("If 'pathways' is a data.frame it must have columns 'pathway_id' and 'gene_id'.", call. = FALSE)
#     }
#     if (!("pathway_name" %in% names(pathways))) pathways$pathway_name <- pathways$pathway_id

#     p_list  <- split(pathways$gene_id, pathways$pathway_id)
#     p_names <- tapply(pathways$pathway_name, pathways$pathway_id, FUN = function(x) x[1])
#   } else if (is.list(pathways)) {
#     p_list  <- pathways
#     p_names <- names(pathways)
#     if (is.null(p_names)) {
#       p_names <- paste0("PWY_", seq_along(pathways))
#       names(p_list) <- p_names
#     }
#   } else {
#     stop("'pathways' must be either a list or a data.frame.", call. = FALSE)
#   }
#   p_list <- lapply(p_list, function(g) tolower(as.character(g)))

#   ## -------- helpers ----------
#   .fix_w <- function(w) {
#     w <- as.numeric(w)
#     bad <- !is.finite(w) | is.na(w) | w <= 0
#     if (any(bad)) w[bad] <- min_abs_w
#     w
#   }

#   .stouffer_from_z <- function(z_i, w_i = NULL) {
#     z_i <- as.numeric(z_i)
#     z_i <- z_i[is.finite(z_i) & !is.na(z_i)]
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

#   ## -------- global pools for permutations ----------
#   ok_pool <- which(is.finite(z_all) & !is.na(z_all) & !is.na(genes_norm) & genes_norm != "")
#   if (!length(ok_pool)) stop("No valid gene-level Z statistics found.", call. = FALSE)

#   z_pool <- as.numeric(z_all[ok_pool])

#   if (!is.null(seed)) set.seed(seed)

#   ## -------- MVN prep (once) ----------
#   genes_out_tab <- NULL
#   cor_pairs     <- NULL
#   if (B_perm > 0L && perm_mode == "mvn") {
#     if (is.null(magma_genes_out) || !file.exists(magma_genes_out)) {
#       stop("perm_mode='mvn' requires magma_genes_out = path to merged MAGMA *.genes.out", call. = FALSE)
#     }
#     genes_out_tab <- magma_read_genes_out(magma_genes_out, gene_col = gene_col, chr_col = "CHR")

#     if (!is.null(magma_cor_file)) {
#       if (!file.exists(magma_cor_file)) stop("magma_cor_file does not exist: ", magma_cor_file, call. = FALSE)
#       cor_pairs <- magma_read_gene_cor_pairs(magma_cor_file, gene1_col = 1, gene2_col = 2, r_col = 3)
#     }
#   }

#   ## -------- result container ----------
#   n_pw <- length(p_list)
#   res  <- data.frame(
#     pathway_id         = names(p_list),
#     pathway_name       = if (is.null(p_names)) names(p_list) else unname(p_names),
#     n_genes            = NA_integer_,
#     gene_names         = NA_character_,
#     gene_zvals         = NA_character_,
#     stouffer_Z         = NA_real_,
#     stouffer_p_analytic= NA_real_,
#     stouffer_p_perm    = NA_real_,
#     stringsAsFactors   = FALSE
#   )

#   ## -------- main loop ----------
#   for (i in seq_len(n_pw)) {
#     genes_i <- p_list[[i]]

#     z_i <- z_vec[genes_i]
#     keep <- is.finite(z_i) & !is.na(z_i)
#     z_i  <- as.numeric(z_i[keep])
#     genes_used <- genes_i[keep]

#     d <- length(z_i)
#     res$n_genes[i] <- d
#     if (!d) next

#     canon_ids <- unname(gene_map[genes_used])
#     canon_ids <- unique(canon_ids[!is.na(canon_ids)])
#     res$gene_names[i] <- paste(canon_ids, collapse = ";")
#     res$gene_zvals[i] <- paste(z_i, collapse = ";")

#     w_i <- NULL
#     if (!is.null(w_vec)) {
#       w_i <- as.numeric(w_vec[genes_used])
#     }

#     st <- .stouffer_from_z(z_i, w_i)
#     res$stouffer_Z[i]          <- st$Z
#     res$stouffer_p_analytic[i] <- st$p

#     ## ---- permutations (two-sided; compare |Z|) ----
#     if (B_perm > 0L && is.finite(st$Z)) {

#       Z_obs <- abs(st$Z)
#       Z_null <- rep(NA_real_, B_perm)

#       if (perm_mode == "resample_global") {

#         for (b in seq_len(B_perm)) {
#           z_perm <- sample(z_pool, size = d, replace = FALSE)
#           st_b <- .stouffer_from_z(z_perm, if (is.null(w_i)) NULL else w_i)
#           Z_null[b] <- abs(st_b$Z)
#         }

#         res$stouffer_p_perm[i] <- (1 + sum(Z_null >= Z_obs, na.rm = TRUE)) / (B_perm + 1)

#       } else if (perm_mode == "mvn") {

#         R_S <- magma_build_R_for_pathway(
#           genes_S       = genes_used,
#           genes_out_tab = genes_out_tab,
#           cor_pairs     = cor_pairs,
#           gene_col      = gene_col,
#           chr_col       = "CHR",
#           make_PD       = make_PD
#         )
#         if (is.null(R_S)) next

#         # simulate correlated Z under null
#         sim <- magma_simulate_null_Zp(
#           genes_S       = rownames(R_S),
#           genes_out_tab = genes_out_tab,
#           R_S           = R_S,
#           B             = B_perm,
#           gene_col      = gene_col,
#           chr_col       = "CHR"
#         )

#         # IMPORTANT: magma_simulate_null_Zp() in your code returns both Z and P.
#         # We prefer Z here if available; otherwise invert P to |Z|.
#         Zmat <- sim$Z
#         if (is.null(Zmat)) {
#           # fallback from P to |Z| (loses sign)
#           Pmat <- sim$P
#           if (is.null(Pmat)) next
#           for (b in seq_len(B_perm)) {
#             p_perm <- as.numeric(Pmat[b, ])
#             p_perm <- p_perm[is.finite(p_perm) & !is.na(p_perm)]
#             if (!length(p_perm)) next
#             z_abs <- stats::qnorm(p_perm / 2, lower.tail = FALSE)
#             st_b <- .stouffer_from_z(z_abs, if (is.null(w_i)) NULL else w_i)
#             Z_null[b] <- abs(st_b$Z)
#           }
#         } else {
#           for (b in seq_len(B_perm)) {
#             z_perm <- as.numeric(Zmat[b, ])
#             z_perm <- z_perm[is.finite(z_perm) & !is.na(z_perm)]
#             if (!length(z_perm)) next
#             st_b <- .stouffer_from_z(z_perm, if (is.null(w_i)) NULL else w_i)
#             Z_null[b] <- abs(st_b$Z)
#           }
#         }

#         res$stouffer_p_perm[i] <- (1 + sum(Z_null >= Z_obs, na.rm = TRUE)) / (B_perm + 1)
#       }
#     }
#   }

#   ## -------- sort ----------
#   if (B_perm > 0L) {
#     ord <- order(res$stouffer_p_perm, decreasing = FALSE, na.last = TRUE)
#   } else {
#     ord <- order(res$stouffer_p_analytic, decreasing = FALSE, na.last = TRUE)
#   }
#   res <- res[ord, , drop = FALSE]

#   ## -------- optional CSV output ----------
#   if (output) {
#     if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#     species_tag <- if (is.null(species)) "custom" else species
#     out_path <- file.path(out_dir, paste0("magcat_stoufferZ_pathways_", species_tag, ".csv"))
#     utils::write.csv(res, out_path, row.names = FALSE)
#     attr(res, "file") <- out_path
#   }

#   res
# } 
