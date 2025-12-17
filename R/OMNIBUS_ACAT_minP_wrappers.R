#' Omnibus pathway test (ACAT-O or minP) using existing MAGCAT wrappers
#'
#' This function:
#'   1) Calls your existing wrappers:
#'        - magcat_acat_pathways()
#'        - magcat_wfisher_pathways()
#'        - magcat_tfisher_pathways()
#'        - magcat_stouffer_pathways()
#'        - magcat_minp_pathways()
#'   2) Extracts analytic p-values from each.
#'   3) Optionally drops pathways with n_genes < 2.
#'   4) Combines method p's by:
#'        - omnibus = "ACAT" -> ACAT-O (sumFREGAT::ACATO)
#'        - omnibus = "minP" -> minP (metap::minimump or Sidák)
#'   5) Optionally does omni-level gene-set permutations to give omni_perm_p.
#'
#' @export
omni_pathways_2 <- function(gene_results,
                                 pathways     = NULL,
                                 species      = NULL,
                                 pmn_gene_col = NULL,
                                 gene_col     = "GENE",
                                 p_col        = "P",
                                 effect_col   = "ZSTAT",
                                 weight_col   = NULL,
                                 is_onetail   = FALSE,
                                 ptrunc       = 0.05,
                                 min_p        = 1e-15,
                                 do_fix       = TRUE,
                                 B_perm       = 0L,        # <-- omni permutations
                                 perm_mode       = c("resample_global", "mvn"),      # NEW
                                 magma_genes_out  = NULL,                  # NEW: path to merged .genes.out
                                 magma_cor_file   = NULL,                  # NEW: optional file with gene-pair correlations
                                 make_PD          = TRUE,                   # NEW
                                 magma_raw_dir = NULL,
                                 seed         = NULL,
                                 omnibus      = c("ACAT", "minP"),
                                 remove_singletons = TRUE,
                                 output       = FALSE,
                                 out_dir      = "magcat_omni") {

  omnibus <- match.arg(omnibus)
  
  perm_mode <- match.arg(perm_mode)

  ## ---------- helpers -------------------------------------------------
  .p_adjust_BH <- function(p) {
    p <- as.numeric(p)
    out <- rep(NA_real_, length(p))
    idx <- which(!is.na(p))
    if (length(idx) > 0L) {
      out[idx] <- stats::p.adjust(p[idx], method = "BH")
    }
    out
  }

  .qvalue_vec <- function(p) {
    p <- as.numeric(p)
    out <- rep(NA_real_, length(p))
    idx <- which(!is.na(p) & p > 0 & p <= 1)
    if (length(idx) > 1L && requireNamespace("qvalue", quietly = TRUE)) {
      qq <- tryCatch(
        qvalue::qvalue(p[idx])$qvalues,
        error = function(e) rep(NA_real_, length(idx))
      )
      out[idx] <- qq
    }
    out
  }

  .combine_minP <- function(p_vec) {
    p_vec <- as.numeric(p_vec)
    p_vec <- p_vec[!is.na(p_vec) & p_vec > 0 & p_vec <= 1]
    if (!length(p_vec)) return(NA_real_)
    if (requireNamespace("metap", quietly = TRUE)) {
      out <- tryCatch(
        metap::minimump(p_vec)$p,
        error = function(e) NA_real_
      )
      return(out)
    } else {
      k     <- length(p_vec)
      p_min <- min(p_vec)
      return(1 - (1 - p_min)^k)  # Sidák
    }
  }

  .combine_omni <- function(p_vec) {
    p_vec <- as.numeric(p_vec)
    p_vec <- p_vec[!is.na(p_vec) & p_vec > 0 & p_vec <= 1]
    if (!length(p_vec)) return(NA_real_)

    if (omnibus == "ACAT") {
      if (!requireNamespace("sumFREGAT", quietly = TRUE)) {
        stop("omni_pathways(): omnibus='ACAT' requires sumFREGAT::ACATO().",
             call. = FALSE)
      }
      p_vec <- fix_p_for_acat(p_vec, min_p = min_p)
      out <- tryCatch(
        as.numeric(sumFREGAT::ACATO(p_vec)),
        error = function(e) NA_real_
      )
      return(out)
    } else {  # "minP"
      return(.combine_minP(p_vec))
    }
  }

  .find_pcol <- function(tab, method) {
    nm <- names(tab)
    patt <- switch(
      method,
      "acat"     = "acat.*p|p.*acat",
      "wfisher"  = "wfisher.*p|p.*wfisher",
      "tfisher"  = "tpm.*p|p.*tpm|tfisher.*p",
      "stouffer" = "stouffer.*p",
      "minp"     = "minp.*p",
      stop("Unknown method in .find_pcol: ", method)
    )
    cand <- grep(patt, nm, value = TRUE, ignore.case = TRUE)
    cand <- cand[!grepl("perm|BH|q|stat|z", cand, ignore.case = TRUE)]
    if (!length(cand)) return(NULL)
    cand[1]
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  ## ---------- 1) CALL YOUR WRAPPERS (analytic p's) --------------------

  acat_tab <- magcat_acat_pathways(
    gene_results = gene_results,
    pathways     = pathways,
    species      = species,
    pmn_gene_col = pmn_gene_col,
    gene_col     = gene_col,
    p_col        = p_col,
    min_p        = min_p,
    do_fix       = do_fix,
    B            = 0L,      # no per-method perms here
    seed         = seed,
    output       = FALSE
  )

  wf_tab <- magcat_wfisher_pathways(
    gene_results = gene_results,
    pathways     = pathways,
    species      = species,
    pmn_gene_col = pmn_gene_col,
    gene_col     = gene_col,
    p_col        = p_col,
    effect_col   = effect_col,
    weight_col   = weight_col,
    is_onetail   = is_onetail,
    min_p        = min_p,
    do_fix       = do_fix,
    output       = FALSE
  )

  tf_tab <- magcat_tfisher_pathways(
    gene_results     = gene_results,
    pathways         = pathways,
    species          = species,
    pmn_gene_col     = pmn_gene_col,
    gene_col         = gene_col,
    p_col            = p_col,
    ptrunc           = ptrunc,
    min_p            = min_p,
    do_fix           = do_fix,
    B_perm           = 0L,
    seed             = seed,
    analytic_logical = TRUE,
    output           = FALSE
  )

  st_tab <- magcat_stouffer_pathways(
    gene_results     = gene_results,
    pathways         = pathways,
    species          = species,
    pmn_gene_col     = pmn_gene_col,
    gene_col         = gene_col,
    p_col            = p_col,
    weight_col       = weight_col,
    min_p            = min_p,
    do_fix           = do_fix,
    B_perm           = 0L,
    seed             = seed,
    analytic_logical = TRUE,
    output           = FALSE
  )

  minp_tab <- magcat_minp_pathways(
    gene_results     = gene_results,
    pathways         = pathways,
    species          = species,
    pmn_gene_col     = pmn_gene_col,
    gene_col         = gene_col,
    p_col            = p_col,
    min_p            = min_p,
    do_fix           = do_fix,
    B_perm           = 0L,
    seed             = seed,
    analytic_logical = TRUE,
    output           = FALSE
  )

  ## ---------- 2) BUILD BASE RESULT TABLE -----------------------------

  if (!all(c("pathway_id", "pathway_name", "n_genes") %in% names(acat_tab))) {
    stop("omni_pathways(): acat table must have pathway_id, pathway_name, n_genes.",
         call. = FALSE)
  }

  res <- acat_tab[, c("pathway_id", "pathway_name", "n_genes"), drop = FALSE]

  pcol_acat     <- .find_pcol(acat_tab,    "acat")
  pcol_wfisher  <- .find_pcol(wf_tab,      "wfisher")
  pcol_tfisher  <- .find_pcol(tf_tab,      "tfisher")
  pcol_stouffer <- .find_pcol(st_tab,      "stouffer")
  pcol_minp     <- .find_pcol(minp_tab,    "minp")

  if (is.null(pcol_acat))     stop("Could not find ACAT p-column in acat_tab.",     call. = FALSE)
  if (is.null(pcol_wfisher))  stop("Could not find wFisher p-column in wf_tab.",    call. = FALSE)
  if (is.null(pcol_tfisher))  stop("Could not find TFisher/TPM p-column in tf_tab.",call. = FALSE)
  if (is.null(pcol_stouffer)) stop("Could not find Stouffer p-column in st_tab.",   call. = FALSE)
  if (is.null(pcol_minp))     stop("Could not find minP p-column in minp_tab.",     call. = FALSE)

  res$acat_p      <- acat_tab[ match(res$pathway_id, acat_tab$pathway_id),  pcol_acat     ]
  res$wfisher_p   <- wf_tab[  match(res$pathway_id, wf_tab$pathway_id),     pcol_wfisher  ]
  res$tpm_p       <- tf_tab[  match(res$pathway_id, tf_tab$pathway_id),     pcol_tfisher  ]
  res$stouffer_p  <- st_tab[  match(res$pathway_id, st_tab$pathway_id),     pcol_stouffer ]
  res$minp_gene_p <- minp_tab[match(res$pathway_id, minp_tab$pathway_id),   pcol_minp     ]

  ## ---------- 3) DROP n_genes < 2 (your request) ---------------------

  if (remove_singletons) {
    keep <- !is.na(res$n_genes) & res$n_genes >= 2L
    res  <- res[keep, , drop = FALSE]
  }

  ## ---------- 4) OMNI (analytic) ------------------------------------

  omni_p <- rep(NA_real_, nrow(res))
  for (i in seq_len(nrow(res))) {
    pv <- c(
      res$acat_p[i],
      res$wfisher_p[i],
      res$tpm_p[i],
      res$stouffer_p[i],
      res$minp_gene_p[i]
    )
    omni_p[i] <- .combine_omni(pv)
  }
  res$omni_p    <- omni_p
  res$omni_p_BH <- .p_adjust_BH(res$omni_p)
  #res$omni_p_q  <- .qvalue_vec(res$omni_p)

  ## ---------- 5) OMNI-LEVEL PERMUTATIONS -----------------------------

  res$omni_perm_p    <- NA_real_
  res$omni_perm_p_BH <- NA_real_
  #res$omni_perm_p_q  <- NA_real_

  if (B_perm > 0L) {

    if (perm_mode == "resample_global") {

      gr <- gene_results
      genes_all <- tolower(as.character(gr[[gene_col]]))
      p_all     <- as.numeric(gr[[p_col]])

      ok <- which(!is.na(genes_all) & genes_all != "" &
                  !is.na(p_all) & p_all > 0 & p_all <= 1)

      genes_all <- genes_all[ok]
      p_all     <- p_all[ok]

      # pathway list (as you already build)
      if (!is.null(pathways)) {
        if (is.data.frame(pathways)) {
          if (!all(c("pathway_id","gene_id") %in% names(pathways))) {
            stop("pathways data.frame must have pathway_id and gene_id", call. = FALSE)
          }
          p_list <- split(pathways$gene_id, pathways$pathway_id)
        } else if (is.list(pathways)) {
          p_list <- pathways
        } else stop("'pathways' must be list or data.frame.", call. = FALSE)
      } else if (!is.null(species)) {
        pw <- magcat_load_pathways(
          species  = species,
          gene_col = if (is.null(pmn_gene_col)) NULL else pmn_gene_col
        )
        p_list <- split(pw$gene_id, pw$pathway_id)
      } else {
        stop("Need 'pathways' or 'species' for permutations.", call. = FALSE)
      }

      p_list <- lapply(p_list, function(g) tolower(as.character(g)))
      p_list <- p_list[res$pathway_id]

      # Map gene -> p (observed)
      p_map <- setNames(p_all, genes_all)

      omni_perm_p <- rep(NA_real_, nrow(res))

      for (i in seq_len(nrow(res))) {

        genes_S <- p_list[[i]]
        if (is.null(genes_S) || length(genes_S) < 2L) next
        if (is.na(res$omni_p[i])) next

        # keep only genes that exist in gene_results
        genes_S <- genes_S[genes_S %in% names(p_map)]
        d <- length(genes_S)
        if (d < 2L) next

        omni_null <- numeric(B_perm)

        for (b in seq_len(B_perm)) {

          # Permute gene resample_global (shuffle p-values across gene IDs)
          perm_idx <- sample.int(length(p_all))
          p_map_b  <- setNames(p_all[perm_idx], genes_all)

          p_i <- as.numeric(p_map_b[genes_S])
          p_i <- fix_p_for_acat(p_i, min_p = min_p)

          # ACAT
          p_acat <- ACAT::ACAT(Pvals = p_i)

          # TPM
          p_tpm <- magcat_tpm_stat(p_i, ptrunc = ptrunc)$pval

          # Stouffer Z from P under null (two-sided)
          z_i  <- stats::qnorm(p_i/2, lower.tail = FALSE)
          z_st <- sum(z_i) / sqrt(length(z_i))
          p_st <- 2 * stats::pnorm(-abs(z_st))

          # wFisher (one-tail under null)
          p_wf <- wFisher(p = p_i, weight = rep(1, length(p_i)), is.onetail = TRUE)$p

          # minP gene
          p_min_gene <- .combine_minP(p_i)

          omni_null[b] <- .combine_omni(c(p_acat, p_wf, p_tpm, p_st, p_min_gene))
        }

        omni_obs <- res$omni_p[i]
        omni_perm_p[i] <- (1 + sum(omni_null <= omni_obs, na.rm = TRUE)) / (B_perm + 1)
      }

      res$omni_perm_p    <- omni_perm_p
      res$omni_perm_p_BH <- .p_adjust_BH(res$omni_perm_p)
      #res$omni_perm_p_q  <- .qvalue_vec(res$omni_perm_p)
    }



    if (perm_mode == "mvn") {

      if (is.null(magma_genes_out) || !file.exists(magma_genes_out)) {
        stop("perm_mode='mvn' requires magma_genes_out = path to merged MAGMA *.genes.out",
            call. = FALSE)
      }

      genes_out_tab <- magma_read_genes_out(magma_genes_out, gene_col = gene_col, chr_col = "CHR")

      cor_pairs <- NULL
      if (!is.null(magma_cor_file)) {
        if (!file.exists(magma_cor_file)) {
          stop("magma_cor_file does not exist: ", magma_cor_file, call. = FALSE)
        }
        # assumes 3-column: gene1 gene2 r
        cor_pairs <- magma_read_gene_cor_pairs(magma_cor_file, gene1_col = 1, gene2_col = 2, r_col = 3)
      }

      # build pathway list (same as you already do)
      if (!is.null(pathways)) {
        if (is.data.frame(pathways)) {
          if (!all(c("pathway_id","gene_id") %in% names(pathways))) {
            stop("If 'pathways' is a data.frame it must have pathway_id and gene_id.", call. = FALSE)
          }
          p_list <- split(pathways$gene_id, pathways$pathway_id)
        } else if (is.list(pathways)) {
          p_list <- pathways
        } else {
          stop("'pathways' must be list or data.frame.", call. = FALSE)
        }
      } else if (!is.null(species)) {
        pw <- magcat_load_pathways(
          species  = species,
          gene_col = if (is.null(pmn_gene_col)) NULL else pmn_gene_col
        )
        p_list <- split(pw$gene_id, pw$pathway_id)
      } else {
        stop("Need 'pathways' or 'species' for perms.", call. = FALSE)
      }

      p_list <- lapply(p_list, function(g) tolower(as.character(g)))
      p_list <- p_list[res$pathway_id]  # same order

      # ---- compute observed omni_p already done above ----

      omni_perm_p <- rep(NA_real_, nrow(res))

      for (i in seq_len(nrow(res))) {
        genes_S <- p_list[[i]]
        if (is.null(genes_S) || length(genes_S) < 2L) next
        if (is.na(res$omni_p[i])) next

        # Build R_S (LD-aware within local neighborhoods; distant assumed 0)
        R_S <- magma_build_R_for_pathway(
          genes_S        = genes_S,
          genes_out_tab  = genes_out_tab,
          cor_pairs      = cor_pairs,
          gene_col       = gene_col,
          chr_col        = "CHR",
          make_PD        = make_PD
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

        # Now recompute component tests from simulated P (and Z for Stouffer)
        # We do "onetail" null (no direction) unless you *really* need signed WFisher.
        omni_null <- numeric(B_perm)

        for (b in seq_len(B_perm)) {
          p_i <- sim$P[b, ]

          # ACAT
          p_acat <- ACAT::ACAT(Pvals = fix_p_for_acat(p_i, min_p = min_p))

          # Fisher / TPM / minP gene-level
          # (use your existing helper for TPM)
          tpm   <- magcat_tpm_stat(p_i, ptrunc = ptrunc)
          p_tpm <- tpm$pval

          # Stouffer from Z directly
          z_i <- sim$Z[b, ]
          z_st <- sum(z_i) / sqrt(length(z_i))
          p_st <- 2 * stats::pnorm(-abs(z_st))

          # wFisher: under null, make it onetail (directionless) to avoid needing signs
          wf <- wFisher(p = p_i, weight = rep(1, length(p_i)), is.onetail = TRUE)
          p_wf <- wf$p

          p_min_gene <- .combine_minP(p_i)

          omni_null[b] <- .combine_omni(c(p_acat, p_wf, p_tpm, p_st, p_min_gene))
        }

        omni_obs <- res$omni_p[i]
        omni_perm_p[i] <- (1 + sum(omni_null <= omni_obs, na.rm = TRUE)) / (B_perm + 1)
      }

      res$omni_perm_p    <- omni_perm_p
      res$omni_perm_p_BH <- .p_adjust_BH(res$omni_perm_p)
      #res$omni_perm_p_q  <- .qvalue_vec(res$omni_perm_p)
    }
  }


  ## ---------- sort by analytic omni_p --------------------------------

  # after permutations block finishes (or even regardless)
  res$omni_final_p <- if (B_perm > 0L) res$omni_perm_p else res$omni_p
  res$omni_final_BH <- .p_adjust_BH(res$omni_final_p)

  # optional: drop qvalue entirely (recommended)
  #res$omni_p_q <- NULL
  #res$omni_perm_p_q <- NULL

  # sort by final (not analytic)
  ord <- order(res$omni_final_p, decreasing = FALSE, na.last = TRUE)
  res <- res[ord, , drop = FALSE]

  ## ---------- optional CSV -------------------------------------------

  if (output) {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    species_tag <- if (is.null(species)) "custom" else species
    out_tag     <- if (omnibus == "ACAT") "acato" else "minp"
    out_path <- file.path(
      out_dir,
      paste0("omni_pathways_", species_tag, "_", out_tag, ".csv")
    )
    utils::write.csv(res, out_path, row.names = FALSE)
    attr(res, "file") <- out_path
  }

  res
}











#' Omnibus pathway test (ACAT-O or minP) using existing MAGCAT wrappers
#'
#' Uses your wrapper functions for BOTH:
#'   (1) analytic method p-values
#'   (2) omni-level permutations (resample_global or mvn)
#'
#' NOTE:
#'   - Fisher is classic Fisher via magcat_fisher_pathways() (no weights/direction).
#'   - TFisher uses magcat_soft_tfisher_pathways(); ptrunc is passed as tau1.
#'
#' @export
omni_pathways55 <- function(gene_results,
                          pathways     = NULL,
                          species      = NULL,
                          pmn_gene_col = NULL,
                          gene_col     = "GENE",
                          p_col        = "P",
                          effect_col   = "ZSTAT",   # kept for backward compat (unused now)
                          weight_col   = NULL,
                          is_onetail   = FALSE,     # kept for backward compat (unused now)
                          ptrunc       = 0.05,
                          min_p        = 1e-15,
                          do_fix       = TRUE,
                          B_perm       = 0L,
                          perm_mode       = c("resample_global", "mvn"),
                          magma_genes_out  = NULL,
                          magma_cor_file   = NULL,
                          make_PD          = TRUE,
                          magma_raw_dir = NULL,
                          seed         = NULL,
                          omnibus      = c("ACAT", "minP"),
                          remove_singletons = TRUE,
                          output       = FALSE,
                          out_dir      = "magcat_omni") {

  omnibus   <- match.arg(omnibus)
  perm_mode <- match.arg(perm_mode)

  ## ---------- sanity: require your wrappers ----------
  req_funs <- c("magcat_acat_pathways",
                "magcat_fisher_pathways",
                "magcat_soft_tfisher_pathways",
                "magcat_stouffer_pathways",
                "magcat_minp_pathways")
  missing_f <- req_funs[!vapply(req_funs, exists, logical(1), mode = "function")]
  if (length(missing_f)) {
    stop("omni_pathways(): missing required MAGCAT functions: ",
         paste(missing_f, collapse = ", "),
         call. = FALSE)
  }

  ## ---------- helpers ----------
  .p_adjust_BH <- function(p) {
    p <- as.numeric(p)
    out <- rep(NA_real_, length(p))
    idx <- which(!is.na(p))
    if (length(idx) > 0L) out[idx] <- stats::p.adjust(p[idx], method = "BH")
    out
  }

  .combine_minP <- function(p_vec) {
    p_vec <- as.numeric(p_vec)
    p_vec <- p_vec[!is.na(p_vec) & p_vec > 0 & p_vec <= 1]
    if (!length(p_vec)) return(NA_real_)
    if (requireNamespace("metap", quietly = TRUE)) {
      out <- tryCatch(metap::minimump(p_vec)$p, error = function(e) NA_real_)
      return(out)
    } else {
      k <- length(p_vec)
      p_min <- min(p_vec)
      return(1 - (1 - p_min)^k)  # Sidák
    }
  }

  .combine_omni <- function(p_vec) {
    p_vec <- as.numeric(p_vec)
    p_vec <- p_vec[!is.na(p_vec) & p_vec > 0 & p_vec <= 1]
    if (!length(p_vec)) return(NA_real_)

    if (omnibus == "ACAT") {
      if (!requireNamespace("sumFREGAT", quietly = TRUE)) {
        stop("omni_pathways(): omnibus='ACAT' requires sumFREGAT::ACATO().",
             call. = FALSE)
      }
      p_vec <- fix_p_for_acat(p_vec, min_p = min_p)
      out <- tryCatch(as.numeric(sumFREGAT::ACATO(p_vec)), error = function(e) NA_real_)
      return(out)
    } else {
      return(.combine_minP(p_vec))
    }
  }

  .find_pcol <- function(tab, method) {
    nm <- names(tab)
    patt <- switch(
      method,
      "acat"     = "acat.*p|p.*acat",
      "fisher"   = "fisher.*p|p.*fisher",
      "tfisher"  = "tfisher.*p|p.*tfisher",
      "stouffer" = "stouffer.*p",
      "minp"     = "minp.*p",
      stop("Unknown method in .find_pcol: ", method)
    )
    cand <- grep(patt, nm, value = TRUE, ignore.case = TRUE)
    cand <- cand[!grepl("perm|BH|q|stat|z", cand, ignore.case = TRUE)]
    if (!length(cand)) return(NULL)
    cand[1]
  }

  ## Use your wrappers to compute component p’s for ONE pathway (id = "TMP")
  .component_pvals_onepath <- function(gene_results_one, genes_S_norm) {
    pw1 <- list(TMP = genes_S_norm)

    ac <- magcat_acat_pathways(
      gene_results = gene_results_one,
      pathways     = pw1,
      gene_col     = gene_col,
      p_col        = p_col,
      min_p        = min_p,
      do_fix       = do_fix,
      B            = 0L,
      seed         = NULL,
      output       = FALSE
    )
    p_acat <- ac$acat_p[1]

    fi <- magcat_fisher_pathways(
      gene_results = gene_results_one,
      pathways     = pw1,
      gene_col     = gene_col,
      p_col        = p_col,
      min_p        = min_p,
      do_fix       = do_fix,
      output       = FALSE
    )
    p_fisher <- fi$fisher_p[1]

    tf <- magcat_soft_tfisher_pathways(
      gene_results     = gene_results_one,
      pathways         = pw1,
      gene_col         = gene_col,
      p_col            = p_col,
      tau1             = ptrunc,     # <--- your ptrunc maps to tau1
      min_p            = min_p,
      do_fix           = do_fix,
      B_perm           = 0L,
      seed             = NULL,
      analytic_logical = TRUE,
      output           = FALSE
    )
    p_tfisher <- tf$tfisher_p_analytic[1]

    st <- magcat_stouffer_pathways(
      gene_results     = gene_results_one,
      pathways         = pw1,
      gene_col         = gene_col,
      p_col            = p_col,
      weight_col       = weight_col,
      min_p            = min_p,
      do_fix           = do_fix,
      B_perm           = 0L,
      seed             = NULL,
      analytic_logical = TRUE,
      output           = FALSE
    )
    p_stouffer <- st$stouffer_p_analytic[1]

    mp <- magcat_minp_pathways(
      gene_results     = gene_results_one,
      pathways         = pw1,
      gene_col         = gene_col,
      p_col            = p_col,
      min_p            = min_p,
      do_fix           = do_fix,
      B_perm           = 0L,
      seed             = NULL,
      analytic_logical = TRUE,
      output           = FALSE
    )
    p_minp <- mp$minp_p_analytic[1]

    c(p_acat, p_fisher, p_tfisher, p_stouffer, p_minp)
  }

  if (!is.null(seed)) set.seed(seed)

  ## ---------- 1) analytic p-values via your wrappers ----------
  acat_tab <- magcat_acat_pathways(
    gene_results = gene_results,
    pathways     = pathways,
    species      = species,
    pmn_gene_col = pmn_gene_col,
    gene_col     = gene_col,
    p_col        = p_col,
    min_p        = min_p,
    do_fix       = do_fix,
    B            = 0L,
    seed         = seed,
    output       = FALSE
  )

  fisher_tab <- magcat_fisher_pathways(
    gene_results = gene_results,
    pathways     = pathways,
    species      = species,
    pmn_gene_col = pmn_gene_col,
    gene_col     = gene_col,
    p_col        = p_col,
    min_p        = min_p,
    do_fix       = do_fix,
    output       = FALSE
  )

  tf_tab <- magcat_soft_tfisher_pathways(
    gene_results     = gene_results,
    pathways         = pathways,
    species          = species,
    pmn_gene_col     = pmn_gene_col,
    gene_col         = gene_col,
    p_col            = p_col,
    tau1             = ptrunc,   # ptrunc -> tau1
    min_p            = min_p,
    do_fix           = do_fix,
    B_perm           = 0L,
    seed             = seed,
    analytic_logical = TRUE,
    output           = FALSE
  )

  st_tab <- magcat_stouffer_pathways(
    gene_results     = gene_results,
    pathways         = pathways,
    species          = species,
    pmn_gene_col     = pmn_gene_col,
    gene_col         = gene_col,
    p_col            = p_col,
    weight_col       = weight_col,
    min_p            = min_p,
    do_fix           = do_fix,
    B_perm           = 0L,
    seed             = seed,
    analytic_logical = TRUE,
    output           = FALSE
  )

  minp_tab <- magcat_minp_pathways(
    gene_results     = gene_results,
    pathways         = pathways,
    species          = species,
    pmn_gene_col     = pmn_gene_col,
    gene_col         = gene_col,
    p_col            = p_col,
    min_p            = min_p,
    do_fix           = do_fix,
    B_perm           = 0L,
    seed             = seed,
    analytic_logical = TRUE,
    output           = FALSE
  )

  ## ---------- 2) base result table ----------
  if (!all(c("pathway_id", "pathway_name", "n_genes") %in% names(acat_tab))) {
    stop("omni_pathways(): acat table must have pathway_id, pathway_name, n_genes.",
         call. = FALSE)
  }

  res <- acat_tab[, c("pathway_id", "pathway_name", "n_genes"), drop = FALSE]

  pcol_acat     <- .find_pcol(acat_tab,   "acat")
  pcol_fisher   <- .find_pcol(fisher_tab, "fisher")
  pcol_tfisher  <- .find_pcol(tf_tab,     "tfisher")
  pcol_stouffer <- .find_pcol(st_tab,     "stouffer")
  pcol_minp     <- .find_pcol(minp_tab,   "minp")

  if (is.null(pcol_acat))     stop("Could not find ACAT p-column in acat_tab.",      call. = FALSE)
  if (is.null(pcol_fisher))   stop("Could not find Fisher p-column in fisher_tab.",  call. = FALSE)
  if (is.null(pcol_tfisher))  stop("Could not find TFisher p-column in tf_tab.",     call. = FALSE)
  if (is.null(pcol_stouffer)) stop("Could not find Stouffer p-column in st_tab.",    call. = FALSE)
  if (is.null(pcol_minp))     stop("Could not find minP p-column in minp_tab.",      call. = FALSE)

  res$acat_p     <- acat_tab[  match(res$pathway_id, acat_tab$pathway_id),      pcol_acat     ]
  res$fisher_p   <- fisher_tab[match(res$pathway_id, fisher_tab$pathway_id),    pcol_fisher   ]
  res$tfisher_p  <- tf_tab[    match(res$pathway_id, tf_tab$pathway_id),        pcol_tfisher  ]
  res$stouffer_p <- st_tab[    match(res$pathway_id, st_tab$pathway_id),        pcol_stouffer ]
  res$minp_gene_p<- minp_tab[  match(res$pathway_id, minp_tab$pathway_id),      pcol_minp     ]

  # Backward-compat aliases (if your downstream code expects old names)
  res$wfisher_p <- res$fisher_p
  res$tpm_p     <- res$tfisher_p

  ## ---------- 3) drop singletons ----------
  if (remove_singletons) {
    keep <- !is.na(res$n_genes) & res$n_genes >= 2L
    res  <- res[keep, , drop = FALSE]
  }

  ## ---------- 4) analytic omni ----------
  omni_p <- rep(NA_real_, nrow(res))
  for (i in seq_len(nrow(res))) {
    pv <- c(res$acat_p[i], res$fisher_p[i], res$tfisher_p[i], res$stouffer_p[i], res$minp_gene_p[i])
    omni_p[i] <- .combine_omni(pv)
  }
  res$omni_p    <- omni_p
  res$omni_p_BH <- .p_adjust_BH(res$omni_p)

  ## ---------- 5) omni-level permutations (USING YOUR WRAPPERS) ----------
  res$omni_perm_p    <- NA_real_
  res$omni_perm_p_BH <- NA_real_

  if (B_perm > 0L) {

    # prep global maps once
    gr0 <- gene_results
    genes_all <- tolower(as.character(gr0[[gene_col]]))
    p_all     <- as.numeric(gr0[[p_col]])
    ok <- which(!is.na(genes_all) & genes_all != "" & !is.na(p_all) & p_all > 0 & p_all <= 1)
    genes_all <- genes_all[ok]
    p_all     <- p_all[ok]
    p_map0    <- stats::setNames(p_all, genes_all)

    # weights map if requested (for stouffer wrapper)
    w_map0 <- NULL
    if (!is.null(weight_col)) {
      if (!weight_col %in% names(gr0)) {
        stop("weight_col = '", weight_col, "' not found in gene_results.", call. = FALSE)
      }
      w_all <- as.numeric(gr0[[weight_col]])
      w_all <- w_all[ok]
      w_map0 <- stats::setNames(w_all, genes_all)
    }

    # build pathway list ONCE (same as you had)
    if (!is.null(pathways)) {
      if (is.data.frame(pathways)) {
        if (!all(c("pathway_id","gene_id") %in% names(pathways))) {
          stop("pathways data.frame must have pathway_id and gene_id", call. = FALSE)
        }
        p_list <- split(pathways$gene_id, pathways$pathway_id)
      } else if (is.list(pathways)) {
        p_list <- pathways
      } else stop("'pathways' must be list or data.frame.", call. = FALSE)
    } else if (!is.null(species)) {
      pw <- magcat_load_pathways(
        species  = species,
        gene_col = if (is.null(pmn_gene_col)) NULL else pmn_gene_col
      )
      p_list <- split(pw$gene_id, pw$pathway_id)
    } else {
      stop("Need 'pathways' or 'species' for permutations.", call. = FALSE)
    }
    p_list <- lapply(p_list, function(g) tolower(as.character(g)))
    p_list <- p_list[res$pathway_id]  # keep same order as res

    omni_perm_p <- rep(NA_real_, nrow(res))

    if (perm_mode == "resample_global") {

      for (i in seq_len(nrow(res))) {

        genes_S <- p_list[[i]]
        if (is.null(genes_S) || length(genes_S) < 2L) next
        if (is.na(res$omni_p[i])) next

        genes_S <- genes_S[genes_S %in% names(p_map0)]
        d <- length(genes_S)
        if (d < 2L) next

        omni_null <- numeric(B_perm)

        for (b in seq_len(B_perm)) {

          # shuffle p-values across gene IDs (global resample)
          perm_idx <- sample.int(length(p_all))
          p_map_b  <- stats::setNames(p_all[perm_idx], genes_all)

          p_i <- as.numeric(p_map_b[genes_S])

          gene_results_b <- data.frame(
            tmp_gene = genes_S,
            tmp_p    = p_i,
            stringsAsFactors = FALSE
          )
          names(gene_results_b)[1:2] <- c(gene_col, p_col)

          if (!is.null(weight_col)) {
            gene_results_b[[weight_col]] <- as.numeric(w_map0[genes_S])
          }

          pv5 <- .component_pvals_onepath(gene_results_b, genes_S)

          omni_null[b] <- .combine_omni(pv5)
        }

        omni_obs <- res$omni_p[i]
        omni_perm_p[i] <- (1 + sum(omni_null <= omni_obs, na.rm = TRUE)) / (B_perm + 1)
      }

      res$omni_perm_p    <- omni_perm_p
      res$omni_perm_p_BH <- .p_adjust_BH(res$omni_perm_p)
    }

    if (perm_mode == "mvn") {

      if (is.null(magma_genes_out) || !file.exists(magma_genes_out)) {
        stop("perm_mode='mvn' requires magma_genes_out = path to merged MAGMA *.genes.out",
             call. = FALSE)
      }

      genes_out_tab <- magma_read_genes_out(magma_genes_out, gene_col = gene_col, chr_col = "CHR")

      cor_pairs <- NULL
      if (!is.null(magma_cor_file)) {
        if (!file.exists(magma_cor_file)) {
          stop("magma_cor_file does not exist: ", magma_cor_file, call. = FALSE)
        }
        cor_pairs <- magma_read_gene_cor_pairs(magma_cor_file, gene1_col = 1, gene2_col = 2, r_col = 3)
      }

      for (i in seq_len(nrow(res))) {

        genes_S <- p_list[[i]]
        if (is.null(genes_S) || length(genes_S) < 2L) next
        if (is.na(res$omni_p[i])) next

        R_S <- magma_build_R_for_pathway(
          genes_S        = genes_S,
          genes_out_tab  = genes_out_tab,
          cor_pairs      = cor_pairs,
          gene_col       = gene_col,
          chr_col        = "CHR",
          make_PD        = make_PD
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

        genes_sim <- colnames(sim$P)
        if (is.null(genes_sim)) genes_sim <- rownames(R_S)

        omni_null <- numeric(B_perm)

        for (b in seq_len(B_perm)) {

          p_i <- as.numeric(sim$P[b, ])

          gene_results_b <- data.frame(
            tmp_gene = genes_sim,
            tmp_p    = p_i,
            stringsAsFactors = FALSE
          )
          names(gene_results_b)[1:2] <- c(gene_col, p_col)

          if (!is.null(weight_col)) {
            gene_results_b[[weight_col]] <- as.numeric(w_map0[genes_sim])
          }

          pv5 <- .component_pvals_onepath(gene_results_b, tolower(genes_sim))

          omni_null[b] <- .combine_omni(pv5)
        }

        omni_obs <- res$omni_p[i]
        omni_perm_p[i] <- (1 + sum(omni_null <= omni_obs, na.rm = TRUE)) / (B_perm + 1)
      }

      res$omni_perm_p    <- omni_perm_p
      res$omni_perm_p_BH <- .p_adjust_BH(res$omni_perm_p)
    }
  }

  ## ---------- final p + sort ----------
  res$omni_final_p  <- if (B_perm > 0L) res$omni_perm_p else res$omni_p
  res$omni_final_BH <- .p_adjust_BH(res$omni_final_p)

  ord <- order(res$omni_final_p, decreasing = FALSE, na.last = TRUE)
  res <- res[ord, , drop = FALSE]

  ## ---------- optional CSV ----------
  if (output) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    species_tag <- if (is.null(species)) "custom" else species
    out_tag     <- if (omnibus == "ACAT") "acato" else "minp"
    out_path <- file.path(out_dir, paste0("omni_pathways_", species_tag, "_", out_tag, ".csv"))
    utils::write.csv(res, out_path, row.names = FALSE)
    attr(res, "file") <- out_path
  }

  res
}










#' Omnibus pathway test (ACAT-O or minP) using MAGCAT *prepared* objects
#'
#' Uses your wrapper pipeline for:
#'   (1) analytic method p-values  (via *_prepare() + *_run_prepared())
#'   (2) omni-level permutations   (resample_global or mvn)
#'
#' Methods included (ONLY): ACAT, Fisher, soft-TFisher, minP, Stouffer
#'

## ============================================================
##  IMPORTANT: patch your magcat_acat_prepared() to store genes_all_u
##  (so gene_names uses original case, not only lowercase)
## ============================================================
## Inside your magcat_acat_prepared(), add this into the returned list:
##   genes_all_u  = genes_all_u,
## (right next to genes_norm_u / p_all_u)


## ============================================================
##  OMNI (prepared analytic + fast permutations)
##  Methods: ACAT, Fisher, soft-TFisher, Stouffer, minP
##  Perm modes: resample_global, mvn
## ============================================================

#' @export
omni_pathways <- function(gene_results,
                          pathways     = NULL,
                          species      = NULL,
                          pmn_gene_col = NULL,
                          gene_col     = "GENE",
                          p_col        = "P",
                          weight_col   = NULL,     # for Stouffer weights (optional)
                          ptrunc       = 0.05,     # tau1 for soft TFisher
                          min_p        = 1e-15,
                          do_fix       = TRUE,
                          B_perm       = 0L,
                          perm_mode    = c("resample_global", "mvn"),
                          magma_genes_out = NULL,
                          magma_cor_file  = NULL,
                          make_PD      = TRUE,
                          seed         = NULL,
                          omnibus      = c("ACAT", "minP"),
                          remove_singletons = TRUE,
                          output       = FALSE,
                          out_dir      = "magcat_omni") {

  omnibus   <- match.arg(omnibus)
  perm_mode <- match.arg(perm_mode)

  if (!is.null(seed)) set.seed(seed)

  ## ---- require prepared API (NOT the old wrappers) ----
  req_funs <- c(
    "magcat_acat_prepare", "magcat_acat_run_prepared",
    "magcat_fisher_prepare", "magcat_fisher_run_prepared",
    "magcat_soft_tfisher_prepare", "magcat_soft_tfisher_run_prepared",
    "magcat_stouffer_prepare", "magcat_stouffer_run_prepared",
    "magcat_minp_prepare", "magcat_minp_run_prepared"
  )
  missing_f <- req_funs[!vapply(req_funs, exists, logical(1), mode = "function")]
  if (length(missing_f)) {
    stop("omni_pathways(): missing required MAGCAT functions: ",
         paste(missing_f, collapse = ", "),
         call. = FALSE)
  }

  ## ---- helpers ----
  .p_adjust_BH <- function(p) {
    p <- as.numeric(p)
    out <- rep(NA_real_, length(p))
    idx <- which(!is.na(p))
    if (length(idx)) out[idx] <- stats::p.adjust(p[idx], method = "BH")
    out
  }

  .sidak_minp <- function(p_vec) {
    p_vec <- as.numeric(p_vec)
    p_vec <- p_vec[!is.na(p_vec) & p_vec > 0 & p_vec < 1]
    if (!length(p_vec)) return(NA_real_)
    k <- length(p_vec)
    pmin <- min(p_vec)
    1 - (1 - pmin)^k
  }

  .combine_omni <- function(p_vec) {
    p_vec <- as.numeric(p_vec)
    p_vec <- p_vec[!is.na(p_vec) & p_vec > 0 & p_vec < 1]
    if (!length(p_vec)) return(NA_real_)

    if (omnibus == "ACAT") {
      if (!requireNamespace("sumFREGAT", quietly = TRUE)) {
        stop("omnibus='ACAT' requires sumFREGAT::ACATO().", call. = FALSE)
      }
      p_vec <- fix_p_for_acat(p_vec, min_p = min_p)
      return(tryCatch(as.numeric(sumFREGAT::ACATO(p_vec)), error = function(e) NA_real_))
    } else {
      return(.sidak_minp(p_vec))
    }
  }

  .stouffer_p_from_p <- function(p_i, w_i = NULL) {
    p_i <- fix_p_for_acat(p_i, min_p = min_p)
    z_i <- stats::qnorm(p_i / 2, lower.tail = FALSE)  # two-sided -> positive Z
    if (is.null(w_i)) w_i <- rep(1, length(z_i))
    w_i <- as.numeric(w_i)

    bad <- is.na(w_i) | w_i <= 0
    if (any(bad)) {
      w_pos <- w_i[!bad]
      repl <- if (length(w_pos)) stats::median(w_pos) else 1
      w_i[bad] <- repl
    }

    z <- sum(w_i * z_i) / sqrt(sum(w_i^2))
    2 * stats::pnorm(-abs(z))
  }

  ## ============================================================
  ## 1) PREP ONCE
  ## ============================================================

  acat_prep <- magcat_acat_prepare(
    gene_results = gene_results,
    pathways     = pathways,
    species      = species,
    pmn_gene_col = pmn_gene_col,
    gene_col     = gene_col,
    p_col        = p_col,
    min_p        = min_p,
    do_fix       = do_fix
  )

  fisher_prep <- magcat_fisher_prepare(
    gene_results = gene_results,
    pathways     = pathways,
    species      = species,
    pmn_gene_col = pmn_gene_col,
    gene_col     = gene_col,
    p_col        = p_col,
    min_p        = min_p,
    do_fix       = do_fix
  )

  tf_prep <- magcat_soft_tfisher_prepare(
    gene_results = gene_results,
    pathways     = pathways,
    species      = species,
    pmn_gene_col = pmn_gene_col,
    gene_col     = gene_col,
    p_col        = p_col,
    tau1         = ptrunc,
    min_p        = min_p,
    do_fix       = do_fix
  )

  st_prep <- magcat_stouffer_prepare(
    gene_results = gene_results,
    pathways     = pathways,
    species      = species,
    pmn_gene_col = pmn_gene_col,
    gene_col     = gene_col,
    p_col        = p_col,
    weight_col   = weight_col,
    min_p        = min_p,
    do_fix       = do_fix
  )

  minp_prep <- magcat_minp_prepare(
    gene_results = gene_results,
    pathways     = pathways,
    species      = species,
    pmn_gene_col = pmn_gene_col,
    gene_col     = gene_col,
    p_col        = p_col,
    min_p        = min_p,
    do_fix       = do_fix
  )

  ## ============================================================
  ## 2) ANALYTIC TABLES (FAST)
  ## ============================================================

  acat_tab   <- magcat_acat_run_prepared(acat_prep,   output = FALSE)
  fisher_tab <- magcat_fisher_run_prepared(fisher_prep, output = FALSE)
  tf_tab     <- magcat_soft_tfisher_run_prepared(tf_prep, output = FALSE)
  st_tab     <- magcat_stouffer_run_prepared(st_prep, output = FALSE)
  minp_tab   <- magcat_minp_run_prepared(minp_prep, output = FALSE)

  res <- data.frame(
    pathway_id   = acat_tab$pathway_id,
    pathway_name = acat_tab$pathway_name,
    n_genes      = acat_tab$n_genes,
    acat_p       = acat_tab$acat_p,
    fisher_p     = fisher_tab$fisher_p,
    tfisher_p    = tf_tab$tfisher_p_analytic,
    stouffer_p   = st_tab$stouffer_p_analytic,
    minp_gene_p  = minp_tab$minp_p_analytic,
    stringsAsFactors = FALSE
  )

  if (remove_singletons) {
    keep <- !is.na(res$n_genes) & res$n_genes >= 2L
    res  <- res[keep, , drop = FALSE]
  }

  ## ============================================================
  ## 3) OMNI (analytic)
  ## ============================================================

  res$omni_p <- NA_real_
  for (i in seq_len(nrow(res))) {
    pv <- c(res$acat_p[i], res$fisher_p[i], res$tfisher_p[i], res$stouffer_p[i], res$minp_gene_p[i])
    res$omni_p[i] <- .combine_omni(pv)
  }
  res$omni_p_BH <- .p_adjust_BH(res$omni_p)

  ## ============================================================
  ## 4) OMNI permutations (FAST; no wrapper calls inside)
  ## ============================================================

  res$omni_perm_p    <- NA_real_
  res$omni_perm_p_BH <- NA_real_

  if (B_perm > 0L) {

    # pathway gene sets from core (all of your prepare objects agree on p_list order)
    p_list <- fisher_prep$p_list
    p_list <- p_list[res$pathway_id]  # align with res rows

    # weights per pathway (fixed gene weights, like your shuffle model)
    w_norm <- st_prep$w_norm
    w_list <- vector("list", length(p_list))
    names(w_list) <- names(p_list)
    for (i in seq_along(p_list)) {
      g <- p_list[[i]]
      if (is.null(g)) { w_list[[i]] <- NULL; next }
      if (is.null(w_norm)) {
        w_list[[i]] <- rep(1, length(g))
      } else {
        w_list[[i]] <- as.numeric(w_norm[g])
      }
    }

    # pool of p-values for resample_global (already filtered to finite)
    pool_p <- as.numeric(fisher_prep$p_all)
    pool_ok <- which(is.finite(pool_p) & pool_p > 0 & pool_p < 1)
    pool_p <- pool_p[pool_ok]
    pool_p <- fix_p_for_acat(pool_p, min_p = min_p)

    omni_perm_p <- rep(NA_real_, nrow(res))

    if (perm_mode == "resample_global") {

      for (i in seq_len(nrow(res))) {

        genes_S <- p_list[[i]]
        if (is.null(genes_S) || length(genes_S) < 2L) next
        if (is.na(res$omni_p[i])) next

        # d = number of genes actually used (match observed n_genes)
        d <- res$n_genes[i]
        if (is.na(d) || d < 2L) next

        omni_null <- numeric(B_perm)

        w_i <- w_list[[i]]
        if (!is.null(w_i)) {
          # ensure weights length matches d
          if (length(w_i) != length(genes_S)) {
            # if your pathway had missing genes filtered in run_prepared,
            # just fall back to equal weights for permutations
            w_i <- rep(1, d)
          } else {
            w_i <- w_i[seq_len(min(length(w_i), d))]
          }
        }

        for (b in seq_len(B_perm)) {

          # sample p-values without replacement (exchangeable null)
          p_i <- sample(pool_p, size = d, replace = FALSE)

          # component p’s
          p_acat   <- ACAT::ACAT(Pvals = fix_p_for_acat(p_i, min_p = min_p))
          p_fisher <- stats::pchisq(-2 * sum(log(fix_p_for_acat(p_i, min_p=min_p))),
                                    df = 2 * d, lower.tail = FALSE)

          stat_tf <- TFisher::stat.soft(p = fix_p_for_acat(p_i, min_p=min_p), tau1 = ptrunc)
          p_tf    <- 1 - as.numeric(TFisher::p.soft(q = stat_tf, n = d, tau1 = ptrunc, M = NULL))

          p_st    <- .stouffer_p_from_p(p_i, w_i)

          p_min   <- .sidak_minp(p_i)

          omni_null[b] <- .combine_omni(c(p_acat, p_fisher, p_tf, p_st, p_min))
        }

        omni_obs <- res$omni_p[i]
        omni_perm_p[i] <- (1 + sum(omni_null <= omni_obs, na.rm = TRUE)) / (B_perm + 1)
      }

      res$omni_perm_p    <- omni_perm_p
      res$omni_perm_p_BH <- .p_adjust_BH(res$omni_perm_p)
    }

    if (perm_mode == "mvn") {

      if (is.null(magma_genes_out) || !file.exists(magma_genes_out)) {
        stop("perm_mode='mvn' requires magma_genes_out = path to merged MAGMA *.genes.out",
             call. = FALSE)
      }

      genes_out_tab <- magma_read_genes_out(magma_genes_out, gene_col = gene_col, chr_col = "CHR")

      cor_pairs <- NULL
      if (!is.null(magma_cor_file)) {
        if (!file.exists(magma_cor_file)) stop("magma_cor_file does not exist: ", magma_cor_file, call. = FALSE)
        cor_pairs <- magma_read_gene_cor_pairs(magma_cor_file, gene1_col = 1, gene2_col = 2, r_col = 3)
      }

      for (i in seq_len(nrow(res))) {

        genes_S <- p_list[[i]]
        if (is.null(genes_S) || length(genes_S) < 2L) next
        if (is.na(res$omni_p[i])) next

        R_S <- magma_build_R_for_pathway(
          genes_S       = genes_S,
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

        # weights fixed to pathway genes (in sim order)
        genes_sim <- colnames(sim$P)
        if (is.null(genes_sim)) genes_sim <- rownames(R_S)
        w_i <- if (is.null(st_prep$w_norm)) rep(1, length(genes_sim)) else as.numeric(st_prep$w_norm[tolower(genes_sim)])

        omni_null <- numeric(B_perm)

        for (b in seq_len(B_perm)) {

          p_i <- as.numeric(sim$P[b, ])
          d   <- length(p_i)

          p_acat   <- ACAT::ACAT(Pvals = fix_p_for_acat(p_i, min_p = min_p))
          p_fisher <- stats::pchisq(-2 * sum(log(fix_p_for_acat(p_i, min_p=min_p))),
                                    df = 2 * d, lower.tail = FALSE)

          stat_tf <- TFisher::stat.soft(p = fix_p_for_acat(p_i, min_p=min_p), tau1 = ptrunc)
          p_tf    <- 1 - as.numeric(TFisher::p.soft(q = stat_tf, n = d, tau1 = ptrunc, M = NULL))

          p_st    <- .stouffer_p_from_p(p_i, w_i)
          p_min   <- .sidak_minp(p_i)

          omni_null[b] <- .combine_omni(c(p_acat, p_fisher, p_tf, p_st, p_min))
        }

        omni_obs <- res$omni_p[i]
        omni_perm_p[i] <- (1 + sum(omni_null <= omni_obs, na.rm = TRUE)) / (B_perm + 1)
      }

      res$omni_perm_p    <- omni_perm_p
      res$omni_perm_p_BH <- .p_adjust_BH(res$omni_perm_p)
    }
  }

  ## final p
  res$omni_final_p  <- if (B_perm > 0L) res$omni_perm_p else res$omni_p
  res$omni_final_BH <- .p_adjust_BH(res$omni_final_p)

  ord <- order(res$omni_final_p, decreasing = FALSE, na.last = TRUE)
  res <- res[ord, , drop = FALSE]

  if (output) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    species_tag <- if (is.null(species)) "custom" else species
    out_tag     <- if (omnibus == "ACAT") "acato" else "minp"
    out_path <- file.path(out_dir, paste0("omni_pathways_", species_tag, "_", out_tag, ".csv"))
    utils::write.csv(res, out_path, row.names = FALSE)
    attr(res, "file") <- out_path
  }

  res
}
