# ==============================================================================
# Null Calibration Diagnostics for CATFISH Pathway Tests
#   - Fixes TFisher mislabeling (tfisher contains "fisher")
#   - Fixes "only omnibus" plots by (i) strict MVN checking and (ii) safe fallback
#   - Avoids omni_p_final by default in diagnostics (uses omni_p_mvn* if available)
# ==============================================================================

devtools::load_all(".")
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(data.table)

# ==============================================================================
# USER PARAMETERS
# ==============================================================================
N_NULL <- 100L
SEED   <- 42L
B_PERM <- 1000L   # permutations per pathway for MVN resampling

TAU_GRID <- c(0.1, 0.01,0.05)

PERM_MODE_ANALYTIC <- "none"
PERM_MODE_MVN      <- "mvn_global"

# If TRUE: MVN diagnostics require MVN component p-values to be non-NA.
# If FALSE: will fall back to analytic component p-values when MVN columns are NA.
STRICT_MVN_COMPONENTS <- TRUE

# ---------- Correlation files (REQUIRED for MVN mode) ----------
COR_BASE <- paste0(
  "/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/",
  "magma_multi_snp_wise_genes_by_chr_N_maize/"
)
COR_FILE_ARABIDOPSIS <- paste0(COR_BASE, "magma_gene_cor_pairs_MLM_arabidopsis.txt")
COR_FILE_FLY_MALE    <- paste0(COR_BASE, "magma_gene_cor_pairs_MLM_Fly_male.txt")
COR_FILE_FLY_FEMALE  <- paste0(COR_BASE, "magma_gene_cor_pairs_MLM_Fly_female.txt")

# ==============================================================================
# THEME
# ==============================================================================
plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.title     = element_text(
      size   = 24,
      face   = "bold",
      hjust  = 0.5,
      margin = margin(b = 10)
    ),
    axis.title.x   = element_text(
      size = 24,
      face = "bold"
    ),
    axis.title.y   = element_text(
      size = 24,
      face = "bold"
    ),
    axis.text.x    = element_text(
      size = 24,
      color = "black"
    ),
    axis.text.y    = element_text(
      size = 24,
      color = "black"
    ),
    axis.line      = element_line(color = "black"),

    ## <-- CHANGED (was element_blank())
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.25),

    legend.position = "top",
    legend.title    = element_blank(),
    legend.text     = element_text(size = 24)
  )
# ==============================================================================
# 1) GLOBAL NULL GENE GENERATOR
# ==============================================================================
generate_null_genes <- function(gene_template, min_p = 1e-10, max_p = 1 - 1e-10) {
  out <- gene_template
  n <- nrow(out)
  # Generate uniform p-values, clamped to avoid Inf in qnorm
  out$P <- runif(n)
  out$P <- pmax(pmin(out$P, max_p), min_p)
  out$ZSTAT <- qnorm(1 - out$P)  # N(0,1) under null
  # Double-check no Inf/NA
  out$ZSTAT[!is.finite(out$ZSTAT)] <- 0
  out
}

# ==============================================================================
# 2) SAFE CALL WRAPPER (prevents "unused argument" errors)
# ==============================================================================
safe_do_call <- function(fun, args) {
  fml <- names(formals(fun))
  args <- args[names(args) %in% fml]
  do.call(fun, args)
}

catfish_call <- function(...) {
  args <- list(...)
  safe_do_call(catfish_omni2_pathways, args)
}

# ==============================================================================
# 3) COLUMN UTILITIES
# ==============================================================================
non_na_n <- function(x) sum(!is.na(x))

# IMPORTANT: tfisher must be checked BEFORE fisher (tfisher contains "fisher")
component_from_col <- function(col) {
  if (grepl("tfisher", col, ignore.case = TRUE)) return("TFisher")
  if (grepl("stouffer", col, ignore.case = TRUE)) return("Stouffer")
  if (grepl("minp", col, ignore.case = TRUE)) return("minP")
  if (grepl("acat", col, ignore.case = TRUE)) return("ACAT")
  if (grepl("fisher", col, ignore.case = TRUE)) return("Fisher")
  if (grepl("^omni|omni_", col, ignore.case = TRUE)) return("Omnibus")
  return(col)
}

pick_first_usable <- function(res, candidates) {
  for (cn in candidates) {
    if (cn %in% names(res) && non_na_n(res[[cn]]) > 0) return(cn)
  }
  NA_character_
}

# Choose which p-value column to use for each component
# NOTE: For diagnostics we prefer omni_p_mvn* over omni_p_final.
get_pcol_map <- function(res, mode = c("analytic", "mvn"), strict_mvn = TRUE) {
  mode <- match.arg(mode)

  map_analytic <- list(
    ACAT     = c("acat_p"),
    Fisher   = c("fisher_p"),
    TFisher  = c("tfisher_p_analytic"),
    minP     = c("minp_p_analytic"),
    Stouffer = c("stouffer_p_analytic"),
    Omnibus  = c("omni_p_analytic")
  )

  map_mvn <- list(
    ACAT     = c("acat_p_mvn_cal", "acat_p_mvn"),
    Fisher   = c("fisher_p_mvn_cal", "fisher_p_mvn"),
    TFisher  = c("tfisher_p_mvn_cal", "tfisher_p_mvn"),
    minP     = c("minp_p_mvn_cal", "minp_p_mvn"),
    Stouffer = c("stouffer_p_mvn_cal", "stouffer_p_mvn"),
    Omnibus  = c("omni_p_mvn_compcal", "omni_p_mvn")  # DO NOT prefer omni_p_final
  )

  if (mode == "analytic") {
    out <- tibble(
      component = names(map_analytic),
      col = vapply(map_analytic, function(cands) pick_first_usable(res, cands), character(1))
    ) %>% filter(!is.na(col))
    return(out)
  }

  # MVN mode
  out_mvn <- tibble(
    component = names(map_mvn),
    col = vapply(map_mvn, function(cands) pick_first_usable(res, cands), character(1))
  )

  # If omnibus MVN is missing, last-resort fallback (but warn)
  if (is.na(out_mvn$col[out_mvn$component == "Omnibus"])) {
    if ("omni_p_final" %in% names(res) && non_na_n(res$omni_p_final) > 0) {
      message("WARNING: omni_p_mvn* not found/usable; falling back to omni_p_final for diagnostics (not ideal).")
      out_mvn$col[out_mvn$component == "Omnibus"] <- "omni_p_final"
    }
  }

  # Check MVN component availability
  missing_components <- out_mvn %>%
    filter(component != "Omnibus" & is.na(col)) %>%
    pull(component)

  if (length(missing_components) > 0) {
    msg <- paste0(
      "MVN component p-values missing/NA for: ",
      paste(missing_components, collapse = ", "),
      ".\nThis is exactly why you only see Omnibus.\n",
      "Either CATFISH isn't producing component MVN calibration in this call, ",
      "or those columns are being returned as all-NA.\n"
    )
    if (strict_mvn) stop(msg)
    message("WARNING: ", msg, "Falling back to ANALYTIC columns for missing components.")
    # fallback: fill missing with analytic
    out_ana <- get_pcol_map(res, mode = "analytic", strict_mvn = FALSE)
    out_mvn <- out_mvn %>%
      left_join(out_ana, by = "component", suffix = c("_mvn", "_ana")) %>%
      mutate(col = ifelse(is.na(col_mvn), col_ana, col_mvn)) %>%
      select(component, col)
  }

  out_mvn %>% filter(!is.na(col))
}

# ==============================================================================
# 4) ONE NULL RUN
# ==============================================================================
run_null_once <- function(gene_template, pathways,
                          gene_col, pmn_gene_col,
                          perm_mode = "none",
                          magma_cor_file = NULL,
                          B_perm = 1000L,
                          seed = NULL,
                          mvn_marginal = "uniform",
                          mvn_calibrate_components = TRUE,
                          mvn_min_p = 1e-15,
                          make_PD = TRUE,
                          ...) {
  if (!is.null(seed)) set.seed(seed)

  null_genes <- generate_null_genes(gene_template)

  # Check for any Inf/NA in generated null genes
  if (any(!is.finite(null_genes$P))) {
    cat("WARNING: null_genes$P has non-finite values\n")
    null_genes$P[!is.finite(null_genes$P)] <- 0.5
  }
  if (any(!is.finite(null_genes$ZSTAT))) {
    cat("WARNING: null_genes$ZSTAT has non-finite values\n")
    null_genes$ZSTAT[!is.finite(null_genes$ZSTAT)] <- 0
  }

  catfish_call(
    gene_results = null_genes,
    pathways     = pathways,
    species      = NULL,
    gene_col     = gene_col,
    pmn_gene_col = pmn_gene_col,
    p_raw_col    = "P",
    z_col        = "ZSTAT",
    tau_grid     = TAU_GRID,
    min_p        = mvn_min_p,
    do_fix       = TRUE,
    stouffer_min_abs_w   = 1e-8,
    stouffer_alternative = "greater",
    include_magma_in_omni = FALSE,
    include_magma_in_perm = FALSE,
    omnibus    = "ACAT",
    perm_mode  = perm_mode,
    B_perm     = B_perm,
    magma_cor_file = magma_cor_file,
    make_PD    = make_PD,
    mvn_marginal = mvn_marginal,
    mvn_calibrate_components = mvn_calibrate_components,
    output     = FALSE,
    ...
  )
}

# ==============================================================================
# 5) POOL N_NULL RUNS
# ==============================================================================
run_null_diagnostics <- function(gene_template, pathways,
                                 n_null = 100L, seed = 42L,
                                 label = "Dataset",
                                 perm_mode = c("none", "mvn", "mvn_global"),
                                 magma_cor_file = NULL,
                                 B_perm = 1000L,
                                 strict_mvn = TRUE,
                                 gene_col = "GENE",
                                 pmn_gene_col = "Gene-name",
                                 ...) {
  perm_mode <- match.arg(perm_mode)

  # Infer mode from perm_mode: mvn/mvn_global -> "mvn", none -> "analytic"
  mode <- if (perm_mode %in% c("mvn", "mvn_global")) "mvn" else "analytic"

  # MVN/mvn_global mode requires a correlation file
  if (perm_mode %in% c("mvn", "mvn_global") && is.null(magma_cor_file)) {
    stop("MVN/mvn_global mode requires magma_cor_file. ",
         "Either provide the path or use perm_mode='none' for analytic.")
  }

  set.seed(seed)
  cat("\n=========================================\n")
  cat("Null diagnostics:", label, "\n")
  cat("  perm_mode =", perm_mode, "(mode inferred:", mode, ")\n")
  cat("  B_perm    =", B_perm, "\n")
  cat("  N_NULL    =", n_null, "\n")
  cat("  Genes     =", nrow(gene_template), "\n")
  cat("=========================================\n")

  results_list <- vector("list", n_null)

  for (i in seq_len(n_null)) {
    if (i == 1 || i %% 10 == 0) cat("  Iteration", i, "of", n_null, "\n")

    res <- tryCatch(
      run_null_once(
        gene_template, pathways,
        gene_col = gene_col,
        pmn_gene_col = pmn_gene_col,
        perm_mode = perm_mode,
        magma_cor_file = magma_cor_file,
        B_perm = B_perm,
        seed = seed + i,
        ...
      ),
      error = function(e) {
        cat("    ERROR in iteration", i, ":", e$message, "\n")
        NULL
      }
    )
    if (is.null(res)) next

    # Show non-NA counts once (first iteration) so you can SEE what's populated
    if (i == 1) {
      pcols_to_check <- c(
        "acat_p_mvn_cal","fisher_p_mvn_cal","tfisher_p_mvn_cal",
        "minp_p_mvn_cal","stouffer_p_mvn_cal",
        "omni_p_mvn","omni_p_mvn_compcal","omni_p_final",
        "acat_p","fisher_p","tfisher_p_analytic","minp_p_analytic","stouffer_p_analytic","omni_p_analytic"
      )
      pcols_to_check <- intersect(pcols_to_check, names(res))
      cat("\n  Non-NA counts (iteration 1):\n")
      print(setNames(lapply(pcols_to_check, function(cc) non_na_n(res[[cc]])), pcols_to_check))
      cat("\n")
    }

    pmap <- get_pcol_map(res, mode = mode, strict_mvn = strict_mvn)
    keep_cols <- unique(pmap$col)

    if (!("pathway_id" %in% names(res))) stop("Result lacks pathway_id column.")

    results_list[[i]] <- res %>%
      select(pathway_id, all_of(keep_cols)) %>%
      mutate(null_iter = i)
  }

  null_df <- bind_rows(results_list) %>%
    mutate(dataset = label, mode = mode, perm_mode = perm_mode)

  n_ok <- sum(!sapply(results_list, is.null))
  n_pw <- length(unique(null_df$pathway_id))
  cat("  Done:", n_ok, "successful iterations |",
      n_pw, "pathways |", nrow(null_df), "rows\n")

  null_df
}

# ==============================================================================
# 6) METRICS
# ==============================================================================
compute_lambda <- function(p_values) {
  p <- p_values[!is.na(p_values) & p_values > 0 & p_values < 1]
  if (length(p) < 10) return(NA_real_)
  median(qchisq(1 - p, df = 1)) / qchisq(0.5, df = 1)
}

compute_type1_error <- function(p_values, alpha = 0.05) {
  p <- p_values[!is.na(p_values)]
  if (length(p) == 0) return(NA_real_)
  mean(p < alpha)
}

null_summary_table <- function(null_df) {
  cols <- setdiff(names(null_df), c("pathway_id", "null_iter", "dataset", "mode", "perm_mode"))
  if (length(cols) == 0) stop("No p-value columns found in null_df.")

  tibble(col = cols) %>%
    mutate(component = vapply(col, component_from_col, character(1))) %>%
    rowwise() %>%
    mutate(
      n_pvals        = sum(!is.na(null_df[[col]])),
      lambda         = compute_lambda(null_df[[col]]),
      type1_alpha005 = compute_type1_error(null_df[[col]], 0.05),
      type1_alpha001 = compute_type1_error(null_df[[col]], 0.01),
      mean_p         = mean(null_df[[col]], na.rm = TRUE),
      median_p       = median(null_df[[col]], na.rm = TRUE)
    ) %>%
    ungroup() %>%
    arrange(factor(component, levels = c("ACAT","Fisher","TFisher","minP","Stouffer","Omnibus")))
}

# ==============================================================================
# 7) PLOTS
# ==============================================================================
qq_plot_uniform <- function(p_values, title = "", lambda_val = NULL) {
  p <- sort(p_values[!is.na(p_values) & p_values > 0 & p_values < 1])
  n <- length(p)
  if (n == 0) return(ggplot() + ggtitle(title))

  df <- data.frame(
    expected = sort(-log10(ppoints(n))),
    observed = -log10(p)
  )

  sub_txt <- if (!is.null(lambda_val)) sprintf("lambda = %.3f", lambda_val) else ""

  ggplot(df, aes(x = expected, y = observed)) +
    geom_abline(slope = 1, intercept = 0,
                color = "red", linetype = "dashed", linewidth = 0.8) +
    geom_point(alpha = 0.25, size = 0.8, color = "steelblue") +
    labs(
      title    = title,
      subtitle = sub_txt,
      x = expression(Expected ~ ~-log[10](p)),
      y = expression(Observed ~ ~-log[10](p))
    ) +
    # Allow panels to fill the available width in the patchwork layout
    plot_theme
}

qq_panel <- function(null_df, dataset_label) {
  summ <- null_summary_table(null_df) %>% filter(n_pvals > 0)

  plots <- lapply(seq_len(nrow(summ)), function(i) {
    qq_plot_uniform(null_df[[summ$col[i]]],
                    title = summ$component[i],
                    lambda_val = summ$lambda[i])
  })

  wrap_plots(plots, ncol = 3) +
    plot_annotation(
      title = paste0("QQ Plots Under Global Null — ", dataset_label),
      subtitle = paste0(
        "mode=", unique(null_df$mode),
        " | perm_mode=", unique(null_df$perm_mode),
        " | N_null=", length(unique(null_df$null_iter))
      ),
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5)
      )
    )
}

lambda_bar_plot <- function(summary_tbl, dataset_label) {
  df <- summary_tbl %>% filter(n_pvals > 0)
  df$component <- factor(df$component, levels = unique(df$component))

  ggplot(df, aes(x = component, y = lambda)) +
    geom_col(fill = "steelblue", alpha = 0.8, width = 0.7) +
    geom_hline(yintercept = 1, linetype = "dashed",
               color = "red", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.3f", lambda)),
              vjust = -0.5, size = 4.2) +
    labs(
      title = paste0("Genomic Control lambda — ", dataset_label),
      subtitle = "Expected: lambda ~ 1 under null",
      x = "", y = expression(lambda)
    ) +
    ylim(0, max(c(df$lambda, 1.1), na.rm = TRUE) * 1.15) +
    plot_theme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

type1_bar_plot <- function(summary_tbl, dataset_label) {
  df <- summary_tbl %>% filter(n_pvals > 0)

  t1_long <- df %>%
    select(component, type1_alpha005, type1_alpha001) %>%
    pivot_longer(cols = starts_with("type1"),
                 names_to = "alpha_level",
                 values_to = "error_rate") %>%
    mutate(alpha_level = ifelse(alpha_level == "type1_alpha005",
                                "alpha = 0.05", "alpha = 0.01"))

  ref_lines <- data.frame(alpha_level = c("alpha = 0.05", "alpha = 0.01"),
                          expected = c(0.05, 0.01))

  t1_long$component <- factor(t1_long$component, levels = unique(df$component))

  ggplot(t1_long, aes(x = component, y = error_rate, fill = alpha_level)) +
    geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
    geom_hline(data = ref_lines,
               aes(yintercept = expected, color = alpha_level),
               linetype = "dashed", linewidth = 0.7, show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.3f", error_rate)),
              position = position_dodge(0.8),
              vjust = -0.5, size = 3.5) +
    scale_fill_manual(values = c("alpha = 0.05" = "steelblue",
                                 "alpha = 0.01" = "coral")) +
    scale_color_manual(values = c("alpha = 0.05" = "steelblue",
                                  "alpha = 0.01" = "coral")) +
    labs(
      title = paste0("Empirical Type I Error — ", dataset_label),
      subtitle = "Dashed lines = nominal alpha",
      x = "", y = "Rejection Rate", fill = ""
    ) +
    ylim(0, max(t1_long$error_rate, na.rm = TRUE) * 1.35) +
    plot_theme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

run_and_save <- function(null_df, dataset_label, prefix) {
  cat("\n--- Generating diagnostics:", dataset_label, "---\n")

  summ <- null_summary_table(null_df)
  print(summ, digits = 4)

  write.csv(summ, paste0(prefix, "_summary.csv"), row.names = FALSE)

  p1 <- qq_panel(null_df, dataset_label)
  ggsave(paste0(prefix, "_qq.png"), p1, width = 14, height = 10, dpi = 300, bg = "white")

  p2 <- lambda_bar_plot(summ, dataset_label)
  ggsave(paste0(prefix, "_lambda.png"), p2, width = 10, height = 7, dpi = 300, bg = "white")

  p3 <- type1_bar_plot(summ, dataset_label)
  ggsave(paste0(prefix, "_type1.png"), p3, width = 10, height = 7, dpi = 300, bg = "white")

  cat("Saved:", prefix, "_{summary,qq,lambda,type1}\n")
  invisible(summ)
}

# ==============================================================================
# 8) DATA LOADING HELPERS (your original)
# ==============================================================================
load_gene_template <- function(magma_dir, magma_pattern,
                               gene_len_file,
                               gene_id_prefix = "") {
  files <- list.files(magma_dir, magma_pattern, full.names = TRUE)
  cat("  MAGMA files:", length(files), "\n")

  genes_raw <- do.call(rbind, lapply(files, function(f) {
    read.table(f, header = TRUE, stringsAsFactors = FALSE, comment.char = "#")
  }))

  if (!"P" %in% names(genes_raw) && ncol(genes_raw) >= 9) {
    colnames(genes_raw)[9] <- "P"
  }

  genes_raw <- genes_raw[order(genes_raw$GENE, genes_raw$P), ]
  genes_raw <- genes_raw[!duplicated(genes_raw$GENE), ]

  gene_len <- read.delim(gene_len_file)
  if (nchar(gene_id_prefix) > 0) {
    gene_len$gene_id <- sub(paste0("^", gene_id_prefix), "", gene_len$gene_id)
  }

  adj <- catfish_adjust_gene_p(
    gene_results = genes_raw,
    gene_lengths = gene_len,
    gene_col     = "GENE",
    nsnp_col     = "NSNPS",
    p_col        = "P",
    z_col        = "ZSTAT",
    len_gene_col = "gene_id",
    len_col      = "length"
  )

  template <- data.frame(
    GENE  = adj[[1]],
    ZSTAT = adj[[2]],
    P     = adj[[3]],
    stringsAsFactors = FALSE
  )
  cat("  Template genes:", nrow(template), "\n")
  template
}

# ==============================================================================
# 9) RUN: All Datasets
# ==============================================================================

# Nice names for plotting
NICE_NAMES <- c("ACAT", "Fisher", "TFisher", "minP", "Stouffer", "Omnibus")



# ---------- ARABIDOPSIS ----------
arab_pw=read.delim("inst/extdata/pathway/aracyc_pathways.20230103")

arab_pw <- arab_pw %>%
  transmute(
    pathway_id   = Pathway.id,
    pathway_name = Pathway.name,
    gene_id      = Gene.id
  ) %>%
  mutate(
    gene_id = sub("^gene:", "", gene_id)  # safe even if no prefix
  ) %>%
  filter(!is.na(gene_id), gene_id != "") %>%
  distinct()
head(arab_pw)

arab_template <- load_gene_template(
  magma_dir = paste0(
    "/Users/nirwantandukar/Documents/Research/results/Arabidopsis/MAGMA/AT_cold_by_chr"
  ),
  magma_pattern = "\\.genes\\.out$",
  gene_len_file = "inst/extdata/Arabidopsis_gene_lengths.tsv",
  gene_id_prefix = "gene:"
)


null_arab_mvn <- run_null_diagnostics(
  gene_template  = arab_template,
  pathways       = arab_pw,
  n_null         = N_NULL,
  seed           = SEED,
  label          = "Arabidopsis (BIO6)",
  perm_mode      = PERM_MODE_MVN,
  magma_cor_file = COR_FILE_ARABIDOPSIS,
  B_perm         = B_PERM,
  strict_mvn     = STRICT_MVN_COMPONENTS,
  gene_col       = "GENE",
  pmn_gene_col   = "Gene-name"
)

summ_at <- run_and_save(null_arab_mvn, "Arabidopsis — MVN", "null_diag_arabidopsis_mvn")

# ---------- FLY MALE ----------
fly_pw <- read.delim("inst/extdata/pathway/Fly_Cyc.tsv")
head(fly_pw)

fly_m_template <- load_gene_template(
  magma_dir = paste0(
    "/Users/nirwantandukar/Documents/Research/",
    "results/DGRP/MAGMA/Fly_magma_genes_by_chr_male"
  ),
  magma_pattern = "\\.genes\\.out$",
  gene_len_file = "inst/extdata/dmel.flybase.fbgn.genes.loc.tsv"
)

null_fly_m_mvn <- run_null_diagnostics(
  gene_template  = fly_m_template,
  pathways       = fly_pw,
  n_null         = N_NULL,
  seed           = SEED + 1000L,
  label          = "Fly Male (Starvation)",
  perm_mode      = PERM_MODE_MVN,
  magma_cor_file = COR_FILE_FLY_MALE,
  B_perm         = B_PERM,
  strict_mvn     = STRICT_MVN_COMPONENTS,
  gene_col       = "GENE",
  pmn_gene_col   = "Gene-name"
)

summ_fly_m <- run_and_save(null_fly_m_mvn, "Fly Male — MVN", "null_diag_fly_male_mvn")

# ---------- FLY FEMALE ----------
fly_f_template <- load_gene_template(
  magma_dir = paste0(
    "/Users/nirwantandukar/Documents/Research/",
    "results/DGRP/MAGMA/Fly_magma_genes_by_chr_female"
  ),
  magma_pattern = "\\.genes\\.out$",
  gene_len_file = "inst/extdata/dmel.flybase.fbgn.genes.loc.tsv"
)

# ---- MVN diagnostics (STRICT by default)
null_fly_f_mvn <- run_null_diagnostics(
  gene_template  = fly_f_template,
  pathways       = fly_pw,
  n_null         = N_NULL,
  seed           = SEED + 5000L,
  label          = "Fly Female (Starvation)",
  perm_mode      = PERM_MODE_MVN,
  make_PD        = TRUE,
  magma_cor_file = COR_FILE_FLY_FEMALE,
  B_perm         = B_PERM,
  strict_mvn     = STRICT_MVN_COMPONENTS,
  gene_col       = "GENE",
  pmn_gene_col   = "Gene-name"
)
saveRDS(null_fly_f_mvn,"null_fly_f_mvn.RDS")
summ_fly_f <- run_and_save(null_fly_f_mvn, "Fly Female — MVN", "null_diag_fly_female_mvn")

# ---- If you want a fallback run that *always plots components* (uses analytic if MVN NA):
# STRICT_MVN_COMPONENTS <- FALSE
# (rerun the block above)







# ==============================================================================
# 10. COMBINED COMPARISON ACROSS DATASETS
# ==============================================================================

cat("\n\n##############################\n")
cat("## COMBINED COMPARISON\n")
cat("##############################\n")

# Load saved null diagnostics results
null_arab_mvn <- readRDS("null_diagnostics/null_arab_mvn.RDS")
null_fly_f_mvn <- readRDS("null_diagnostics/null_fly_f_mvn.RDS")

# Compute summary statistics for each dataset
summ_at <- null_summary_table(null_arab_mvn) %>% mutate(dataset = "Arabidopsis")
summ_fly_f <- null_summary_table(null_fly_f_mvn) %>% mutate(dataset = "Fly Female")

# Merge summaries
combined_summ <- bind_rows(summ_at, summ_fly_f)
combined_summ$component <- factor(
  combined_summ$component,
  levels = c("ACAT", "Fisher", "TFisher", "minP", "Stouffer", "Omnibus")
)

write.csv(combined_summ, "null_diag_combined_summary.csv", row.names = FALSE)

cat("\nCombined summary:\n")
print(combined_summ, digits = 4)

# ==============================================================================
# Build QQ data for combined plotting
# ==============================================================================
build_qq_data <- function(null_df, dataset_label) {
  summ <- null_summary_table(null_df) %>% filter(n_pvals > 0)

  qq_list <- lapply(seq_len(nrow(summ)), function(i) {
    p <- sort(null_df[[summ$col[i]]][!is.na(null_df[[summ$col[i]]]) &
                                       null_df[[summ$col[i]]] > 0 &
                                       null_df[[summ$col[i]]] < 1])
    n <- length(p)
    if (n == 0) return(NULL)

    data.frame(
      expected  = sort(-log10(ppoints(n))),
      observed  = -log10(p),
      component = summ$component[i],
      lambda    = summ$lambda[i],
      dataset   = dataset_label,
      stringsAsFactors = FALSE
    )
  })
  bind_rows(qq_list)
}

qq_combined_df <- bind_rows(
  build_qq_data(null_arab_mvn, "Arabidopsis"),
  build_qq_data(null_fly_f_mvn, "Fly Female")
)

qq_combined_df$component <- factor(
  qq_combined_df$component,
  levels = c("ACAT", "Fisher", "TFisher", "minP", "Stouffer", "Omnibus")
)
qq_combined_df$dataset <- factor(
  qq_combined_df$dataset,
  levels = c("Arabidopsis", "Fly Female")
)

# Dataset colors
dataset_colors <- c("Arabidopsis" = "#2E7D32", "Fly Female" = "#D84315")

# ==============================================================================
# Create individual QQ plots for each component (A panels - 6 plots)
# ==============================================================================
component_names <- c("ACAT", "Fisher", "TFisher", "minP", "Stouffer", "Omnibus")

qq_plots <- lapply(component_names, function(comp) {
  df_comp <- qq_combined_df %>% filter(component == comp)

  # Get lambda values for subtitle
  lambda_vals <- combined_summ %>%
    filter(component == comp) %>%
    arrange(dataset)
  lambda_text <- paste(
    sapply(seq_len(nrow(lambda_vals)), function(j) {
      sprintf("%s: %.2f", lambda_vals$dataset[j], lambda_vals$lambda[j])
    }),
    collapse = " | "
  )

  ggplot(df_comp, aes(x = expected, y = observed, color = dataset)) +
    geom_abline(slope = 1, intercept = 0, color = "black",
                linetype = "dashed", linewidth = 0.5) +
    geom_point(alpha = 0.4, size = 0.8) +
    scale_color_manual(values = dataset_colors) +
    labs(
      title = comp,
      subtitle = bquote(lambda ~ ": " ~ .(lambda_text)),
      x = expression(-log[10](expected)),
      y = expression(-log[10](observed))
    ) +
    # Avoid fixed aspect ratio so the QQ panel fills the row width
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 8, hjust = 0.5),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    ) + plot_theme
})

# Combine 6 QQ plots in 3 rows x 2 columns
qq_panel <- wrap_plots(qq_plots, ncol = 2, nrow = 3)

# ==============================================================================
# Lambda bar plot (B panel)
# ==============================================================================
p_lambda <- ggplot(
  combined_summ %>% filter(!is.na(lambda)),
  aes(x = component, y = lambda, fill = dataset)
) +
  geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.85) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.7) +
  geom_text(
    aes(label = sprintf("%.2f", lambda)),
    position = position_dodge(0.8),
    vjust = -0.3, size = 2.5
  ) +
  scale_fill_manual(values = dataset_colors) +
  labs(
    title = "Genomic Control (lambda)",
    x = "", y = expression(lambda), fill = ""
  ) +
  ylim(0, max(combined_summ$lambda, na.rm = TRUE) * 1.25) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  ) + plot_theme

# ==============================================================================
# Type I error bar plot (C panel)
# ==============================================================================
p_type1 <- ggplot(
  combined_summ %>% filter(!is.na(type1_alpha005)),
  aes(x = component, y = type1_alpha005, fill = dataset)
) +
  geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.85) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.7) +
  geom_text(
    aes(label = sprintf("%.3f", type1_alpha005)),
    position = position_dodge(0.8),
    vjust = -0.3, size = 2.5
  ) +
  scale_fill_manual(values = dataset_colors) +
  labs(
    title = "Type I Error (alpha = 0.05)",
    x = "", y = "Rejection Rate", fill = ""
  ) +
  ylim(0, max(combined_summ$type1_alpha005, na.rm = TRUE) * 1.35) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  ) + plot_theme

# ==============================================================================
# Combine all panels: 3 rows QQ (A) + 1 row with Lambda (B) and Type I (C)
# ==============================================================================
# Bottom row: lambda and type I error side by side
bottom_row <- p_lambda + p_type1 +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

# Full combined figure
# Heights: QQ panel (3 rows) vs bottom row - give QQ plots more vertical space
combined_fig <- qq_panel / bottom_row +
  plot_layout(heights = c(5, 1.5)) +
  plot_annotation(
    title = "Null Calibration Diagnostics",
    subtitle = paste0("MVN resampling | N_null = ", N_NULL, " iterations"),
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5)
    )
  ) 

# Save
quartz()
print(combined_fig)
ggsave(
  "Figures/null_calibration_combined.png", combined_fig,
  width = 24, height = 24, dpi = 300, bg = "white"
)
ggsave(
  "Figures/null_calibration_combined.pdf", combined_fig,
  width = 8, height = 10, bg = "white"
)

cat("\nSaved: Figures/null_calibration_combined.png\n")
cat("Saved: Figures/null_calibration_combined.pdf\n")

# ==============================================================================
# 11. FINAL REPORT
# ==============================================================================

cat("\n\n")
cat(strrep("=", 60), "\n")
cat("NULL CALIBRATION DIAGNOSTICS — COMPLETE\n")
cat(strrep("=", 60), "\n\n")

cat("Iterations per dataset:", N_NULL, "\n\n")

for (ds in c("Arabidopsis", "Fly Female")) {
  s <- combined_summ %>% filter(dataset == ds)
  cat("---", ds, "---\n")
  for (i in seq_len(nrow(s))) {
    cat(sprintf(
      "  %-12s  lambda = %.3f  |  T1@0.05 = %.3f  |  T1@0.01 = %.3f\n",
      s$component[i], s$lambda[i],
      s$type1_alpha005[i], s$type1_alpha001[i]
    ))
  }
  cat("\n")
}

cat("Interpretation guide:\n")
cat("  lambda ~ 1.0   : well-calibrated\n")
cat("  lambda > 1.0   : anti-conservative (inflated p)\n")
cat("  lambda < 1.0   : conservative (deflated p)\n")
cat("  T1@0.05 ~ 0.05 : correct type I error control\n")
cat("  T1@0.05 > 0.05 : anti-conservative\n")
cat("  T1@0.05 < 0.05 : conservative\n")

cat("\nOutput files:\n")
cat("  Per-dataset: null_diag_{arabidopsis,fly_female}_mvn_",
    "{summary.csv, qq.png, lambda.png, type1.png}\n")
cat("  Combined:    null_diag_combined_summary.csv\n")
cat("               Figures/null_calibration_combined.png\n")
cat("               Figures/null_calibration_combined.pdf\n")
cat("\n========== DONE ==========\n")
