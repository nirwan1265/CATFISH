## =============================================================================
## MAGCAT/CATFISH Validation & Diagnostics
## =============================================================================
##
## This script provides validation analyses to demonstrate that the MAGCAT
## pipeline is working correctly. It includes:
##
##   1. Gene-level adjustment validation (lm diagnostics)
##   2. P-value calibration (QQ plots)
##   3. Method comparison (correlation between methods)
##   4. Pathway enrichment visualization
##   5. MVN resampling validation
##
## =============================================================================

library(MAGCAT)
library(ggplot2)
library(dplyr)

## =============================================================================
## SECTION 1: Gene-Level Adjustment Validation
## =============================================================================
##
## The magcat_adjust_gene_p() function regresses out gene length and SNP count
## confounding from MAGMA gene-level Z-scores. This section validates that
## the adjustment successfully removes these biases.

#' Validate gene-level p-value adjustment
#'
#' Creates diagnostic plots showing that gene length and SNP count biases
#' are removed after adjustment.
#'
#' @param adj_result Output from magcat_adjust_gene_p() containing $fit
#' @param output_dir Directory to save plots (NULL for no saving)
#' @return List of ggplot objects and correlation statistics
validate_gene_adjustment <- function(adj_result, output_dir = NULL) {

 if (is.null(adj_result$fit)) {
   stop("adj_result must contain $fit (lm object from magcat_adjust_gene_p)")
 }

 lm_fit <- adj_result$fit

 # Build diagnostic data frame
 df_diag <- cbind(
   model.frame(lm_fit),
   Z_fitted = fitted(lm_fit),
   Z_resid  = resid(lm_fit)
 )

 # Calculate correlations before/after adjustment
 cor_before_len  <- cor(df_diag$Z_raw, df_diag$log_gene_length, use = "complete.obs")
 cor_after_len   <- cor(df_diag$Z_resid, df_diag$log_gene_length, use = "complete.obs")
 cor_before_nsnp <- cor(df_diag$Z_raw, df_diag$log_nsnp, use = "complete.obs")
 cor_after_nsnp  <- cor(df_diag$Z_resid, df_diag$log_nsnp, use = "complete.obs")

 correlations <- data.frame(
   variable = c("Gene Length", "Gene Length", "SNP Count", "SNP Count"),
   stage    = c("Before", "After", "Before", "After"),
   correlation = c(cor_before_len, cor_after_len, cor_before_nsnp, cor_after_nsnp)
 )

 # Plot 1: Raw Z vs log(gene length)
 p1 <- ggplot(df_diag, aes(x = log_gene_length, y = Z_raw)) +
   geom_hex(bins = 50) +
   geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1) +
   scale_fill_viridis_c() +
   labs(
     x = "log(Gene Length)",
     y = "Raw MAGMA Z-score",
     title = "Before Adjustment: Z vs Gene Length",
     subtitle = sprintf("Correlation: r = %.3f", cor_before_len)
   ) +
   theme_minimal() +
   theme(plot.title = element_text(face = "bold"))

 # Plot 2: Adjusted Z vs log(gene length)
 p2 <- ggplot(df_diag, aes(x = log_gene_length, y = Z_resid)) +
   geom_hex(bins = 50) +
   geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1) +
   scale_fill_viridis_c() +
   labs(
     x = "log(Gene Length)",
     y = "Adjusted Z-score (Residual)",
     title = "After Adjustment: Z vs Gene Length",
     subtitle = sprintf("Correlation: r = %.3f (bias removed)", cor_after_len)
   ) +
   theme_minimal() +
   theme(plot.title = element_text(face = "bold"))

 # Plot 3: Raw Z vs log(#SNPs)
 p3 <- ggplot(df_diag, aes(x = log_nsnp, y = Z_raw)) +
   geom_hex(bins = 50) +
   geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1) +
   scale_fill_viridis_c() +
   labs(
     x = "log(# SNPs per Gene)",
     y = "Raw MAGMA Z-score",
     title = "Before Adjustment: Z vs SNP Count",
     subtitle = sprintf("Correlation: r = %.3f", cor_before_nsnp)
   ) +
   theme_minimal() +
   theme(plot.title = element_text(face = "bold"))

 # Plot 4: Adjusted Z vs log(#SNPs)
 p4 <- ggplot(df_diag, aes(x = log_nsnp, y = Z_resid)) +
   geom_hex(bins = 50) +
   geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1) +
   scale_fill_viridis_c() +
   labs(
     x = "log(# SNPs per Gene)",
     y = "Adjusted Z-score (Residual)",
     title = "After Adjustment: Z vs SNP Count",
     subtitle = sprintf("Correlation: r = %.3f (bias removed)", cor_after_nsnp)
   ) +
   theme_minimal() +
   theme(plot.title = element_text(face = "bold"))

 # Save plots if output_dir provided
 if (!is.null(output_dir)) {
   dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
   ggsave(file.path(output_dir, "adjustment_gene_length_before.png"), p1, width = 8, height = 6)
   ggsave(file.path(output_dir, "adjustment_gene_length_after.png"), p2, width = 8, height = 6)
   ggsave(file.path(output_dir, "adjustment_snp_count_before.png"), p3, width = 8, height = 6)
   ggsave(file.path(output_dir, "adjustment_snp_count_after.png"), p4, width = 8, height = 6)
   cat("Plots saved to:", output_dir, "\n")
 }

 # Print summary
 cat("\n=== Gene-Level Adjustment Validation ===\n")
 cat("\nCorrelations with Gene Length:\n")
 cat(sprintf("  Before adjustment: r = %.4f\n", cor_before_len))
 cat(sprintf("  After adjustment:  r = %.4f\n", cor_after_len))
 cat("\nCorrelations with SNP Count:\n")
 cat(sprintf("  Before adjustment: r = %.4f\n", cor_before_nsnp))
 cat(sprintf("  After adjustment:  r = %.4f\n", cor_after_nsnp))
 cat("\nInterpretation: Correlations near 0 after adjustment indicate\n")
 cat("successful removal of gene size and SNP density confounding.\n")

 list(
   plots = list(
     gene_length_before = p1,
     gene_length_after  = p2,
     snp_count_before   = p3,
     snp_count_after    = p4
   ),
   correlations = correlations,
   lm_summary   = summary(lm_fit),
   df_diag      = df_diag
 )
}


## =============================================================================
## SECTION 2: P-Value Calibration (QQ Plots)
## =============================================================================
##
## QQ plots compare observed p-values to the expected uniform distribution
## under the null hypothesis. Well-calibrated p-values should follow the
## diagonal, with deviation at the tail indicating true signal.

#' Create QQ plot for p-values
#'
#' @param pvals Numeric vector of p-values
#' @param title Plot title
#' @param show_lambda If TRUE, show genomic inflation factor
#' @return ggplot object
qq_plot <- function(pvals, title = "QQ Plot", show_lambda = TRUE) {

 pvals <- pvals[!is.na(pvals) & pvals > 0 & pvals <= 1]
 n <- length(pvals)

 observed <- -log10(sort(pvals))
 expected <- -log10(ppoints(n))

 # Calculate genomic inflation factor (lambda)
 chisq <- qchisq(1 - pvals, 1)
 lambda <- median(chisq) / qchisq(0.5, 1)

 df <- data.frame(expected = expected, observed = observed)

 p <- ggplot(df, aes(x = expected, y = observed)) +
   geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
   geom_point(alpha = 0.5, size = 1) +
   labs(
     x = expression(Expected ~ -log[10](p)),
     y = expression(Observed ~ -log[10](p)),
     title = title
   ) +
   theme_minimal() +
   theme(plot.title = element_text(face = "bold"))

 if (show_lambda) {
   p <- p + annotate(
     "text", x = max(expected) * 0.1, y = max(observed) * 0.9,
     label = sprintf("lambda == %.3f", lambda),
     parse = TRUE, hjust = 0, size = 4
   )
 }

 p
}


#' Create QQ plots for pathway p-values from omnibus results
#'
#' @param omni_results Output from magcat_omni2_pathways()
#' @param output_dir Directory to save plots (NULL for no saving)
#' @return List of ggplot objects
validate_pvalue_calibration <- function(omni_results, output_dir = NULL) {

 plots <- list()

 # Component methods
 if ("acat_p" %in% names(omni_results)) {
   plots$acat <- qq_plot(omni_results$acat_p, "ACAT P-values")
 }
 if ("fisher_p" %in% names(omni_results)) {
   plots$fisher <- qq_plot(omni_results$fisher_p, "Fisher P-values")
 }
 if ("minp_p" %in% names(omni_results)) {
   plots$minp <- qq_plot(omni_results$minp_p, "MinP P-values")
 }
 if ("stouffer_p_analytic" %in% names(omni_results)) {
   plots$stouffer <- qq_plot(omni_results$stouffer_p_analytic, "Stouffer P-values")
 }

 # Omnibus p-values
 if ("omni_p_analytic" %in% names(omni_results)) {
   plots$omni_analytic <- qq_plot(omni_results$omni_p_analytic, "Omnibus (Analytic)")
 }
 if ("omni_p_global" %in% names(omni_results)) {
   plots$omni_global <- qq_plot(omni_results$omni_p_global, "Omnibus (Global Perm)")
 }
 if ("omni_p_mvn" %in% names(omni_results)) {
   plots$omni_mvn <- qq_plot(omni_results$omni_p_mvn, "Omnibus (MVN Calibrated)")
 }

 # Save plots
 if (!is.null(output_dir)) {
   dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
   for (name in names(plots)) {
     ggsave(
       file.path(output_dir, paste0("qq_", name, ".png")),
       plots[[name]], width = 6, height = 6
     )
   }
   cat("QQ plots saved to:", output_dir, "\n")
 }

 plots
}


## =============================================================================
## SECTION 3: Method Comparison
## =============================================================================
##
## Compare p-values across different methods to understand their correlation
## and identify pathways that are consistently significant.

#' Compare p-values across methods
#'
#' @param omni_results Output from magcat_omni2_pathways()
#' @param output_dir Directory to save plots (NULL for no saving)
#' @return List containing correlation matrix and plots
compare_methods <- function(omni_results, output_dir = NULL) {

 # Extract available p-value columns
 p_cols <- c("acat_p", "fisher_p", "minp_p", "stouffer_p_analytic",
             "tfisher_p_analytic", "omni_p_analytic", "omni_p_mvn")
 p_cols <- p_cols[p_cols %in% names(omni_results)]

 if (length(p_cols) < 2) {
   warning("Need at least 2 p-value columns to compare methods")
   return(NULL)
 }

 # Create -log10 transformed matrix
 p_mat <- as.matrix(omni_results[, p_cols])
 logp_mat <- -log10(p_mat)
 logp_mat[!is.finite(logp_mat)] <- NA

 # Correlation matrix
 cor_mat <- cor(logp_mat, use = "pairwise.complete.obs")

 # Nicer column names for display
 nice_names <- c(
   acat_p = "ACAT",
   fisher_p = "Fisher",
   minp_p = "MinP",
   stouffer_p_analytic = "Stouffer",
   tfisher_p_analytic = "TFisher",
   omni_p_analytic = "Omnibus",
   omni_p_mvn = "Omnibus (MVN)"
 )
 rownames(cor_mat) <- nice_names[rownames(cor_mat)]
 colnames(cor_mat) <- nice_names[colnames(cor_mat)]

 # Heatmap
 cor_df <- as.data.frame(as.table(cor_mat))
 names(cor_df) <- c("Method1", "Method2", "Correlation")

 p_heatmap <- ggplot(cor_df, aes(x = Method1, y = Method2, fill = Correlation)) +
   geom_tile() +
   geom_text(aes(label = sprintf("%.2f", Correlation)), size = 3) +
   scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5) +
   labs(title = "Correlation Between Methods (-log10 p-values)") +
   theme_minimal() +
   theme(
     axis.text.x = element_text(angle = 45, hjust = 1),
     plot.title = element_text(face = "bold")
   )

 # Scatter plot: ACAT vs Fisher (example)
 if (all(c("acat_p", "fisher_p") %in% names(omni_results))) {
   p_scatter <- ggplot(omni_results, aes(x = -log10(acat_p), y = -log10(fisher_p))) +
     geom_point(alpha = 0.5) +
     geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
     geom_smooth(method = "lm", se = FALSE, color = "blue") +
     labs(
       x = expression(-log[10](ACAT ~ p)),
       y = expression(-log[10](Fisher ~ p)),
       title = "ACAT vs Fisher P-values"
     ) +
     theme_minimal()
 } else {
   p_scatter <- NULL
 }

 # Save plots
 if (!is.null(output_dir)) {
   dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
   ggsave(file.path(output_dir, "method_correlation_heatmap.png"), p_heatmap, width = 8, height = 7)
   if (!is.null(p_scatter)) {
     ggsave(file.path(output_dir, "acat_vs_fisher_scatter.png"), p_scatter, width = 7, height = 6)
   }
   cat("Method comparison plots saved to:", output_dir, "\n")
 }

 cat("\n=== Method Correlation Matrix ===\n")
 print(round(cor_mat, 3))

 list(
   correlation_matrix = cor_mat,
   heatmap = p_heatmap,
   scatter = p_scatter
 )
}


## =============================================================================
## SECTION 4: Pathway Enrichment Visualization
## =============================================================================

#' Create volcano-style plot for pathway results
#'
#' @param omni_results Output from magcat_omni2_pathways()
#' @param p_col Which p-value column to use
#' @param fdr_threshold FDR threshold for significance
#' @param output_dir Directory to save plot
#' @return ggplot object
pathway_volcano <- function(omni_results,
                           p_col = "omni_p_mvn",
                           fdr_threshold = 0.05,
                           output_dir = NULL) {

 if (!p_col %in% names(omni_results)) {
   stop("Column '", p_col, "' not found in results")
 }

 df <- omni_results
 df$neglog10p <- -log10(df[[p_col]])
 df$fdr <- p.adjust(df[[p_col]], method = "BH")
 df$significant <- df$fdr < fdr_threshold

 # Label top pathways
 df <- df[order(df[[p_col]]), ]
 df$label <- ifelse(rank(df[[p_col]]) <= 10, df$pathway_name, NA)

 p <- ggplot(df, aes(x = n_genes, y = neglog10p)) +
   geom_point(aes(color = significant), alpha = 0.6) +
   geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
   scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "red")) +
   labs(
     x = "Number of Genes in Pathway",
     y = expression(-log[10](p)),
     title = "Pathway Enrichment Results",
     subtitle = sprintf("FDR < %.2f shown in red", fdr_threshold),
     color = "Significant"
   ) +
   theme_minimal() +
   theme(plot.title = element_text(face = "bold"))

 # Add labels for top pathways if ggrepel is available
 if (requireNamespace("ggrepel", quietly = TRUE)) {
   p <- p + ggrepel::geom_text_repel(
     aes(label = label),
     size = 2.5, max.overlaps = 15,
     na.rm = TRUE
   )
 }

 if (!is.null(output_dir)) {
   dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
   ggsave(file.path(output_dir, "pathway_volcano.png"), p, width = 10, height = 8)
   cat("Volcano plot saved to:", output_dir, "\n")
 }

 p
}


#' Create bar plot of top pathways
#'
#' @param omni_results Output from magcat_omni2_pathways()
#' @param p_col Which p-value column to use
#' @param top_n Number of top pathways to show
#' @param output_dir Directory to save plot
#' @return ggplot object
top_pathways_barplot <- function(omni_results,
                                 p_col = "omni_p_mvn",
                                 top_n = 20,
                                 output_dir = NULL) {

 df <- omni_results[order(omni_results[[p_col]]), ]
 df <- head(df, top_n)
 df$neglog10p <- -log10(df[[p_col]])

 # Truncate long pathway names
 df$pathway_label <- ifelse(
   nchar(df$pathway_name) > 40,
   paste0(substr(df$pathway_name, 1, 37), "..."),
   df$pathway_name
 )

 # Order by p-value
 df$pathway_label <- factor(df$pathway_label, levels = rev(df$pathway_label))

 p <- ggplot(df, aes(x = pathway_label, y = neglog10p)) +
   geom_col(fill = "steelblue") +
   geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
   coord_flip() +
   labs(
     x = NULL,
     y = expression(-log[10](p)),
     title = paste("Top", top_n, "Enriched Pathways"),
     subtitle = "Red dashed line: p = 0.05"
   ) +
   theme_minimal() +
   theme(
     plot.title = element_text(face = "bold"),
     axis.text.y = element_text(size = 8)
   )

 if (!is.null(output_dir)) {
   dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
   ggsave(file.path(output_dir, "top_pathways_barplot.png"), p, width = 10, height = 8)
   cat("Bar plot saved to:", output_dir, "\n")
 }

 p
}


## =============================================================================
## SECTION 5: MVN Resampling Validation
## =============================================================================
##
## Compare analytic p-values to MVN-calibrated p-values to assess the
## impact of accounting for gene-gene correlations.

#' Compare analytic vs MVN-calibrated p-values
#'
#' @param omni_results Output from magcat_omni2_pathways() with MVN enabled
#' @param output_dir Directory to save plots
#' @return List of plots and statistics
validate_mvn_calibration <- function(omni_results, output_dir = NULL) {

 if (!"omni_p_analytic" %in% names(omni_results) ||
     !"omni_p_mvn" %in% names(omni_results)) {
   stop("Results must contain both omni_p_analytic and omni_p_mvn columns")
 }

 df <- omni_results[!is.na(omni_results$omni_p_analytic) &
                    !is.na(omni_results$omni_p_mvn), ]

 df$logp_analytic <- -log10(df$omni_p_analytic)
 df$logp_mvn <- -log10(df$omni_p_mvn)

 # Scatter plot
 p_scatter <- ggplot(df, aes(x = logp_analytic, y = logp_mvn)) +
   geom_point(alpha = 0.5) +
   geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
   geom_smooth(method = "lm", se = FALSE, color = "blue") +
   labs(
     x = expression(-log[10](p[analytic])),
     y = expression(-log[10](p[MVN])),
     title = "Analytic vs MVN-Calibrated P-values",
     subtitle = "Points below diagonal: MVN is more conservative"
   ) +
   theme_minimal() +
   theme(plot.title = element_text(face = "bold"))

 # Histogram of log-ratio
 df$log_ratio <- log10(df$omni_p_mvn / df$omni_p_analytic)

 p_hist <- ggplot(df, aes(x = log_ratio)) +
   geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
   geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
   labs(
     x = expression(log[10](p[MVN] / p[analytic])),
     y = "Count",
     title = "Distribution of P-value Ratio",
     subtitle = "Positive = MVN is more conservative"
   ) +
   theme_minimal() +
   theme(plot.title = element_text(face = "bold"))

 # Statistics
 median_ratio <- median(df$log_ratio, na.rm = TRUE)
 pct_conservative <- mean(df$omni_p_mvn > df$omni_p_analytic, na.rm = TRUE) * 100
 correlation <- cor(df$logp_analytic, df$logp_mvn, use = "complete.obs")

 cat("\n=== MVN Calibration Summary ===\n")
 cat(sprintf("Correlation (analytic vs MVN): r = %.3f\n", correlation))
 cat(sprintf("Median log10 ratio (MVN/analytic): %.3f\n", median_ratio))
 cat(sprintf("Pathways where MVN is more conservative: %.1f%%\n", pct_conservative))
 cat("\nInterpretation: MVN calibration accounts for gene-gene LD,\n")
 cat("typically resulting in more conservative (larger) p-values.\n")

 # Save plots
 if (!is.null(output_dir)) {
   dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
   ggsave(file.path(output_dir, "mvn_vs_analytic_scatter.png"), p_scatter, width = 7, height = 6)
   ggsave(file.path(output_dir, "mvn_ratio_histogram.png"), p_hist, width = 7, height = 5)
   cat("\nMVN validation plots saved to:", output_dir, "\n")
 }

 list(
   scatter = p_scatter,
   histogram = p_hist,
   stats = list(
     correlation = correlation,
     median_log_ratio = median_ratio,
     pct_conservative = pct_conservative
   )
 )
}


## =============================================================================
## SECTION 6: Full Validation Report
## =============================================================================

#' Run all validation analyses
#'
#' @param gene_results Gene-level results from MAGMA
#' @param gene_lengths Gene length data
#' @param omni_results Output from magcat_omni2_pathways()
#' @param output_dir Directory to save all plots and report
#' @return List of all validation results
run_full_validation <- function(gene_results = NULL,
                                gene_lengths = NULL,
                                omni_results,
                                output_dir = "validation_results") {

 results <- list()

 cat("=============================================================\n")
 cat("           MAGCAT Validation Report\n")
 cat("=============================================================\n\n")

 # 1. Gene adjustment validation (if inputs provided)
 if (!is.null(gene_results) && !is.null(gene_lengths)) {
   cat("Running gene adjustment validation...\n")

   adj_out <- magcat_adjust_gene_p(
     gene_results = gene_results,
     gene_lengths = gene_lengths,
     gene_col     = "GENE",
     nsnp_col     = "NSNPS",
     p_col        = "P",
     z_col        = "ZSTAT"
   )

   results$adjustment <- validate_gene_adjustment(
     adj_out,
     output_dir = file.path(output_dir, "gene_adjustment")
   )
 }

 # 2. P-value calibration (QQ plots)
 cat("\nRunning p-value calibration checks...\n")
 results$qq_plots <- validate_pvalue_calibration(
   omni_results,
   output_dir = file.path(output_dir, "qq_plots")
 )

 # 3. Method comparison
 cat("\nComparing methods...\n")
 results$method_comparison <- compare_methods(
   omni_results,
   output_dir = file.path(output_dir, "method_comparison")
 )

 # 4. Pathway visualization
 cat("\nCreating pathway visualizations...\n")
 results$volcano <- pathway_volcano(
   omni_results,
   output_dir = file.path(output_dir, "pathway_plots")
 )
 results$top_pathways <- top_pathways_barplot(
   omni_results,
   output_dir = file.path(output_dir, "pathway_plots")
 )

 # 5. MVN validation (if available)
 if ("omni_p_mvn" %in% names(omni_results)) {
   cat("\nValidating MVN calibration...\n")
   results$mvn <- validate_mvn_calibration(
     omni_results,
     output_dir = file.path(output_dir, "mvn_validation")
   )
 }

 cat("\n=============================================================\n")
 cat("Validation complete! Results saved to:", output_dir, "\n")
 cat("=============================================================\n")

 invisible(results)
}


## =============================================================================
## Example Usage
## =============================================================================
##
## # After running the main workflow (usage2.R), validate results:
##
## # 1. Validate gene-level adjustment
## adj_out <- magcat_adjust_gene_p(
##   gene_results = genes_all,
##   gene_lengths = gene_lengths,
##   gene_col = "GENE", nsnp_col = "NSNPS",
##   p_col = "P", z_col = "ZSTAT"
## )
## adj_validation <- validate_gene_adjustment(adj_out, output_dir = "validation")
##
## # 2. Run full validation on omnibus results
## validation <- run_full_validation(
##   gene_results = genes_all,
##   gene_lengths = gene_lengths,
##   omni_results = omni_results,
##   output_dir   = "validation_results"
## )
##
## # 3. View individual plots
## print(validation$adjustment$plots$gene_length_before)
## print(validation$adjustment$plots$gene_length_after)
## print(validation$qq_plots$omni_mvn)
## print(validation$volcano)
