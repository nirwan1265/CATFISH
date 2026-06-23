# ==============================================================================
# Compare Fly Male vs Female: Component Tests and OMNIBUS
# Shows: Individual tests differ between sexes, but OMNIBUS converges
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

if (!requireNamespace("ggVennDiagram", quietly = TRUE)) {
  install.packages("ggVennDiagram")
}
library(ggVennDiagram)

if (!requireNamespace("UpSetR", quietly = TRUE)) {
  install.packages("UpSetR")
}
library(UpSetR)

if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}
library(patchwork)

# ==============================================================================
# USER PARAMETERS
# ==============================================================================

# ---------- Significance threshold mode ----------
# Options: "top_k", "top_pct", "bonferroni", "fdr", "nominal"
SIG_MODE <- "top_k"

# Parameters for each mode
TOP_K_PATHWAYS   <- 20      # for SIG_MODE = "top_k"
TOP_PCT_PATHWAYS <- 5       # for SIG_MODE = "top_pct"
FDR_THRESHOLD    <- 0.1     # for SIG_MODE = "fdr"
NOMINAL_ALPHA    <- 0.05    # for SIG_MODE = "nominal"
# Note: "bonferroni" uses 0.05 / n_pathways automatically

# ==============================================================================
# Plot Theme
# ==============================================================================

plot_theme <- theme(
  plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
  plot.subtitle = element_text(size = 11, hjust = 0.5),
  axis.title = element_text(size = 12),
  axis.text = element_text(size = 10),
  legend.title = element_text(size = 11),
  legend.text = element_text(size = 10),
  strip.text = element_text(size = 11, face = "bold")
)

# ==============================================================================
# 1. Load Results
# ==============================================================================

# Update these paths if needed
omni_female <- readRDS(
  "/Users/nirwantandukar/Documents/Github/MAGCAT/omni_Fly_female.RDS"
)
omni_male <- readRDS(
  "/Users/nirwantandukar/Documents/Github/MAGCAT/omni_Fly_male.RDS"
)

cat("Female pathways:", nrow(omni_female), "\n")
cat("Male pathways:", nrow(omni_male), "\n")

# Check column names match
cat("\nFemale columns:\n")
print(names(omni_female))

# ==============================================================================
# 2. Define P-value Columns
# ==============================================================================

# Check for calibrated columns
calibrated_cols <- c(
  "ACAT"     = "acat_p_mvn_cal",
  "Fisher"   = "fisher_p_mvn_cal",
  "TFisher"  = "tfisher_p_mvn_cal",
  "minP"     = "minp_p_mvn_cal",
  "Stouffer" = "stouffer_p_mvn_cal"
)

analytic_cols <- c(
  "ACAT"     = "acat_p",
  "Fisher"   = "fisher_p",
  "TFisher"  = "tfisher_p_analytic",
  "minP"     = "minp_p_analytic",
  "Stouffer" = "stouffer_p_analytic"
)

# Use calibrated if available
calibrated_available <- all(calibrated_cols %in% names(omni_female))

if (calibrated_available) {
  comp_cols <- calibrated_cols
  cat("\nUsing MVN-calibrated component p-values\n")
} else {
  comp_cols <- analytic_cols
  cat("\nUsing analytic component p-values\n")
}

# OMNIBUS column
omni_col <- "omni_p_final"
if (!omni_col %in% names(omni_female)) {
  omni_col <- "omni_p_mvn"
  if (!omni_col %in% names(omni_female)) {
    omni_col <- "omni_p_analytic"
  }
}
cat("Using OMNIBUS column:", omni_col, "\n")

# ==============================================================================
# 3. Find Common Pathways Between Male and Female
# ==============================================================================

common_pathways <- intersect(omni_female$pathway_id, omni_male$pathway_id)
cat("\nCommon pathways between male and female:", length(common_pathways), "\n")

# Filter to common pathways for fair comparison
omni_female_common <- omni_female[omni_female$pathway_id %in% common_pathways, ]
omni_male_common <- omni_male[omni_male$pathway_id %in% common_pathways, ]

# Align row order
omni_female_common <- omni_female_common[
  order(omni_female_common$pathway_id), ]
omni_male_common <- omni_male_common[
  order(omni_male_common$pathway_id), ]

# ==============================================================================
# 4. Determine Significance Threshold
# ==============================================================================

n_pathways <- length(common_pathways)
bonferroni_alpha <- 0.05 / n_pathways

cat("\nBonferroni threshold:", format(bonferroni_alpha, digits = 3), "\n")

# Set threshold based on user parameter SIG_MODE
if (SIG_MODE == "bonferroni") {
  sig_threshold <- bonferroni_alpha
  threshold_label <- paste0("Bonferroni (", format(bonferroni_alpha, digits = 3), ")")
  top_n_fallback <- NULL
} else if (SIG_MODE == "fdr") {
  sig_threshold <- FDR_THRESHOLD
  threshold_label <- paste0("FDR < ", FDR_THRESHOLD)
  top_n_fallback <- NULL
} else if (SIG_MODE == "nominal") {
  sig_threshold <- NOMINAL_ALPHA
  threshold_label <- paste0("Nominal (p < ", NOMINAL_ALPHA, ")")
  top_n_fallback <- NULL
} else if (SIG_MODE == "top_pct") {
  # Calculate top X% as number of pathways
  top_n_fallback <- ceiling(n_pathways * TOP_PCT_PATHWAYS / 100)
  sig_threshold <- NULL
  threshold_label <- paste0("Top ", TOP_PCT_PATHWAYS, "% (n=", top_n_fallback, ")")
} else {
 # Default: top_k
  top_n_fallback <- TOP_K_PATHWAYS
  sig_threshold <- NULL
  threshold_label <- paste0("Top ", TOP_K_PATHWAYS)
}

cat("Using threshold mode:", SIG_MODE, "\n")
cat("Threshold label:", threshold_label, "\n")

# ==============================================================================
# 5. Helper Function to Get Significant Pathways
# ==============================================================================

get_sig_pathways <- function(df, col, threshold = NULL, top_n = 20) {
  if (!col %in% names(df)) return(character(0))
  pvals <- df[[col]]
  pathway_ids <- df$pathway_id
  valid_idx <- !is.na(pvals)
  pvals <- pvals[valid_idx]
  pathway_ids <- pathway_ids[valid_idx]

  if (is.null(threshold)) {
    ord <- order(pvals)
    n_take <- min(top_n, length(pvals))
    return(pathway_ids[ord[seq_len(n_take)]])
  } else {
    return(pathway_ids[pvals < threshold])
  }
}

# ==============================================================================
# 6. Get Significant Pathways for Each Test x Sex
# ==============================================================================

# Use top_n_fallback if sig_threshold is NULL (top_k or top_pct mode)
top_n_to_use <- if (is.null(sig_threshold)) top_n_fallback else 20

# Component tests
comp_sig_female <- lapply(names(comp_cols), function(m) {
  get_sig_pathways(omni_female_common, comp_cols[m], sig_threshold, top_n_to_use)
})
names(comp_sig_female) <- names(comp_cols)

comp_sig_male <- lapply(names(comp_cols), function(m) {
  get_sig_pathways(omni_male_common, comp_cols[m], sig_threshold, top_n_to_use)
})
names(comp_sig_male) <- names(comp_cols)

# OMNIBUS
omni_sig_female <- get_sig_pathways(
  omni_female_common, omni_col, sig_threshold, top_n_to_use
)
omni_sig_male <- get_sig_pathways(
  omni_male_common, omni_col, sig_threshold, top_n_to_use
)

cat("\nSignificant pathways per method:\n")
cat("\n  Female:\n")
for (m in names(comp_cols)) {
  cat("    ", m, ":", length(comp_sig_female[[m]]), "\n")
}
cat("    OMNIBUS:", length(omni_sig_female), "\n")

cat("\n  Male:\n")
for (m in names(comp_cols)) {
  cat("    ", m, ":", length(comp_sig_male[[m]]), "\n")
}
cat("    OMNIBUS:", length(omni_sig_male), "\n")

# ==============================================================================
# 7. Jaccard Similarity: Male vs Female for Each Method
# ==============================================================================

jaccard <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(1)
  inter <- length(intersect(a, b))
  union_size <- length(union(a, b))
  if (union_size == 0) return(0)
  inter / union_size
}

# Calculate Jaccard for each component and OMNIBUS
jaccard_results <- data.frame(
  Method = c(names(comp_cols), "OMNIBUS"),
  Jaccard_Male_vs_Female = NA_real_,
  Overlap = NA_integer_,
  Union = NA_integer_,
  stringsAsFactors = FALSE
)

for (i in seq_along(comp_cols)) {
  m <- names(comp_cols)[i]
  j <- jaccard(comp_sig_female[[m]], comp_sig_male[[m]])
  jaccard_results$Jaccard_Male_vs_Female[i] <- j
  jaccard_results$Overlap[i] <- length(
    intersect(comp_sig_female[[m]], comp_sig_male[[m]])
  )
  jaccard_results$Union[i] <- length(
    union(comp_sig_female[[m]], comp_sig_male[[m]])
  )
}

# OMNIBUS
j_omni <- jaccard(omni_sig_female, omni_sig_male)
jaccard_results$Jaccard_Male_vs_Female[6] <- j_omni
jaccard_results$Overlap[6] <- length(intersect(omni_sig_female, omni_sig_male))
jaccard_results$Union[6] <- length(union(omni_sig_female, omni_sig_male))

cat("\n\n========================================\n")
cat("Male vs Female Agreement (Jaccard Similarity)\n")
cat("========================================\n\n")
print(jaccard_results)

# Bar plot of Jaccard similarity
jaccard_results$Method <- factor(
  jaccard_results$Method,
  levels = c("ACAT", "Fisher", "TFisher", "minP", "Stouffer", "OMNIBUS")
)

jaccard_plot <- ggplot(
  jaccard_results,
  aes(x = Method, y = Jaccard_Male_vs_Female, fill = Method)
) +
  geom_col(alpha = 0.8) +
  geom_text(
    aes(label = round(Jaccard_Male_vs_Female, 3)),
    vjust = -0.5, size = 4
  ) +
  scale_fill_manual(values = c(
    "ACAT" = "steelblue",
    "Fisher" = "steelblue",
    "TFisher" = "steelblue",
    "minP" = "steelblue",
    "Stouffer" = "steelblue",
    "OMNIBUS" = "darkred"
  )) +
  labs(
    title = "Male vs Female Agreement: Component Tests vs OMNIBUS",
    subtitle = paste0(
      "Jaccard similarity of significant pathways (", threshold_label, ")"
    ),
    x = "",
    y = "Jaccard Similarity"
  ) +
  ylim(0, max(jaccard_results$Jaccard_Male_vs_Female, na.rm = TRUE) * 1.15) +
  theme_minimal() +
  plot_theme +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    legend.position = "none"
  )

quartz()
print(jaccard_plot)
ggsave("fly_male_female_jaccard.png", jaccard_plot,
       width = 8, height = 8, dpi = 300, bg = "white")

# ==============================================================================
# 8. Venn Diagrams: Male vs Female for Each Method
# ==============================================================================

venn_plots <- list()

# Component tests
for (m in names(comp_cols)) {
  venn_list <- list(
    "Female" = comp_sig_female[[m]],
    "Male" = comp_sig_male[[m]]
  )

  p <- ggVennDiagram(venn_list, label = "count", label_alpha = 0) +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(title = m) +
    theme(plot.title = element_text(hjust = 0.5))

  venn_plots[[m]] <- p
}

# OMNIBUS
venn_omni <- list(
  "Female" = omni_sig_female,
  "Male" = omni_sig_male
)

p_omni <- ggVennDiagram(venn_omni, label = "count", label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "darkred") +
  labs(title = "OMNIBUS") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

venn_plots[["OMNIBUS"]] <- p_omni

# Combine all venns
combined_venn <- (
  venn_plots[["ACAT"]] + venn_plots[["Fisher"]] + venn_plots[["TFisher"]]
) / (
  venn_plots[["minP"]] + venn_plots[["Stouffer"]] + venn_plots[["OMNIBUS"]]
) +
  plot_annotation(
    title = "Male vs Female: Component Tests and OMNIBUS",
    subtitle = paste0("Overlap of significant pathways (", threshold_label, ")"),
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11)
    )
  )

quartz()
print(combined_venn)
ggsave("fly_male_female_venn_all.png", combined_venn,
       width = 14, height = 10, dpi = 300, bg = "white")

# ==============================================================================
# 9. Spearman Rank Correlation: Male vs Female
# ==============================================================================

spearman_results <- data.frame(
  Method = c(names(comp_cols), "OMNIBUS"),
  Spearman_Cor = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(comp_cols)) {
  m <- names(comp_cols)[i]
  col <- comp_cols[m]
  if (col %in% names(omni_female_common) && col %in% names(omni_male_common)) {
    cor_val <- cor(
      omni_female_common[[col]],
      omni_male_common[[col]],
      use = "pairwise.complete.obs",
      method = "spearman"
    )
    spearman_results$Spearman_Cor[i] <- cor_val
  }
}

# OMNIBUS
cor_omni <- cor(
  omni_female_common[[omni_col]],
  omni_male_common[[omni_col]],
  use = "pairwise.complete.obs",
  method = "spearman"
)
spearman_results$Spearman_Cor[6] <- cor_omni

cat("\n\n========================================\n")
cat("Spearman Rank Correlation (Male vs Female)\n")
cat("========================================\n\n")
print(spearman_results)

# Bar plot
spearman_results$Method <- factor(
  spearman_results$Method,
  levels = c("ACAT", "Fisher", "TFisher", "minP", "Stouffer", "OMNIBUS")
)

spearman_plot <- ggplot(
  spearman_results,
  aes(x = Method, y = Spearman_Cor, fill = Method)
) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = round(Spearman_Cor, 3)), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c(
    "ACAT" = "steelblue",
    "Fisher" = "steelblue",
    "TFisher" = "steelblue",
    "minP" = "steelblue",
    "Stouffer" = "steelblue",
    "OMNIBUS" = "darkred"
  )) +
  labs(
    title = "Male vs Female Rank Correlation: Component Tests vs OMNIBUS",
    subtitle = "Spearman correlation of pathway rankings (all pathways)",
    x = "",
    y = "Spearman Correlation"
  ) +
  ylim(0, 1) +
  theme_minimal() +
  plot_theme +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    legend.position = "none"
  )

quartz()
print(spearman_plot)
ggsave("fly_male_female_spearman.png", spearman_plot,
       width = 8, height = 8, dpi = 300, bg = "white")

# ==============================================================================
# 10. Combined Metrics Plot
# ==============================================================================

metrics_df <- data.frame(
  Method = rep(
    c("ACAT", "Fisher", "TFisher", "minP", "Stouffer", "OMNIBUS"),
    2
  ),
  Metric = rep(c("Jaccard", "Spearman"), each = 6),
  Value = c(
    jaccard_results$Jaccard_Male_vs_Female,
    spearman_results$Spearman_Cor
  )
)

metrics_df$Method <- factor(
  metrics_df$Method,
  levels = c("ACAT", "Fisher", "TFisher", "minP", "Stouffer", "OMNIBUS")
)

metrics_plot <- ggplot(
  metrics_df,
  aes(x = Method, y = Value, fill = Method)
) +
  geom_col(alpha = 0.8, position = "dodge") +
  geom_text(
    aes(label = round(Value, 2)),
    vjust = -0.5, size = 3, position = position_dodge(width = 0.9)
  ) +
  facet_wrap(~Metric, scales = "free_y") +
  scale_fill_manual(values = c(
    "ACAT" = "steelblue",
    "Fisher" = "steelblue",
    "TFisher" = "steelblue",
    "minP" = "steelblue",
    "Stouffer" = "steelblue",
    "OMNIBUS" = "darkred"
  )) +
  labs(
    title = "Male vs Female Agreement: Component Tests vs OMNIBUS",
    subtitle = "Higher values = more agreement between sexes",
    x = "",
    y = "Value"
  ) +
  theme_minimal() +
  plot_theme +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    legend.position = "none"
  )

quartz()
print(metrics_plot)
ggsave("fly_male_female_metrics.png", metrics_plot,
       width = 12, height = 6, dpi = 300, bg = "white")

# ==============================================================================
# 11. Consensus Pathways: Found by OMNIBUS in BOTH Male and Female
# ==============================================================================

omni_consensus <- intersect(omni_sig_female, omni_sig_male)

cat("\n\n========================================\n")
cat("OMNIBUS Consensus Pathways (Both Sexes)\n")
cat("========================================\n")
cat("Pathways significant in BOTH male and female by OMNIBUS:\n")
cat("n =", length(omni_consensus), "\n\n")

if (length(omni_consensus) > 0) {
  # Get row indices for consensus pathways
 female_idx <- match(omni_consensus, omni_female_common$pathway_id)
  male_idx <- match(omni_consensus, omni_male_common$pathway_id)

  # Build comprehensive consensus dataframe
  consensus_df <- data.frame(
    pathway_id = omni_consensus,
    pathway_name = omni_female_common$pathway_name[female_idx],
    stringsAsFactors = FALSE
  )

  # Add gene counts if available
  if ("n_genes" %in% names(omni_female_common)) {
    consensus_df$n_genes_female <- omni_female_common$n_genes[female_idx]
    consensus_df$n_genes_male <- omni_male_common$n_genes[male_idx]
  }

  # Add gene names if available
  if ("gene_names" %in% names(omni_female_common)) {
    consensus_df$genes_female <- omni_female_common$gene_names[female_idx]
    consensus_df$genes_male <- omni_male_common$gene_names[male_idx]
  }

  # Add OMNIBUS p-values
  consensus_df$omni_p_female <- omni_female_common[[omni_col]][female_idx]
  consensus_df$omni_p_male <- omni_male_common[[omni_col]][male_idx]

  # Add component test p-values for both sexes
  for (m in names(comp_cols)) {
    col <- comp_cols[m]
    if (col %in% names(omni_female_common)) {
      consensus_df[[paste0(m, "_p_female")]] <-
        omni_female_common[[col]][female_idx]
      consensus_df[[paste0(m, "_p_male")]] <-
        omni_male_common[[col]][male_idx]
    }
  }

  # Sort by average OMNIBUS p-value
  consensus_df$avg_omni_p <- (consensus_df$omni_p_female +
                               consensus_df$omni_p_male) / 2
  consensus_df <- consensus_df[order(consensus_df$avg_omni_p), ]
  consensus_df$avg_omni_p <- NULL  # Remove helper column

  # Print summary
  for (i in seq_len(nrow(consensus_df))) {
    cat(
      i, ". ", consensus_df$pathway_id[i], "\n",
      "   ", consensus_df$pathway_name[i], "\n",
      "   Female p: ", format(consensus_df$omni_p_female[i], digits = 3),
      " | Male p: ", format(consensus_df$omni_p_male[i], digits = 3), "\n",
      sep = ""
    )
    if ("n_genes_female" %in% names(consensus_df)) {
      cat("   Genes (F/M): ", consensus_df$n_genes_female[i], "/",
          consensus_df$n_genes_male[i], "\n", sep = "")
    }
    cat("\n")
  }

  write.csv(consensus_df, "fly_omnibus_consensus_pathways.csv",
            row.names = FALSE)
  cat("Saved to fly_omnibus_consensus_pathways.csv\n")

  # ---------- Create figure for consensus pathways ----------

  # Prepare data for plotting - OMNIBUS comparison
  consensus_plot_df <- consensus_df %>%
    mutate(
      pathway_label = ifelse(
        nchar(pathway_name) > 40,
        paste0(substr(pathway_name, 1, 37), "..."),
        pathway_name
      )
    ) %>%
    select(pathway_id, pathway_label, omni_p_female, omni_p_male) %>%
    pivot_longer(
      cols = c(omni_p_female, omni_p_male),
      names_to = "Sex",
      values_to = "pvalue"
    ) %>%
    mutate(
      Sex = ifelse(Sex == "omni_p_female", "Female", "Male"),
      neg_log10_p = -log10(pvalue)
    )

  # Order pathways by average -log10(p)
  pathway_order <- consensus_plot_df %>%
    group_by(pathway_label) %>%
    summarize(avg_p = mean(neg_log10_p), .groups = "drop") %>%
    arrange(desc(avg_p)) %>%
    pull(pathway_label)

  consensus_plot_df$pathway_label <- factor(
    consensus_plot_df$pathway_label,
    levels = rev(pathway_order)
  )

  # Dot plot comparing male vs female OMNIBUS p-values
  consensus_dot_plot <- ggplot(
    consensus_plot_df,
    aes(x = neg_log10_p, y = pathway_label, color = Sex, shape = Sex)
  ) +
    geom_point(size = 4, alpha = 0.8) +
    geom_line(
      aes(group = pathway_label),
      color = "gray60", linewidth = 0.5, alpha = 0.5
    ) +
    geom_vline(
      xintercept = -log10(0.05),
      linetype = "dashed", color = "gray50"
    ) +
    scale_color_manual(values = c("Female" = "hotpink", "Male" = "steelblue")) +
    scale_shape_manual(values = c("Female" = 16, "Male" = 17)) +
    labs(
      title = "Consensus Pathways: Male vs Female OMNIBUS P-values",
      subtitle = paste0(
        "Pathways significant in both sexes (", threshold_label, ")"
      ),
      x = expression(-log[10](OMNIBUS~p)),
      y = ""
    ) +
    theme_minimal() +
    plot_theme +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(size = 9)
    )

  quartz()
  print(consensus_dot_plot)
  ggsave("fly_consensus_pathways_comparison.png", consensus_dot_plot,
         width = 10, height = max(6, nrow(consensus_df) * 0.4),
         dpi = 300, bg = "white")

  # ---------- Heatmap of all methods for consensus pathways ----------

  # Reshape for heatmap - include all component tests
  consensus_heatmap_df <- consensus_df %>%
    mutate(
      pathway_label = ifelse(
        nchar(pathway_name) > 35,
        paste0(substr(pathway_name, 1, 32), "..."),
        pathway_name
      )
    ) %>%
    select(pathway_id, pathway_label,
           omni_p_female, omni_p_male,
           starts_with("ACAT_"), starts_with("Fisher_"),
           starts_with("TFisher_"), starts_with("minP_"),
           starts_with("Stouffer_")) %>%
    pivot_longer(
      cols = -c(pathway_id, pathway_label),
      names_to = "Method_Sex",
      values_to = "pvalue"
    ) %>%
    mutate(
      Sex = ifelse(grepl("female", Method_Sex), "Female", "Male"),
      Method = gsub("_p_(female|male)", "", Method_Sex),
      Method = gsub("omni", "OMNIBUS", Method),
      neg_log10_p = -log10(pvalue)
    )

  # Order pathways
  consensus_heatmap_df$pathway_label <- factor(
    consensus_heatmap_df$pathway_label,
    levels = rev(pathway_order)
  )

  consensus_heatmap_df$Method <- factor(
    consensus_heatmap_df$Method,
    levels = c("OMNIBUS", "ACAT", "Fisher", "TFisher", "minP", "Stouffer")
  )

  consensus_heatmap <- ggplot(
    consensus_heatmap_df,
    aes(x = interaction(Method, Sex, sep = "\n"),
        y = pathway_label, fill = neg_log10_p)
  ) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(
      low = "white", mid = "steelblue", high = "darkred",
      midpoint = -log10(0.05),
      name = expression(-log[10](p))
    ) +
    labs(
      title = "Consensus Pathways: All Methods by Sex",
      subtitle = paste0(
        "Pathways significant in both sexes (", threshold_label, ")"
      ),
      x = "",
      y = ""
    ) +
    theme_minimal() +
    plot_theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 9)
    )

  quartz()
  print(consensus_heatmap)
  ggsave("fly_consensus_pathways_heatmap.png", consensus_heatmap,
         width = 12, height = max(6, nrow(consensus_df) * 0.5),
         dpi = 300, bg = "white")
}

# ==============================================================================
# 12. Sex-Specific Pathways (OMNIBUS)
# ==============================================================================

omni_female_only <- setdiff(omni_sig_female, omni_sig_male)
omni_male_only <- setdiff(omni_sig_male, omni_sig_female)

cat("\n\n========================================\n")
cat("Sex-Specific OMNIBUS Pathways\n")
cat("========================================\n")

cat("\nFemale-specific (n =", length(omni_female_only), "):\n")
if (length(omni_female_only) > 0) {
  f_only_df <- omni_female_common[
    omni_female_common$pathway_id %in% omni_female_only,
    c("pathway_id", "pathway_name")
  ]
  for (i in seq_len(nrow(f_only_df))) {
    cat("  - ", f_only_df$pathway_id[i], ": ",
        f_only_df$pathway_name[i], "\n", sep = "")
  }
}

cat("\nMale-specific (n =", length(omni_male_only), "):\n")
if (length(omni_male_only) > 0) {
  m_only_df <- omni_male_common[
    omni_male_common$pathway_id %in% omni_male_only,
    c("pathway_id", "pathway_name")
  ]
  for (i in seq_len(nrow(m_only_df))) {
    cat("  - ", m_only_df$pathway_id[i], ": ",
        m_only_df$pathway_name[i], "\n", sep = "")
  }
}

# ==============================================================================
# 13. Heatmap: Compare Top Pathways Across Both Sexes
# ==============================================================================

# Get top pathways from either sex (union of top hits)
top_union <- union(
  head(omni_female_common$pathway_id[
    order(omni_female_common[[omni_col]])], 15),
  head(omni_male_common$pathway_id[
    order(omni_male_common[[omni_col]])], 15)
)

# Build heatmap data
heatmap_data <- data.frame(
  pathway_id = top_union,
  pathway_name = omni_female_common$pathway_name[
    match(top_union, omni_female_common$pathway_id)
  ]
)

# Add OMNIBUS p-values
heatmap_data$Female_OMNIBUS <- omni_female_common[[omni_col]][
  match(top_union, omni_female_common$pathway_id)
]
heatmap_data$Male_OMNIBUS <- omni_male_common[[omni_col]][
  match(top_union, omni_male_common$pathway_id)
]

# Add component p-values
for (m in names(comp_cols)) {
  col <- comp_cols[m]
  heatmap_data[[paste0("Female_", m)]] <- omni_female_common[[col]][
    match(top_union, omni_female_common$pathway_id)
  ]
  heatmap_data[[paste0("Male_", m)]] <- omni_male_common[[col]][
    match(top_union, omni_male_common$pathway_id)
  ]
}

# Reshape for plotting
heatmap_long <- heatmap_data %>%
  pivot_longer(
    cols = -c(pathway_id, pathway_name),
    names_to = "Method_Sex",
    values_to = "pvalue"
  ) %>%
  mutate(
    Sex = ifelse(grepl("^Female", Method_Sex), "Female", "Male"),
    Method = gsub("^(Female|Male)_", "", Method_Sex),
    neg_log10_p = -log10(pvalue),
    pathway_label = ifelse(
      nchar(pathway_name) > 30,
      paste0(substr(pathway_name, 1, 27), "..."),
      pathway_name
    )
  )

# Order pathways by average OMNIBUS p
pathway_order <- heatmap_data %>%
  mutate(avg_omni = (Female_OMNIBUS + Male_OMNIBUS) / 2) %>%
  arrange(avg_omni) %>%
  pull(pathway_id)

heatmap_long$pathway_id <- factor(
  heatmap_long$pathway_id,
  levels = rev(pathway_order)
)

heatmap_long$Method <- factor(
  heatmap_long$Method,
  levels = c("OMNIBUS", "ACAT", "Fisher", "TFisher", "minP", "Stouffer")
)

# Create heatmap
heatmap_plot <- ggplot(
  heatmap_long,
  aes(x = interaction(Method, Sex), y = pathway_id, fill = neg_log10_p)
) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "white", mid = "steelblue", high = "darkred",
    midpoint = -log10(0.05),
    name = expression(-log[10](p))
  ) +
  scale_x_discrete(
    labels = function(x) {
      parts <- strsplit(x, "\\.")
      sapply(parts, function(p) paste0(p[1], "\n(", p[2], ")"))
    }
  ) +
  labs(
    title = "P-value Comparison: Male vs Female Across All Methods",
    subtitle = paste0("Top pathways from either sex (", threshold_label, ")"),
    x = "Method (Sex)",
    y = "Pathway"
  ) +
  theme_minimal() +
  plot_theme +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8),
    axis.text.y = element_text(size = 8)
  )

quartz()
print(heatmap_plot)
ggsave("fly_male_female_heatmap.png", heatmap_plot,
       width = 16, height = 10, dpi = 300, bg = "white")

# ==============================================================================
# 14. Summary Statistics
# ==============================================================================

cat("\n\n")
cat(strrep("=", 60), "\n")
cat("SUMMARY: Male vs Female Comparison\n")
cat(strrep("=", 60), "\n\n")

cat("Threshold used:", threshold_label, "\n\n")

cat("Jaccard Similarity (Male vs Female):\n")
cat("  Component tests average:",
    round(mean(jaccard_results$Jaccard_Male_vs_Female[1:5], na.rm = TRUE), 3),
    "\n")
cat("  OMNIBUS:",
    round(jaccard_results$Jaccard_Male_vs_Female[6], 3), "\n")

omni_vs_comp_jaccard <- jaccard_results$Jaccard_Male_vs_Female[6] -
  mean(jaccard_results$Jaccard_Male_vs_Female[1:5], na.rm = TRUE)
cat("  OMNIBUS improvement: +",
    round(omni_vs_comp_jaccard, 3), "\n")

cat("\nSpearman Rank Correlation (Male vs Female):\n")
cat("  Component tests average:",
    round(mean(spearman_results$Spearman_Cor[1:5], na.rm = TRUE), 3), "\n")
cat("  OMNIBUS:",
    round(spearman_results$Spearman_Cor[6], 3), "\n")

omni_vs_comp_spearman <- spearman_results$Spearman_Cor[6] -
  mean(spearman_results$Spearman_Cor[1:5], na.rm = TRUE)
cat("  OMNIBUS improvement: +",
    round(omni_vs_comp_spearman, 3), "\n")

cat("\nOverlap counts:\n")
cat("  OMNIBUS consensus (both sexes):", length(omni_consensus), "\n")
cat("  Female-specific:", length(omni_female_only), "\n")
cat("  Male-specific:", length(omni_male_only), "\n")

cat("\n** Key Insight **\n")
if (omni_vs_comp_jaccard > 0) {
  cat("OMNIBUS shows HIGHER agreement between male and female\n")
  cat("than individual component tests, demonstrating that the\n")
  cat("holistic combination identifies more robust biological signals\n")
  cat("that transcend sex-specific noise in individual methods.\n")
} else {
  cat("Component tests and OMNIBUS show similar male-female agreement.\n")
  cat("The biological signal may be sex-specific for this phenotype.\n")
}

cat("\n========== DONE ==========\n")
cat("Output files:\n")
cat("  - fly_male_female_jaccard.png\n")
cat("  - fly_male_female_venn_all.png\n")
cat("  - fly_male_female_spearman.png\n")
cat("  - fly_male_female_metrics.png\n")
cat("  - fly_male_female_heatmap.png\n")
cat("  - fly_omnibus_consensus_pathways.csv\n")
cat("  - fly_consensus_pathways_comparison.png\n")
cat("  - fly_consensus_pathways_heatmap.png\n")
