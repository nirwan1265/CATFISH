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
nominal_alpha <- 0.1

cat("\nBonferroni threshold:", format(bonferroni_alpha, digits = 3), "\n")

# Check significance counts for OMNIBUS
omni_f_bonf <- sum(omni_female_common[[omni_col]] < bonferroni_alpha,
                   na.rm = TRUE)
omni_m_bonf <- sum(omni_male_common[[omni_col]] < bonferroni_alpha,
                   na.rm = TRUE)
omni_f_nom <- sum(omni_female_common[[omni_col]] < nominal_alpha, na.rm = TRUE)
omni_m_nom <- sum(omni_male_common[[omni_col]] < nominal_alpha, na.rm = TRUE)

cat("\nOMNIBUS significant (Bonferroni):\n")
cat("  Female:", omni_f_bonf, "\n")
cat("  Male:", omni_m_bonf, "\n")

cat("\nOMNIBUS significant (nominal 0.1):\n")
cat("  Female:", omni_f_nom, "\n")
cat("  Male:", omni_m_nom, "\n")

# Decide threshold
if (omni_f_bonf > 0 || omni_m_bonf > 0) {
  sig_threshold <- bonferroni_alpha
  threshold_label <- "Bonferroni"
} else if (omni_f_nom >= 5 || omni_m_nom >= 5) {
  sig_threshold <- nominal_alpha
  threshold_label <- "Nominal (0.1)"
} else {
  sig_threshold <- NULL
  threshold_label <- "Top 20"
}

cat("Using threshold:", threshold_label, "\n")

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

# Component tests
comp_sig_female <- lapply(names(comp_cols), function(m) {
  get_sig_pathways(omni_female_common, comp_cols[m], sig_threshold, 20)
})
names(comp_sig_female) <- names(comp_cols)

comp_sig_male <- lapply(names(comp_cols), function(m) {
  get_sig_pathways(omni_male_common, comp_cols[m], sig_threshold, 20)
})
names(comp_sig_male) <- names(comp_cols)

# OMNIBUS
omni_sig_female <- get_sig_pathways(
  omni_female_common, omni_col, sig_threshold, 20
)
omni_sig_male <- get_sig_pathways(
  omni_male_common, omni_col, sig_threshold, 20
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
  aes(x = Method, y = Jaccard_Male_vs_Female,
      fill = Method == "OMNIBUS")
) +
  geom_col(alpha = 0.8) +
  geom_text(
    aes(label = round(Jaccard_Male_vs_Female, 3)),
    vjust = -0.5, size = 4
  ) +
  scale_fill_manual(
    values = c("TRUE" = "darkred", "FALSE" = "steelblue"),
    guide = "none"
  ) +
  labs(
    title = "Male vs Female Agreement: Component Tests vs OMNIBUS",
    subtitle = paste0(
      "Jaccard similarity of significant pathways (", threshold_label, ")"
    ),
    x = "Method",
    y = "Jaccard Similarity (Male vs Female)"
  ) +
  ylim(0, max(jaccard_results$Jaccard_Male_vs_Female, na.rm = TRUE) * 1.15) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

quartz()
print(jaccard_plot)
ggsave("fly_male_female_jaccard.pdf", jaccard_plot, width = 10, height = 6)
ggsave("fly_male_female_jaccard.png", jaccard_plot, width = 10, height = 6,
       dpi = 300)

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
    subtitle = paste0(
      "Overlap of significant pathways (", threshold_label, ")"
    ),
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5)
    )
  )

quartz()
print(combined_venn)
ggsave("fly_male_female_venn_all.pdf", combined_venn, width = 14, height = 10)
ggsave("fly_male_female_venn_all.png", combined_venn, width = 14, height = 10,
       dpi = 300)

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
  aes(x = Method, y = Spearman_Cor, fill = Method == "OMNIBUS")
) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = round(Spearman_Cor, 3)), vjust = -0.5, size = 4) +
  scale_fill_manual(
    values = c("TRUE" = "darkred", "FALSE" = "steelblue"),
    guide = "none"
  ) +
  labs(
    title = "Male vs Female Rank Correlation: Component Tests vs OMNIBUS",
    subtitle = "Spearman correlation of pathway rankings (all pathways)",
    x = "Method",
    y = "Spearman Correlation"
  ) +
  ylim(0, 1) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

quartz()
print(spearman_plot)
ggsave("fly_male_female_spearman.pdf", spearman_plot, width = 10, height = 6)
ggsave("fly_male_female_spearman.png", spearman_plot, width = 10, height = 6,
       dpi = 300)

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
  aes(x = Method, y = Value, fill = Method == "OMNIBUS")
) +
  geom_col(alpha = 0.8, position = "dodge") +
  geom_text(
    aes(label = round(Value, 2)),
    vjust = -0.5, size = 3, position = position_dodge(width = 0.9)
  ) +
  facet_wrap(~Metric, scales = "free_y") +
  scale_fill_manual(
    values = c("TRUE" = "darkred", "FALSE" = "steelblue"),
    name = "",
    labels = c("Component Test", "OMNIBUS")
  ) +
  labs(
    title = "Male vs Female Agreement: Component Tests vs OMNIBUS",
    subtitle = "Higher values = more agreement between sexes",
    x = "Method",
    y = "Value"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

quartz()
print(metrics_plot)
ggsave("fly_male_female_metrics.pdf", metrics_plot, width = 12, height = 6)
ggsave("fly_male_female_metrics.png", metrics_plot, width = 12, height = 6,
       dpi = 300)

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
  consensus_df <- omni_female_common[
    omni_female_common$pathway_id %in% omni_consensus,
    c("pathway_id", "pathway_name")
  ]

  # Add p-values from both
  consensus_df$omni_p_female <- omni_female_common[[omni_col]][
    match(consensus_df$pathway_id, omni_female_common$pathway_id)
  ]
  consensus_df$omni_p_male <- omni_male_common[[omni_col]][
    match(consensus_df$pathway_id, omni_male_common$pathway_id)
  ]

  consensus_df <- consensus_df[order(consensus_df$omni_p_female), ]

  for (i in seq_len(nrow(consensus_df))) {
    cat(
      i, ". ", consensus_df$pathway_id[i], "\n",
      "   ", consensus_df$pathway_name[i], "\n",
      "   Female p: ", format(consensus_df$omni_p_female[i], digits = 3),
      " | Male p: ", format(consensus_df$omni_p_male[i], digits = 3), "\n\n",
      sep = ""
    )
  }

  write.csv(consensus_df, "fly_omnibus_consensus_pathways.csv",
            row.names = FALSE)
  cat("Saved to fly_omnibus_consensus_pathways.csv\n")
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
    subtitle = "Top pathways from either sex",
    x = "Method (Sex)",
    y = "Pathway"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 7),
    axis.text.y = element_text(size = 7),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

quartz()
print(heatmap_plot)
ggsave("fly_male_female_heatmap.pdf", heatmap_plot, width = 16, height = 10)
ggsave("fly_male_female_heatmap.png", heatmap_plot, width = 16, height = 10,
       dpi = 300)

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
cat("  - fly_male_female_jaccard.pdf/png\n")
cat("  - fly_male_female_venn_all.pdf/png\n")
cat("  - fly_male_female_spearman.pdf/png\n")
cat("  - fly_male_female_metrics.pdf/png\n")
cat("  - fly_male_female_heatmap.pdf/png\n")
cat("  - fly_omnibus_consensus_pathways.csv\n")
