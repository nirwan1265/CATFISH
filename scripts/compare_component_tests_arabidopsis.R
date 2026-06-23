# ==============================================================================
# Compare Component Tests from CATFISH Omnibus Results
# Shows that ACAT, Fisher, TFisher, minP, and Stouffer give different results
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

# For Venn diagrams
if (!requireNamespace("ggVennDiagram", quietly = TRUE)) {
  install.packages("ggVennDiagram")
}
library(ggVennDiagram)

# For upset plots (alternative to Venn for >3 sets)
if (!requireNamespace("UpSetR", quietly = TRUE)) {
  install.packages("UpSetR")
}
library(UpSetR)

# ==============================================================================
# 1. Plot Theme
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
# 1. Load Results
# ==============================================================================

omni_results <- readRDS("/Users/nirwantandukar/Documents/Github/MAGCAT/omni_results_arabidopsis_ACAT_B10000.RDS")

head(omni_results)
str(omni_results)

# Check available columns
names(omni_results)

# ==============================================================================
# 2. Define Component P-value Columns
# ==============================================================================

# Analytic (uncalibrated) component p-values
analytic_cols <- c(
  "ACAT"     = "acat_p",
  "Fisher"   = "fisher_p",
  "TFisher"  = "tfisher_p_analytic",
  "minP"     = "minp_p_analytic",
  "Stouffer" = "stouffer_p_analytic"
)

# MVN-calibrated component p-values (if available)
calibrated_cols <- c(
  "ACAT"     = "acat_p_mvn_cal",
  "Fisher"   = "fisher_p_mvn_cal",
  "TFisher"  = "tfisher_p_mvn_cal",
  "minP"     = "minp_p_mvn_cal",
  "Stouffer" = "stouffer_p_mvn_cal"
)

# Check which calibrated columns exist
calibrated_available <- all(calibrated_cols %in% names(omni_results))
cat("MVN-calibrated component p-values available:", calibrated_available, "\n")

# Choose which p-values to use for comparison
# Prefer calibrated if available
if (calibrated_available) {
  use_cols <- calibrated_cols
  cat("Using MVN-calibrated component p-values for comparison\n")
} else {
  use_cols <- analytic_cols
  cat("Using analytic (uncalibrated) component p-values for comparison\n")
}

# ==============================================================================
# 3. Determine Significance Threshold
# ==============================================================================

# Bonferroni correction
n_pathways <- nrow(omni_results)
bonferroni_alpha <- 0.005 / n_pathways
cat("Number of pathways:", n_pathways, "\n")
cat("Bonferroni threshold (alpha = 0.05):", bonferroni_alpha, "\n")

# Nominal threshold
nominal_alpha <- 0.0005

# Check how many significant under each method at Bonferroni level
bonf_sig_counts <- sapply(use_cols, function(col) {
  if (col %in% names(omni_results)) {
    sum(omni_results[[col]] < bonferroni_alpha, na.rm = TRUE)
  } else {
    NA
  }
})
cat("\nSignificant pathways at Bonferroni level:\n")
print(bonf_sig_counts)

# Check at nominal level
nom_sig_counts <- sapply(use_cols, function(col) {
  if (col %in% names(omni_results)) {
    sum(omni_results[[col]] < nominal_alpha, na.rm = TRUE)
  } else {
    NA
  }
})
cat("\nSignificant pathways at nominal (0.05) level:\n")
print(nom_sig_counts)

# Decide which threshold to use
# If any method has Bonferroni-significant results, use Bonferroni
# Otherwise, use top 20 pathways per method
use_bonferroni <- any(bonf_sig_counts > 0, na.rm = TRUE)

if (use_bonferroni) {
  sig_threshold <- bonferroni_alpha
  threshold_label <- "Bonferroni"
  cat("\nUsing Bonferroni threshold for significance\n")
} else {
  cat("\nNo Bonferroni-significant pathways found.\n")
  cat("Will use top 20 pathways per method for comparison.\n")
  sig_threshold <- NULL
  threshold_label <- "Top 20"
}

# ==============================================================================
# 4. Create Significant Pathway Sets
# ==============================================================================

get_significant_pathways <- function(df, col, threshold = NULL, top_n = 20) {
  if (!col %in% names(df)) return(character(0))

  pvals <- df[[col]]
  pathway_ids <- df$pathway_id

  # Remove NAs
  valid_idx <- !is.na(pvals)
  pvals <- pvals[valid_idx]
  pathway_ids <- pathway_ids[valid_idx]

  if (is.null(threshold)) {
    # Use top N
    if (length(pvals) == 0) return(character(0))
    ord <- order(pvals)
    n_take <- min(top_n, length(pvals))
    return(pathway_ids[ord[1:n_take]])
  } else {
    # Use threshold
    return(pathway_ids[pvals < threshold])
  }
}

# Get significant/top pathways for each method
sig_pathways <- lapply(names(use_cols), function(method) {
  col <- use_cols[method]
  get_significant_pathways(omni_results, col, threshold = sig_threshold, top_n = 20)
})
names(sig_pathways) <- names(use_cols)

# Print counts
cat("\nPathways selected per method:\n")
sapply(sig_pathways, length)

# ==============================================================================
# 5. Venn Diagram Comparison
# ==============================================================================

# ggVennDiagram can handle up to 7 sets, but best with 2-5
# We have 5 methods, which works

# Create Venn diagram
## IMPORTANT: ggVennDiagram/ggplot text "size" is NOT in pt.
## It's in mm-ish ggplot units. To get true “24 pt”, use: 24 / ggplot2::.pt

pt <- ggplot2::.pt  # conversion helper

venn_plot <- ggVennDiagram(
  sig_pathways,
  label = "count",
  label_alpha = 0,
  edge_size = 0.6,
  set_label_size = 24 / pt,  # <-- ACAT / Fisher / TFisher / Stouffer / minP (24 pt)
  label_size     = 20 / pt   # <-- region counts (tweak if you want)
) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    title = "Overlap of Significant Pathways Across Component Tests",
    subtitle = paste0(
      "Threshold: ", threshold_label,
      if (!is.null(sig_threshold)) {
        paste0(" (p < ", format(sig_threshold, digits = 3), ")")
      } else {
        " per method"
      }
    )
  ) +
  plot_theme +
  theme(
    ## this figure: no grid lines (even though plot_theme has them)
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),

    ## no axes (so no weird X/Y leftovers)
    axis.title   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    axis.line    = element_blank(),

    ## keep legend + make it readable
    legend.position  = "top",
    legend.text      = element_text(size = 24),
    legend.title     = element_blank(),
    legend.key.width = unit(4.5, "cm"),   # make colorbar longer
    legend.key.height= unit(0.6, "cm")    # and taller
  ) +
  guides(
    fill = guide_colourbar(
      barwidth  = unit(6, "cm"),
      barheight = unit(0.6, "cm"),
      label.position = "bottom",
      title.position = "top"
    )
  )

## When saving: force a white background (your PNG looks like it saved transparent)
ggsave(
  "venn_component_tests_arabidopsis.png",
  venn_plot,
  width = 12, height = 12, dpi = 300,
  bg = "white"
)

tune_ggvenn <- function(p, set_pt = 24, count_pt = 18, legend_pt = 14) {
  pt <- ggplot2::.pt

  # --- A) figure-specific legend override (your plot_theme has legend.text=52)
  p <- p + theme(
    legend.position  = "top",
    legend.text      = element_text(size = legend_pt),
    legend.title     = element_blank(),
    legend.key.width = unit(2.5, "cm"),
    legend.key.height= unit(0.4, "cm")
  )

  # --- B) force-set sizes on the actual ggVennDiagram text layers
  for (i in seq_along(p$layers)) {
    lyr <- p$layers[[i]]
    if (inherits(lyr$geom, "GeomText") || inherits(lyr$geom, "GeomLabel")) {

      # If this layer labels set names (ACAT/Fisher/...)
      if (!is.null(lyr$mapping$label) &&
          deparse(lyr$mapping$label) %in% c("name", "set", "sets")) {
        p$layers[[i]]$aes_params$size <- set_pt / pt
      }

      # If this layer labels region counts
      if (!is.null(lyr$mapping$label) &&
          deparse(lyr$mapping$label) %in% c("count", "Counts", "n")) {
        p$layers[[i]]$aes_params$size <- count_pt / pt
      }
    }
  }

  p
}

## Build your plot EXACTLY how you already do it (NO ggsave here)
venn_plot <- ggVennDiagram(
  sig_pathways,
  label = "count",
  label_alpha = 0,
  edge_size = 0.6
) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    title = "Overlap of Significant Pathways Across Component Tests",
    subtitle = paste0("Threshold: ", threshold_label, " per method")
  ) +
  plot_theme +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.line  = element_blank()
  )

## Now force the sizes you want:
venn_plot <- tune_ggvenn(
  venn_plot,
  set_pt    = 24,  # ACAT/minP/etc
  count_pt  = 20,  # counts inside regions
  legend_pt = 14   # legend text (make this smaller/bigger)
)
venn_plot <- venn_plot +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.title   = element_blank(),  # extra safety

    axis.text.x  = element_blank(),
    axis.text.y  = element_blank(),
    axis.text    = element_blank(),  # extra safety

    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks   = element_blank(),  # extra safety

    axis.line    = element_blank()
  ) +
  labs(x = NULL, y = NULL) 

quartz()
print(venn_plot)

#ggsave("venn_component_tests_arabidopsis.pdf", venn_plot, width = 10, height = 8)
ggsave("venn_component_tests_arabidopsis.png", venn_plot, width = 8, height = 8, dpi = 300)


cat("\nVenn diagram saved to venn_component_tests_arabidopsis.pdf/png\n")

# ==============================================================================
# 6. UpSet Plot (Better for 5+ Sets)
# ==============================================================================

# Convert to binary matrix for UpSetR
all_pathways <- unique(unlist(sig_pathways))

upset_matrix <- sapply(names(sig_pathways), function(method) {
  as.integer(all_pathways %in% sig_pathways[[method]])
})
rownames(upset_matrix) <- all_pathways
upset_df <- as.data.frame(upset_matrix)

png("upset_component_tests_arabidopsis.png", width = 1600, height = 1600, res = 300)
upset(
  upset_df,
  sets = names(sig_pathways),
  order.by = "freq",
  decreasing = TRUE,
  mainbar.y.label = "Intersection Size",
  sets.x.label = "Set Size",
  text.scale = 2.2,
  point.size = 4.5,
  line.size = 1.5,
  mb.ratio = c(0.6, 0.4),
  sets.bar.color = "steelblue"
)
dev.off()

cat("UpSet plot saved to upset_component_tests_arabidopsis.pdf\n")

# ==============================================================================
# 7. Pairwise Jaccard Similarity
# ==============================================================================

jaccard <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(1)
  inter <- length(intersect(a, b))
  union <- length(union(a, b))
  if (union == 0) return(0)
  inter / union
}

methods <- names(sig_pathways)
n_methods <- length(methods)

jaccard_matrix <- matrix(NA, nrow = n_methods, ncol = n_methods,
                         dimnames = list(methods, methods))

for (i in 1:n_methods) {
  for (j in 1:n_methods) {
    jaccard_matrix[i, j] <- jaccard(sig_pathways[[i]], sig_pathways[[j]])
  }
}

cat("\nJaccard Similarity Matrix:\n")
print(round(jaccard_matrix, 3))

# Heatmap of Jaccard similarity
jaccard_df <- as.data.frame(as.table(jaccard_matrix))
names(jaccard_df) <- c("Method1", "Method2", "Jaccard")

jaccard_heatmap <- ggplot(jaccard_df, aes(x = Method1, y = Method2, fill = Jaccard)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Jaccard, 2)), color = "black", size = 6) +
  scale_fill_gradient2(low = "white", mid = "steelblue", high = "darkblue",
                       midpoint = 0.5, limits = c(0, 1)) +
  labs(title = "Pairwise Jaccard Similarity Between Component Tests",
       subtitle = paste0("Based on ", threshold_label, " significant pathways"),
       x = "", y = "") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "right") 


quartz()
print(jaccard_heatmap)
#ggsave("jaccard_heatmap_component_tests.pdf", jaccard_heatmap, width = 8, height = 7)
ggsave("jaccard_heatmap_component_tests.png", jaccard_heatmap, width = 8, height = 8, dpi = 300)

cat("Jaccard heatmap saved to jaccard_heatmap_component_tests.pdf/png\n")

# ==============================================================================
# 8. Rank Correlation Analysis
# ==============================================================================

# Compare rankings across methods
rank_df <- data.frame(pathway_id = omni_results$pathway_id)

for (method in names(use_cols)) {
  col <- use_cols[method]
  if (col %in% names(omni_results)) {
    rank_df[[paste0(method, "_rank")]] <- rank(omni_results[[col]], ties.method = "average", na.last = "keep")
  }
}

# Spearman correlation of ranks
rank_cols <- grep("_rank$", names(rank_df), value = TRUE)
rank_matrix <- as.matrix(rank_df[, rank_cols])
colnames(rank_matrix) <- sub("_rank$", "", colnames(rank_matrix))

spearman_cor <- cor(rank_matrix, use = "pairwise.complete.obs", method = "spearman")

cat("\nSpearman Rank Correlation Matrix:\n")
print(round(spearman_cor, 3))

# Heatmap of Spearman correlations
spearman_df <- as.data.frame(as.table(spearman_cor))
names(spearman_df) <- c("Method1", "Method2", "Correlation")

spearman_heatmap <- ggplot(spearman_df, aes(x = Method1, y = Method2, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Correlation, 2)), color = "black", size = 6) +
  scale_fill_gradient2(low = "red", mid = "white", high = "darkblue",
                       midpoint = 0.5, limits = c(0, 1)) +
  labs(title = "Spearman Rank Correlation Between Component Tests",
       subtitle = "Correlation of pathway rankings across all pathways",
       x = "", y = "") +
  theme_minimal(base_size=18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

quartz()
print(spearman_heatmap)
#ggsave("spearman_heatmap_component_tests.pdf", spearman_heatmap, width = 8, height = 7)
ggsave("spearman_heatmap_component_tests.png", spearman_heatmap, width = 8, height = 8, dpi = 300)

cat("Spearman heatmap saved to spearman_heatmap_component_tests.pdf/png\n")

# ==============================================================================
# 9. Summary Table: Top Pathways Per Method
# ==============================================================================

top_n <- 10

top_pathways_list <- lapply(names(use_cols), function(method) {
  col <- use_cols[method]
  if (!col %in% names(omni_results)) return(NULL)

  df <- omni_results[, c("pathway_id", "pathway_name", col)]
  names(df)[3] <- "pvalue"
  df <- df[order(df$pvalue), ]
  df <- head(df, top_n)
  df$method <- method
  df$rank <- 1:nrow(df)
  df
})
names(top_pathways_list) <- names(use_cols)

top_pathways_combined <- do.call(rbind, top_pathways_list)
rownames(top_pathways_combined) <- NULL

cat("\n\nTop", top_n, "Pathways per Component Test:\n")
cat("=========================================\n")

for (method in names(use_cols)) {
  cat("\n", method, ":\n", sep = "")
  subset_df <- top_pathways_list[[method]]
  if (!is.null(subset_df)) {
    print(subset_df[, c("rank", "pathway_id", "pathway_name", "pvalue")], row.names = FALSE)
  }
}

# Save combined table
write.csv(top_pathways_combined, "top_pathways_per_method.csv", row.names = FALSE)
cat("\nTop pathways table saved to top_pathways_per_method.csv\n")

# ==============================================================================
# 10. Unique Pathways Per Method (Found ONLY by that method)
# ==============================================================================

cat("\n\nUnique Pathways (found ONLY by that method in top/significant set):\n")
cat("====================================================================\n")

for (method in names(sig_pathways)) {
  other_methods <- setdiff(names(sig_pathways), method)
  other_pathways <- unique(unlist(sig_pathways[other_methods]))
  unique_to_method <- setdiff(sig_pathways[[method]], other_pathways)

  cat("\n", method, " (n = ", length(unique_to_method), "):\n", sep = "")

  if (length(unique_to_method) > 0) {
    # Get pathway names
    unique_info <- omni_results[omni_results$pathway_id %in% unique_to_method,
                                c("pathway_id", "pathway_name")]
    for (i in 1:nrow(unique_info)) {
      cat("  - ", unique_info$pathway_id[i], ": ", unique_info$pathway_name[i], "\n", sep = "")
    }
  } else {
    cat("  (none)\n")
  }
}

# ==============================================================================
# 11. P-value Distribution Comparison
# ==============================================================================

# Reshape for plotting
pval_long <- omni_results %>%
  select(pathway_id, all_of(unname(use_cols[use_cols %in% names(omni_results)]))) %>%
  pivot_longer(
    cols = -pathway_id,
    names_to = "method",
    values_to = "pvalue"
  ) %>%
  mutate(
    method = case_when(
      method == "acat_p" | method == "acat_p_mvn_cal" ~ "ACAT",
      method == "fisher_p" | method == "fisher_p_mvn_cal" ~ "Fisher",
      method == "tfisher_p_analytic" | method == "tfisher_p_mvn_cal" ~ "TFisher",
      method == "minp_p_analytic" | method == "minp_p_mvn_cal" ~ "minP",
      method == "stouffer_p_analytic" | method == "stouffer_p_mvn_cal" ~ "Stouffer",
      TRUE ~ method
    ),
    log10_pval = -log10(pvalue)
  )

# Density plot
density_plot <- ggplot(pval_long, aes(x = log10_pval, fill = method, color = method)) +
  geom_density(alpha = 0.3) +
  labs(
    title = "P-value Distribution Across Component Tests",
    x = expression(-log[10](p)),
    y = "Density"
  ) +
  theme_minimal() +
  theme(legend.position = "right") + plot_theme
quartz()
print(density_plot)
#ggsave("pvalue_density_component_tests.pdf", density_plot, width = 10, height = 6)
ggsave("pvalue_density_component_tests.png", density_plot, width = 8, height = 8, dpi = 300)

# Boxplot
boxplot_plot <- ggplot(pval_long, aes(x = method, y = log10_pval, fill = method)) +
  geom_boxplot(alpha = 0.7) +
  labs(
    title = "P-value Distribution Across Component Tests",
    x = "Method",
    y = expression(-log[10](p))
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) + plot_theme

print(boxplot_plot)
#ggsave("pvalue_boxplot_component_tests.pdf", boxplot_plot, width = 8, height = 6)
ggsave("pvalue_boxplot_component_tests.png", boxplot_plot, width = 8, height = 8, dpi = 300)

cat("\nP-value distribution plots saved\n")

# ==============================================================================
# 12. Pairwise Scatterplots of P-values
# ==============================================================================

# Create wide format for scatterplots
pval_wide <- omni_results %>%
  select(pathway_id, all_of(unname(use_cols[use_cols %in% names(omni_results)])))

# Rename columns for clarity
col_rename <- setNames(names(use_cols), use_cols)
for (old_name in names(col_rename)) {
  if (old_name %in% names(pval_wide)) {
    names(pval_wide)[names(pval_wide) == old_name] <- col_rename[old_name]
  }
}

# Pairs plot (using GGally if available)
if (requireNamespace("GGally", quietly = TRUE)) {
  library(GGally)

  pairs_plot <- ggpairs(
    pval_wide %>% select(-pathway_id) %>% mutate(across(everything(), ~-log10(.))),
    title = "Pairwise Comparison of -log10(p) Across Methods",
    upper = list(continuous = wrap("cor", size = 4)),
    lower = list(continuous = wrap("points", alpha = 0.3, size = 0.5)),
    diag = list(continuous = wrap("densityDiag", alpha = 0.5))
  ) +
    theme_minimal()

  ggsave("pairwise_pvalue_scatter.pdf", pairs_plot, width = 12, height = 12)
  ggsave("pairwise_pvalue_scatter.png", pairs_plot, width = 12, height = 12, dpi = 300)
  cat("Pairwise scatterplot saved to pairwise_pvalue_scatter.pdf/png\n")
} else {
  cat("Install GGally for pairwise scatterplot: install.packages('GGally')\n")
}

# ==============================================================================
# 13. Summary Statistics
# ==============================================================================

cat("\n\n========== SUMMARY ==========\n")
cat("Total pathways analyzed:", n_pathways, "\n")
cat("Bonferroni threshold:", format(bonferroni_alpha, digits = 3), "\n")
cat("Threshold used for comparison:", threshold_label, "\n\n")

cat("Component test agreement summary:\n")
cat("- Average Jaccard similarity:", round(mean(jaccard_matrix[lower.tri(jaccard_matrix)]), 3), "\n")
cat("- Average Spearman correlation:", round(mean(spearman_cor[lower.tri(spearman_cor)]), 3), "\n")

# Count pathways found by exactly 1, 2, 3, 4, 5 methods
pathway_counts <- table(rowSums(upset_df))
cat("\nPathways found by N methods:\n")
for (n in 1:5) {
  if (as.character(n) %in% names(pathway_counts)) {
    cat("  ", n, " method(s):", pathway_counts[as.character(n)], "\n")
  }
}

# Core pathways (found by all methods)
core_pathways <- all_pathways[rowSums(upset_df) == length(sig_pathways)]
cat("\nCore pathways (found by ALL methods): n =", length(core_pathways), "\n")
if (length(core_pathways) > 0 && length(core_pathways) <= 20) {
  core_info <- omni_results[omni_results$pathway_id %in% core_pathways, c("pathway_id", "pathway_name")]
  for (i in 1:nrow(core_info)) {
    cat("  - ", core_info$pathway_id[i], ": ", core_info$pathway_name[i], "\n", sep = "")
  }
}

cat("\n========== DONE ==========\n")
cat("All output files saved in:", getwd(), "\n")
