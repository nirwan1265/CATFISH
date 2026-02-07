# ==============================================================================
# Compare OMNIBUS Test vs Individual Component Tests
# Shows that OMNIBUS is a holistic combination, not just picking from one test
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggvenn)

# For Venn diagrams
if (!requireNamespace("ggVennDiagram", quietly = TRUE)) {
  install.packages("ggVennDiagram")
}
library(ggVennDiagram)

# For upset plots
if (!requireNamespace("UpSetR", quietly = TRUE)) {
  install.packages("UpSetR")
}
library(UpSetR)

# ==============================================================================
# 0. Plot Theme
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

omni_results <- readRDS(
  "/Users/nirwantandukar/Documents/Github/MAGCAT/omni_results_arabidopsis_ACAT_B10000.RDS"
)

head(omni_results)
names(omni_results)
write.csv(omni_results,"omni_results_ACAT.csv",row.names = F)
# ==============================================================================
# 2. Define P-value Columns
# ==============================================================================

# Component p-values (use calibrated if available)
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

# Check if calibrated available
calibrated_available <- all(calibrated_cols %in% names(omni_results))

if (calibrated_available) {
  comp_cols <- calibrated_cols
  cat("Using MVN-calibrated component p-values\n")
} else {
  comp_cols <- analytic_cols
  cat("Using analytic component p-values\n")
}

# OMNIBUS p-value (use final which is MVN > global > analytic)
omni_col <- "omni_p_final"
if (!omni_col %in% names(omni_results)) {
  omni_col <- "omni_p_mvn"
  if (!omni_col %in% names(omni_results)) {
    omni_col <- "omni_p_analytic"
  }
}
cat("Using OMNIBUS column:", omni_col, "\n")

# ==============================================================================
# 3. Determine Significance Threshold
# ==============================================================================

n_pathways <- nrow(omni_results)
bonferroni_alpha <- 0.05 / n_pathways
#nominal_alpha <- 0.1335 # changed to get 20 pathways
nominal_alpha <- 0.1

cat("\nNumber of pathways:", n_pathways, "\n")
cat("Bonferroni threshold:", format(bonferroni_alpha, digits = 3), "\n")

# Check omnibus significant count
omni_bonf_sig <- sum(omni_results[[omni_col]] < bonferroni_alpha, na.rm = TRUE)
omni_nom_sig <- sum(omni_results[[omni_col]] < nominal_alpha, na.rm = TRUE)

cat("\nOMNIBUS significant pathways:\n")
cat("  Bonferroni:", omni_bonf_sig, "\n")
cat("  Nominal (0.1):", omni_nom_sig, "\n")

# Use Bonferroni if we have hits, otherwise nominal, otherwise top 20
if (omni_bonf_sig > 0) {
  sig_threshold <- bonferroni_alpha
  threshold_label <- "Bonferroni"
} else if (omni_nom_sig >= 5) {
  sig_threshold <- nominal_alpha
  threshold_label <- "Nominal (0.1)"
} else {
  sig_threshold <- NULL
  threshold_label <- "Top 20"
}

cat("Using threshold:", threshold_label, "\n")

# ==============================================================================
# 4. Get Significant Pathways
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

# Get OMNIBUS significant pathways
omni_sig <- get_sig_pathways(omni_results, omni_col, sig_threshold, top_n = 40)
cat("\nOMNIBUS selected pathways:", length(omni_sig), "\n")

# Get component significant pathways
comp_sig <- lapply(names(comp_cols), function(m) {
  get_sig_pathways(omni_results, comp_cols[m], sig_threshold, top_n = 40)
})
names(comp_sig) <- names(comp_cols)

cat("Component test selected pathways:\n")
for (m in names(comp_sig)) {
  cat("  ", m, ":", length(comp_sig[[m]]), "\n")
}

# ==============================================================================
# 5. Venn: OMNIBUS vs All Components Combined
# ==============================================================================

# Union of all component-significant pathways
all_comp_union <- unique(unlist(comp_sig))

venn_omni_vs_union <- list(
  "OMNIBUS" = omni_sig,
  "Component tests" = all_comp_union
)

omni <- venn_omni_vs_union$OMNIBUS
comp <- venn_omni_vs_union$`Component tests`

n_omni_only  <- length(setdiff(omni, comp))      # 0
n_both       <- length(intersect(omni, comp))    # 15
n_comp_only  <- length(setdiff(comp, omni))      # 24

# make sure names match what you want printed
venn_list <- list(
  OMNIBUS = venn_omni_vs_union$OMNIBUS,
  `Any Component` = venn_omni_vs_union$`Component tests`
)
library(ggVennDiagram)
library(ggplot2)

# Union of all component-significant pathways
all_comp_union <- unique(unlist(comp_sig))

venn_omni_vs_union <- list(
  "OMNIBUS" = omni_sig,
  "Component tests" = all_comp_union
)

venn_plot1 <- ggVennDiagram(
  venn_omni_vs_union,
  label = "count",
  label_alpha = 0,
  edge_size = 0.8
) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    title = "OMNIBUS vs Union of All Component Tests",
    subtitle = paste0("Threshold: ", threshold_label),
    fill = "count"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 24),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(size = 12)
  )

venn_plot1 <- ggVennDiagram(
  venn_omni_vs_union,
  label       = "count",
  label_alpha = 0,
  edge_size   = 0.8,
  label_size  = 10,          # <-- increase the numbers inside regions
  set_size    = 0            # <-- remove the set names (OMNIBUS / Component tests)
) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    title    = "OMNIBUS vs Union of All Component Tests",
    subtitle = paste0("Threshold: ", threshold_label),
    fill     = "count"
  ) +
  theme_void() +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 24),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(size = 14)
  )
quartz()
print(venn_plot1)



# ggsave(
#   "venn_omnibus_vs_components_union.pdf",
#   venn_plot1, width = 8, height = 8
# )
ggsave(
  "venn_omnibus_vs_components_union.png",
  venn_plot1, width = 8, height = 8, dpi = 300
)

# ==============================================================================
# 6. Venn: OMNIBUS vs Each Component (5 separate plots)
# ==============================================================================

# venn_plots <- list()

# for (method in names(comp_sig)) {
#   venn_list <- list(
#     "OMNIBUS" = omni_sig,
#     method = comp_sig[[method]]
#   )
#   names(venn_list)[2] <- method

#   p <- ggVennDiagram(
#     venn_list,
#     label = "count",
#     label_alpha = 0
#   ) +
#     scale_fill_gradient(low = "white", high = "coral") +
#     labs(title = paste("OMNIBUS vs", method)) +
#     theme(plot.title = element_text(hjust = 0.5))

#   venn_plots[[method]] <- p
# }

venn_plots <- list()

for (method in names(comp_sig)) {

  venn_list <- list(
    "OMNIBUS" = omni_sig,
    method   = comp_sig[[method]]
  )
  names(venn_list)[2] <- method

  p <- ggVennDiagram(
    venn_list,
    label       = "count",
    label_alpha = 0,
    label_size  = 8,   # bigger numbers
    set_size    = 0    # remove set labels
  ) +
    scale_fill_gradient(low = "white", high = "coral") +
    labs(
      title = paste("OMNIBUS vs", method),
      fill  = "count"
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position = "top"   # legend on top
    )

  venn_plots[[method]] <- p
}


# Combine into one figure
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  combined_venn <- (venn_plots[[1]] + venn_plots[[2]] + venn_plots[[3]]) /
    (venn_plots[[4]] + venn_plots[[5]] + plot_spacer()) +
    plot_annotation(
      title = "OMNIBUS vs Individual Component Tests",
      subtitle = paste0("Threshold: ", threshold_label)
    )
  quartz()
  print(combined_venn)
  # ggsave(
  #   "venn_omnibus_vs_each_component.pdf",
  #   combined_venn, width = 14, height = 10
  # )
  ggsave(
    "venn_omnibus_vs_each_component.png",
    combined_venn, width = 16, height = 8, dpi = 300
  )
} else {
  cat("Install patchwork for combined plot: install.packages('patchwork')\n")
}

# ==============================================================================
# 7. KEY ANALYSIS: How Many Components Support Each OMNIBUS Pathway?
# ==============================================================================

# For each OMNIBUS-significant pathway, count how many component tests
# also found it significant

omni_sig_df <- omni_results[omni_results$pathway_id %in% omni_sig, ]

# Count component support for each omnibus pathway
omni_sig_df$n_components_support <- sapply(
  omni_sig_df$pathway_id,
  function(pid) {
    sum(sapply(comp_sig, function(cs) pid %in% cs))
  }
)

# Which specific components support it?
omni_sig_df$supporting_components <- sapply(
  omni_sig_df$pathway_id,
  function(pid) {
    supporters <- names(comp_sig)[sapply(comp_sig, function(cs) pid %in% cs)]
    if (length(supporters) == 0) return("None")
    paste(supporters, collapse = ", ")
  }
)

# Summary
cat("\n\n========================================\n")
cat("OMNIBUS Pathways: Component Support Analysis\n")
cat("========================================\n\n")

support_table <- table(omni_sig_df$n_components_support)
cat("OMNIBUS pathways supported by N component tests:\n")
for (n in sort(unique(omni_sig_df$n_components_support))) {
  pct <- round(100 * support_table[as.character(n)] / length(omni_sig), 1)
  cat("  ", n, " components:", support_table[as.character(n)],
      "(", pct, "%)\n")
}

# Bar plot of support
support_bar <- ggplot(
  omni_sig_df,
  aes(x = factor(n_components_support))
) +
  geom_bar(fill = "steelblue", alpha = 0.8) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  labs(
    title = "Component Test Support for OMNIBUS-Significant Pathways",
    subtitle = paste0(
      "How many individual tests also identify each OMNIBUS pathway? (",
      threshold_label, ")"
    ),
    x = "Number of Supporting Component Tests (out of 5)",
    y = "Number of Pathways"
  ) +
  scale_x_discrete(
    labels = c("0" = "0\n(OMNIBUS unique)",
               "1" = "1", "2" = "2", "3" = "3", "4" = "4",
               "5" = "5\n(All agree)")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  ) + plot_theme

quartz()
print(support_bar)

#ggsave("omnibus_component_support.pdf", support_bar, width = 10, height = 6)
ggsave("omnibus_component_support.png", support_bar, width = 8, height = 8,
       dpi = 300)

# ==============================================================================
# 8. Detailed Table: OMNIBUS Pathways with Support Info
# ==============================================================================

# Prepare nice output table
output_table <- omni_sig_df %>%
  select(
    pathway_id, pathway_name, n_genes,
    all_of(omni_col),
    n_components_support, supporting_components
  ) %>%
  arrange(!!sym(omni_col))

names(output_table)[4] <- "omnibus_p"

cat("\n\nOMNIBUS-Significant Pathways with Component Support:\n")
cat(strrep("=", 80), "\n")

print(
  output_table %>%
    select(pathway_id, pathway_name, omnibus_p,
           n_components_support, supporting_components),
  n = 30
)

write.csv(output_table, "omnibus_pathways_with_support.csv", row.names = FALSE)
cat("\nTable saved to omnibus_pathways_with_support.csv\n")

# ==============================================================================
# 9. Pathways UNIQUE to OMNIBUS (not found by ANY individual component)
# ==============================================================================

omni_unique <- omni_sig_df$pathway_id[omni_sig_df$n_components_support == 0]

cat("\n\n========================================\n")
cat("Pathways Found ONLY by OMNIBUS (n =", length(omni_unique), ")\n")
cat("========================================\n")
cat("These are pathways where OMNIBUS detected signal that NO single\n")
cat("component test found - demonstrating the holistic combination!\n\n")

if (length(omni_unique) > 0) {
  unique_df <- omni_sig_df[omni_sig_df$pathway_id %in% omni_unique, ]
  for (i in seq_len(nrow(unique_df))) {
    cat(
      i, ". ", unique_df$pathway_id[i], "\n",
      "   Name: ", unique_df$pathway_name[i], "\n",
      "   OMNIBUS p: ", format(unique_df[[omni_col]][i], digits = 3), "\n",
      "   Component p-values:\n", sep = ""
    )
    for (m in names(comp_cols)) {
      col <- comp_cols[m]
      if (col %in% names(unique_df)) {
        pval <- unique_df[[col]][i]
        cat("     ", m, ": ", format(pval, digits = 3), "\n", sep = "")
      }
    }
    cat("\n")
  }
} else {
  cat("No pathways found exclusively by OMNIBUS at this threshold.\n")
  cat("This may indicate strong agreement between methods,\n")
  cat("or a less stringent threshold might reveal OMNIBUS-unique hits.\n")
}

# ==============================================================================
# 10. Pathways Found by ALL Components (Consensus)
# ==============================================================================

consensus <- omni_sig_df$pathway_id[omni_sig_df$n_components_support == 5]

cat("\n\n========================================\n")
cat("Consensus Pathways (Found by ALL 5 Components, n =",
    length(consensus), ")\n")
cat("========================================\n")

if (length(consensus) > 0) {
  consensus_df <- omni_sig_df[omni_sig_df$pathway_id %in% consensus, ]
  for (i in seq_len(nrow(consensus_df))) {
    cat(
      i, ". ", consensus_df$pathway_id[i], ": ",
      consensus_df$pathway_name[i], "\n", sep = ""
    )
  }
}

# ==============================================================================
# 11. UpSet Plot: OMNIBUS + All Components
# ==============================================================================

# Combine omnibus with components for upset
all_sets <- c(list("OMNIBUS" = omni_sig), comp_sig)
all_pathways <- unique(unlist(all_sets))

upset_matrix <- sapply(names(all_sets), function(s) {
  as.integer(all_pathways %in% all_sets[[s]])
})
rownames(upset_matrix) <- all_pathways
upset_df <- as.data.frame(upset_matrix)

png("upset_omnibus_and_components.png", width = 1600, height = 1600, res = 300)
upset(
  upset_df,
  sets = c("OMNIBUS", names(comp_cols)),
  keep.order = TRUE,
  order.by = "freq",
  decreasing = TRUE,
  mainbar.y.label = "Intersection Size",
  sets.x.label = "Set Size",
  text.scale = 2.2,
  point.size = 4.5,
  line.size = 1.5,
  mb.ratio = c(0.6, 0.4),
  sets.bar.color = c("darkred", rep("steelblue", 5)),
  main.bar.color = "gray30"
)
dev.off()
cat("\nUpSet plot saved to upset_omnibus_and_components.pdf\n")

# ==============================================================================
# 12. Heatmap: P-value Comparison for OMNIBUS Pathways
# ==============================================================================

# Get component p-values for omnibus-significant pathways
heatmap_df <- omni_sig_df %>%
  select(pathway_id, pathway_name, all_of(c(omni_col, unname(comp_cols)))) %>%
  mutate(across(-c(pathway_id, pathway_name), ~ -log10(.))) %>%
  arrange(desc(!!sym(omni_col)))

# Reshape for heatmap
heatmap_long <- heatmap_df %>%
  pivot_longer(
    cols = -c(pathway_id, pathway_name),
    names_to = "method",
    values_to = "neg_log10_p"
  ) %>%
  mutate(
    method = case_when(
      method == omni_col ~ "OMNIBUS",
      method == "acat_p" | method == "acat_p_mvn_cal" ~ "ACAT",
      method == "fisher_p" | method == "fisher_p_mvn_cal" ~ "Fisher",
      method == "tfisher_p_analytic" |
        method == "tfisher_p_mvn_cal" ~ "TFisher",
      method == "minp_p_analytic" | method == "minp_p_mvn_cal" ~ "minP",
      method == "stouffer_p_analytic" |
        method == "stouffer_p_mvn_cal" ~ "Stouffer",
      TRUE ~ method
    ),
    method = factor(
      method,
      levels = c("OMNIBUS", "ACAT", "Fisher", "TFisher", "minP", "Stouffer")
    )
  )

# Use pathway_name for y-axis, truncate if too long
heatmap_long <- heatmap_long %>%
  mutate(
    pathway_label = ifelse(
      nchar(pathway_name) > 40,
      paste0(substr(pathway_name, 1, 37), "..."),
      pathway_name
    ),
    pathway_label = factor(
      pathway_label,
      levels = rev(unique(pathway_label))
    )
  )

# Significance line
if (!is.null(sig_threshold)) {
  sig_line <- -log10(sig_threshold)
} else {
  sig_line <- -log10(0.05)  # nominal for reference
}

heatmap_plot <- ggplot(
  heatmap_long,
  aes(x = method, y = pathway_label, fill = neg_log10_p)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(
    low = "white", mid = "steelblue", high = "darkred",
    midpoint = sig_line,
    name = expression(-log[10](p))
  ) +
  labs(
    title = "P-value Comparison: OMNIBUS vs Component Tests",
    subtitle = paste0(
      "For OMNIBUS-significant pathways (", threshold_label, ")"
    ),
    x = "Method",
    y = "Pathway"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  ) + plot_theme

quartz()
print(heatmap_plot)

# Adjust height based on number of pathways
plot_height <- max(8, length(omni_sig) * 0.3)
ggsave(
  "heatmap_omnibus_vs_components.png",
  heatmap_plot, width = 16, height = plot_height, dpi = 300
)

cat("Heatmap saved to heatmap_omnibus_vs_components.pdf/png\n")

# ==============================================================================
# 13. Summary Statistics
# ==============================================================================

cat("\n\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n")
cat("SUMMARY: OMNIBUS as a Holistic Test\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n\n")

cat("Total OMNIBUS-significant pathways:", length(omni_sig), "\n\n")

cat("Breakdown by component support:\n")
for (n in 0:5) {
  count <- sum(omni_sig_df$n_components_support == n)
  pct <- round(100 * count / length(omni_sig), 1)
  label <- case_when(
    n == 0 ~ "OMNIBUS-unique (holistic signal)",
    n == 1 ~ "Supported by 1 component",
    n == 5 ~ "Consensus (all 5 agree)",
    TRUE ~ paste0("Supported by ", n, " components")
  )
  cat("  ", n, " components: ", count, " (", pct, "%) - ", label, "\n",
      sep = "")
}

# Jaccard of OMNIBUS vs each component
cat("\nJaccard similarity (OMNIBUS vs each component):\n")
for (m in names(comp_sig)) {
  inter <- length(intersect(omni_sig, comp_sig[[m]]))
  union <- length(union(omni_sig, comp_sig[[m]]))
  jacc <- if (union > 0) round(inter / union, 3) else 0
  cat("  OMNIBUS vs ", m, ": ", jacc, "\n", sep = "")
}

# Key insight
holistic_pct <- round(
  100 * sum(omni_sig_df$n_components_support < 5) / length(omni_sig),
  1
)
cat("\n** Key Insight **\n")
cat(holistic_pct, "% of OMNIBUS pathways are NOT found by all 5 components,\n")
cat("demonstrating that OMNIBUS synthesizes complementary signals\n")
cat("rather than just picking the most obvious hits.\n")


# ==============================================================================
# 14. Cumulative Support: Fraction of OMNIBUS Pathways Supported by >= k Components
# ==============================================================================

n_omni <- nrow(omni_sig_df)
cum_support <- data.frame(
  k = 0:5,
  count = sapply(0:5, function(k) sum(omni_sig_df$n_components_support >= k)),
  stringsAsFactors = FALSE
)
cum_support$fraction <- cum_support$count / n_omni
cum_support$label <- paste0(round(100 * cum_support$fraction, 1), "%")

cum_bar <- ggplot(cum_support, aes(x = factor(k), y = fraction)) +
  geom_col(fill = "steelblue", alpha = 0.85, width = 0.65) +
  geom_text(aes(label = label), vjust = -0.4, size = 6, fontface = "bold") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1.08),
    expand = c(0, 0)
  ) +
  labs(
    title = "Omnibus Behaves as a Broad Union, Not a Narrow Intersection",
    subtitle = paste0(
      "Fraction of OMNIBUS pathways also found by \u2265 k component tests (",
      threshold_label, ")"
    ),
    x = "Minimum Number of Supporting Components (k)",
    y = "Fraction of OMNIBUS Pathways"
  ) +
  plot_theme +
  theme(
    plot.title    = element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 16)
  )

quartz()
print(cum_bar)
ggsave("cumulative_support_omnibus.png", cum_bar,
       width = 8, height = 8, dpi = 300)
cat("Cumulative support bar saved to cumulative_support_omnibus.png\n")

# ==============================================================================
# 15. Jaccard & Overlap Proportion: OMNIBUS vs Each Component
# ==============================================================================

jacc_df <- data.frame(
  Component = names(comp_sig),
  Jaccard = sapply(comp_sig, function(cs) {
    inter <- length(intersect(omni_sig, cs))
    uni   <- length(union(omni_sig, cs))
    if (uni > 0) inter / uni else 0
  }),
  Overlap = sapply(comp_sig, function(cs) {
    inter <- length(intersect(omni_sig, cs))
    denom <- min(length(omni_sig), length(cs))
    if (denom > 0) inter / denom else 0
  }),
  stringsAsFactors = FALSE
)
rownames(jacc_df) <- NULL

jacc_long <- jacc_df %>%
  pivot_longer(cols = c(Jaccard, Overlap),
               names_to = "Metric", values_to = "Value")

jacc_plot <- ggplot(jacc_long,
                    aes(x = Component, y = Value, fill = Metric)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_text(
    aes(label = sprintf("%.2f", Value)),
    position = position_dodge(width = 0.7),
    vjust = -0.4, size = 5, fontface = "bold"
  ) +
  scale_fill_manual(values = c("Jaccard" = "steelblue", "Overlap" = "coral")) +
  scale_y_continuous(limits = c(0, 1.08), expand = c(0, 0)) +
  labs(
    title = "Set Similarity: OMNIBUS vs Each Component Test",
    subtitle = paste0(
      "Jaccard = |A\u2229B|/|A\u222AB|,  Overlap = |A\u2229B|/min(|A|,|B|)  (",
      threshold_label, ")"
    ),
    x = "Component Test",
    y = "Similarity"
  ) +
  plot_theme +
  theme(
    plot.title    = element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    legend.position = "top"
  )

quartz()
print(jacc_plot)
ggsave("jaccard_overlap_omnibus_vs_components.png", jacc_plot,
       width = 8, height = 8, dpi = 300)
cat("Jaccard/Overlap bar saved to jaccard_overlap_omnibus_vs_components.png\n")

# ==============================================================================
# 16. Heatmap with "Winner" Strip: Which Component Drives Each OMNIBUS Pathway?
# ==============================================================================

if (!requireNamespace("ggnewscale", quietly = TRUE)) {
  install.packages("ggnewscale")
}
library(ggnewscale)

# For each omnibus pathway find the component with the smallest p-value
winner_df <- omni_sig_df %>%
  rowwise() %>%
  mutate(
    best_comp_p = min(c_across(all_of(unname(comp_cols))), na.rm = TRUE),
    winner = {
      pvals <- c_across(all_of(unname(comp_cols)))
      names(comp_cols)[which.min(pvals)]
    }
  ) %>%
  ungroup() %>%
  select(pathway_id, pathway_name, all_of(c(omni_col, unname(comp_cols))), winner)

# Build heatmap data
heat_df <- winner_df %>%
  mutate(across(
    all_of(c(omni_col, unname(comp_cols))),
    ~ -log10(.)
  )) %>%
  arrange(desc(!!sym(omni_col)))

# Truncate names, then make unique to avoid duplicate levels
pathway_order <- ifelse(
  nchar(heat_df$pathway_name) > 40,
  paste0(substr(heat_df$pathway_name, 1, 37), "..."),
  heat_df$pathway_name
)
pathway_order <- make.unique(pathway_order, sep = " ")

# Map pathway_id -> unique label (preserves row order)
id_to_label <- setNames(pathway_order, heat_df$pathway_id)

heat_long <- heat_df %>%
  mutate(pathway_label = id_to_label[pathway_id]) %>%
  pivot_longer(
    cols = all_of(c(omni_col, unname(comp_cols))),
    names_to = "method",
    values_to = "neg_log10_p"
  ) %>%
  mutate(
    method = case_when(
      method == omni_col ~ "OMNIBUS",
      method == "acat_p" | method == "acat_p_mvn_cal" ~ "ACAT",
      method == "fisher_p" | method == "fisher_p_mvn_cal" ~ "Fisher",
      method == "tfisher_p_analytic" |
        method == "tfisher_p_mvn_cal" ~ "TFisher",
      method == "minp_p_analytic" |
        method == "minp_p_mvn_cal" ~ "minP",
      method == "stouffer_p_analytic" |
        method == "stouffer_p_mvn_cal" ~ "Stouffer",
      TRUE ~ method
    ),
    method = factor(
      method,
      levels = c("OMNIBUS", "ACAT", "Fisher",
                  "TFisher", "minP", "Stouffer")
    ),
    pathway_label = factor(
      pathway_label, levels = rev(pathway_order)
    )
  )

# Component colors (fixed per component)
comp_colors <- c(
  "ACAT"     = "#E41A1C",
  "Fisher"   = "#377EB8",
  "TFisher"  = "#4DAF4A",
  "minP"     = "#984EA3",
  "Stouffer" = "#FF7F00"
)

# For each pathway, find the order of components (1st best, 5th worst)
comp_names <- c("ACAT", "Fisher", "TFisher", "minP", "Stouffer")
comp_col_vals <- unname(comp_cols)

# Build rank strip: which component is in each position (1st, 2nd, 3rd, 4th, 5th)
rank_strip_list <- lapply(seq_len(nrow(winner_df)), function(i) {
  pvals <- as.numeric(winner_df[i, comp_col_vals])
  ord <- order(pvals)  # indices sorted by p-value (smallest first)
  data.frame(
    pathway_id = winner_df$pathway_id[i],
    pos_1 = comp_names[ord[1]],
    pos_2 = comp_names[ord[2]],
    pos_3 = comp_names[ord[3]],
    pos_4 = comp_names[ord[4]],
    pos_5 = comp_names[ord[5]],
    stringsAsFactors = FALSE
  )
})
rank_strip_df <- do.call(rbind, rank_strip_list)
rank_strip_df$pathway_label <- id_to_label[rank_strip_df$pathway_id]

# Pivot to long format
rank_long <- rank_strip_df %>%
  pivot_longer(
    cols = starts_with("pos_"),
    names_to = "position",
    values_to = "component"
  ) %>%
  mutate(
    pathway_label = factor(pathway_label, levels = rev(pathway_order)),
    position = factor(position, levels = c("pos_1", "pos_2", "pos_3", "pos_4", "pos_5")),
    component = factor(component, levels = comp_names)
  )

p_heat2 <- ggplot() +
  # Main heatmap
  geom_tile(
    data = heat_long,
    aes(x = method, y = pathway_label, fill = neg_log10_p),
    color = "white", linewidth = 0.5
  ) +
  scale_fill_gradient2(
    low = "white", mid = "steelblue", high = "darkred",
    midpoint = sig_line,
    name = expression(-log[10](p))
  ) +
  new_scale_fill() +
  # Rank strip: thin tiles colored by component
  geom_tile(
    data = rank_long,
    aes(x = position, y = pathway_label, fill = component),
    width = 0.4, color = "white", linewidth = 0.2
  ) +
  scale_fill_manual(values = comp_colors, name = "Component") +
  scale_x_discrete(
    limits = c("OMNIBUS", "ACAT", "Fisher", "TFisher", "minP", "Stouffer",
               "pos_1", "pos_2", "pos_3", "pos_4", "pos_5"),
    labels = c("OMNIBUS", "ACAT", "Fisher", "TFisher", "minP", "Stouffer",
               "1", "2", "3", "4", "5")
  ) +
  labs(
    title = "P-value Heatmap with Ranking Strip (Arabidopsis)",
    subtitle = paste0("Right strip: component rank order (1=best, 5=worst) | ", threshold_label),
    x = "Method",
    y = "Pathway"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "right"
  ) + plot_theme

quartz()
print(p_heat2)

plot_height2 <- max(8, length(omni_sig) * 0.3)
ggsave("heatmap_with_winner_strip.png", p_heat2,
       width = 18, height = plot_height2, dpi = 300, bg = "white")
cat("Heatmap with winner strip saved to heatmap_with_winner_strip.png\n")

# Winner distribution summary
cat("\nBest-component distribution across OMNIBUS pathways:\n")
winner_tab <- table(winner_df$winner)
for (w in sort(names(winner_tab))) {
  cat("  ", w, ": ", winner_tab[w], " (",
      round(100 * winner_tab[w] / sum(winner_tab), 1), "%)\n", sep = "")
}

cat("\n========== DONE ==========\n")
cat("Output files:\n")
cat("  - venn_omnibus_vs_components_union.pdf/png\n")
cat("  - venn_omnibus_vs_each_component.pdf/png\n")
cat("  - omnibus_component_support.pdf/png\n")
cat("  - omnibus_pathways_with_support.csv\n")
cat("  - upset_omnibus_and_components.pdf\n")
cat("  - heatmap_omnibus_vs_components.pdf/png\n")
cat("  - cumulative_support_omnibus.png\n")
cat("  - jaccard_overlap_omnibus_vs_components.png\n")
cat("  - heatmap_with_winner_strip.png\n")

