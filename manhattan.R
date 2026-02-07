# ==============================================================================
# Manhattan Plots: GWAS (SNP-level) and MAGMA (Gene-level)
# Arabidopsis, Fly Male, Fly Female
# ==============================================================================

library(data.table)
library(dplyr)
library(ggplot2)

# Install qqman if needed
if (!requireNamespace("qqman", quietly = TRUE)) {
  install.packages("qqman")
}
library(qqman)

# ==============================================================================
# Helper function for ggplot2 Manhattan (more customizable)
# ==============================================================================

gg_manhattan <- function(df, chr_col = "CHR", bp_col = "BP", p_col = "P",

                         title = "Manhattan Plot",
                                                  suggestive = 1e-5, genome_wide = NULL,
                         point_size = 1, alpha = 0.8) {

  # Prepare data
  df <- df %>%
    rename(CHR = !!sym(chr_col), BP = !!sym(bp_col), P = !!sym(p_col)) %>%
    filter(!is.na(P), !is.na(CHR), !is.na(BP)) %>%
    mutate(
      CHR = as.character(CHR),
      CHR_num = as.numeric(factor(CHR, levels = unique(CHR))),
      logP = -log10(P)
    )

  # Calculate cumulative position
  chr_info <- df %>%
    group_by(CHR, CHR_num) %>%
    summarise(chr_len = max(BP), .groups = "drop") %>%
    arrange(CHR_num) %>%
    mutate(
      cum_start = lag(cumsum(chr_len), default = 0),
      chr_center = cum_start + chr_len / 2
    )

  df <- df %>%
    left_join(chr_info %>% select(CHR, cum_start), by = "CHR") %>%
    mutate(BP_cum = BP + cum_start)

  # Color palette
  colors <- rep(c("steelblue", "coral"), length.out = length(unique(df$CHR)))

  # Bonferroni threshold
  if (is.null(genome_wide)) {
    genome_wide <- 0.05 / nrow(df)
  }

  p <- ggplot(df, aes(x = BP_cum, y = logP, color = factor(CHR_num))) +
    geom_point(size = point_size, alpha = alpha) +
    scale_color_manual(values = colors, guide = "none") +
    scale_x_continuous(
      breaks = chr_info$chr_center,
      labels = chr_info$CHR,
      expand = c(0.01, 0.01)
    ) +
    geom_hline(yintercept = -log10(suggestive),
               linetype = "dashed", color = "blue", linewidth = 0.5) +
    geom_hline(yintercept = -log10(genome_wide),
               linetype = "solid", color = "red", linewidth = 0.5) +
    labs(
      title = title,
      x = "Chromosome",
      y = expression(-log[10](p))
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    annotate("text", x = max(df$BP_cum) * 0.98,
             y = -log10(suggestive) + 0.3,
             label = paste0("Suggestive (", suggestive, ")"),
             hjust = 1, size = 3, color = "blue") +
    annotate("text", x = max(df$BP_cum) * 0.98,
             y = -log10(genome_wide) + 0.3,
             label = paste0("Bonferroni (", format(genome_wide, digits = 2), ")"),
             hjust = 1, size = 3, color = "red")

  return(p)
}

# ==============================================================================
# 1. GWAS Manhattan Plots
# ==============================================================================

cat("\n========================================\n")
cat("GWAS Manhattan Plots (SNP-level)\n")
cat("========================================\n")

# --- 1a. Arabidopsis GWAS ---
cat("\nLoading Arabidopsis GWAS...\n")

gwas_at <- fread(
  paste0(
    "/Users/nirwantandukar/Documents/Research/data/Arabidopsis/",
    "raw_gwas/Bio6_coldest_temp.csv"
  )
)

names(gwas_at)[names(gwas_at) == "P-Value"] <- "P"
names(gwas_at)[names(gwas_at) == "Positions"] <- "BP"

cat("  Total SNPs:", nrow(gwas_at), "\n")

# qqman Manhattan
png("manhattan_gwas_arabidopsis.png", width = 1600, height = 1200, res = 300)
manhattan(
  gwas_at,
  chr = "CHR", bp = "BP", p = "P", snp = "SNP",
  main = "Arabidopsis GWAS - Cold Temperature (BIO6)",
  col = c("steelblue", "coral"),
  suggestiveline = -log10(1e-5),
  genomewideline = -log10(0.05 / nrow(gwas_at)),
  cex = 0.6
)
dev.off()
cat("  Saved: manhattan_gwas_arabidopsis.png\n")

# QQ plot
png("qq_gwas_arabidopsis.png", width = 1600, height = 1200, res = 300)
qq(gwas_at$P, main = "Arabidopsis GWAS QQ Plot")
dev.off()
cat("  Saved: qq_gwas_arabidopsis.png\n")

# --- 1b. Fly Male GWAS ---
cat("\nLoading Fly Male GWAS...\n")

gwas_fly_m <- fread(
  paste0(
    "/Users/nirwantandukar/Documents/Research/data/DGRP/",
    "Starvation_stress/raw_gwas/raw_GWAS_Starvation_stress_male_DGRP.csv"
  )
)

names(gwas_fly_m)[names(gwas_fly_m) == "P-Value"] <- "P"
names(gwas_fly_m)[names(gwas_fly_m) == "Positions"] <- "BP"

cat("  Total SNPs:", nrow(gwas_fly_m), "\n")

# Convert chromosome names to numeric for qqman
gwas_fly_m_plot <- gwas_fly_m %>%
  mutate(
    CHR_orig = CHR,
    CHR = case_when(
      CHR == "2L" ~ 1,
      CHR == "2R" ~ 2,
      CHR == "3L" ~ 3,
      CHR == "3R" ~ 4,
      CHR == "4"  ~ 5,
      CHR == "X"  ~ 6,
      CHR == "Y"  ~ 7,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(CHR), !is.na(P), !is.na(BP))

png("manhattan_gwas_fly_male.png", width = 1600, height = 1200, res = 300)
manhattan(
  gwas_fly_m_plot,
  chr = "CHR", bp = "BP", p = "P", snp = "SNP",
  main = "Fly Male GWAS - Starvation Stress",
  col = c("steelblue", "coral"),
  suggestiveline = -log10(1e-5),
  genomewideline = -log10(0.05 / nrow(gwas_fly_m_plot)),
  cex = 0.6,
  chrlabs = c("2L", "2R", "3L", "3R", "4", "X", "Y")
)
dev.off()
cat("  Saved: manhattan_gwas_fly_male.png\n")

png("qq_gwas_fly_male.png", width = 1600, height = 1200, res = 300)
qq(gwas_fly_m_plot$P, main = "Fly Male GWAS QQ Plot")
dev.off()
cat("  Saved: qq_gwas_fly_male.png\n")

# --- 1c. Fly Female GWAS ---
cat("\nLoading Fly Female GWAS...\n")

gwas_fly_f <- fread(
  paste0(
    "/Users/nirwantandukar/Documents/Research/data/DGRP/",
    "Starvation_stress/raw_gwas/raw_GWAS_Starvation_stress_female_DGRP.csv"
  )
)

names(gwas_fly_f)[names(gwas_fly_f) == "P-Value"] <- "P"
names(gwas_fly_f)[names(gwas_fly_f) == "Positions"] <- "BP"

cat("  Total SNPs:", nrow(gwas_fly_f), "\n")

gwas_fly_f_plot <- gwas_fly_f %>%
  mutate(
    CHR_orig = CHR,
    CHR = case_when(
      CHR == "2L" ~ 1,
      CHR == "2R" ~ 2,
      CHR == "3L" ~ 3,
      CHR == "3R" ~ 4,
      CHR == "4"  ~ 5,
      CHR == "X"  ~ 6,
      CHR == "Y"  ~ 7,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(CHR), !is.na(P), !is.na(BP))

png("manhattan_gwas_fly_female.png", width = 1600, height = 1200, res = 300)
manhattan(
  gwas_fly_f_plot,
  chr = "CHR", bp = "BP", p = "P", snp = "SNP",
  main = "Fly Female GWAS - Starvation Stress",
  col = c("steelblue", "coral"),
  suggestiveline = -log10(1e-5),
  genomewideline = -log10(0.05 / nrow(gwas_fly_f_plot)),
  cex = 0.6,
  chrlabs = c("2L", "2R", "3L", "3R", "4", "X", "Y")
)
dev.off()
cat("  Saved: manhattan_gwas_fly_female.png\n")

png("qq_gwas_fly_female.png", width = 1600, height = 1200, res = 300)
qq(gwas_fly_f_plot$P, main = "Fly Female GWAS QQ Plot")
dev.off()
cat("  Saved: qq_gwas_fly_female.png\n")


# ==============================================================================
# ==============================================================================
# 2. MAGMA Gene Manhattan Plots
# ==============================================================================
# ==============================================================================


cat("\n========================================\n")
cat("MAGMA Manhattan Plots (Gene-level)\n")
cat("========================================\n")

# --- 2a. Arabidopsis MAGMA ---
cat("\nLoading Arabidopsis MAGMA...\n")

magma_at_files <- list.files(
  path = paste0(
    "/Users/nirwantandukar/Documents/Research/results/Arabidopsis/",
    "MAGMA/AT_cold_by_chr"
  ),
  pattern = "\\.genes\\.out$",
  full.names = TRUE
)

magma_at_list <- lapply(magma_at_files, function(f) {
  read.table(f, header = TRUE, stringsAsFactors = FALSE, comment.char = "#")
})
magma_at <- do.call(rbind, magma_at_list)

# Use P_MULTI or P column
if ("P_MULTI" %in% names(magma_at)) {
  magma_at$P <- magma_at$P_MULTI
} else if (!"P" %in% names(magma_at)) {
  magma_at$P <- magma_at[[9]]
}

# Use gene midpoint as BP
magma_at$BP <- (magma_at$START + magma_at$STOP) / 2

# Remove duplicates
magma_at <- magma_at %>%
  group_by(GENE) %>%
  slice_min(P, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("  Total genes:", nrow(magma_at), "\n")

png("manhattan_magma_arabidopsis.png", width = 1600, height = 1200, res = 300)
manhattan(
  magma_at,
  chr = "CHR", bp = "BP", p = "P", snp = "GENE",
  main = "Arabidopsis MAGMA Gene Analysis - Cold Temperature",
  col = c("steelblue", "coral"),
  suggestiveline = -log10(1e-4),
  genomewideline = -log10(0.05 / nrow(magma_at)),
  cex = 0.8,
  cex.axis = 0.8
)
dev.off()
cat("  Saved: manhattan_magma_arabidopsis.png\n")

png("qq_magma_arabidopsis.png", width = 1600, height = 1200, res = 300)
qq(magma_at$P, main = "Arabidopsis MAGMA QQ Plot")
dev.off()
cat("  Saved: qq_magma_arabidopsis.png\n")

# --- 2b. Fly Male MAGMA ---
cat("\nLoading Fly Male MAGMA...\n")

magma_fly_m_files <- list.files(
  path = paste0(
    "/Users/nirwantandukar/Documents/Research/results/DGRP/MAGMA/",
    "Fly_magma_genes_by_chr_male"
  ),
  pattern = "\\.genes\\.out$",
  full.names = TRUE
)

magma_fly_m_list <- lapply(magma_fly_m_files, function(f) {
  read.table(f, header = TRUE, stringsAsFactors = FALSE, comment.char = "#")
})
magma_fly_m <- do.call(rbind, magma_fly_m_list)

if ("P_MULTI" %in% names(magma_fly_m)) {
  magma_fly_m$P <- magma_fly_m$P_MULTI
} else if (!"P" %in% names(magma_fly_m)) {
  magma_fly_m$P <- magma_fly_m[[9]]
}

magma_fly_m$BP <- (magma_fly_m$START + magma_fly_m$STOP) / 2

magma_fly_m <- magma_fly_m %>%
  group_by(GENE) %>%
  slice_min(P, n = 1, with_ties = FALSE) %>%
  ungroup()

magma_fly_m_plot = magma_fly_m

cat("  Total genes:", nrow(magma_fly_m), "\n")

# Convert chr names
# magma_fly_m_plot <- magma_fly_m %>%
#   mutate(
#     CHR_orig = CHR,
#     CHR = case_when(
#       CHR == "2L" ~ 1,
#       CHR == "2R" ~ 2,
#       CHR == "3L" ~ 3,
#       CHR == "3R" ~ 4,
#       CHR == "4"  ~ 5,
#       CHR == "X"  ~ 6,
#       CHR == "Y"  ~ 7,
#       TRUE ~ NA_real_
#     )
#   ) %>%
#   filter(!is.na(CHR), !is.na(P), !is.na(BP))

png("manhattan_magma_fly_male.png", width = 1600, height = 1200, res = 300)
manhattan(
  magma_fly_m_plot,
  chr = "CHR", bp = "BP", p = "P", snp = "GENE",
  main = "Fly Male MAGMA Gene Analysis - Starvation Stress",
  col = c("steelblue", "coral"),
  suggestiveline = -log10(1e-4),
  genomewideline = -log10(0.05 / nrow(magma_fly_m_plot)),
  cex = 0.8,
  chrlabs = c("2L", "2R", "3L", "3R", "4", "X", "Y")
)
dev.off()
cat("  Saved: manhattan_magma_fly_male.png\n")

png("qq_magma_fly_male.png", width = 1600, height = 1200, res = 300)
qq(magma_fly_m_plot$P, main = "Fly Male MAGMA QQ Plot")
dev.off()
cat("  Saved: qq_magma_fly_male.png\n")

# --- 2c. Fly Female MAGMA ---
cat("\nLoading Fly Female MAGMA...\n")

magma_fly_f_files <- list.files(
  path = paste0(
    "/Users/nirwantandukar/Documents/Research/results/DGRP/MAGMA/",
    "Fly_magma_genes_by_chr_female"
  ),
  pattern = "\\.genes\\.out$",
  full.names = TRUE
)

magma_fly_f_list <- lapply(magma_fly_f_files, function(f) {
  read.table(f, header = TRUE, stringsAsFactors = FALSE, comment.char = "#")
})
magma_fly_f <- do.call(rbind, magma_fly_f_list)

if ("P_MULTI" %in% names(magma_fly_f)) {
  magma_fly_f$P <- magma_fly_f$P_MULTI
} else if (!"P" %in% names(magma_fly_f)) {
  magma_fly_f$P <- magma_fly_f[[9]]
}

magma_fly_f$BP <- (magma_fly_f$START + magma_fly_f$STOP) / 2

magma_fly_f <- magma_fly_f %>%
  group_by(GENE) %>%
  slice_min(P, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("  Total genes:", nrow(magma_fly_f), "\n")

magma_fly_f_plot = magma_fly_f

# magma_fly_f_plot <- magma_fly_f %>%
#   mutate(
#     CHR_orig = CHR,
#     CHR = case_when(
#       CHR == "2L" ~ 1,
#       CHR == "2R" ~ 2,
#       CHR == "3L" ~ 3,
#       CHR == "3R" ~ 4,
#       CHR == "4"  ~ 5,
#       CHR == "X"  ~ 6,
#       CHR == "Y"  ~ 7,
#       TRUE ~ NA_real_
#     )
#   ) %>%
#   filter(!is.na(CHR), !is.na(P), !is.na(BP))

png("manhattan_magma_fly_female.png", width = 1600, height = 1200, res = 300)
manhattan(
  magma_fly_f_plot,
  chr = "CHR", bp = "BP", p = "P", snp = "GENE",
  main = "Fly Female MAGMA Gene Analysis - Starvation Stress",
  col = c("steelblue", "coral"),
  suggestiveline = -log10(1e-4),
  genomewideline = -log10(0.05 / nrow(magma_fly_f_plot)),
  cex = 0.8,
  chrlabs = c("2L", "2R", "3L", "3R", "4", "X", "Y")
)
dev.off()
cat("  Saved: manhattan_magma_fly_female.png\n")

png("qq_magma_fly_female.png", width = 1600, height = 1200, res = 300)
qq(magma_fly_f_plot$P, main = "Fly Female MAGMA QQ Plot")
dev.off()
cat("  Saved: qq_magma_fly_female.png\n")

# ==============================================================================
# 3. Summary Statistics
# ==============================================================================

cat("\n========================================\n")
cat("Summary Statistics\n")
cat("========================================\n")

cat("\nGWAS SNP Counts:\n")
cat("  Arabidopsis:", nrow(gwas_at), "SNPs\n")
cat("  Fly Male:", nrow(gwas_fly_m), "SNPs\n")
cat("  Fly Female:", nrow(gwas_fly_f), "SNPs\n")

cat("\nMAGMA Gene Counts:\n")
cat("  Arabidopsis:", nrow(magma_at), "genes\n")
cat("  Fly Male:", nrow(magma_fly_m), "genes\n")
cat("  Fly Female:", nrow(magma_fly_f), "genes\n")

cat("\nGWAS Bonferroni Thresholds (p = 0.05/n):\n")
cat("  Arabidopsis:", format(0.05 / nrow(gwas_at), digits = 3), "\n")
cat("  Fly Male:", format(0.05 / nrow(gwas_fly_m), digits = 3), "\n")
cat("  Fly Female:", format(0.05 / nrow(gwas_fly_f), digits = 3), "\n")

cat("\nMAGMA Bonferroni Thresholds (p = 0.05/n):\n")
cat("  Arabidopsis:", format(0.05 / nrow(magma_at), digits = 3), "\n")
cat("  Fly Male:", format(0.05 / nrow(magma_fly_m), digits = 3), "\n")
cat("  Fly Female:", format(0.05 / nrow(magma_fly_f), digits = 3), "\n")

cat("\nSignificant SNPs (Bonferroni):\n")
cat("  Arabidopsis:", sum(gwas_at$P < 0.05 / nrow(gwas_at), na.rm = TRUE), "\n")
cat("  Fly Male:", sum(gwas_fly_m$P < 0.05 / nrow(gwas_fly_m), na.rm = TRUE), "\n")
cat("  Fly Female:", sum(gwas_fly_f$P < 0.05 / nrow(gwas_fly_f), na.rm = TRUE), "\n")

cat("\nSignificant Genes (Bonferroni):\n")
cat("  Arabidopsis:", sum(magma_at$P < 0.05 / nrow(magma_at), na.rm = TRUE), "\n")
cat("  Fly Male:",
    sum(magma_fly_m$P < 0.05 / nrow(magma_fly_m), na.rm = TRUE), "\n")
cat("  Fly Female:",
    sum(magma_fly_f$P < 0.05 / nrow(magma_fly_f), na.rm = TRUE), "\n")

cat("\nSuggestive SNPs (p < 1e-5):\n")
cat("  Arabidopsis:", sum(gwas_at$P < 1e-5, na.rm = TRUE), "\n")
cat("  Fly Male:", sum(gwas_fly_m$P < 1e-5, na.rm = TRUE), "\n")
cat("  Fly Female:", sum(gwas_fly_f$P < 1e-5, na.rm = TRUE), "\n")

cat("\nSuggestive Genes (p < 1e-4):\n")
cat("  Arabidopsis:", sum(magma_at$P < 1e-4, na.rm = TRUE), "\n")
cat("  Fly Male:", sum(magma_fly_m$P < 1e-4, na.rm = TRUE), "\n")
cat("  Fly Female:", sum(magma_fly_f$P < 1e-4, na.rm = TRUE), "\n")

cat("\n========== DONE ==========\n")
cat("Output files:\n")
cat("  GWAS Manhattan:\n")
cat("    - manhattan_gwas_arabidopsis.png\n")
cat("    - manhattan_gwas_fly_male.png\n")
cat("    - manhattan_gwas_fly_female.png\n")
cat("  GWAS QQ:\n")
cat("    - qq_gwas_arabidopsis.png\n")
cat("    - qq_gwas_fly_male.png\n")
cat("    - qq_gwas_fly_female.png\n")
cat("  MAGMA Manhattan:\n")
cat("    - manhattan_magma_arabidopsis.png\n")
cat("    - manhattan_magma_fly_male.png\n")
cat("    - manhattan_magma_fly_female.png\n")
cat("  MAGMA QQ:\n")
cat("    - qq_magma_arabidopsis.png\n")
cat("    - qq_magma_fly_male.png\n")
cat("    - qq_magma_fly_female.png\n")
