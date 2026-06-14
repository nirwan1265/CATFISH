suppressPackageStartupMessages({
  source("/Users/nirwantandukar/Documents/Github/MAGCAT/all_figs.R")
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
  library(stringr)
  library(scales)
})

# ------------------------------------------------------------------------------
# Final Figure 6 script
# Dry tons per acre: GWAS + MAGMA + CATFISH pathway integration
# Layout:
#   A. GWAS Manhattan plot
#   B. MAGMA gene Manhattan plot
#   C. CATFISH pathway heatmap
# ------------------------------------------------------------------------------

final_fig_dir <- "/Users/nirwantandukar/Documents/Github/MAGCAT/Final/final_main_fig"
dir.create(final_fig_dir, recursive = TRUE, showWarnings = FALSE)

gwas_dir <- "/Users/nirwantandukar/Documents/Research/results/GWAS/MLM/BAP/Dry_tons_per_acre"
magma_dir <- "/Users/nirwantandukar/Documents/Research/results/CATFISH/MAGMA/Dry_tons_per_acre/CATFISH_permutation_B1000000_mvn_GPD_strict_tau"
magma_file <- file.path(magma_dir, "Dry_tons_per_acre_combined_genes.tsv")
pathway_file <- file.path(magma_dir, "Dry_tons_per_acre_CATFISH_ACAT_mvn_B1000000_GPD_strict_tau.csv")

if (!dir.exists(gwas_dir)) stop("GWAS dir not found: ", gwas_dir, call. = FALSE)
if (!file.exists(magma_file)) stop("MAGMA file not found: ", magma_file, call. = FALSE)
if (!file.exists(pathway_file)) stop("Pathway file not found: ", pathway_file, call. = FALSE)

panel_theme <- plot_theme +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0),
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    plot.margin = margin(6, 8, 6, 8)
  )

odd_chr_col <- "#2F3E46"
even_chr_col <- "#9AA5B1"
sig_col <- "#D62728"
subsig_col <- "#2A6FDB"

build_chr_positions <- function(df, chr_col, pos_col) {
  df %>%
    group_by(.data[[chr_col]]) %>%
    summarise(chr_len = max(.data[[pos_col]], na.rm = TRUE), .groups = "drop") %>%
    arrange(.data[[chr_col]]) %>%
    mutate(
      cum_start = lag(cumsum(chr_len), default = 0),
      chr_mid = cum_start + chr_len / 2
    )
}

read_all_gwas <- function(dir_path) {
  files <- list.files(dir_path, pattern = "\\.assoc\\.txt$", full.names = TRUE)
  chr_num <- as.integer(sub(".*Chr([0-9]+)\\.maf01\\.assoc\\.txt$", "\\1", basename(files)))
  files <- files[order(chr_num)]

  pieces <- lapply(files, function(f) {
    x <- fread(f, header = TRUE)
    setnames(x, old = c("chr", "rs", "ps", "p_wald"), new = c("CHR", "SNP_ID", "POS", "P"), skip_absent = TRUE)
    x
  })

  rbindlist(pieces, use.names = TRUE, fill = TRUE)
}

# ------------------------------------------------------------------------------
# Panel A: GWAS Manhattan
# ------------------------------------------------------------------------------
gwas_raw <- read_all_gwas(gwas_dir)

gwas_plot_data <- gwas_raw %>%
  filter(!is.na(P), P > 0, !is.na(POS), !is.na(CHR)) %>%
  mutate(
    CHR = as.integer(CHR),
    POS = as.numeric(POS),
    logP = -log10(P)
  ) %>%
  filter(is.finite(logP), is.finite(POS), is.finite(CHR))

gwas_chr <- build_chr_positions(gwas_plot_data, "CHR", "POS")
gwas_plot_data <- gwas_plot_data %>%
  left_join(gwas_chr %>% dplyr::select(CHR, cum_start), by = "CHR") %>%
  mutate(
    cum_pos = POS + cum_start,
    chr_group = if_else(CHR %% 2L == 0L, "even", "odd")
  )

gwas_plot_sample <- gwas_plot_data

gwas_threshold <- 7
suggestive_threshold <- 5

# ------------------------------------------------------------------------------
# Panel B: MAGMA Manhattan
# ------------------------------------------------------------------------------
magma_gene <- read.delim(magma_file, stringsAsFactors = FALSE, check.names = FALSE) %>%
  mutate(
    CHR = as.integer(CHR),
    START = as.numeric(START),
    STOP = as.numeric(STOP),
    gene_mid = (START + STOP) / 2,
    logP = -log10(P)
  ) %>%
  filter(is.finite(CHR), is.finite(gene_mid), is.finite(logP))

magma_chr <- build_chr_positions(magma_gene, "CHR", "STOP")
magma_plot_data <- magma_gene %>%
  left_join(magma_chr %>% dplyr::select(CHR, cum_start), by = "CHR") %>%
  mutate(
    cum_pos = gene_mid + cum_start,
    chr_group = if_else(CHR %% 2L == 0L, "even", "odd")
  )

magma_bonf_threshold <- -log10(0.05 / nrow(magma_plot_data))
magma_fdr <- magma_plot_data %>%
  mutate(fdr = p.adjust(P, method = "BH"))
magma_fdr_threshold <- magma_fdr %>%
  filter(fdr < 0.05) %>%
  summarise(th = max(-log10(P), na.rm = TRUE)) %>%
  pull(th)
if (!length(magma_fdr_threshold) || !is.finite(magma_fdr_threshold)) {
  magma_fdr_threshold <- NA_real_
}

common_y_max <- ceiling(max(c(gwas_plot_sample$logP, magma_plot_data$logP), na.rm = TRUE))
common_y_max <- max(common_y_max, 12)

panel_a <- ggplot(gwas_plot_sample, aes(x = cum_pos, y = logP, color = chr_group)) +
  geom_point(alpha = 0.72, size = 0.75, stroke = 0) +
  geom_hline(yintercept = gwas_threshold, color = sig_col, linetype = "dashed", linewidth = 0.45) +
  geom_hline(yintercept = suggestive_threshold, color = subsig_col, linetype = "dashed", linewidth = 0.45) +
  scale_color_manual(values = c(odd = odd_chr_col, even = even_chr_col)) +
  scale_x_continuous(
    breaks = gwas_chr$chr_mid,
    labels = gwas_chr$CHR,
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(
    limits = c(0, common_y_max),
    breaks = pretty_breaks(n = 5)
  ) +
  labs(
    title = "GWAS SNP-Level Association",
    x = "Chromosome",
    y = expression(-log[10](italic(p)))
  ) +
  panel_theme +
  theme(legend.position = "none")

panel_b <- ggplot(magma_plot_data, aes(x = cum_pos, y = logP, color = chr_group)) +
  geom_point(alpha = 0.8, size = 0.85, stroke = 0) +
  geom_hline(yintercept = magma_bonf_threshold, color = sig_col, linetype = "dashed", linewidth = 0.45) +
  {if (is.finite(magma_fdr_threshold)) geom_hline(yintercept = magma_fdr_threshold, color = subsig_col, linetype = "dashed", linewidth = 0.45)} +
  scale_color_manual(values = c(odd = odd_chr_col, even = even_chr_col)) +
  scale_x_continuous(
    breaks = magma_chr$chr_mid,
    labels = magma_chr$CHR,
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(
    limits = c(0, common_y_max),
    breaks = pretty_breaks(n = 5)
  ) +
  labs(
    title = "MAGMA Gene-Level Association",
    x = "Chromosome",
    y = expression(-log[10](italic(p)))
  ) +
  panel_theme +
  theme(legend.position = "none")

# ------------------------------------------------------------------------------
# Panel C: CATFISH pathway heatmap
# ------------------------------------------------------------------------------
pathway_results <- fread(pathway_file, data.table = FALSE)

pick_component_col <- function(df, candidates) {
  for (col in candidates) {
    if (!col %in% names(df)) next
    vals <- suppressWarnings(as.numeric(df[[col]]))
    if (any(is.finite(vals) & !is.na(vals))) return(col)
  }
  stop("Could not find a populated component column among: ", paste(candidates, collapse = ", "), call. = FALSE)
}

acat_col <- pick_component_col(pathway_results, c("acat_p_mvn_cal", "acat_p"))
fisher_col <- pick_component_col(pathway_results, c("fisher_p_mvn_cal", "fisher_p"))
stouffer_col <- pick_component_col(pathway_results, c("stouffer_p_mvn_cal", "stouffer_p_analytic"))
minp_col <- pick_component_col(pathway_results, c("minp_p_mvn_cal", "minp_p_analytic"))
tfisher_col <- pick_component_col(pathway_results, c("tfisher_p_mvn_cal", "tfisher_p_analytic"))
omni_col <- pick_component_col(pathway_results, c("omni_p_final", "omni_p_mvn", "omni_p_analytic"))

top_pathways <- pathway_results %>%
  arrange(.data[[omni_col]]) %>%
  slice_head(n = 18) %>%
  mutate(
    pathway_label = paste0(str_wrap(pathway_name, width = 38), " [", n_genes, "]"),
    pathway_label = factor(pathway_label, levels = rev(pathway_label))
  )

heatmap_df <- top_pathways %>%
  transmute(
    pathway_label,
    ACAT = .data[[acat_col]],
    Fisher = .data[[fisher_col]],
    Stouffer = .data[[stouffer_col]],
    minP = .data[[minp_col]],
    TFisher = .data[[tfisher_col]],
    Omnibus = .data[[omni_col]]
  ) %>%
  pivot_longer(
    cols = c("ACAT", "Fisher", "Stouffer", "minP", "TFisher", "Omnibus"),
    names_to = "Method",
    values_to = "p_value"
  ) %>%
  mutate(
    logP = -log10(pmax(p_value, 1e-50)),
    Method = factor(Method, levels = c("ACAT", "Fisher", "Stouffer", "minP", "TFisher", "Omnibus"))
  )

heat_cap <- max(12, ceiling(stats::quantile(heatmap_df$logP, probs = 0.98, na.rm = TRUE)))
heatmap_df$logP_cap <- pmin(heatmap_df$logP, heat_cap)

panel_c <- ggplot(heatmap_df, aes(x = Method, y = pathway_label, fill = logP_cap)) +
  geom_tile(color = "white", linewidth = 0.45) +
  scale_fill_gradientn(
    colours = c("#EEF3F8", "#8FB3D9", "#C65A46", "#7F0000"),
    values = rescale(c(0, heat_cap * 0.25, heat_cap * 0.65, heat_cap)),
    limits = c(0, heat_cap),
    name = expression(-log[10](italic(p)))
  ) +
  labs(
    title = "CATFISH Pathway Enrichment",
    x = NULL,
    y = NULL
  ) +
  panel_theme +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, face = "bold"),
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9)
  )

combined_plot <-
  ((panel_a | panel_b) / panel_c) +
  plot_layout(heights = c(0.92, 1.28)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 18))

ggsave(
  filename = file.path(final_fig_dir, "Fig6.png"),
  plot = combined_plot,
  width = 14,
  height = 12,
  dpi = 300,
  bg = "white"
)

message("Fig6 output written to Final/final_main_fig/Fig6.png")
