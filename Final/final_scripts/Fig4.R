suppressPackageStartupMessages({
  source("/Users/nirwantandukar/Documents/Github/MAGCAT/all_figs.R")
  library(patchwork)
})

# ------------------------------------------------------------------------------
# Final Figure 4 script
# Combined null-calibration and missing-correlation calibration summary
# Panels:
#   A. Null calibration lambda under LD
#   B. Null calibration type I error under LD
#   C. Missing-correlation lambda
#   D. Missing-correlation type I error
# ------------------------------------------------------------------------------

results_dir <- "/Users/nirwantandukar/Documents/Github/MAGCAT/simulation_results"
final_fig_dir <- "/Users/nirwantandukar/Documents/Github/MAGCAT/Final/final_main_fig"
dir.create(final_fig_dir, recursive = TRUE, showWarnings = FALSE)

env <- get_legacy_sim_env()

block_b <- readRDS(file.path(results_dir, "block_b.rds"))
plots_b <- env$plot_block_a_null(block_b)
plots_b$lambda <- retitle_block_plot(plots_b$lambda, "Block A", "Block B")
plots_b$type1  <- retitle_block_plot(plots_b$type1, "Block A", "Block B")
plots_b$lambda <- strip_geom_text_layers(plots_b$lambda)
plots_b$type1  <- strip_geom_text_layers(plots_b$type1)
plots_b$lambda <- strip_plot_headers(plots_b$lambda) +
  facet_grid(pathway_size ~ cor_structure,
             scales = "free_y",
             labeller = labeller(pathway_size = function(x) paste0("m=", x))) +
  labs(y = expression(lambda), x = NULL) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
    plot.margin = margin(5, 5, 5, 5)
  )
plots_b$type1 <- strip_plot_headers(plots_b$type1) +
  facet_grid(pathway_size ~ cor_structure,
             scales = "free_y",
             labeller = labeller(pathway_size = function(x) paste0("m=", x))) +
  labs(y = "Type I error", x = "Method") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
    plot.margin = margin(5, 5, 5, 5)
  )

block_c <- readRDS(file.path(results_dir, "block_c.rds"))
plots_c <- env$plot_block_c(block_c)
plots_c$lambda <- strip_plot_headers(plots_c$lambda) +
  labs(y = expression(lambda), x = NULL) +
  theme(
    legend.position = "top",
    plot.margin = margin(5, 5, 5, 5)
  )
plots_c$type1 <- strip_plot_headers(plots_c$type1) +
  labs(y = "Type I error", x = "Fraction missing") +
  theme(
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  )

combined_plot <-
  (plots_b$lambda / plots_b$type1 / plots_c$lambda / plots_c$type1) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 18))

ggsave(
  filename = file.path(final_fig_dir, "Fig4.png"),
  plot = combined_plot,
  width = 14,
  height = 22,
  dpi = 300,
  bg = "white"
)

message("Fig4 output written to Final/final_main_fig/Fig4.png")
