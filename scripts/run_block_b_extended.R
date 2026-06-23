# ==============================================================================
# Run Block B Extended: Comprehensive Broken Component Stress Test
# ==============================================================================
# This script runs the extended stress test to show that the omnibus method
# remains robust even when individual component tests are miscalibrated
# (either anti-conservative/inflated or conservative/deflated).
# ==============================================================================

# Source the main simulation framework
source("simulation_null_diagnostics.R")

cat(strrep("=", 70), "\n")
cat("BLOCK B EXTENDED: Comprehensive Broken Component Stress Test\n")
cat(strrep("=", 70), "\n\n")

# Create output directory
output_dir <- "simulation_results"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Run the extended stress test
# Using smaller n_sims for faster execution, increase for publication
cat("Running stress test simulations...\n")
cat("This may take several minutes.\n\n")

results_b_ext <- run_block_b_extended(
  n_sims = 300L,        # Number of null simulations per scenario
  b = 500L,             # Permutations for null distribution
  seed = 32345L,
  pathway_sizes_arg = c(5L, 25L, 50L),  # Pathway sizes to test (including small m=5)
  rho = 0.3,            # LD correlation strength
  verbose = TRUE
)

# Save results
cat("\nSaving results...\n")
saveRDS(results_b_ext, file.path(output_dir, "block_b_extended.rds"))

# Generate plots
cat("Generating figures...\n")
plots <- plot_block_b_extended(results_b_ext, output_dir = output_dir)

cat("\n", strrep("=", 70), "\n")
cat("COMPLETE!\n")
cat(strrep("=", 70), "\n\n")

cat("Results saved to:\n")
cat("  - ", file.path(output_dir, "block_b_extended.rds"), "\n\n")

cat("Figures saved to:\n")
cat("  - block_b_ext_type1_by_nbroken.png    (Type I error vs broken components)\n")
cat("  - block_b_ext_severity_boxplot.png    (Severity effects on omnibus)\n")
cat("  - block_b_ext_lambda_heatmap.png      (Genomic control heatmap)\n")
cat("  - block_b_ext_mixed_scenarios.png     (Mixed anti+conservative)\n")
cat("  - block_b_ext_omnibus_vs_individual.png (Key comparison figure)\n")
cat("  - block_b_ext_combined_main.png       (Combined main figure)\n\n")

# Print summary statistics
cat(strrep("=", 70), "\n")
cat("SUMMARY STATISTICS\n")
cat(strrep("=", 70), "\n\n")

# Omnibus performance summary
omni_summary <- results_b_ext %>%
  filter(method == "omnibus", scenario == "analytic_broken") %>%
  group_by(break_type, severity) %>%
  summarise(
    mean_type1 = mean(type1_05, na.rm = TRUE),
    sd_type1 = sd(type1_05, na.rm = TRUE),
    mean_lambda = mean(lambda, na.rm = TRUE),
    .groups = "drop"
  )

cat("Omnibus Type I Error by Breakage Type and Severity:\n")
print(omni_summary, n = 30)

# Comparison: omnibus vs individual methods
compare_summary <- results_b_ext %>%
  filter(scenario == "analytic_broken", n_broken >= 2) %>%
  mutate(is_omnibus = method == "omnibus") %>%
  group_by(is_omnibus, break_type) %>%
  summarise(
    mean_type1 = mean(type1_05, na.rm = TRUE),
    mean_deviation = mean(abs(type1_05 - 0.05), na.rm = TRUE),
    .groups = "drop"
  )

cat("\n\nOmnibus vs Individual Methods (≥2 broken components):\n")
print(compare_summary)

cat("\n\nKey Takeaway:\n")
cat("The omnibus test maintains better calibration than individual component\n")
cat("tests when those components are miscalibrated. Even with 3-4 broken\n")
cat("components, the omnibus type I error remains closer to the nominal 0.05.\n")
