suppressPackageStartupMessages({
  source("/Users/nirwantandukar/Documents/Github/MAGCAT/all_figs.R")
})

# ------------------------------------------------------------------------------
# Final Figure 5 script
# Omnibus power regret across effect sizes and archetypes
# ------------------------------------------------------------------------------

dir.create("/Users/nirwantandukar/Documents/Github/MAGCAT/Final/final_sup/supp_figs",
           recursive = TRUE, showWarnings = FALSE)

file.copy(
  from = "/Users/nirwantandukar/Documents/Github/MAGCAT/simulation_results/block_h_omnibus_regret.png",
  to = "/Users/nirwantandukar/Documents/Github/MAGCAT/Final/final_sup/supp_figs/SuppFig5_block_h_omnibus_regret.png",
  overwrite = TRUE
)

message("SuppFig5 output written to Final/final_sup/supp_figs/SuppFig5_block_h_omnibus_regret.png")
