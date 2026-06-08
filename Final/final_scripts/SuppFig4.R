suppressPackageStartupMessages({
  source("/Users/nirwantandukar/Documents/Github/MAGCAT/all_figs.R")
})

# ------------------------------------------------------------------------------
# Final Figure 4 script
# Effect-size sensitivity of CATFISH component and omnibus tests across archetypes
# ------------------------------------------------------------------------------

dir.create("/Users/nirwantandukar/Documents/Github/MAGCAT/Final/final_sup/supp_figs",
           recursive = TRUE, showWarnings = FALSE)

file.copy(
  from = "/Users/nirwantandukar/Documents/Github/MAGCAT/simulation_results/block_h_power.png",
  to = "/Users/nirwantandukar/Documents/Github/MAGCAT/Final/final_sup/supp_figs/SuppFig4_block_h_power.png",
  overwrite = TRUE
)

message("SuppFig4 output written to Final/final_sup/supp_figs/SuppFig4_block_h_power.png")
