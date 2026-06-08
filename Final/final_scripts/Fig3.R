suppressPackageStartupMessages({
  source("/Users/nirwantandukar/Documents/Github/MAGCAT/all_figs.R")
})

# ------------------------------------------------------------------------------
# Final Figure 3 script
# Pathway-size sensitivity of archetype verification for CATFISH component tests
# Produces:
#   - main figure copy
#   - compact supplementary figure
#   - compact supplementary table
#   - full supplementary value table
# ------------------------------------------------------------------------------

write_block_a_by_m_compact(
  summary_csv = "/Users/nirwantandukar/Documents/Github/MAGCAT/simulation_results/block_a_archetype_summary_by_m.csv",
  results_dir = "/Users/nirwantandukar/Documents/Github/MAGCAT/Final/final_main_fig",
  supp_dir = "/Users/nirwantandukar/Documents/Github/MAGCAT/Final/final_sup"
)

file.copy(
  from = "/Users/nirwantandukar/Documents/Github/MAGCAT/Final/final_main_fig/block_a_archetype_recovery_by_m_compact.png",
  to = "/Users/nirwantandukar/Documents/Github/MAGCAT/Final/final_main_fig/Fig3.png",
  overwrite = TRUE
)

write_block_a_by_m_tables(
  summary_csv = "/Users/nirwantandukar/Documents/Github/MAGCAT/simulation_results/block_a_archetype_summary_by_m.csv",
  tables_dir = "/Users/nirwantandukar/Documents/Github/MAGCAT/Final/final_sup/tables"
)

message("Fig3 outputs written to Final/final_main_fig and Final/final_sup/tables")
