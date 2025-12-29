# ## ============================================================
# ## Build Drosophila pathway long tables (FlyCyc + FlyBase)
# ## Output format: pathway_id, pathway_name, gene
# ## ============================================================

# suppressPackageStartupMessages({
#   library(readr)
#   library(dplyr)
#   library(tidyr)
#   library(stringr)
# })


# ### GET PATHWAY INFORMATION FROM BIOCYC
# ## ---- inputs (your files) -----------------------------------
# flycyc_csv <- "inst/extdata/pathway/fly_cyc_unedited.csv"
# flybase_met_tsv <- "inst/extdata/pathway/metabolic_pathway_group_data_fb_2025_05.tsv"
# flybase_sig_tsv <- "inst/extdata/pathway/signaling_pathway_group_data_fb_2025_05.tsv"


# x <- read_csv(flycyc_csv, show_col_types = FALSE)
# # Robustly pick columns:
# # - gene column must be exactly "Genes" OR contain "Gene"
# # - pathway id must contain "ID" (e.g., Pathway_ID)
# # - pathway name is the remaining pathway-ish column (e.g., Pathway_name)
# gene_col <- names(x)[tolower(names(x)) %in% c("gene-name", "gene", "genes")]
# if (length(gene_col) == 0) {
#   gene_col <- names(x)[stringr::str_detect(tolower(names(x)), "gene")]
# }

# pid_col  <- names(x)[stringr::str_detect(tolower(names(x)), "pathway") &
#                      stringr::str_detect(tolower(names(x)), "id")]

# pname_col <- setdiff(
#   names(x)[stringr::str_detect(tolower(names(x)), "pathway")],
#   pid_col
# )

# stopifnot(length(gene_col)  == 1)
# stopifnot(length(pid_col)   == 1)
# stopifnot(length(pname_col) == 1)

# flycyc_long <- x %>%
#   transmute(
#     pathway_name = .data[[pname_col]],
#     pathway_id   = .data[[pid_col]],
#     gene         = .data[[gene_col]]
#   ) %>%
#   tidyr::fill(pathway_id, pathway_name, .direction = "down") %>%
#   mutate(
#     pathway_id   = paste0("FLYCYC:", str_trim(as.character(pathway_id))),
#     pathway_name = str_trim(as.character(pathway_name)),
#     gene         = str_trim(as.character(gene))
#   ) %>%
#   filter(!is.na(gene), gene != "") %>%
#   distinct()
# flycyc_long

# write_tsv(flycyc_long, "Fly_Cyc.tsv")


# getwd()
