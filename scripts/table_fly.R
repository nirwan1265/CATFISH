# library(dplyr)

# # Get the files
# files_male <- list.files(path = "/Users/nirwantandukar/Documents/Research/results/DGRP/MAGMA/Fly_magma_genes_by_chr_male",
#         pattern = "^Male_starvation_fly_.*\\.genes\\.out$",
#         full.names = TRUE)
# files_female <- list.files(path = "/Users/nirwantandukar/Documents/Research/results/DGRP/MAGMA/Fly_magma_genes_by_chr_female",
#         pattern = "^Female_starvation_fly_.*\\.genes\\.out$",
#         full.names = TRUE)


# # Read all files
# gene_list_male <- lapply(files_male, function(f) {
#   if (!file.exists(f)) stop("File not found: ", f)
#   utils::read.table(f, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
# })

# gene_list_female <- lapply(files_female, function(f) {
#   if (!file.exists(f)) stop("File not found: ", f)
#   utils::read.table(f, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
# })


# # Stack into one big data.frame
# genes_all_raw_male <- do.call(rbind, gene_list_male)
# genes_all_raw_female <- do.call(rbind, gene_list_female)

# colnames(genes_all_raw_male)[9]="P"
# colnames(genes_all_raw_female)[9]="P"

# # ---- (Optional but recommended) deduplicate by GENE: keep smallest P ----
# if (!"GENE" %in% names(genes_all_raw_male) || !"P" %in% names(genes_all_raw_male)) {
#   stop("Combined MAGMA file must have columns 'GENE' and 'P'.")
# }

# if (!"GENE" %in% names(genes_all_raw_female) || !"P" %in% names(genes_all_raw_female)) {
#   stop("Combined MAGMA file must have columns 'GENE' and 'P'.")
# }

# head(genes_all_raw_male)
# head(genes_all_raw_female)

# # Order
# genes_all_raw_male <- genes_all_raw_male %>%
#   mutate(P = as.numeric(P)) %>%
#   arrange(P)

# genes_all_raw_female <- genes_all_raw_female %>%
#   mutate(P = as.numeric(P)) %>%
#   arrange(P)

# # check
# head(genes_all_raw_male)
# head(genes_all_raw_female)

# write.csv(genes_all_raw_male,"Tables/MAGMA_genes_all_raw_male_fly.csv", row.names=F)
# write.csv(genes_all_raw_female,"Tables/MAGMA_genes_all_raw_female_fly.csv", row.names=F)




# # Fig VENN

# ### TOP 20
# library(tidyverse)
# library(ggVennDiagram)

# # -------------------------
# # 1) Venn of ALL "sig" genes (P < 0.05)
# # -------------------------
# alpha <- 0.001

# male_sig <- genes_all_raw_male %>%
#   mutate(P = as.numeric(P)) %>%
#   filter(!is.na(P), P < alpha) %>%
#   pull(GENE) %>%
#   unique()

# female_sig <- genes_all_raw_female %>%
#   mutate(P = as.numeric(P)) %>%
#   filter(!is.na(P), P < alpha) %>%
#   pull(GENE) %>%
#   unique()

# p_venn_sig <- ggVennDiagram(list(Male = male_sig, Female = female_sig), label = "count") +
#   ggtitle(paste0("Sig genes overlap (P < ", alpha, ")"))

# # -------------------------
# # 2) Venn of TOP 20 genes (by smallest P)
# # -------------------------
# male_top20 <- genes_all_raw_male %>%
#   mutate(P = as.numeric(P)) %>%
#   filter(!is.na(P)) %>%
#   arrange(P) %>%
#   slice_head(n = 50) %>%
#   pull(GENE) %>%
#   unique()

# female_top20 <- genes_all_raw_female %>%
#   mutate(P = as.numeric(P)) %>%
#   filter(!is.na(P)) %>%
#   arrange(P) %>%
#   slice_head(n = 50) %>%
#   pull(GENE) %>%
#   unique()

# p_venn_top20 <- ggVennDiagram(list(Male = male_top20, Female = female_top20), label = "count") +
#   ggtitle("Top 20 genes overlap (smallest P)")

# # -------------------------
# # Show
# # -------------------------
# quartz()
# p_venn_sig

# quartz()
# p_venn_top20

# # Optional: print counts
# cat("\n--- P < 0.05 ---\n")
# cat("Male:", length(male_sig), " Female:", length(female_sig),
#     " Shared:", length(intersect(male_sig, female_sig)), "\n")

# cat("\n--- Top 20 ---\n")
# cat("Male:", length(male_top20), " Female:", length(female_top20),
#     " Shared:", length(intersect(male_top20, female_top20)), "\n")



