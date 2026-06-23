#!/usr/bin/env Rscript
# Extract top 20 genes by GWAS rank, MAGMA rank, and CATFISH score
# from the full all-genes scoring files, then write a Word table.

suppressPackageStartupMessages({
  for (pkg in c("data.table","officer","flextable")) {
    if (!requireNamespace(pkg, quietly=TRUE))
      install.packages(pkg, repos="https://cloud.r-project.org", quiet=TRUE)
    library(pkg, character.only=TRUE)
  }
})

F_FILE <- "/Users/nirwantandukar/Documents/Research/results/CATFISH/DGRP/Fly_starvation_female/candidate_gene_scoring_B100000_GPD/candidate_genes_all_Starvation_stress_female_B100000_GPD.csv"
M_FILE <- "/Users/nirwantandukar/Documents/Research/results/CATFISH/DGRP/Fly_starvation_male/candidate_gene_scoring_B100000_GPD/candidate_genes_all_Starvation_stress_male_B100000_GPD.csv"
OUT    <- "/Users/nirwantandukar/Documents/Github/MAGCAT/Starvation_Candidate_Genes_Top20.docx"

message("Reading female file...")
f <- fread(F_FILE)
message("Reading male file...")
m <- fread(M_FILE)

top20 <- function(dt, col, decreasing = FALSE) {
  dt_sorted <- dt[order(dt[[col]], decreasing = decreasing)]
  dt_sorted <- dt_sorted[!is.na(dt_sorted[[col]])]
  genes <- head(dt_sorted$GENE, 20)
  layers <- head(dt_sorted$support_layers, 20)
  pathway <- head(dt_sorted$best_pathway_name, 20)
  # Add ★ for 3-layer genes
  label <- ifelse(layers == 3,
                  paste0(genes, " \u2605"),
                  genes)
  label
}

# Female
fg <- top20(f, "gwas_rank",   decreasing = FALSE)
fm <- top20(f, "magma_rank",  decreasing = FALSE)
fc <- top20(f, "score",       decreasing = TRUE)

# Male
mg <- top20(m, "gwas_rank",   decreasing = FALSE)
mm <- top20(m, "magma_rank",  decreasing = FALSE)
mc <- top20(m, "score",       decreasing = TRUE)

# Build data frame
tbl <- data.frame(
  Rank            = 1:20,
  Female_GWAS     = fg,
  Female_MAGMA    = fm,
  Female_CATFISH  = fc,
  Male_GWAS       = mg,
  Male_MAGMA      = mm,
  Male_CATFISH    = mc,
  stringsAsFactors = FALSE
)

message("Top 20 table:")
print(tbl)

# ── Build Word doc with flextable ──────────────────────────────────────────────
ft <- flextable(tbl) |>
  set_header_labels(
    Rank           = "Rank",
    Female_GWAS    = "Female\nGWAS",
    Female_MAGMA   = "Female\nMAGMA",
    Female_CATFISH = "Female\nCATFISH",
    Male_GWAS      = "Male\nGWAS",
    Male_MAGMA     = "Male\nMAGMA",
    Male_CATFISH   = "Male\nCATFISH"
  ) |>
  bold(part = "header") |>
  bg(part = "header", bg = "#D6E4F0") |>
  color(part = "header", color = "#1F4E79") |>
  fontsize(size = 8, part = "all") |>
  font(fontname = "Arial", part = "all") |>
  align(align = "center", part = "all") |>
  # Highlight 3-layer genes (contain ★) in green
  bg(j = "Female_CATFISH",
     i = ~ grepl("\u2605", Female_CATFISH), bg = "#E2EFDA") |>
  bg(j = "Male_CATFISH",
     i = ~ grepl("\u2605", Male_CATFISH),   bg = "#E2EFDA") |>
  autofit() |>
  set_caption("Top 20 candidate genes for starvation resistance by method. \u2605 = 3-layer pathway-elevated gene.")

doc <- read_docx() |>
  body_add_par("Starvation GWAS — Top 20 Candidate Genes by Method",
               style = "heading 1") |>
  body_add_par("Ranked independently: GWAS (min SNP p-value per gene), MAGMA (gene-level p-value), CATFISH (multilayer score). \u2605 = 3-layer pathway-elevated gene.",
               style = "Normal") |>
  body_add_par("", style = "Normal") |>
  body_add_flextable(ft) |>
  body_add_par("", style = "Normal") |>
  body_add_par("Note: GWAS and MAGMA rankings are fully independent of CATFISH — all genes in the genome are ranked, not just CATFISH candidates.",
               style = "Normal")

print(doc, target = OUT)
message("Written: ", OUT)
