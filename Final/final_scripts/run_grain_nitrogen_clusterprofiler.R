suppressPackageStartupMessages({
  library(biomaRt)
  library(clusterProfiler)
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tidyr)
})

candidate_file <- "/Users/nirwantandukar/Documents/Github/Grain_Nitrogen/results/all_candidate_genes.csv"
RELAX_LANDRACE_GWAS <- TRUE
RELAXED_LANDRACE_P_CUTOFF <- 1e-6

out_dir <- if (RELAX_LANDRACE_GWAS) {
  file.path(
    "/Users/nirwantandukar/Documents/Github/MAGCAT/Final/final_supp_tables",
    paste0(
      "grain_nitrogen_clusterprofiler_relaxed_landrace_gwas_p",
      format(RELAXED_LANDRACE_P_CUTOFF, scientific = TRUE)
    )
  )
} else {
  "/Users/nirwantandukar/Documents/Github/MAGCAT/Final/final_supp_tables/grain_nitrogen_clusterprofiler"
}
gene_list_dir <- file.path(out_dir, "gene_lists")
enrichment_dir <- file.path(out_dir, "enrichment_tables")
biomart_cache_file <- "/Users/nirwantandukar/Documents/Github/MAGCAT/Final/final_supp_tables/grain_nitrogen_clusterprofiler/maize_biomart_go_background.csv"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gene_list_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(enrichment_dir, recursive = TRUE, showWarnings = FALSE)

clean_set_name <- function(x) {
  x |>
    str_replace("^\ufeff", "") |>
    str_replace_all("\\s+", "_") |>
    str_replace_all("[^A-Za-z0-9_]+", "_") |>
    str_replace_all("_+", "_") |>
    str_replace("^_|_$", "")
}

load_relaxed_landrace_gwas_sets <- function(p_cutoff) {
  sources <- tribble(
    ~gene_set, ~path,
    "GWAS_GBS_landraces_N", "/Users/nirwantandukar/Documents/Github/Grain_Nitrogen/results/GWAS_results/annotation_maize_N_5PC_farmcpu_full_maize_MLM.csv",
    "GWAS_PHG_landraces_N", "/Users/nirwantandukar/Documents/Github/Grain_Nitrogen/results/GWAS_results/annotation_maize_NHx_5PC_farmcpu_full_maize_MLM.csv"
  )

  sources |>
    mutate(
      genes = map(path, function(path_i) {
        read_csv(path_i, show_col_types = FALSE) |>
          filter(!is.na(PValue), PValue <= p_cutoff, !is.na(GeneID), GeneID != "") |>
          pull(GeneID) |>
          unique() |>
          sort()
      })
    ) |>
    select(gene_set, genes)
}

fetch_maize_go_background <- function(cache_file) {
  if (file.exists(cache_file)) {
    return(read_csv(cache_file, show_col_types = FALSE))
  }

  mart <- useEnsemblGenomes(
    biomart = "plants_mart",
    dataset = "zmays_eg_gene"
  )

  biomart_df <- getBM(
    attributes = c("ensembl_gene_id", "go_id", "name_1006", "namespace_1003"),
    mart = mart
  ) |>
    as_tibble() |>
    rename(
      GeneID = ensembl_gene_id,
      GO_ID = go_id,
      GO_term_name = name_1006,
      ontology = namespace_1003
    ) |>
    filter(
      !is.na(GeneID), GeneID != "",
      !is.na(GO_ID), GO_ID != "",
      !is.na(GO_term_name), GO_term_name != "",
      !is.na(ontology), ontology != ""
    ) |>
    mutate(
      ontology = recode(
        ontology,
        biological_process = "BP",
        molecular_function = "MF",
        cellular_component = "CC",
        .default = ontology
      )
    ) |>
    distinct(GeneID, ontology, GO_ID, GO_term_name)

  write_csv(biomart_df, cache_file)
  biomart_df
}

candidate_df <- read_csv(candidate_file, show_col_types = FALSE, name_repair = "minimal")
names(candidate_df) <- names(candidate_df) |> vapply(clean_set_name, character(1))
candidate_df <- candidate_df[, names(candidate_df) != "", drop = FALSE]

candidate_sets <- lapply(names(candidate_df), function(set_name) {
  genes <- candidate_df[[set_name]] |>
    as.character() |>
    trimws()
  genes <- genes[!is.na(genes) & genes != ""]
  tibble(gene_set = set_name, GeneID = unique(genes))
}) |>
  bind_rows()

if (RELAX_LANDRACE_GWAS) {
  relaxed_landrace_sets <- load_relaxed_landrace_gwas_sets(RELAXED_LANDRACE_P_CUTOFF)

  candidate_sets <- candidate_sets |>
    filter(!gene_set %in% relaxed_landrace_sets$gene_set) |>
    bind_rows(
      relaxed_landrace_sets |>
        mutate(tbl = map2(gene_set, genes, ~ tibble(gene_set = .x, GeneID = .y))) |>
        pull(tbl) |>
        bind_rows()
    )

  write_csv(
    relaxed_landrace_sets |>
      transmute(
        gene_set,
        relaxed_p_cutoff = RELAXED_LANDRACE_P_CUTOFF,
        n_candidate_genes = map_int(genes, length)
      ),
    file.path(out_dir, "relaxed_landrace_gwas_summary.csv")
  )
}

go_map <- fetch_maize_go_background(biomart_cache_file)

term2gene <- go_map |>
  select(GO_ID, GeneID) |>
  distinct()

term2name <- go_map |>
  select(GO_ID, GO_term_name, ontology) |>
  distinct()

universe_genes <- sort(unique(term2gene$GeneID))

gene_set_summary <- candidate_sets |>
  group_by(gene_set) |>
  summarise(
    n_candidate_genes = n_distinct(GeneID),
    n_go_mapped_genes = n_distinct(GeneID[GeneID %in% universe_genes]),
    pct_go_mapped = round(100 * n_go_mapped_genes / n_candidate_genes, 1),
    .groups = "drop"
  ) |>
  arrange(desc(n_candidate_genes), gene_set)

write_csv(gene_set_summary, file.path(out_dir, "candidate_gene_set_summary.csv"))
write_csv(candidate_sets, file.path(out_dir, "all_candidate_gene_lists_long.csv"))

split(candidate_sets$GeneID, candidate_sets$gene_set) |>
  iwalk(function(genes, set_name) {
    gene_tbl <- tibble(GeneID = sort(unique(genes)))
    write_csv(gene_tbl, file.path(gene_list_dir, paste0(set_name, "_genes.csv")))
  })

go_background_summary <- term2gene |>
  left_join(term2name, by = "GO_ID") |>
  group_by(ontology, GO_ID, GO_term_name) |>
  summarise(
    term_size = n_distinct(GeneID),
    .groups = "drop"
  ) |>
  arrange(ontology, desc(term_size), GO_ID)

write_csv(go_background_summary, file.path(out_dir, "go_background_terms.csv"))

run_enrichment_for_set <- function(genes, set_name) {
  genes <- sort(unique(genes[genes %in% universe_genes]))

  if (length(genes) < 3L) {
    return(tibble(
      gene_set = set_name,
      ontology = character(),
      ID = character(),
      Description = character(),
      GeneRatio = character(),
      BgRatio = character(),
      pvalue = numeric(),
      p.adjust = numeric(),
      qvalue = numeric(),
      geneID = character(),
      Count = integer(),
      n_input_genes = integer(),
      n_mapped_genes = integer()
    ))
  }

  enrich_res <- enricher(
    gene = genes,
    universe = universe_genes,
    TERM2GENE = term2gene,
    TERM2NAME = term2name |> select(GO_ID, GO_term_name),
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 3,
    maxGSSize = 2000
  )

  if (is.null(enrich_res) || nrow(as.data.frame(enrich_res)) == 0L) {
    return(tibble(
      gene_set = set_name,
      ontology = character(),
      ID = character(),
      Description = character(),
      GeneRatio = character(),
      BgRatio = character(),
      pvalue = numeric(),
      p.adjust = numeric(),
      qvalue = numeric(),
      geneID = character(),
      Count = integer(),
      n_input_genes = integer(),
      n_mapped_genes = integer()
    ))
  }

  as_tibble(enrich_res@result) |>
    left_join(term2name, by = c("ID" = "GO_ID", "Description" = "GO_term_name")) |>
    mutate(
      gene_set = set_name,
      n_input_genes = sum(candidate_sets$gene_set == set_name),
      n_mapped_genes = length(genes),
      .before = 1
    ) |>
    relocate(ontology, .after = gene_set) |>
    arrange(p.adjust, pvalue, desc(Count), ID)
}

enrichment_results <- split(candidate_sets$GeneID, candidate_sets$gene_set) |>
  imap(run_enrichment_for_set) |>
  bind_rows()

if (nrow(enrichment_results) > 0L) {
  split(enrichment_results, enrichment_results$gene_set) |>
    iwalk(function(tbl, set_name) {
      write_csv(tbl, file.path(enrichment_dir, paste0(set_name, "_clusterProfiler_GO.csv")))
    })
}

write_csv(enrichment_results, file.path(out_dir, "clusterProfiler_GO_all_sets.csv"))

sig_results <- enrichment_results |>
  filter(!is.na(p.adjust), p.adjust < 0.05)

write_csv(sig_results, file.path(out_dir, "clusterProfiler_GO_BH05_all_sets.csv"))

shared_terms <- sig_results |>
  group_by(ontology, ID, Description) |>
  summarise(
    n_gene_sets = n_distinct(gene_set),
    gene_sets = paste(sort(unique(gene_set)), collapse = "; "),
    best_p_adjust = min(p.adjust, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(n_gene_sets), best_p_adjust, ontology, Description)

write_csv(shared_terms, file.path(out_dir, "shared_enriched_GO_terms_across_sets.csv"))

message("Wrote outputs to: ", out_dir)
