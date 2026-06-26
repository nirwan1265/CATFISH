#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(dplyr); library(data.table); library(readr) })

GWAS_MODE <- "top_pct"; GWAS_TOP_PCT <- 1; GWAS_P_THRESHOLD <- 1e-5
MAGMA_MODE <- "top_pct"; MAGMA_TOP_PCT <- 5; MAGMA_FDR_THRESHOLD <- 0.30
WINDOW_SIZE <- 25000
TRAIT <- "First"
CHROMS <- c("2L","2R","3L","3R","X")
MAGMA_BASE <- "/Users/nirwantandukar/Documents/Research/results/CATFISH/MAGMA/First"
BASE_DIR   <- "/Users/nirwantandukar/Documents/Research/results/CATFISH/DGRP/First"
CATFISH_SUBDIR <- "CATFISH_permutation_B1000000_mvn_GPD"
GWAS_FILE  <- "/Users/nirwantandukar/Documents/Research/results/DGRP/DGRP_chilipeppers/GWAS_BLUEs/First.assoc.txt"
GENE_LOC_FILE <- "/Users/nirwantandukar/Documents/Github/MAGCAT/inst/extdata/fly.genes.loc"
OUT_DIR <- file.path(BASE_DIR, "candidate_gene_scoring")
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)
CHR_MAP <- c("2L"=1L,"2R"=2L,"3L"=3L,"3R"=4L,"4"=5L,"X"=6L,"Y"=7L)

pick_completed_catfish_csv <- function(base_dir, subdir, trait) {
  csvs <- c(file.path(base_dir, subdir, paste0(trait, "_CATFISH_ACAT_mvn_B1000000_GPD.csv")),
            file.path(base_dir, subdir, "omni_acat_mvn.csv"))
  hit <- csvs[file.exists(csvs)]
  if (!length(hit)) stop("No CATFISH CSV found. Run run_catfish_First.R first.", call.=FALSE)
  message("Using: ", hit[[1]]); hit[[1]]
}
pick_completed_gene_table <- function(base_dir, subdir, trait) {
  hit <- file.path(base_dir, subdir, paste0(trait, "_combined_genes.tsv"))
  if (!file.exists(hit)) stop("No gene table found.", call.=FALSE); hit
}
read_gene_loc_clean <- function(path) {
  x <- read.delim(path, stringsAsFactors=FALSE, check.names=FALSE)
  x <- x[, c("GENE","CHR","START","STOP"), drop=FALSE]
  x$GENE <- as.character(x$GENE); x$CHR <- as.integer(x$CHR)
  x$START <- as.integer(x$START); x$STOP <- as.integer(x$STOP)
  x <- x[!duplicated(x$GENE),]
  x$START_EXT <- pmax(0L, x$START - WINDOW_SIZE); x$STOP_EXT <- x$STOP + WINDOW_SIZE; x
}
read_assoc_gwas <- function(gwas_file, chr_map) {
  x <- fread(gwas_file, header=TRUE)
  col_map <- c("p_wald"="P","P.value"="P","p.value"="P","P-Value"="P","pvalue"="P",
               "ps"="POS","pos"="POS","BP"="POS","Pos"="POS",
               "chr"="CHR","Chr"="CHR","chromosome"="CHR","rs"="SNP","SNP"="SNP")
  for (old in names(col_map))
    if (old %in% names(x) && !col_map[[old]] %in% names(x)) setnames(x, old, col_map[[old]])
  if (is.character(x$CHR) || any(x$CHR %in% names(chr_map))) x[, CHR := chr_map[as.character(CHR)]]
  x <- x[!is.na(CHR) & !is.na(POS)]; x[, CHR := as.integer(CHR)]; x[, POS := as.integer(POS)]; x
}

CATFISH_CSV <- pick_completed_catfish_csv(BASE_DIR, CATFISH_SUBDIR, TRAIT)
MAGMA_GENE_FILE <- pick_completed_gene_table(BASE_DIR, CATFISH_SUBDIR, TRAIT)
cat("CATFISH CSV:      ", CATFISH_CSV, "\n", sep="")
cat("MAGMA gene table: ", MAGMA_GENE_FILE, "\n", sep="")
cat("GWAS file:        ", GWAS_FILE, "\n", sep="")
cat("Output dir:       ", OUT_DIR, "\n\n", sep="")

gene_loc <- read_gene_loc_clean(GENE_LOC_FILE)
cat("Genes in annotation: ", nrow(gene_loc), "\n", sep="")

gwas_dt <- read_assoc_gwas(GWAS_FILE, CHR_MAP)
cat("Total GWAS SNPs: ", nrow(gwas_dt), "\n", sep="")
gene_dt <- as.data.table(gene_loc); setkey(gene_dt, CHR, START_EXT, STOP_EXT)
gwas_dt[, c("START_EXT","STOP_EXT") := .(POS, POS)]; setkey(gwas_dt, CHR, START_EXT, STOP_EXT)
snp_gene_map <- foverlaps(gwas_dt, gene_dt, by.x=c("CHR","START_EXT","STOP_EXT"),
                           by.y=c("CHR","START_EXT","STOP_EXT"), type="within", nomatch=NULL)
gwas_gene <- snp_gene_map %>% as.data.frame() %>% group_by(GENE) %>%
  summarise(gwas_min_p=min(P,na.rm=TRUE), gwas_n_snps=n(), .groups="drop") %>%
  arrange(gwas_min_p, GENE) %>% mutate(gwas_rank=dplyr::row_number())
gwas_eff_threshold <- sort(gwas_gene$gwas_min_p)[ceiling(nrow(gwas_gene)*GWAS_TOP_PCT/100)]
cat("Genes with GWAS evidence: ", nrow(gwas_gene), " | threshold: ",
    format(gwas_eff_threshold, digits=4), "\n\n", sep="")

magma_gene <- read.delim(MAGMA_GENE_FILE, stringsAsFactors=FALSE, check.names=FALSE)
if (!"P" %in% names(magma_gene) && "P_MULTI" %in% names(magma_gene)) magma_gene$P <- magma_gene$P_MULTI
if (!"ZSTAT" %in% names(magma_gene)) magma_gene$ZSTAT <- NA_real_
magma_gene <- magma_gene %>% rename(magma_p=P, magma_z=ZSTAT) %>%
  mutate(magma_fdr=p.adjust(magma_p,"BH"), magma_bonf=p.adjust(magma_p,"bonferroni")) %>%
  arrange(magma_p, GENE) %>% mutate(magma_rank=dplyr::row_number())
magma_eff_threshold <- sort(magma_gene$magma_p)[ceiling(nrow(magma_gene)*MAGMA_TOP_PCT/100)]
cat("Genes with MAGMA results: ", nrow(magma_gene), " | threshold: ",
    format(magma_eff_threshold, digits=4), "\n\n", sep="")

omni_results <- read.csv(CATFISH_CSV, stringsAsFactors=FALSE)
omni_col <- if ("omni_p_final" %in% names(omni_results)) "omni_p_final" else
            if ("omni_p_mvn" %in% names(omni_results)) "omni_p_mvn" else "omni_p_analytic"
if (!"pathway_name" %in% names(omni_results)) omni_results$pathway_name <- omni_results$pathway_id
omni_results <- omni_results %>%
  arrange(.data[[omni_col]], pathway_name, pathway_id) %>% mutate(pathway_rank=dplyr::row_number())

pathway_genes_list <- strsplit(omni_results$genes_used, ";")
pathway_genes <- data.frame(
  pathway_id=rep(omni_results$pathway_id, lengths(pathway_genes_list)),
  pathway_name=rep(omni_results$pathway_name, lengths(pathway_genes_list)),
  pathway_p=rep(omni_results[[omni_col]], lengths(pathway_genes_list)),
  pathway_rank=rep(omni_results$pathway_rank, lengths(pathway_genes_list)),
  GENE=trimws(unlist(pathway_genes_list)), stringsAsFactors=FALSE
) %>% filter(!is.na(GENE), nzchar(GENE))

pathway_best <- pathway_genes %>% arrange(pathway_rank, pathway_name, pathway_id) %>%
  group_by(GENE) %>% slice_head(n=1) %>% ungroup() %>%
  transmute(GENE, best_pathway_id=pathway_id, best_pathway_name=pathway_name,
            best_pathway_p=pathway_p, pathway_rank=pathway_rank)

pathway_gene_support <- pathway_genes %>% group_by(GENE) %>%
  summarise(n_top_pathways=n_distinct(pathway_id),
            mean_pathway_mlog10p=mean(-log10(pathway_p), na.rm=TRUE),
            pathways=paste(pathway_id, collapse="; "),
            pathway_names=paste(pathway_name, collapse="; "), .groups="drop") %>%
  left_join(pathway_best, by="GENE") %>% arrange(pathway_rank, GENE)

n_pathways_total <- nrow(omni_results); n_genes_total <- nrow(magma_gene)
cat("Pathways: ", n_pathways_total, " | Genes with pathway membership: ",
    nrow(pathway_gene_support), "\n\n", sep="")

gene_evidence <- magma_gene %>%
  left_join(gwas_gene, by="GENE") %>%
  left_join(pathway_gene_support, by="GENE") %>%
  mutate(
    hit_gwas       = !is.na(gwas_min_p) & gwas_min_p <= gwas_eff_threshold,
    hit_magma      = !is.na(magma_p) & magma_p <= magma_eff_threshold,
    hit_pathway    = !is.na(n_top_pathways) & n_top_pathways >= 1,
    support_layers = hit_gwas + hit_magma + hit_pathway,
    gwas_rank_score    = ifelse(!is.na(gwas_rank),    -log10(gwas_rank    / n_genes_total),    0),
    magma_rank_score   = ifelse(!is.na(magma_rank),   -log10(magma_rank   / n_genes_total),    0),
    pathway_rank_score = ifelse(!is.na(pathway_rank), -log10(pathway_rank / n_pathways_total), 0),
    multi_pathway_bonus = ifelse(!is.na(n_top_pathways), log10(n_top_pathways + 1), 0),
    score = gwas_rank_score + magma_rank_score + 0.5*pathway_rank_score + 0.3*multi_pathway_bonus
  ) %>%
  arrange(desc(score), gwas_rank, magma_rank, pathway_rank, GENE)

select_cols <- intersect(c("GENE","score","support_layers",
  "gwas_rank_score","magma_rank_score","pathway_rank_score","multi_pathway_bonus",
  "magma_p","magma_fdr","magma_rank","magma_z","gwas_min_p","gwas_n_snps","gwas_rank",
  "n_top_pathways","best_pathway_p","pathway_rank","best_pathway_id","best_pathway_name",
  "pathways","pathway_names","hit_gwas","hit_magma","hit_pathway"), names(gene_evidence))

top50  <- gene_evidence %>% select(all_of(select_cols)) %>% slice_head(n=50)
top200 <- gene_evidence %>% select(all_of(select_cols)) %>% slice_head(n=200)
write.csv(gene_evidence %>% select(all_of(select_cols)),
          file.path(OUT_DIR, paste0("candidate_genes_all_",   TRAIT, ".csv")), row.names=FALSE)
write.csv(top50,  file.path(OUT_DIR, paste0("candidate_genes_top50_",  TRAIT, ".csv")), row.names=FALSE)
write.csv(top200, file.path(OUT_DIR, paste0("candidate_genes_top200_", TRAIT, ".csv")), row.names=FALSE)
cat("Saved to:", OUT_DIR, "\n\nTop 10 genes by score:\n")
print(top50 %>% slice_head(n=10))
