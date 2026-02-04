devtools::document()
library(CATFISH)
library(metapro)
library(sumFREGAT)
library(metap)
library(rtracklayer)
library(ggplot2)
library(TFisher)
library(data.table)
library(dplyr)



# Fly
gff_path <- "/Users/nirwantandukar/Documents/Research/data/Pathway/Drosophila_melanogaster.BDGP6.54.115.gff3"
# Fly
gene_loc_out <- "/Users/nirwantandukar/Downloads/fly.genes.loc"


# Fly
files <- list.files(path = "/Users/nirwantandukar/Documents/Research/results/DGRP/MAGMA/Fly_magma_genes_by_chr_female",
        pattern = "^Female_starvation_fly_.*\\.genes\\.out$",
        full.names = TRUE)


# Read all files
gene_list <- lapply(files, function(f) {
  if (!file.exists(f)) stop("File not found: ", f)
  utils::read.table(f, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
})

# Stack into one big data.frame
genes_all_raw <- do.call(rbind, gene_list)

head(genes_all_raw)
colnames(genes_all_raw)[9] = "P"

# Order by p-value, then drop duplicates so we keep the best P per gene
o <- order(genes_all_raw$GENE, genes_all_raw$P)
genes_all <- genes_all_raw[o, ]
genes_all <- genes_all[!duplicated(genes_all$GENE), ]

# Nice to have: sort by chromosome/start for sanity
if (all(c("CHR", "START") %in% names(genes_all))) {
  genes_all <- genes_all[order(genes_all$CHR, genes_all$START), ]
}

# Inspect
head(genes_all)
nrow(genes_all)

length(unique(genes_all$GENE))

# Fly
gff3_path <- "/Users/nirwantandukar/Documents/Research/data/GFF3/"

# For flies gene length:
maize_gene_len=read.delim("inst/extdata/dmel.flybase.fbgn.genes.loc.tsv")
head(maize_gene_len)


### Adjust pvalues based on Gene length and NSPS
adj_out <- catfish_adjust_gene_p(
  gene_results = genes_all,
  gene_lengths = maize_gene_len,
  gene_col     = "GENE",
  nsnp_col     = "NSNPS",
  p_col        = "P",
  z_col        = "ZSTAT",       # or NULL to reconstruct from P
  len_gene_col = "gene_id",
  len_col      = "length"
)
head(adj_out)
genes_adj <- adj_out[,c(1,2,3)]
colnames(genes_adj)=c("GENE", "ZSTAT","P")
head(genes_adj)
head(genes_all)


# Fly
fly_pw=read.delim("inst/extdata/pathway/Fly_Cyc.tsv")
head(fly_pw)


omni_results <- catfish_omni2_pathways(
  gene_results   = genes_adj,
  species        = NULL,                     # load PMN maize pathways automatically
  pathways       = fly_pw,
  pmn_gene_col   = "Gene-name",                 # column in PMN file
  gene_col       = "GENE",                      # column in your gene results
  p_raw_col      = "P",                         # use MAGMA P column
  z_col          = "ZSTAT",                     # use MAGMA ZSTAT column for Stouffer
  weight_col     = NULL,                        # optional if you have custom weights
  tau_grid       = c(0.1, 0.05, 0.02, 0.01, 0.005, 0.001),
  min_p          = 1e-15,
  do_fix         = TRUE,
  stouffer_min_abs_w = 1e-8,
  stouffer_alternative = "greater",
  #magma_out      = out,              # MAGMA competitive p-values
  include_magma_in_omni = FALSE,
  include_magma_in_perm = FALSE,                # only for analytic omnibus, no MAGMA in permutations
  omnibus        = "ACAT",                      # "ACAT" or "minP"
  B_perm         = 10000L,                        # number of permutations for omnibus
  perm_mode      = "mvn",                       # "mvn"or "global", "both", "none"
  magma_cor_file = "/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_multi_snp_wise_genes_by_chr_N_maize/magma_gene_cor_pairs_MLM_Fly_female.txt",  # 3-column file gene1 gene2 r
  #magma_cor_file = "/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_multi_snp_wise_genes_by_chr_N_maize/magma_gene_cor_pairs_MLM_Fly_male.txt",  # 3-column file gene1 gene2 r
  make_PD        = TRUE,
  mvn_marginal = "uniform", # “uniform”, “empirical”
  mvn_calibrate_components = TRUE, # NEW
  seed           = 123,
  output         = TRUE,
  out_dir        = "catfish_omnibus_results_Fly_female_B10000"
)

#saveRDS(omni_results,"omni_Fly_male_ACAT_B10000.RDS")
saveRDS(omni_results,"omni_Fly_female_ACAT_B10000.RDS")
