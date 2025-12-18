devtools::document()
devtools::load_all()


library(MAGCAT)
library(metapro)
library(sumFREGAT)
library(metap)
library(rtracklayer)
library(ggplot2)
library(TFisher)
library(data.table)
# 1. See where it's installed
.libPaths()






# gff3_to_geneloc
gff_path <- "/Users/nirwantandukar/Downloads/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3"

gene_loc_out <- "/Users/nirwantandukar/Downloads/maize.genes.loc"

gff3_to_geneloc(
  gff        = "/Users/nirwantandukar/Downloads/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3",
  out        = "inst/extdata/maize.genes.loc",
  chr_prefix = "chr"   # default; strips "chr"
)

# Gene-loc function
# If it has
magma_annotate(
  stats_file = "/Users/nirwantandukar/Documents/Research/results/GWAS/GAPIT/raw_GWAS_MLM_3PC_N.txt",
  rename_columns = c(
    CHR    = "Chr",
    SNP    = "SNP",
    POS    = "Pos",
    PVALUE = "P.value"   # not used here but keeps things consistent
  ),
  species    = "maize",        # uses the built-in maize.genes.loc  
  out_prefix = "N_maize_MLM",
  out_dir    = "annot",
  window     = c(25, 25)
)

# Run MAGMA on gene-level
# If it has NMISS
magma_gene(
  bfile      = "/Users/nirwantandukar/Documents/Research/data/MAGMA/maize/bed_bim_fam_file/all_maize2",
  gene_annot = "annot/N_maize_MLM.genes.annot",
  stats_file = "/Users/nirwantandukar/Documents/Research/results/GWAS/MLM/nitrogen/GWAS_results/nitrogen_0-5cm_maize_LMM.txt",
  n_total    = 3107,
  rename_columns = c(
    CHR    = "chr",
    SNP    = "rs",
    POS    = "ps",
    PVALUE = "p_wald",
    NMISS  = "n_miss"
  ),
  out_prefix = "N_maize_MLM",
  out_dir    = "magma_genes",
  gene_model = c("multi=snp-wise"),
)

# Run MAGMA on gene-level per chromosome
magma_gene(
  bfile      = "/Users/nirwantandukar/Documents/Research/data/MAGMA/maize/bed_bim_fam_file/all_maize2",
  gene_annot = "annot/N_maize_MLM.genes.annot",
  stats_file = "/Users/nirwantandukar/Documents/Research/results/GWAS/GAPIT/raw_GWAS_MLM_3PC_N.txt",
  n_total    = 3107,
  rename_columns = c(
    CHR    = "Chr",
    SNP    = "SNP",
    POS    = "Pos",
    PVALUE = "P.value",
    NMISS  = "n_miss"
  ),
  out_prefix = "N_maize_MLM",
  out_dir    = "magma_genes_by_chr",
  gene_model = c("multi=snp-wise"),
  chroms     = 1:10,
  n_threads  = 10       # use up to 10 parallel workers
)


# Run MAGMA on gene-level per chromosome
# If it has total number of observations
magma_gene(
  bfile      = "/Users/nirwantandukar/Documents/Research/data/MAGMA/maize/bed_bim_fam_file/all_maize2",
  gene_annot = "annot/N_maize_MLM.genes.annot",
  stats_file = "/Users/nirwantandukar/Documents/Research/results/GWAS/GAPIT/raw_GWAS_MLM_3PC_N.txt",
  # n_total    = 3107, this is ignored.
  rename_columns = c(
    CHR    = "Chr",
    SNP    = "SNP",
    POS    = "Pos",
    PVALUE = "P.value",
    NOBS   = "nobs"
  ),
  out_prefix = "N_maize_MLM",
  out_dir    = "magma_multi_snp_wise_genes_by_chr_N_maize",
  gene_model = c("multi=snp-wise"),
  chroms     = 1:10,
  n_threads  = 10
)



# Combine all the chromosome
## Vector of MAGMA gene files (chr1–chr10)
files <- sprintf("/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_multi_snp_wise_genes_by_chr_N_maize/N_maize_MLM_chr%d.multi_snp_wise.genes.out", 1:10)

# Optional: if they’re in a subdir, prepend the path, e.g.:
# files <- file.path("magma_genes", sprintf("N_maize_MLM_chr%d.snp_wise_top.genes.out", 1:10))

# Read all files
gene_list <- lapply(files, function(f) {
  if (!file.exists(f)) stop("File not found: ", f)
  utils::read.table(f, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
})

# Stack into one big data.frame
genes_all_raw <- do.call(rbind, gene_list)

# ---- (Optional but recommended) deduplicate by GENE: keep smallest P ----
if (!"GENE" %in% names(genes_all_raw) || !"P" %in% names(genes_all_raw)) {
  stop("Combined MAGMA file must have columns 'GENE' and 'P'.")
}

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
write.table(
  genes_all,
  file      = "/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_multi_snp_wise_genes_by_chr_N_maize/magma_N_maize.txt",
  sep       = "\t",      # tab-separated
  quote     = FALSE,     # don't quote strings
  row.names = FALSE      # don't write row numbers
)


# Get gene length
gff3_path <- "/Users/nirwantandukar/Documents/Research/data/GFF3/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.62.chr.gff3"

maize_gene_len <- get_gene_lengths(
  gff3_file  = gff3_path,
  output     = TRUE,
  output_dir = "inst/extdata",
  file_name  = "Zea_mays_gene_lengths.tsv"
)



# remove "gene:""
maize_gene_len=read.delim("inst/extdata/Zea_mays_gene_lengths.tsv")
head(maize_gene_len)

head(genes_all)

### Adjust pvalues based on Gene length and NSPS
adj_out <- magcat_adjust_gene_p(
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
colnames(genes_adj)=c("GENE", "Z_adj","P_adj")
head(genes_adj)
head(genes_all)

write.csv(genes_adj,"genes_adj.csv", row.names=F)

lm_fit    <- adj_out$fit   # if you want to look at coefficients
summary(lm_fit)

df_diag <- cbind(
  model.frame(lm_fit),                    # Z_raw, log_gene_length, log_nsnp
  Z_fitted = fitted(lm_fit),              # predicted Z
  Z_resid  = resid(lm_fit)                # residual Z (size-adjusted)
)

str(df_diag)
head(df_diag)

library(ggplot2)

# Raw Z vs log gene length
p1 <- ggplot(df_diag, aes(x = log_gene_length, y = Z_raw)) +
  geom_hex(bins = 60) +                     # use geom_point if you want
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "log(gene length)",
    y = "Raw MAGMA Z-score",
    title = "Dependence of gene-level Z on gene length (raw)"
  )

# Residual Z vs log gene length
p2 <- ggplot(df_diag, aes(x = log_gene_length, y = Z_resid)) +
  geom_hex(bins = 60) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "log(gene length)",
    y = "Size-adjusted Z-score (residual)",
    title = "Dependence largely removed after adjustment"
  )
quartz()
p1
quartz()
p2


# Raw Z vs log(#SNPs)
p3 <- ggplot(df_diag, aes(x = log_nsnp, y = Z_raw)) +
  geom_hex(bins = 60) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "log(# SNPs per gene)",
    y = "Raw MAGMA Z-score",
    title = "Dependence of gene-level Z on SNP density (raw)"
  )

# Residual Z vs log(#SNPs)
p4 <- ggplot(df_diag, aes(x = log_nsnp, y = Z_resid)) +
  geom_hex(bins = 60) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "log(# SNPs per gene)",
    y = "Size-adjusted Z-score",
    title = "Dependence on SNP density removed"
  )
quartz()
p3
quartz()
p4


cor_before_len  <- cor(df_diag$Z_raw,  df_diag$log_gene_length, use = "complete.obs")
cor_after_len   <- cor(df_diag$Z_resid, df_diag$log_gene_length, use = "complete.obs")

cor_before_nsnp <- cor(df_diag$Z_raw,  df_diag$log_nsnp, use = "complete.obs")
cor_after_nsnp  <- cor(df_diag$Z_resid, df_diag$log_nsnp, use = "complete.obs")

c(
  cor_before_len,  cor_after_len,
  cor_before_nsnp, cor_after_nsnp
)

#“Raw gene Z-scores showed a weak but detectable correlation with log(gene length) and log(#SNPs) (r = 0.0X). After regression, correlations were reduced to ≈0, indicating effective removal of these biases.”


cor_p <- cor(genes_adj$P, genes_adj$P_adj, use = "complete.obs")
cor_logp <- cor(-log10(genes_adj$P), -log10(genes_adj$P_adj), use = "complete.obs")

cor_p
cor_logp

library(ggplot2)
quartz()
ggplot(genes_adj, aes(x = -log10(P), y = -log10(P_adj))) +
  geom_point(alpha = 0.2, size = 0.6) +
  geom_abline(slope = 1, intercept = 0, colour = "red") +
  labs(
    x = expression(-log[10](italic(P)[raw])),
    y = expression(-log[10](italic(P)[adj])),
    title = "Gene-level p-values before vs after size/SNP adjustment"
  )


ggplot(genes_adj, aes(x = P, y = P_adj)) +
  geom_point(alpha = 0.2, size = 0.6) +
  geom_abline(slope = 1, intercept = 0, colour = "red") +
  labs(
    x = expression(-log[10](italic(P)[raw])),
    y = expression(-log[10](italic(P)[adj])),
    title = "Gene-level p-values before vs after size/SNP adjustment"
  )


# ACAT
# 1) MAGMA gene-level results (your merged file, with GENE / P / CHR etc.)
genes_all <- read.table("/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_multi_snp_wise_genes_by_chr_N_maize/magma_N_maize.txt",
                        header = TRUE, stringsAsFactors = FALSE)
head(genes_all)

# 2) Load maize pathways from CornCyc (inst/extdata/pathway/corncyc_...)
maize_pw <- magcat_load_pathways("maize", gene_col = "Gene-name")
head(maize_pw)

# 3) Run ACAT per pathway (no permutations first):
pw_res <- magcat_acat_pathways(
  gene_results = genes_adj,
  species      = "maize",
  gene_col     = "GENE",
  p_col        = "P_adj",
  output       = TRUE,
  out_dir      = "acat_results"
)

head(pw_res)

# 4) If you want permutation-calibrated:
pw_res <- magcat_acat_pathways(
  gene_results = genes_all,
  species      = "maize",
  gene_col     = "GENE",
  #pathways     = maize_pw,
  p_col        = "P",
  B            = 1000,
  seed         = 42,
  output       = TRUE,
  out_dir      = "acat_results"
)

head(pw_res)
# Check the path
attr(pw_res, "file")

# 5.) ordmeta
ord_res <- magcat_ordmeta_pathways(
  gene_results = genes_adj,
  species      = "maize",
  gene_col     = "GENE",
  p_col        = "P_adj",
  effect_col   = "ZSTAT",
  is_onetail   = FALSE
)
write.csv(ord_res,"ordmeta.csv")

# wFisher, e.g. weight by NSNPS
wf_res <- magcat_wfisher_pathways(
  gene_results = genes_all,
  species      = "maize",
  gene_col     = "GENE",
  p_col        = "P",
  effect_col   = "ZSTAT",
  weight_col   = "NSNPS",
  is_onetail   = FALSE
)

write.csv(wf_res,"weighted_fisher.csv")

wf_res <- magcat_wfisher_pathways(
  gene_results = genes_adj,
  species      = "maize",
  gene_col     = "GENE",
  p_col        = "P_adj",
  effect_col   = "ZSTAT",
  #weight_col   = "NSNPS",
  is_onetail   = FALSE
)
write.csv(wf_res,"weighted_fisher_adj_gene_length.csv")

# TFisher
tf_res <- magcat_tfisher_pathways(
  gene_results    = genes_all,
  species         = "maize",
  gene_col        = "GENE",
  p_col           = "P",
  ptrunc          = 0.05,    # truncation threshold
  B_perm          = 1000L,   # permutations for empirical p
  seed            = 123,
  analytic_logical= TRUE,    # chi-square approx (optional)
  output          = TRUE,
  out_dir         = "magcat_tfisher"
)

head(tf_res)
write.csv(tf_res,"tf_res2.csv")

tf_res <- magcat_tfisher_pathways(
  gene_results    = genes_adj,
  species         = "maize",
  gene_col        = "GENE",
  p_col           = "P_adj",
  ptrunc          = 0.05,    # truncation threshold
  B_perm          = 10000L,   # permutations for empirical p
  seed            = 123,
  analytic_logical= TRUE,    # chi-square approx (optional)
  output          = TRUE,
  out_dir         = "magcat_tfisher"
)

write.csv(tf_res,"tf_res_adj.csv")


# Soft Tfisher
tf_res <- magcat_soft_tfisher_pathways(
  gene_results    = genes_adj,
  species         = "maize",
  gene_col        = "GENE",
  p_col           = "P_adj",
  tau1            = 0.05,
  B_perm          = 10000L,
  seed            = 123,
  analytic_logical= TRUE,
  output          = TRUE,
  out_dir         = "magcat_tfisher_soft"
)

tf_res
write.csv(tf_res,"soft_tf_res_adj_gene_length.csv")


# Stouffer
stouf_res <- magcat_stoufferZ_pathways(
  gene_results = genes_adj,
  species      = "maize",
  gene_col     = "GENE",
  z_col        = "Z_adj",
  weight_col   = NULL,   # equal weights
  B_perm       = 0L,     # <- NO permutations
  seed         = NULL,
  output       = TRUE,
  out_dir      = "magcat_stouffer"
)

head(stouf_res)
attr(stouf_res, "file")


# MinP
minp_res <- magcat_minp_pathways(
  gene_results = genes_adj,  # your .genes.out as data.frame
  species      = "maize",      # or "sorghum"/"arabidopsis"/"plant"
  gene_col     = "GENE",
  p_col        = "P_adj",
  B_perm       = 0L,           # no permutations, analytic only
  min_p        = 1e-15,
  do_fix       = TRUE,
  output       = TRUE,
  out_dir      = "magcat_minp_maize"
)



############################################################
## 2) Adaptive soft-TFisher (analytic only, NO permutations)
##    (Assuming your adaptive function is named:
##     magcat_soft_tfisher_adaptive_pathways)
############################################################

tf_adapt_res <- magcat_soft_tfisher_adaptive_pathways(
  gene_results     = genes_adj,
  species          = "maize",
  gene_col         = "GENE",
  p_col            = "P_adj",
  tau_grid         = c(0.001, 0.005, 0.01, 0.02, 0.05, 0.1),
  #min_p            = 1e-15,
  do_fix           = TRUE,
  B_perm           = 0L,        # <- NO permutations (must be integer)
  perm_mode        = "resample_global",  # ignored when B_perm=0L
  seed             = NULL,
  analytic_logical = TRUE,
  output           = TRUE,
  out_dir          = "magcat_tfisher_adaptive"
)

head(tf_adapt_res)
attr(tf_adapt_res, "file")
args(MAGCAT::magcat_soft_tfisher_adaptive_pathways)


### MAGMA-competitive
## First combine the raw files my guy
## Combine per-chromosome MAGMA gene results (.genes.raw-style) into one file
## (works for your custom v1.10 format with "# VERSION" and correlation tail)

in_dir  <- "/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_multi_snp_wise_genes_by_chr_N_maize"
pattern <- "^N_maize_MLM_chr[0-9]+\\.multi_snp_wise.genes.raw"  
out_file <- file.path(in_dir, "N_maize_MLM_ALLCHR.multi_snp_wise.genes")

files <- list.files(in_dir, pattern = pattern, full.names = TRUE)
stopifnot(length(files) > 0)

# order chr1, chr2, ... chr10, ...
chr_num <- function(x) as.integer(sub(".*_chr([0-9]+)\\..*", "\\1", basename(x)))
files <- files[order(chr_num(files))]

# helper: extract header lines + data lines
read_genes_file <- function(f) {
  lines <- readLines(f, warn = FALSE)
  hdr <- lines[grepl("^#", lines)]
  dat <- lines[!grepl("^#", lines) & nzchar(lines)]
  list(hdr = hdr, dat = dat)
}

x1 <- read_genes_file(files[1])

# keep header from the first file only
header_keep <- x1$hdr

# drop any "# COVAR" lines from subsequent files (if present) and keep only data
data_all <- x1$dat

if (length(files) > 1) {
  for (f in files[-1]) {
    xi <- read_genes_file(f)
    data_all <- c(data_all, xi$dat)
  }
}

# write merged
writeLines(c(header_keep, data_all), con = out_file)

cat("Wrote:", out_file, "\n")
cat("Chromosomes merged:", length(files), "\n")



# EDIT PATHWAY
library(data.table)
head(maize_pw)
# maize_pw: pathway_id, pathway_name, gene_id
dt <- as.data.table(maize_pw)[, .(set_id = as.character(pathway_id),
                                 gene_id = sub("^gene:", "", as.character(gene_id)))]

dt <- dt[!is.na(set_id) & nzchar(set_id) & !is.na(gene_id) & nzchar(gene_id)]
dt <- unique(dt)

# (optional but recommended) keep only genes that exist in your .genes.raw
genes_raw <- "/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_multi_snp_wise_genes_by_chr_N_maize/N_maize_MLM_ALLCHR.multi_snp_wise.genes"
raw_lines <- readLines(genes_raw, warn = FALSE)
raw_dat   <- raw_lines[!grepl("^#", raw_lines) & nzchar(raw_lines)]
genes_ok  <- unique(sub("\\s+.*$", "", raw_dat))
dt <- dt[gene_id %in% genes_ok]

# build wide lines: set_id \t gene1 \t gene2 ...
setlist <- split(dt$gene_id, dt$set_id)
setlist <- lapply(setlist, unique)

set_annot_wide <- "annot/maize_pathways.sets.annot.WIDE"
dir.create(dirname(set_annot_wide), recursive = TRUE, showWarnings = FALSE)

con <- file(set_annot_wide, open = "wt")
for (sid in names(setlist)) {
  genes <- setlist[[sid]]
  if (length(genes) == 0) next
  writeLines(paste(c(sid, genes), collapse = "\t"), con = con)
}
close(con)

set_annot_wide


## Example usage (your workflow becomes ONE call)
# ONE call (but internally it does exactly what you described)
# -------------------------
# Example usage (your case):
# -------------------------
maize_pw <- magcat_load_pathways("maize", gene_col = "Gene-name")


out <- magma_geneset_competitive(
  gene_results_raw = "/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_multi_snp_wise_genes_by_chr_N_maize/N_maize_MLM_ALLCHR.multi_snp_wise.genes",
  #set_annot        = "annot/maize_pathways.sets.annot.WIDE",  # OR set_annot = NULL to auto-build
  set_annot         = NULL,
  out_prefix       = "N_maize_MLM_ALLCHR.PMN_COMP",
  out_dir          = "magma_geneset",
  pathways         = maize_pw,
  genes_all        = genes_all,
  filter_to_gene_results = TRUE
)


# 1) User-friendly species mode (auto-load PMN + auto-build WIDE set-annot)


# -----------------------------------------
# 3) How to store your built-in files (pkg)
# -----------------------------------------
# Put these in your package repo:
#   inst/extdata/MAGMA_pathway/maize_pathway.sets.annot.WIDE
#   inst/extdata/MAGMA_pathway/maize_pathways.long.tsv
#
# NOTE: You said "maize_pathway.txt" — you *can* name it that, but keep a consistent
# suffix so it’s obvious it’s a WIDE set-annot file. I’d keep:
#   *.sets.annot.WIDE

# ----------------------------
# 4) ONE-call usage (final)
# ----------------------------
out <- magma_geneset_competitive(
  gene_results_raw = "/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_multi_snp_wise_genes_by_chr_N_maize/N_maize_MLM_ALLCHR.multi_snp_wise.genes",
  out_prefix       = "N_maize_MLM_ALLCHR.PMN_COMP",
  out_dir          = "magma_geneset",
  species          = "maize",
  genes_all        = genes_all,
  write_tidy       = TRUE
)




head(out)
attr(out, "set_annot_wide")
attr(out, "tidy_file")








# OMNI
omni_minp <- omni_pathways(
  gene_results      = genes_adj,
  species           = "maize",
  gene_col          = "GENE",
  p_col             = "P_adj",
  effect_col        = "Z_adj",
  #weight_col        = "NSNPS",
  is_onetail        = FALSE,
  ptrunc            = 0.05,
  min_p             = 1e-15,
  do_fix            = TRUE,
  omnibus           = "ACAT",      # minP or "ACAT"
  B_perm            = 10000L,
  seed              = 123,
  perm_mode    = "mvn",       # mvn or 
  magma_genes_out = "/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_multi_snp_wise_genes_by_chr_N_maize/magma_N_maize.txt",
  remove_singletons = TRUE,
  output            = TRUE,
  out_dir           = "magcat_omni_full"
)


omni_mvn <- omni_pathways(
  gene_results      = genes_adj,
  species           = "maize",
  gene_col          = "GENE",
  p_col             = "P_adj",
  ptrunc            = 0.05,
  min_p             = 1e-15,
  do_fix            = TRUE,
  omnibus           = "minP",
  B_perm            = 10000L,
  seed              = 123,
  perm_mode         = "mvn",
  magma_genes_out   = "/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_multi_snp_wise_genes_by_chr_N_maize/magma_N_maize.txt",
  magma_cor_file    = NULL,
  remove_singletons = TRUE,
  output            = TRUE,
  out_dir           = "magcat_omni_full"
)



args(MAGCAT::omni_pathways)

m  <- sum(!is.na(omni_mvn$omni_perm_p))
p1 <- min(omni_mvn$omni_perm_p, na.rm = TRUE)

c(m = m, min_perm_p = p1, p1_times_m = p1 * m)
head(sort(omni_mvn$omni_perm_p), 20)
head(sort(omni_mvn$omni_perm_p_BH), 20)



p <- genes_adj$P_adj
c(n = length(p),
  n_na = sum(is.na(p)),
  n_eq1 = sum(!is.na(p) & p == 1),
  frac_eq1 = mean(!is.na(p) & p == 1),
  n_gt0_le1 = sum(!is.na(p) & p > 0 & p <= 1),
  n_gt0_lt1 = sum(!is.na(p) & p > 0 & p < 1))



old <- omni_pathways(
  gene_results=genes_adj, species="maize",
  gene_col="GENE", p_col="P_adj",
  ptrunc=0.05, min_p=1e-15, do_fix=TRUE,
  B_perm=100L, seed=123,
  perm_mode="mvn",
  magma_genes_out="/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_multi_snp_wise_genes_by_chr_N_maize/magma_N_maize.txt",
  remove_singletons=TRUE
)

new <- omni_pathways(
  gene_results=genes_adj, species="maize",
  gene_col="GENE", p_col="P_adj",
  ptrunc=0.05, min_p=1e-15, do_fix=TRUE,
  B_perm=100L, seed=123,
  perm_mode="mvn",
  magma_genes_out="/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_multi_snp_wise_genes_by_chr_N_maize/magma_N_maize.txt",
  remove_singletons=TRUE
)

c(m_old = sum(!is.na(old$omni_perm_p)),
  min_old = min(old$omni_perm_p, na.rm=TRUE),
  m_new = sum(!is.na(new$omni_perm_p)),
  min_new = min(new$omni_perm_p, na.rm=TRUE))

quantile(old$omni_perm_p, c(.01,.05,.1,.2,.5), na.rm=TRUE)
quantile(new$omni_perm_p, c(.01,.05,.1,.2,.5), na.rm=TRUE)
