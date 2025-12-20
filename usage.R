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

# Remove duplicate genes from inst/extdata/maize.genes.loc
# Keeps the first occurrence of each GENE (stable), writes a cleaned file.

remove_duplicate_genes_loc <- function(in_file  = "inst/extdata/maize.genes.loc",
                                      out_file = "inst/extdata/maize.genes.loc",
                                      keep     = c("first","best_span")) {
  keep <- match.arg(keep)

  if (!file.exists(in_file)) stop("File not found: ", in_file, call. = FALSE)

  # MAGMA gene-loc is whitespace-separated (tab or spaces), with header
  x <- utils::read.table(in_file, header = TRUE, sep = "", stringsAsFactors = FALSE, check.names = FALSE)

  req <- c("GENE","CHR","START","STOP")
  miss <- setdiff(req, names(x))
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse=", "), call. = FALSE)

  # Normalize types
  x$GENE  <- as.character(x$GENE)
  x$CHR   <- as.character(x$CHR)
  x$START <- as.integer(x$START)
  x$STOP  <- as.integer(x$STOP)

  # Drop rows with missing gene id
  x <- x[!is.na(x$GENE) & nzchar(x$GENE), , drop = FALSE]

  n_before <- nrow(x)
  dup_genes <- unique(x$GENE[duplicated(x$GENE)])

  if (keep == "first") {
    x2 <- x[!duplicated(x$GENE), , drop = FALSE]
  } else {
    # keep the "best_span": largest interval (STOP-START), tie -> first
    span <- x$STOP - x$START
    ord <- order(x$GENE, -span, seq_len(nrow(x)))
    xs <- x[ord, , drop = FALSE]
    x2 <- xs[!duplicated(xs$GENE), , drop = FALSE]
    # restore original-ish ordering by CHR then START (optional; comment if you want original)
    x2 <- x2[order(as.numeric(gsub("[^0-9]+","", x2$CHR)), x2$START), , drop = FALSE]
  }

  n_after <- nrow(x2)

  # Write back (tab-separated, no quotes) in MAGMA-friendly format
  utils::write.table(
    x2[, req],
    file = out_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )

  invisible(list(
    in_file = in_file,
    out_file = out_file,
    n_before = n_before,
    n_after  = n_after,
    n_removed = n_before - n_after,
    dup_genes_example = head(dup_genes, 20)
  ))
}

# ---- run it (in-place overwrite) ----
res <- remove_duplicate_genes_loc(
  in_file  = "inst/extdata/maize.genes.loc",
  out_file = "inst/extdata/maize.genes.loc",
  keep     = "first"      # or "best_span"
)

print(res)



# Gene-loc function
# If it has
magma_annotate(
  #stats_file = "/Users/nirwantandukar/Documents/Research/results/GWAS/GAPIT/raw_GWAS_MLM_3PC_N.txt",
  stats_file = "/Users/nirwantandukar/Documents/Research/results/GWAS/GLM/nitrogen/gwas_raw/GLM_maize_nitrogen_0_5_stat2.txt",
  rename_columns = c(
    CHR    = "Chr",
    SNP    = "Marker",
    POS    = "Pos",
    PVALUE = "p"   # not used here but keeps things consistent
  ),
  species    = "maize",        # uses the built-in maize.genes.loc  
  out_prefix = "N_maize_GLM",
  out_dir    = "annot",
  window     = c(25, 25)
)

# Run MAGMA on gene-level
# If it has NMISS
magma_gene(
  bfile      = "/Users/nirwantandukar/Documents/Research/data/MAGMA/maize/bed_bim_fam_file/all_maize2",
  gene_annot = "annot/N_maize_GLM.genes.annot",
  stats_file = "/Users/nirwantandukar/Documents/Research/results/GWAS/MLM/nitrogen/GWAS_results/nitrogen_0-5cm_maize_LMM.txt",
  n_total    = 3107,
  rename_columns = c(
    CHR    = "Chr",
    SNP    = "Marker",
    POS    = "Pos",
    PVALUE = "p",
    #NMISS  = "n_miss"
  ),
  out_prefix = "N_maize_MLM",
  out_dir    = "magma_genes",
  gene_model = c("multi=snp-wise"),
)

# Run MAGMA on gene-level per chromosome
magma_gene(
  bfile      = "/Users/nirwantandukar/Documents/Research/data/MAGMA/maize/bed_bim_fam_file/all_maize2",
  gene_annot = "annot/N_maize_MLM.genes.annot",
  #stats_file = "/Users/nirwantandukar/Documents/Research/results/GWAS/GAPIT/raw_GWAS_MLM_3PC_N.txt",
  stats_file = "/Users/nirwantandukar/Documents/Research/results/GWAS/GLM/nitrogen/gwas_raw/GLM_maize_nitrogen_0_5_stat.txt",
  n_total    = 3107,
  rename_columns = c(
    CHR    = "Chr",
    SNP    = "Marker",
    POS    = "Pos",
    PVALUE = "p",
    NMISS  = NULL
  ),
  out_prefix = "N_maize_GLM",
  out_dir    = "magma_genes_by_chr",
  gene_model = c("multi=snp-wise"),
  chroms     = 1:10,
  n_threads  = 10       # use up to 10 parallel workers
)

?magma_gene()

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
#files <- sprintf("/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_genes_by_chr/N_maize_GLM_chr%d.multi_snp_wise.genes.out", 1:10)

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
# write.table(
#   genes_all,
#   file      = "/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_multi_snp_wise_genes_by_chr_N_maize/magma_N__GLM_maize.txt",
#   sep       = "\t",      # tab-separated
#   quote     = FALSE,     # don't quote strings
#   row.names = FALSE      # don't write row numbers
# )


# Get gene length
gff3_path <- "/Users/nirwantandukar/Documents/Research/data/GFF3/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.62.chr.gff3"

# maize_gene_len <- get_gene_lengths(
#   gff3_file  = gff3_path,
#   output     = TRUE,
#   output_dir = "inst/extdata",
#   file_name  = "Zea_mays_gene_lengths.tsv"
# )



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


genes_adj <- adj_out[,c(1,2,3)]
colnames(genes_adj)=c("GENE", "ZSTAT","P")
head(genes_adj)
head(genes_all)

#write.csv(genes_adj,"genes_adj.csv", row.names=F)

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



# MAGMA gene-level results (your merged file, with GENE / P / CHR etc.)
# genes_all <- read.table("/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_multi_snp_wise_genes_by_chr_N_maize/magma_N_maize.txt",
#                         header = TRUE, stringsAsFactors = FALSE)
head(genes_all)

# Load maize pathways from CornCyc (inst/extdata/pathway/corncyc_...)
maize_pw <- magcat_load_pathways("maize", gene_col = "Gene-name")
head(maize_pw)

# Run ACAT per pathway (no permutations first):
pw_res <- magcat_acat_pathways(
  gene_results = genes_adj,
  species      = "maize",
  gene_col     = "GENE",
  min_p = 1e-15,
  p_col        = "P",
  output       = TRUE,
  out_dir      = "acat_results_GLM"
)
args(MAGCAT::magcat_acat_pathways)

# Fisher
f_res= magcat_fisher_pathways( genes_adj,
  species         = "maize",
  gene_col        = "GENE",
  p_col           = "P",
  min_p        = 1e-15,
  output          = TRUE,
  out_dir         = "magcat_fisher"
)
args(MAGCAT::magcat_fisher_pathways)

# Stouffer
stouf_res <- magcat_stoufferZ_pathways(
  gene_results = genes_adj,
  species      = "maize",
  gene_col     = "GENE",
  z_col        = "ZSTAT",
  weight_col   = NULL,   # equal weights
  B_perm       = 0L,     # <- NO permutations
  seed         = NULL,
  alternative = "greater",
  output       = TRUE,
  out_dir      = "magcat_stouffer"
)
args(MAGCAT::magcat_stoufferZ_pathways)

# MinP
minp_res <- magcat_minp_pathways(
  gene_results = genes_adj,  # your .genes.out as data.frame
  species      = "maize",      # or "sorghum"/"arabidopsis"/"plant"
  gene_col     = "GENE",
  p_col        = "P",
  B_perm       = 0L,           # no permutations, analytic only
  min_p        = 1e-15,
  do_fix       = TRUE,
  output       = TRUE,
  out_dir      = "magcat_minp_maize"
)
args(MAGCAT::magcat_minp_pathways)


## Adaptive soft-TFisher (analytic only, NO permutations)
tf_adapt_res <- magcat_soft_tfisher_adaptive_pathways(
  gene_results     = genes_adj,
  species          = "maize",
  gene_col         = "GENE",
  p_col            = "P",
  tau_grid         = c(0.01, 0.05,0.5,1),
  #min_p            = 1e-15,
  do_fix           = TRUE,
  B_perm           = 10000L,        # <- NO permutations (must be integer)
  perm_mode        = "resample_global",  # ignored when B_perm=0L
  seed             = NULL,
  analytic_logical = TRUE,
  output           = TRUE,
  out_dir          = "magcat_tfisher_adaptive"
)
args(MAGCAT::magcat_soft_tfisher_adaptive_pathways)



head(tf_adapt_res)
colnames(tf_adapt_res)
attr(tf_adapt_res, "file")
args(MAGCAT::magcat_soft_tfisher_adaptive_pathways)
?MAGCAT::magma_geneset_competitive()

### MAGMA-competitive
## First combine the raw files my guy
## Combine per-chromosome MAGMA gene results (.genes.raw-style) into one file
## (works for your custom v1.10 format with "# VERSION" and correlation tail)

# in_dir  <- "/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_multi_snp_wise_genes_by_chr_N_maize"
# pattern <- "^N_maize_MLM_chr[0-9]+\\.multi_snp_wise.genes.raw"  
# out_file <- file.path(in_dir, "N_maize_MLM_ALLCHR.multi_snp_wise.genes")

# files <- list.files(in_dir, pattern = pattern, full.names = TRUE)
# stopifnot(length(files) > 0)

# # order chr1, chr2, ... chr10, ...
# chr_num <- function(x) as.integer(sub(".*_chr([0-9]+)\\..*", "\\1", basename(x)))
# files <- files[order(chr_num(files))]

# # helper: extract header lines + data lines
# read_genes_file <- function(f) {
#   lines <- readLines(f, warn = FALSE)
#   hdr <- lines[grepl("^#", lines)]
#   dat <- lines[!grepl("^#", lines) & nzchar(lines)]
#   list(hdr = hdr, dat = dat)
# }

# x1 <- read_genes_file(files[1])

# # keep header from the first file only
# header_keep <- x1$hdr

# # drop any "# COVAR" lines from subsequent files (if present) and keep only data
# data_all <- x1$dat

# if (length(files) > 1) {
#   for (f in files[-1]) {
#     xi <- read_genes_file(f)
#     data_all <- c(data_all, xi$dat)
#   }
# }

# # write merged
# writeLines(c(header_keep, data_all), con = out_file)

# cat("Wrote:", out_file, "\n")
# cat("Chromosomes merged:", length(files), "\n")



# # EDIT PATHWAY
# library(data.table)
# head(maize_pw)
# # maize_pw: pathway_id, pathway_name, gene_id
# dt <- as.data.table(maize_pw)[, .(set_id = as.character(pathway_id),
#                                  gene_id = sub("^gene:", "", as.character(gene_id)))]

# dt <- dt[!is.na(set_id) & nzchar(set_id) & !is.na(gene_id) & nzchar(gene_id)]
# dt <- unique(dt)

# # (optional but recommended) keep only genes that exist in your .genes.raw
# genes_raw <- "/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_multi_snp_wise_genes_by_chr_N_maize/N_maize_MLM_ALLCHR.multi_snp_wise.genes"
# raw_lines <- readLines(genes_raw, warn = FALSE)
# raw_dat   <- raw_lines[!grepl("^#", raw_lines) & nzchar(raw_lines)]
# genes_ok  <- unique(sub("\\s+.*$", "", raw_dat))
# dt <- dt[gene_id %in% genes_ok]

# # build wide lines: set_id \t gene1 \t gene2 ...
# setlist <- split(dt$gene_id, dt$set_id)
# setlist <- lapply(setlist, unique)

# set_annot_wide <- "annot/maize_pathways.sets.annot.WIDE"
# dir.create(dirname(set_annot_wide), recursive = TRUE, showWarnings = FALSE)

# con <- file(set_annot_wide, open = "wt")
# for (sid in names(setlist)) {
#   genes <- setlist[[sid]]
#   if (length(genes) == 0) next
#   writeLines(paste(c(sid, genes), collapse = "\t"), con = con)
# }
# close(con)

# set_annot_wide


# MAGMA
out <- magma_geneset_competitive(
  gene_results_raw = "/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_multi_snp_wise_genes_by_chr_N_maize/N_maize_MLM_ALLCHR.multi_snp_wise.genes",
  out_prefix       = "N_maize_MLM_ALLCHR.PMN_COMP",
  out_dir          = "magma_geneset",
  species          = "maize",
  genes_all        = genes_all,
  write_tidy       = TRUE
)
args(MAGCAT::magma_geneset_competitive)

### Run OMNI

# Combine adjusted and non-adjusted
head(genes_adj)
head(genes_all)
colnames(genes_adj)=c("GENE","Z_adj","P_adj")
genes_adj$P=genes_all$P
genes_adj$ZSTAT=genes_all$ZSTAT
head(genes_adj)

# Run correlation calculations
raw_dir <- "/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_multi_snp_wise_genes_by_chr_N_maize/"
chr_files <- Sys.glob(file.path(raw_dir, "N_maize_MLM_chr*.multi_snp_wise.genes.raw"))

chr_num <- as.integer(sub(".*_chr([0-9]+)\\..*$", "\\1", basename(chr_files)))
chr_files <- chr_files[order(chr_num)]

out_pairs <- file.path(raw_dir, "magma_gene_cor_pairs_MLM.txt")
if (file.exists(out_pairs)) file.remove(out_pairs)

first <- TRUE
for (f in chr_files) {
  tmp <- paste0(tempfile(), ".txt")

  magma_genesraw_to_cor_pairs_banded(
    genes_raw_file = f,
    out_pairs_file = tmp,
    keep_abs_r_ge  = 0,
    overwrite      = TRUE,
    verbose        = FALSE
  )

  x <- readLines(tmp, warn = FALSE)
  if (!length(x)) next
  if (!first) x <- x[-1]  # drop header
  writeLines(x, out_pairs, useBytes = TRUE, sep = "\n", append = !first)
  first <- FALSE
}

head(x)
out_pairs

x <- x[nzchar(x)]          # <-- drop blank lines

if (!first) x <- x[-1]     # drop header after first file

cat(paste(x, collapse = "\n"), "\n",
    file = out_pairs, append = !first)
first <- FALSE
head(x)

writeLines(x, out_pairs, useBytes = TRUE, sep = "\n")


# Run the whole OMNI on all chromosome
mni_omni <- magcat_omni2_pathways(
  gene_results   = genes_adj,
  species        = "maize",                     # load PMN maize pathways automatically
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
  magma_out      = out,              # MAGMA competitive p-values
  include_magma_in_omni = TRUE,
  include_magma_in_perm = FALSE,                # only for analytic omnibus, no MAGMA in permutations
  omnibus        = "ACAT",                      # or "minP"
  B_perm         = 5000L,                        # number of permutations for omnibus
  perm_mode      = "mvn",                       # or "global", "both", "none"
  magma_cor_file = "/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_multi_snp_wise_genes_by_chr_N_maize/magma_gene_cor_pairs_MLM.txt",  # 3-column file gene1 gene2 r
  make_PD        = TRUE,
  seed           = 123,
  output         = TRUE,
  out_dir        = "magcat_omnibus_results"
)

args(MAGCAT::magcat_omni2_pathways)


res <- magcat_omni2_pathways(
  gene_results = genes_adj,
  species      = "maize",
  gene_col     = "GENE",
  p_raw_col    = "P",
  p_adj_col    = "P_adj",
  z_col        = "Z_adj",
  omnibus      = "ACAT",

  B_perm    = 100L,
  perm_mode = "both",
  magma_out = out,       # data.frame(pathway_id, magma_pvalue)
  include_magma_in_omni = TRUE,     # include as a 6th component
  include_magma_in_perm = FALSE,    # keep FALSE (recommended)

  perm_pool      = "obs",
  magma_cor_file="/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_multi_snp_wise_genes_by_chr_N_maize/magma_gene_cor_pairs_MLM.txt",
  make_PD         = TRUE,
  mvn_marginal    = "uniform",
  seed=123, output=TRUE, out_dir="magcat_omni2_results_MLM"
)
