# =============================================================================
# CATFISH Pipeline: Sorghum Stem Volume GWAS
# =============================================================================
# Complete workflow from MAGMA gene-level analysis to pathway enrichment
#
# Pipeline Steps:
#   0. Preprocessing (filter GWAS, convert HMP to PLINK)
#   1. Setup and paths
#   2. MAGMA annotation (SNPs to genes)
#   3. MAGMA gene analysis (gene-level p-values)
#   4. Combine chromosome results
#   5. (Optional) Adjust gene p-values for length/SNP density
#   6. Run CATFISH pathway analysis
#   7. Run omnibus test
# =============================================================================

# -----------------------------------------------------------------------------
# 1. SETUP
# -----------------------------------------------------------------------------

# Load packages
devtools::load_all(".")
library(CATFISH)
library(data.table)
library(dplyr)

# Set tool paths
options(magma.path = "/Users/nirwantandukar/Documents/software/magma_v1.10_mac/magma")
tassel_path <- "/Users/nirwantandukar/Documents/software/tassel-5-standalone/run_pipeline.pl"
plink_path  <- "/Users/nirwantandukar/Documents/software/plink_mac_20231211/plink"

# Verify MAGMA is found
magma_path()

# -----------------------------------------------------------------------------
# USER CONFIGURATION
# -----------------------------------------------------------------------------

# Input files
gwas_file_raw  <- "/Users/nirwantandukar/Documents/Research/results/GWAS/SAP/Stem_diameter/Stem_volume_mod_sub_stem_volume_SAP_bialleles_MAF_0.05_11.assoc.txt"
hmp_dir        <- "/Users/nirwantandukar/Documents/Research/data/SAP/genotype"

# Sample size (SAP panel = 357 individuals)
n_samples      <- 357

# Output settings
out_prefix     <- "sorghum_stem_vol"
out_dir        <- "sorghum_catfish_results"

# Chromosomes to analyze
chroms         <- 1:10  # Sorghum has 10 chromosomes

# Create output directories
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "genotype"), recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# STEP 0: PREPROCESSING
# =============================================================================

# -----------------------------------------------------------------------------
# 0A. Get list of valid SNPs from HMP files
# -----------------------------------------------------------------------------

cat("=== Step 0A: Extracting SNP list from HMP files ===\n")

hmp_files <- list.files(
  path       = hmp_dir,
  pattern    = "^SAP_only_samples_bialleles_MAF_0.05_chr[0-9]+_filtered_80perc_hets_\\.hmp\\.txt$",
  full.names = TRUE
)

cat("Found", length(hmp_files), "HMP files\n")

# Extract SNP IDs from all HMP files (column 1 = rs#)
valid_snps <- c()
for (f in hmp_files) {
  cat("Reading:", basename(f), "\n")
  hmp <- fread(f, header = TRUE, select = 1, nThread = 4)
  valid_snps <- c(valid_snps, hmp[[1]])
}

cat("Total valid SNPs from HMP files:", length(valid_snps), "\n")

# -----------------------------------------------------------------------------
# 0B. Filter GWAS results to only valid SNPs
# -----------------------------------------------------------------------------

cat("\n=== Step 0B: Filtering GWAS results ===\n")

gwas_raw <- fread(gwas_file_raw, header = TRUE)
cat("Original GWAS SNPs:", nrow(gwas_raw), "\n")

gwas_filtered <- gwas_raw[rs %in% valid_snps]
cat("After filtering:", nrow(gwas_filtered), "\n")
cat("Removed:", nrow(gwas_raw) - nrow(gwas_filtered), "SNPs\n")

# Add chromosome prefix to SNP IDs to match genotype file
gwas_filtered[, rs := paste0("chr", chr, "_", rs)]

# Save filtered GWAS
gwas_file <- file.path(out_dir, paste0(out_prefix, "_gwas_filtered.txt"))
fwrite(gwas_filtered, gwas_file, sep = "\t")
cat("Saved filtered GWAS to:", gwas_file, "\n")

# -----------------------------------------------------------------------------
# 0C. Convert HMP to PLINK format (fast vectorized converter)
# -----------------------------------------------------------------------------

cat("\n=== Step 0C: Converting HMP to PLINK format ===\n")

plink_out_dir <- file.path(out_dir, "genotype")

# Fast HMP to PLINK converter using data.table
hmp_to_plink_fast <- function(hmp_file, out_prefix, chr_num) {
  cat("  Converting chr", chr_num, ":", basename(hmp_file), "\n")

  # Read HMP file efficiently
  hmp <- fread(hmp_file, header = TRUE, nThread = 4)

  # Get sample names (columns 12 onwards)
  sample_cols <- 12:ncol(hmp)
  samples <- colnames(hmp)[sample_cols]
  n_samples <- length(samples)
  cat("    Samples:", n_samples, "\n")

  # Get SNP info
  snp_ids <- hmp[["rs#"]]
  alleles <- hmp[["alleles"]]
  positions <- hmp[["pos"]]
  n_snps <- nrow(hmp)

  # Parse alleles (vectorized)
  allele_parts <- tstrsplit(alleles, "/", fixed = TRUE)
  allele1 <- allele_parts[[1]]
  allele2 <- allele_parts[[2]]

  # Make SNP IDs unique by adding chromosome prefix
  snp_ids <- paste0("chr", chr_num, "_", snp_ids)

  # Write MAP file
  map_file <- paste0(out_prefix, ".map")
  map_dt <- data.table(CHR = chr_num, SNP = snp_ids, CM = 0, POS = positions)
  fwrite(map_dt, map_file, sep = "\t", col.names = FALSE)
  cat("    Wrote MAP:", n_snps, "SNPs\n")

  # Convert genotypes to PLINK format (vectorized)
  # Create lookup for each SNP based on its alleles
  # Genotype coding: single letter = homozygous, two letters = het, NA/N = missing

  # Extract genotype columns as matrix for speed
  geno_dt <- hmp[, ..sample_cols]

  # Write PED file sample by sample (memory efficient)
  ped_file <- paste0(out_prefix, ".ped")
  con <- file(ped_file, "w")

  for (i in seq_len(n_samples)) {
    genos <- as.character(geno_dt[[i]])

    # Vectorized conversion
    # Missing data
    is_missing <- genos %in% c("NA", "N", "") | is.na(genos)

    # Single letter = homozygous
    is_single <- nchar(genos) == 1 & !is_missing

    # Two letters = heterozygous
    is_double <- nchar(genos) == 2 & !is_missing

    # Build output
    out_genos <- character(n_snps)
    out_genos[is_missing] <- "0 0"
    out_genos[is_single] <- paste(genos[is_single], genos[is_single])
    out_genos[is_double] <- paste(substr(genos[is_double], 1, 1),
                                   substr(genos[is_double], 2, 2))

    # Any remaining (shouldn't happen)
    remaining <- which(out_genos == "")
    if (length(remaining) > 0) out_genos[remaining] <- "0 0"

    # Write line
    line <- paste(c(samples[i], samples[i], 0, 0, 0, -9, out_genos), collapse = "\t")
    writeLines(line, con)

    if (i %% 50 == 0) cat("\r    Processing sample", i, "/", n_samples)
  }
  close(con)
  cat("\n    Wrote PED:", n_samples, "samples x", n_snps, "SNPs\n")

  # Clean up memory
  rm(hmp, geno_dt)
  gc(verbose = FALSE)

  return(TRUE)
}

# Convert each chromosome
for (chr in chroms) {
  hmp_file <- file.path(hmp_dir,
    paste0("SAP_only_samples_bialleles_MAF_0.05_chr", chr, "_filtered_80perc_hets_.hmp.txt"))

  if (!file.exists(hmp_file)) {
    cat("Warning: HMP file not found for chr", chr, "\n")
    next
  }

  plink_prefix <- file.path(plink_out_dir, paste0("sorghum_sap_chr", chr))
  hmp_to_plink_fast(hmp_file, plink_prefix, chr)
}

# -----------------------------------------------------------------------------
# 0D. Convert each chromosome to binary and merge
# -----------------------------------------------------------------------------

cat("\n=== Step 0D: Converting to binary and merging PLINK files ===\n")

# First convert each chromosome to binary format
for (chr in chroms) {
  ped_file <- file.path(plink_out_dir, paste0("sorghum_sap_chr", chr, ".ped"))
  if (!file.exists(ped_file)) next

  chr_prefix <- file.path(plink_out_dir, paste0("sorghum_sap_chr", chr))

  cat("  Converting chr", chr, "to binary...\n")
  conv_cmd <- paste0(
    plink_path,
    " --file ", chr_prefix,
    " --make-bed",
    " --out ", chr_prefix,
    " --allow-no-sex",
    " --allow-extra-chr"
  )
  system(conv_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
}

# Create merge list file (all chromosomes except first)
merge_list <- file.path(plink_out_dir, "merge_list.txt")
merge_prefixes <- file.path(plink_out_dir, paste0("sorghum_sap_chr", 2:10))
merge_prefixes <- merge_prefixes[file.exists(paste0(merge_prefixes, ".bed"))]

writeLines(merge_prefixes, merge_list)

# Final merged output
bfile <- file.path(plink_out_dir, "sorghum_sap_merged")

# PLINK merge command (now using binary files)
cat("  Merging all chromosomes (first pass)...\n")
plink_cmd <- paste0(
  plink_path,
  " --bfile ", file.path(plink_out_dir, "sorghum_sap_chr1"),
  " --merge-list ", merge_list,
  " --make-bed",
  " --out ", bfile,
  " --allow-no-sex",
  " --allow-extra-chr"
)

result <- system(plink_cmd)

# If merge failed due to multiallelic variants, exclude them and retry
missnp_file <- paste0(bfile, "-merge.missnp")
if (file.exists(missnp_file)) {
  cat("  Found multiallelic variants, excluding and retrying...\n")
  n_missnp <- length(readLines(missnp_file))
  cat("  Excluding", n_missnp, "problematic SNPs\n")

  # Remove problematic SNPs from each chromosome file
  for (chr in chroms) {
    chr_prefix <- file.path(plink_out_dir, paste0("sorghum_sap_chr", chr))
    if (!file.exists(paste0(chr_prefix, ".bed"))) next

    cat("    Filtering chr", chr, "...\n")
    filter_cmd <- paste0(
      plink_path,
      " --bfile ", chr_prefix,
      " --exclude ", missnp_file,
      " --make-bed",
      " --out ", chr_prefix, "_clean",
      " --allow-no-sex",
      " --allow-extra-chr"
    )
    system(filter_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  }

  # Update merge list with cleaned files
  merge_prefixes_clean <- file.path(plink_out_dir, paste0("sorghum_sap_chr", 2:10, "_clean"))
  merge_prefixes_clean <- merge_prefixes_clean[file.exists(paste0(merge_prefixes_clean, ".bed"))]
  writeLines(merge_prefixes_clean, merge_list)

  # Retry merge with cleaned files
  cat("  Merging cleaned files...\n")
  plink_cmd2 <- paste0(
    plink_path,
    " --bfile ", file.path(plink_out_dir, "sorghum_sap_chr1_clean"),
    " --merge-list ", merge_list,
    " --make-bed",
    " --out ", bfile,
    " --allow-no-sex",
    " --allow-extra-chr"
  )
  system(plink_cmd2)
}

cat("\nPLINK files created:", bfile, ".bed/.bim/.fam\n")

# =============================================================================
# CATFISH PIPELINE STARTS HERE
# =============================================================================

# Column names in filtered GWAS file
gwas_columns <- c(
  CHR    = "chr",
  SNP    = "rs",
  POS    = "ps",
  PVALUE = "p_wald"
)

sep <- "\t"

# -----------------------------------------------------------------------------
# 2. GENE LOCATION FILE (Optional - built-in exists for sorghum)
# -----------------------------------------------------------------------------

# The package includes a pre-built sorghum gene location file:
# inst/extdata/sorghum.genes.loc
#
# If you need to create a new one from your own GFF3:
#
# gff3_to_geneloc(
#   gff        = gff3_file,
#   out        = file.path(out_dir, "sorghum.genes.loc"),
#   strict_chr = FALSE
# )

# Use the built-in sorghum gene location file
gene_loc <- system.file("extdata", "sorghum.genes.loc", package = "CATFISH")

# -----------------------------------------------------------------------------
# 3. MAGMA ANNOTATION
# -----------------------------------------------------------------------------

# Annotate SNPs to genes (10kb window upstream/downstream)
ann <- magma_annotate(
  stats_file     = gwas_file,
  rename_columns = gwas_columns,
  gene_loc       = gene_loc,
  out_prefix     = out_prefix,
  out_dir        = file.path(out_dir, "annot"),
  window         = c(10, 10),  # 10kb upstream, 10kb downstream
  sep            = sep,
  nonhuman       = TRUE
)

cat("Annotation complete:", ann$gene_annot, "\n")

# -----------------------------------------------------------------------------
# 4. MAGMA GENE ANALYSIS (by chromosome)
# -----------------------------------------------------------------------------

# Run MAGMA gene analysis per chromosome
magma_gene(
  bfile          = bfile,
  gene_annot     = ann$gene_annot,
  stats_file     = gwas_file,
  sep            = sep,
  n_total        = n_samples,
  rename_columns = gwas_columns,
  out_prefix     = out_prefix,
  out_dir        = file.path(out_dir, "magma_genes"),
  gene_model     = "multi=snp-wise",
  chroms         = chroms,
  n_threads      = min(length(chroms), 6)  # Parallel processing
)

# -----------------------------------------------------------------------------
# 5. COMBINE CHROMOSOME RESULTS
# -----------------------------------------------------------------------------

# Find all chromosome gene output files
magma_out_dir <- file.path(out_dir, "magma_genes")
files <- list.files(
  path       = magma_out_dir,
  pattern    = paste0("^", out_prefix, "_chr.*\\.multi_snp_wise\\.genes\\.out$"),
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No MAGMA output files found. Check if MAGMA ran successfully.")
}

# Read and combine all chromosome files
gene_list <- lapply(files, function(f) {
  utils::read.table(f, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
})
genes_all_raw <- do.call(rbind, gene_list)

# Ensure P column name is correct (MAGMA multi-snp-wise uses P_SNPWISE_TOP1)
if (!"P" %in% names(genes_all_raw)) {
  if ("P_SNPWISE_TOP1" %in% names(genes_all_raw)) {
    genes_all_raw$P <- genes_all_raw$P_SNPWISE_TOP1
  } else if ("P_MULTI" %in% names(genes_all_raw)) {
    genes_all_raw$P <- genes_all_raw$P_MULTI
  } else if ("P.multi" %in% names(genes_all_raw)) {
    genes_all_raw$P <- genes_all_raw$P.multi
  }
}

# Deduplicate by GENE: keep smallest P
genes_all_raw <- genes_all_raw[order(genes_all_raw$GENE, genes_all_raw$P), ]
genes_all <- genes_all_raw[!duplicated(genes_all_raw$GENE), ]

# Sort by chromosome and position
if (all(c("CHR", "START") %in% names(genes_all))) {
  genes_all <- genes_all[order(genes_all$CHR, genes_all$START), ]
}

cat("Combined gene results:", nrow(genes_all), "genes\n")
head(genes_all)

# Save combined results
write.table(
  genes_all,
  file      = file.path(out_dir, paste0(out_prefix, "_genes_combined.txt")),
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 6. (OPTIONAL) ADJUST GENE P-VALUES FOR LENGTH/SNP DENSITY
# -----------------------------------------------------------------------------

# If you have gene length information, you can adjust for confounding:
#
# gene_lengths <- get_gene_lengths(
#   gff3_file  = gff3_file,
#   output     = TRUE,
#   output_dir = out_dir,
#   file_name  = "sorghum_gene_lengths.tsv"
# )
#
# adj_out <- catfish_adjust_gene_p(
#   gene_results = genes_all,
#   gene_lengths = gene_lengths,
#   gene_col     = "GENE",
#   nsnp_col     = "NSNPS",
#   p_col        = "P",
#   z_col        = "ZSTAT"
# )
#
# genes_adj <- adj_out[, c("gene_id", "z_adj", "p_adj")]
# colnames(genes_adj) <- c("GENE", "ZSTAT", "P")

# For now, use unadjusted values
genes_adj <- genes_all

# -----------------------------------------------------------------------------
# 7. RUN CATFISH PATHWAY ANALYSIS
# -----------------------------------------------------------------------------

# Load sorghum pathways
sorghum_pw <- catfish_load_pathways("sorghum", gene_col = "Gene-name")
cat("Loaded", length(unique(sorghum_pw$pathway_id)), "pathways\n")

# --- ACAT ---
cat("\n--- Running ACAT ---\n")
acat_res <- catfish_acat_pathways(
  gene_results = genes_adj,
  species      = "sorghum",
  gene_col     = "GENE",
  p_col        = "P",
  min_p        = 1e-15,
  output       = TRUE,
  out_dir      = file.path(out_dir, "acat")
)

# --- Fisher ---
cat("\n--- Running Fisher ---\n")
fisher_res <- catfish_fisher_pathways(
  gene_results = genes_adj,
  species      = "sorghum",
  gene_col     = "GENE",
  p_col        = "P",
  min_p        = 1e-15,
  output       = TRUE,
  out_dir      = file.path(out_dir, "fisher")
)

# --- Stouffer ---
cat("\n--- Running Stouffer ---\n")
stouffer_res <- catfish_stoufferZ_pathways(
  gene_results = genes_adj,
  species      = "sorghum",
  gene_col     = "GENE",
  z_col        = "ZSTAT",
  alternative  = "greater",
  output       = TRUE,
  out_dir      = file.path(out_dir, "stouffer")
)

# --- MinP ---
cat("\n--- Running MinP ---\n")
minp_res <- catfish_minp_pathways(
  gene_results = genes_adj,
  species      = "sorghum",
  gene_col     = "GENE",
  p_col        = "P",
  min_p        = 1e-15,
  output       = TRUE,
  out_dir      = file.path(out_dir, "minp")
)

# --- Adaptive Soft TFisher ---
cat("\n--- Running Adaptive Soft TFisher ---\n")
tfisher_res <- catfish_soft_tfisher_adaptive_pathways(
  gene_results = genes_adj,
  species      = "sorghum",
  gene_col     = "GENE",
  p_col        = "P",
  tau_grid     = c(0.1, 0.05, 0.02, 0.01, 0.005, 0.001),
  do_fix       = TRUE,
  B_perm       = 0L,  # Analytic only for component test
  output       = TRUE,
  out_dir      = file.path(out_dir, "tfisher")
)

# -----------------------------------------------------------------------------
# 8. EXTRACT GENE CORRELATIONS FOR MVN RESAMPLING
# -----------------------------------------------------------------------------

cat("\n--- Extracting gene correlations for MVN resampling ---\n")

# For LD-aware MVN resampling, extract gene correlations from MAGMA .genes.raw files
raw_files <- list.files(
  path       = magma_out_dir,
  pattern    = paste0("^", out_prefix, "_chr.*\\.multi_snp_wise\\.genes\\.raw$"),
  full.names = TRUE
)

out_pairs <- file.path(out_dir, paste0(out_prefix, "_gene_cor_pairs.txt"))
if (file.exists(out_pairs)) file.remove(out_pairs)

first <- TRUE
for (f in raw_files) {
  cat("Processing:", basename(f), "\n")
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
  if (!first) x <- x[-1]
  writeLines(x, out_pairs, useBytes = TRUE, sep = "\n", append = !first)
  first <- FALSE
}

cat("Gene correlations saved to:", out_pairs, "\n")

# -----------------------------------------------------------------------------
# 9. RUN OMNIBUS TEST
# -----------------------------------------------------------------------------

cat("\n--- Running Omnibus Test ---\n")

omni_results <- catfish_omni2_pathways(
  gene_results   = genes_adj,
  species        = "sorghum",
  pmn_gene_col   = "Gene-name",
  gene_col       = "GENE",
  p_raw_col      = "P",
  z_col          = "ZSTAT",
  weight_col     = NULL,
  tau_grid       = c(0.1, 0.05, 0.02, 0.01, 0.005, 0.001),
  min_p          = 1e-15,
  do_fix         = TRUE,
  stouffer_min_abs_w     = 1e-8,
  stouffer_alternative   = "greater",
  include_magma_in_omni  = FALSE,
  include_magma_in_perm  = FALSE,
  omnibus        = "ACAT",
  B_perm         = 10000L,
  perm_mode      = "mvn",          # MVN resampling (LD-aware)
  magma_cor_file = out_pairs,      # Gene correlation file
  make_PD        = TRUE,
  mvn_marginal   = "uniform",
  mvn_calibrate_components = TRUE,
  seed           = 123,
  output         = TRUE,
  out_dir        = file.path(out_dir, "omnibus")
)

# -----------------------------------------------------------------------------
# 10. RESULTS SUMMARY
# -----------------------------------------------------------------------------

cat("\n=== CATFISH Analysis Complete ===\n")
cat("Results saved to:", out_dir, "\n\n")

# Show top pathways
cat("Top 20 pathways by omnibus p-value:\n")
top_pathways <- head(omni_results[order(omni_results$omni_p_final), ], 20)
print(top_pathways[, c("pathway_id", "pathway_name", "n_genes", "omni_p_final")])

# Add FDR correction
omni_results$fdr <- p.adjust(omni_results$omni_p_final, method = "BH")

# Significant pathways (FDR < 0.05)
sig_pathways <- omni_results[omni_results$fdr < 0.05, ]
cat("\nSignificant pathways (FDR < 0.05):", nrow(sig_pathways), "\n")

if (nrow(sig_pathways) > 0) {
  print(sig_pathways[, c("pathway_id", "pathway_name", "n_genes", "omni_p_final", "fdr")])
}

# Save final results with FDR
write.csv(
  omni_results,
  file      = file.path(out_dir, paste0(out_prefix, "_catfish_omnibus_results.csv")),
  row.names = FALSE
)

cat("\nFinal results saved:", file.path(out_dir, paste0(out_prefix, "_catfish_omnibus_results.csv")), "\n")
