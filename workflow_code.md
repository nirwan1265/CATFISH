---
title: "MAGCAT Workflow: From GWAS to Pathway Enrichment"
author: "MAGCAT Package"
date: "`r Sys.Date()`"
output:
  github_document:
    toc: true
    toc_depth: 3
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  echo = TRUE,
  eval = FALSE,
  message = FALSE,
  warning = FALSE
)
```

## Overview

This document provides a complete workflow for running pathway enrichment analysis using MAGCAT. The pipeline consists of:

1. **Data Preparation** - Convert GFF3 to gene locations
2. **MAGMA Gene Analysis** - Map SNPs to genes and compute gene-level p-values
3. **Pathway Analysis** - Test pathways using multiple methods (ACAT, Fisher, minP, etc.)
4. **Omnibus Testing** - Combine methods with permutation-based calibration

---

## Prerequisites

### Install MAGCAT

```{r install-magcat}
# Install from GitHub
devtools::install_github("nirwantandukar/MAGCAT")

# Load the package
library(MAGCAT)
```

### Install MAGMA

MAGMA must be installed separately. Download from: https://ctg.cncr.nl/software/magma

```{r set-magma-path}
# Option 1: Add MAGMA to your PATH

# Option 2: Set the path in R
options(magma.path = "/path/to/magma")

# Verify MAGMA is found
magma_path()
```

### Required Files

You will need:

- **GWAS summary statistics** (SNP, CHR, BP, P-value)
- **Reference panel** (PLINK format: .bed/.bim/.fam)
- **Gene annotation** (GFF3 file for your species)

---

## Step 1: Prepare Gene Location File

Convert your GFF3 annotation to MAGMA's gene location format.

```{r step1-geneloc}
# Convert GFF3 to MAGMA gene-loc format
gff3_to_geneloc(
  gff = "data/reference_genome.gff3",
  out = "output/genes.loc",
  feature_type = "gene",
  id_fields = c("gene_id", "ID", "Name"),
  drop_scaffolds = TRUE,
  recode_chr = "auto",        # Auto-recode non-numeric chromosomes
  write_chr_map = TRUE        # Save chromosome mapping
)
```

### For non-numeric chromosomes (e.g., Drosophila 2L/2R)

```{r step1-geneloc-fly}
# Specify chromosome order for species with non-numeric chromosomes
gff3_to_geneloc(
  gff = "data/drosophila.gff3",
  out = "output/fly_genes.loc",
  recode_chr = "order",
  chr_order = c("2L", "2R", "3L", "3R", "4", "X", "Y"),
  write_chr_map = TRUE
)
```

---

## Step 2: MAGMA Annotation

Annotate SNPs to genes based on genomic position.
(This step is done via command line or using MAGCAT's internal functions)

```{r step2-annotate, eval=FALSE}
# Command line example (not R code):
# magma --annotate window=25,25 \
#   --snp-loc snp_locations.txt \
#   --gene-loc output/genes.loc \
#   --out output/magma_annotated
```

---

## Step 3: MAGMA Gene Analysis

Run MAGMA gene analysis to get gene-level p-values from SNP-level GWAS results.

```{r step3-magma-gene}
# Run MAGMA gene analysis
magma_gene(
  bfile = "data/reference_panel",           # PLINK prefix

  gene_annot = "output/magma_annotated.genes.annot",
  stats_file = "data/gwas_summary_stats.txt",
  n_total = 10000,                          # Sample size
  rename_columns = c(
    CHR = "chr",
    SNP = "rs",
    POS = "ps",
    PVALUE = "p_wald"
  ),
  out_prefix = "magma_output",
  out_dir = "output/",
  gene_model = "snp-wise=top"
)
```

### Load Gene Results

```{r step3-load-results}
# Read MAGMA gene results
gene_results <- read.delim(
  "output/magma_output.snp_wise_top.genes.out",
  stringsAsFactors = FALSE
)

# View structure
head(gene_results)
#   GENE CHR START  STOP NSNPS NPARAM     N  ZSTAT       P
# 1 Zm001   1  1000 2000    15     10 10000  2.34  0.0096
# 2 Zm002   1  3000 4500    22     15 10000  1.56  0.0594
# ...
```

---

## Step 4: (Optional) Adjust Gene P-values

Adjust for gene length and number of SNPs to reduce bias.

```{r step4-adjust}
# Get gene lengths from GFF3
gene_lengths <- get_gene_lengths(
  gff3_file = "data/reference_genome.gff3",
  feature = "gene",
  id_key = "ID"
)

# Adjust gene p-values
adjusted <- magcat_adjust_gene_p(

  gene_results = gene_results,
  gene_lengths = gene_lengths,
  gene_col = "GENE",
  nsnp_col = "NSNPS",
  p_col = "P",
  log1p_covars = TRUE
)

# Add adjusted p-values back to gene_results
gene_results$P_adj <- adjusted$p_adj[match(gene_results$GENE, adjusted$gene_id)]
```

---

## Step 5: Load Pathway Definitions

### Option A: Use Built-in PMN Pathways

```{r step5-pmn-pathways}
# Load maize pathways from Plant Metabolic Network
pathways <- magcat_load_pathways(species = "maize")

# Available species: "maize", "sorghum", "arabidopsis", "plant"
head(pathways)
#   pathway_id        pathway_name    gene_id
# 1    PWY-7634  Anthocyanin biosyn   Zm00001
# 2    PWY-7634  Anthocyanin biosyn   Zm00045
# ...
```

### Option B: Use Custom Pathways

```{r step5-custom-pathways}
# As a named list
custom_pathways <- list(
  "pathway_A" = c("gene1", "gene2", "gene3", "gene4"),
  "pathway_B" = c("gene5", "gene6", "gene7"),
  "pathway_C" = c("gene8", "gene9", "gene10", "gene11", "gene12")
)

# Or as a data.frame
custom_pathways_df <- data.frame(
  pathway_id = c("pwy1", "pwy1", "pwy1", "pwy2", "pwy2"),
  gene_id = c("gene1", "gene2", "gene3", "gene4", "gene5"),
  pathway_name = c("Pathway 1", "Pathway 1", "Pathway 1", "Pathway 2", "Pathway 2")
)
```

---

## Step 6: Run Individual Pathway Tests

### ACAT (Cauchy Combination)

```{r step6-acat}
acat_results <- magcat_acat_pathways(
 gene_results = gene_results,
  species = "maize",          # Or: pathways = custom_pathways
  gene_col = "GENE",
  p_col = "P"
)

head(acat_results)
```

### Fisher's Method

```{r step6-fisher}
fisher_results <- magcat_fisher_pathways(
  gene_results = gene_results,
  species = "maize",
  gene_col = "GENE",
  p_col = "P"
)

head(fisher_results)
```

### Minimum P-value (Tippett/Wilkinson)

```{r step6-minp}
minp_results <- magcat_minp_pathways(
  gene_results = gene_results,
  species = "maize",
  gene_col = "GENE",
  p_col = "P"
)

head(minp_results)
```

### Stouffer's Z-score Method

```{r step6-stouffer}
stouffer_results <- magcat_stoufferZ_pathways(
  gene_results = gene_results,
  species = "maize",
  gene_col = "GENE",
  z_col = "ZSTAT",
  alternative = "greater"
)

head(stouffer_results)
```

### Adaptive Soft TFisher

```{r step6-tfisher}
tfisher_results <- magcat_soft_tfisher_adaptive_pathways(
  gene_results = gene_results,
  species = "maize",
  gene_col = "GENE",
  p_col = "P",
  tau_grid = c(0.2, 0.1, 0.05, 0.01, 0.005, 0.001)
)

head(tfisher_results)
```

---

## Step 7: Omnibus Test (Combining Methods)

The omnibus test combines multiple methods and optionally uses permutation-based calibration to account for gene-gene correlations.

### Basic Omnibus (Analytic Only)

```{r step7-omnibus-basic}
omni_results <- magcat_omni2_pathways(
  gene_results = gene_results,
  species = "maize",
  gene_col = "GENE",
  p_raw_col = "P",
  z_col = "ZSTAT",
  perm_mode = "none",         # No permutation
  omnibus = "ACAT"            # Combine methods using ACAT
)

head(omni_results)
```

### Omnibus with Global Resampling

```{r step7-omnibus-global}
omni_global <- magcat_omni2_pathways(
  gene_results = gene_results,
  species = "maize",
  gene_col = "GENE",
  p_raw_col = "P",
  z_col = "ZSTAT",
  perm_mode = "global",
  B_perm = 10000,             # Number of permutations
  seed = 42
)

head(omni_global)
```

### Omnibus with MVN Resampling (LD-aware)

For MVN resampling, you need gene-gene correlation data from MAGMA.

```{r step7-omnibus-mvn}
# First, generate gene correlations from MAGMA (command line):
# magma --bfile reference_panel \
#   --gene-annot output/magma_annotated.genes.annot \
#   --gene-cor output/gene_correlations

# Then run omnibus with MVN
omni_mvn <- magcat_omni2_pathways(
  gene_results = gene_results,
  species = "maize",
  gene_col = "GENE",
  p_raw_col = "P",
  z_col = "ZSTAT",
  perm_mode = "mvn",
  B_perm = 10000,
  magma_cor_file = "output/gene_correlations.genes.raw",
  mvn_marginal = "uniform",
  seed = 42
)

head(omni_mvn)
```

### Full Omnibus (Both Methods)

```{r step7-omnibus-both}
omni_full <- magcat_omni2_pathways(
  gene_results = gene_results,
  species = "maize",
  gene_col = "GENE",
  p_raw_col = "P",
  z_col = "ZSTAT",
  perm_mode = "both",         # Run both global and MVN
  B_perm = 10000,
  magma_cor_file = "output/gene_correlations.genes.raw",
  output = TRUE,              # Save results to CSV
  out_dir = "output/magcat_results"
)
```

---

## Step 8: Interpret Results

### Key Output Columns

| Column | Description |
|--------|-------------|
| `pathway_id` | Pathway identifier |
| `pathway_name` | Human-readable name |
| `n_genes` | Number of genes tested |
| `acat_p` | ACAT p-value |
| `fisher_p` | Fisher's method p-value |
| `minp_p` | Minimum p-value statistic |
| `stouffer_p` | Stouffer Z p-value |
| `omni_p_analytic` | Omnibus p-value (analytic) |
| `omni_p_global` | Omnibus p-value (global resampling) |
| `omni_p_mvn` | Omnibus p-value (MVN resampling) |
| `omni_p_final` | Final omnibus p-value (best available) |

### Filter Significant Pathways

```{r step8-filter}
# Get significant pathways (FDR < 0.05)
sig_pathways <- omni_full[omni_full$omni_p_final < 0.05, ]

# Sort by p-value
sig_pathways <- sig_pathways[order(sig_pathways$omni_p_final), ]

# View top hits
head(sig_pathways[, c("pathway_id", "pathway_name", "n_genes", "omni_p_final")])
```

### Multiple Testing Correction

```{r step8-correction}
# Add Benjamini-Hochberg adjusted p-values
omni_full$fdr <- p.adjust(omni_full$omni_p_final, method = "BH")

# Filter by FDR
sig_fdr <- omni_full[omni_full$fdr < 0.05, ]
```

---

## Complete Example Script

Here's a minimal complete workflow:

```{r complete-example}
library(MAGCAT)

# 1. Load gene results (from MAGMA output)
gene_results <- read.delim("magma_output.genes.out")

# 2. Run omnibus pathway analysis
results <- magcat_omni2_pathways(
  gene_results = gene_results,
  species = "maize",
  gene_col = "GENE",
  p_raw_col = "P",
  z_col = "ZSTAT",
  perm_mode = "global",
  B_perm = 10000,
  seed = 42,
  output = TRUE,
  out_dir = "results/"
)

# 3. Get significant pathways
sig <- results[p.adjust(results$omni_p_final, "BH") < 0.05, ]
print(sig[, c("pathway_name", "n_genes", "omni_p_final")])
```

---

## Session Info

```{r session-info, eval=TRUE}
sessionInfo()
```
