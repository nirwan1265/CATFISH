---
title: "CATFISH tutorial"
subtitle: "Pathway analysis and candidate genes from GWAS summary statistics"
---

# 1. What CATFISH does

CATFISH takes GWAS summary statistics, turns them into gene-level results with
MAGMA, and tests pathways. Every pathway is scored by five tests that suit
different signal shapes (ACAT, Fisher, adaptive soft TFisher, minP, Stouffer),
combined into one omnibus with ACAT. The omnibus is calibrated by a
Westfall--Young procedure: each null replicate draws one genome-wide vector of
correlated gene z-scores, chromosome by chromosome, from the MAGMA gene--gene
correlations, and scores every pathway on that same draw. Each pathway gets a
calibrated p-value and a family-wise p-value that controls the error rate over
all pathways tested. Finally, genes are ranked as candidates from SNP, gene and
pathway evidence.

There are two ways to run it:

| | Function | When |
|---|---|---|
| **Everything in one call** | `catfish_run()` | From GWAS summary statistics, or from MAGMA output you already have, to pathway results and candidate genes |
| **Step by step** | `magma_annotate()`, `magma_gene()`, `magma_read_genes()`, `magma_genesraw_to_cor_pairs_banded()`, `catfish_load_pathways()`, the five component tests, `catfish_wy_pathways()`, `catfish_omni2_pathways()`, `gene_prioritize()` | To run one part, change a setting of one step, or split MAGMA into cluster jobs |

`catfish_run()` calls exactly the step functions, so both ways give the same
results with the same settings.

The GWAS itself (GEMMA, GAPIT, PLINK, ...) is run outside CATFISH.

# 2. Installation

```r
install.packages("remotes")
remotes::install_github("nirwan1265/CATFISH")
```

This installs the required packages (ACAT from GitHub, TFisher and data.table
from CRAN). Optional packages:

```r
install.packages(c("metap", "evir", "BiocManager"))
BiocManager::install(c("qvalue", "IHW", "rtracklayer", "GenomicRanges", "S4Vectors"))
```

MAGMA is a separate program (https://ctg.cncr.nl/software/magma). It is needed
only to run MAGMA from R:

```r
library(CATFISH)
options(magma.path = "/full/path/to/magma")
magma_path()     # checks that it is found
```

# 3. Quick start on simulated data

This runs in a few seconds and needs no MAGMA. It simulates 300 genes on three
chromosomes with neighbouring genes correlated, and puts six associated genes
into one pathway.

```r
library(CATFISH)
set.seed(1)
genes <- data.frame(GENE = paste0("g", 1:300), CHR = rep(1:3, each = 100))
R <- diag(100)
R[cbind(1:99, 2:100)] <- R[cbind(2:100, 1:99)] <- 0.3
z <- as.vector(replicate(3, t(chol(R)) %*% rnorm(100)))
z[1:6] <- z[1:6] + 4
genes$ZSTAT <- z
genes$P     <- pnorm(z, lower.tail = FALSE)
cor_pairs <- do.call(rbind, lapply(split(genes$GENE, genes$CHR), function(g)
  data.frame(gene1 = head(g, -1), gene2 = tail(g, -1), r = 0.3)))
pathways <- c(list(PWY_signal = paste0("g", 1:15)),
              setNames(lapply(1:19, function(i) sample(genes$GENE[-(1:15)], 15)),
                       paste0("PWY_", 1:19)))

res <- catfish_wy_pathways(genes, pathways = pathways, cor_pairs = cor_pairs,
                           B = 200, min_genes = 10, seed = 1)
head(res[, c("pathway_id", "n_genes", "z_std", "p_calibrated", "p_fwer_std")])

cand <- gene_prioritize(magma = genes, pathway_results = res, pathways = pathways)
head(cand[, c("rank", "gene", "score", "magma_rank", "pathway_rank", "best_pathway")])
```

Only `PWY_signal` should have `p_fwer_std <= 0.05`, and the six associated genes
should rank first.

# 4. Everything in one call: `catfish_run()`

## 4.1 Starting from GWAS summary statistics

```r
res <- catfish_run(
  gwas_file      = "dry_tons.assoc.txt",          # GEMMA output, one SNP per row
  rename_columns = c(CHR = "chr", SNP = "rs", POS = "ps",
                     PVALUE = "p_wald", NMISS = "n_miss"),
  n_total        = 327,                           # lines in the GWAS
  bfile          = "bap_genotypes",               # PLINK LD reference
  species        = "sorghum",                     # built-in gene locations and pathways
  window         = c(25, 25),                     # kb around each gene
  gene_model     = "snp-wise=mean",
  chroms         = 1:10, n_threads = 4,           # one MAGMA run per chromosome
  B              = 10000,
  min_genes      = 5, max_genes = 200,
  out_dir        = "dry_tons", out_prefix = "dry_tons"
)
```

`rename_columns` maps CATFISH's names to your file's columns: `CHR`, `SNP`,
`POS`, `PVALUE` are required; the sample size is either a per-SNP column (`N`,
`NOBS`, or `NMISS` together with `n_total`) or `n_total` alone.

## 4.2 Starting from MAGMA output

On a cluster, MAGMA is usually run one chromosome per job (section 5.3). Point
`catfish_run()` at the folder of `.genes.out` and `.genes.raw` files; MAGMA is
then not run:

```r
res <- catfish_run(
  magma_dir     = "magma",
  magma_pattern = "^dry_tons_chr.*\\.genes\\.out$",
  species       = "sorghum",
  gwas_file     = "dry_tons.assoc.txt",           # optional: SNP layer of the gene ranking
  rename_columns = c(CHR = "chr", SNP = "rs", POS = "ps", PVALUE = "p_wald"),
  gene_annot    = "magma/dry_tons.genes.annot",   # optional: SNP-to-gene map for that layer
  out_dir       = "dry_tons"
)
```

## 4.3 What it does, and what it returns

| Step | Done by | Notes |
|---|---|---|
| 1 | `magma_annotate()`, `magma_gene()` | only with `gwas_file` and no `magma_dir` |
| 2 | `magma_read_genes()` | combines all chromosomes; prints the median of every gene p column |
| 3 | `magma_genesraw_to_cor_pairs_banded()` | gene--gene correlations of all chromosomes in one file; reused on a rerun |
| 4 | `catfish_wy_pathways()` | Westfall--Young pathway test |
| 5 | `catfish_omni2_pathways()`, analytic only | the five component p-values, added to the pathway table (`components = TRUE`) |
| 6 | `gene_prioritize()` | candidate genes (`prioritize = TRUE`) |

The gene p-value column is `P_SNPWISE_MEAN` when present (a `multi=snp-wise`
run) and `P` otherwise; set `p_col` to choose another.

`res` is a list:

| Element | Content |
|---|---|
| `res$pathways` | one row per pathway: Westfall--Young results plus `acat_p`, `fisher_p`, `tfisher_p`, `minp_p`, `stouffer_p` |
| `res$candidates` | ranked genes |
| `res$genes` | gene-level results used |
| `res$pathway_sets` | pathway memberships used |
| `res$settings` | `p_col`, `gene_regex`, `B`, `tau_grid`, ... |
| `res$files` | everything written to `out_dir` (`<out_prefix>_pathways.csv`, `<out_prefix>_candidate_genes.csv`, the correlation file, MAGMA output) |

# 5. Step by step

## 5.1 Gene locations

```r
gff3_to_geneloc(gff = "reference.gff3", out = "species.genes.loc")
```

Built-in files exist for maize, sorghum and arabidopsis (`species =` in
`magma_annotate()` and `catfish_run()`). Non-numeric chromosome names can be
recoded to numbers (`recode_chr = "auto"`), with a map file written beside the
output.

## 5.2 MAGMA annotation (SNPs to genes)

```r
ann <- magma_annotate(
  stats_file     = "dry_tons.assoc.txt",
  rename_columns = c(CHR = "chr", SNP = "rs", POS = "ps", PVALUE = "p_wald"),
  gene_loc       = "species.genes.loc",
  out_prefix     = "dry_tons",
  out_dir        = "magma",
  window         = c(25, 25),
  nonhuman       = TRUE
)
ann$gene_annot     # the .genes.annot file for the next step
```

## 5.3 MAGMA gene analysis

```r
magma_gene(
  bfile          = "bap_genotypes",
  gene_annot     = ann$gene_annot,
  stats_file     = "dry_tons.assoc.txt",
  rename_columns = c(CHR = "chr", SNP = "rs", PVALUE = "p_wald", NMISS = "n_miss"),
  n_total        = 327,
  out_prefix     = "dry_tons",
  out_dir        = "magma",
  gene_model     = "snp-wise=mean",
  chroms         = 1:10,      # writes dry_tons_chr1.snp_wise_mean.genes.out, ...
  n_threads      = 4          # 1 = one chromosome after another
)
```

For a cluster array job, give the job's chromosome: `chroms = 5`.

## 5.4 Gene results and gene--gene correlations

```r
genes <- magma_read_genes(dir = "magma", pattern = "^dry_tons_chr.*\\.genes\\.out$")

raw <- sub("\\.genes\\.out$", ".genes.raw", attr(genes, "files"))
magma_genesraw_to_cor_pairs_banded(raw, "dry_tons_gene_cor_pairs.txt",
                                   gene_regex = "^SORBI")
```

`magma_read_genes()` prints the median of every gene p-value column; about 0.5
is expected where most genes carry no signal (see section 6). MAGMA computes
gene--gene correlations only within a chromosome, so all chromosome files are
needed. `gene_regex` must match your gene IDs: `"^SORBI"`, `"^Zm"`, `"^AT"`,
`"^FBgn"`, `"^[0-9]+$"` (Entrez).

## 5.5 Pathways

```r
pw <- catfish_load_pathways(species = "sorghum")
# built in: "maize", "sorghum", "arabidopsis", "plant", "fly", "human", "human_reactome"

pw <- catfish_load_pathways(file = "sorghumbicolorcyc_pathways.tsv")   # your PMN file
pw <- catfish_load_gmt("h.all.v2024.1.Hs.entrez.gmt")                  # GMT
pw <- catfish_load_biocyc("humancyc_export.csv", gene_col = "NCBI Gene ID")
```

A pathway set is either a data.frame with `pathway_id`, `gene_id` (and
optionally `pathway_name`) or a named list of gene vectors. Gene IDs are
matched to the gene results without regard to case. `catfish_write_set_annot()`
writes a set as a MAGMA `--set-annot` file.

## 5.6 The five component tests, one at a time

Each takes the gene results and the pathways, and returns one row per pathway.
These p-values are analytic: they treat genes as independent, so with
correlated genes they are anti-conservative. Use them to see which kind of
signal drives a pathway; use the Westfall--Young results (5.7) for inference.

| Test | Function | Sensitive to | Statistic |
|---|---|---|---|
| ACAT | `catfish_acat_pathways()` | a few very strong genes | Cauchy combination of gene p-values |
| Fisher | `catfish_fisher_pathways()` | many moderate genes | -2 sum log p |
| Adaptive soft TFisher | `catfish_soft_tfisher_adaptive_pathways()` | strong drivers plus moderate support | Fisher on p-values below a threshold tau, best tau of a grid |
| minP | `catfish_minp_pathways()` | one dominant gene | smallest gene p, Sidak-corrected |
| Stouffer | `catfish_stoufferZ_pathways()` | a weak shift across the whole pathway | sum of gene z-scores / sqrt(n) |

```r
acat  <- catfish_acat_pathways(genes, pathways = pw, gene_col = "GENE", p_col = "P",
                               min_p = 1e-15)
fish  <- catfish_fisher_pathways(genes, pathways = pw, p_col = "P")
tfish <- catfish_soft_tfisher_adaptive_pathways(genes, pathways = pw, p_col = "P",
                                                tau_grid = c(0.01, 0.05, 0.5))
minp  <- catfish_minp_pathways(genes, pathways = pw, p_col = "P")
stou  <- catfish_stoufferZ_pathways(genes, pathways = pw, z_col = "ZSTAT",
                                    weight_col = NULL, alternative = "greater")
```

Parameters they share:

| Parameter | Meaning |
|---|---|
| `gene_results` | gene table (`GENE`, `P`, ...) |
| `pathways` / `species` | pathway set, or a built-in one |
| `pmn_gene_col` | gene column of a built-in file |
| `gene_col`, `p_col` | gene ID and p-value columns (Stouffer: `z_col`) |
| `min_p` | p-values below this are set to it before combining (default `1e-15`) |
| `do_fix` | clean p-values to lie strictly in (0, 1) (default `TRUE`) |
| `output`, `out_dir` | write a CSV |

Test-specific parameters:

| Function | Parameter | Meaning |
|---|---|---|
| soft TFisher | `tau_grid` | thresholds; a numeric vector, or `"auto"` (quantiles of all gene p-values) |
| Stouffer | `z_col` | gene z column (default `ZSTAT`) |
| Stouffer | `weight_col` | optional per-gene weights (e.g. sqrt(NSNPS)); default equal |
| Stouffer | `min_abs_w` | replaces missing or non-positive weights |
| Stouffer | `alternative` | `"greater"` (enrichment; default), `"two.sided"`, `"less"` |

Output columns: `pathway_id`, `pathway_name`, `n_genes`, `gene_names`, and the
test's own: `acat_p`; `fisher_p`; `tau_hat`, `tfisher_stat_hat`,
`tfisher_p_hat`, `tfisher_p_min`, `tfisher_p_omni`; `minp_stat`, `minp_p`,
`minp_method`; `stouffer_Z`, `stouffer_p`.

## 5.7 Westfall--Young pathway test

```r
res <- catfish_wy_pathways(
  gene_results = genes,
  pathways     = pw,
  gene_col     = "GENE", p_col = "P", chr_col = "CHR",
  cor_file     = "dry_tons_gene_cor_pairs.txt",
  B            = 10000,
  use_tfisher  = TRUE,
  tau_grid     = c(0.01, 0.05, 0.5),
  min_genes    = 5, max_genes = 200,
  seed         = 20260706
)
attr(res, "wy_threshold_05")     # omnibus p needed for FWER 5% on the raw scale
```

The omnibus combines ACAT, Fisher, Stouffer, minP and (with `use_tfisher`) soft
TFisher with ACAT, all in log space, so gene p-values far below 1e-300 keep
their order. `tau_grid` must be fixed in advance: a grid estimated from the
gene p-values would make each pathway's statistic depend on genes outside it and
break the family-wise control. Without `cor_file` (or `cor_pairs`) every gene is
treated as independent, and CATFISH warns.

Before testing, the gene z-scores are centred at their median (gene-level
genomic control; `gc = TRUE` by default). A GWAS can carry a genome-wide shift in its gene
z-scores, shared across chromosomes, that the correlated null does not contain
and that a self-contained test would read as pathway signal. Centring is one
constant per GWAS: the order of genes and the gene--gene correlations are
unchanged, and the shift removed is returned as `gc_shift`. Pass the whole
genome's gene results, because the median is taken over the genes supplied.

Centring cannot tell a technical shift from real polygenic signal. When the
median gene is itself associated (a highly polygenic, well-powered GWAS),
centring asks whether a pathway is more associated than the median gene, which
is closer to a competitive test, and CATFISH prints a note when the median gene
z exceeds 0.5. Use `gc = FALSE` for the strict self-contained test, for gene
results that do not cover the whole genome, for gene p-values already corrected
in another way, or to reproduce uncentred results.

## 5.8 Components and calibration together: `catfish_omni2_pathways()`

The same Westfall--Young calibration, with the component p-values, multiple
testing adjustments and diagnostics in one table:

```r
om <- catfish_omni2_pathways(
  gene_results = genes, pathways = pw,
  gene_col = "GENE", p_raw_col = "P", z_col = "ZSTAT", chr_col = "CHR",
  min_genes = 5, max_genes = 200,
  perm_mode = "wy", B_perm = 10000,
  wy_tau_grid = c(0.01, 0.05, 0.5),
  magma_cor_file = "dry_tons_gene_cor_pairs.txt",
  seed = 20260706
)
```

`perm_mode = "none"` gives the analytic table only; `"global"` resamples genes
from a pooled distribution (a faster heuristic that ignores gene correlation);
`"both"` runs both. With `perm_mode = "wy"`, `omni_p_cal` and
`omni_p_fwer_wy_std` equal `p_calibrated` and `p_fwer_std` of
`catfish_wy_pathways()` for the same seed.

## 5.9 Candidate genes

```r
cand <- gene_prioritize(
  magma           = genes,  magma_p_col = "P",
  pathway_results = res,                       # from 5.7
  pathways        = pw,
  gwas            = "dry_tons.assoc.txt",      # optional SNP layer
  gwas_snp_col    = "rs", gwas_p_col = "p_wald",
  gene_annot      = ann$gene_annot,            # SNP-to-gene map (or gene_loc + window)
  targets         = c("GENE_A", "GENE_B")      # optional known genes
)
attr(cand, "targets")     # where the known genes rank
```

Each gene's score is

    score = w_gwas * S_gwas + w_magma * S_magma + w_pathway * S_pathway
            + w_hub * log10(1 + n_pathways)

with default weights 1, 1, 0.5 and 0.3. Each S is a rank score, -log10(rank /
N): the strongest SNP in or near the gene, the MAGMA gene p-value, and the
strongest tested pathway containing the gene (by `z_std`). SNPs are assigned to
genes by `gene_annot` (the MAGMA annotation) or, failing that, by position
within `window` bp of each gene in `gene_loc`. Give the SNP layer as -log10(p)
(`gwas_neglog10p_col`) when the p column is floored (e.g. at 1e-300), since a
floor ties the strongest genes.

# 6. Choosing the settings

| Setting | Recommendation | Why |
|---|---|---|
| MAGMA gene model | `snp-wise=mean` | The Westfall--Young null assumes gene p-values are uniform where there is no association. The SNP-wise mean model is calibrated when the GWAS is; the SNP-wise top model, and `multi=snp-wise` which includes it, are often inflated in panels with long-range LD (median gene p of 0.11 and 0.18 against 0.45 in a sorghum panel), and inflation is read as pathway signal. Check the medians printed by `magma_read_genes()`. |
| `window` | 25 kb in plants with extended LD; smaller (a few kb) where LD decays fast | SNPs outside genes are assigned to genes within the window. |
| `B` | 10,000 for final results; 200--1,000 to try settings | The smallest calibrated p is about 1/(B+1). |
| `tau_grid` | `c(0.01, 0.05, 0.5)` | Fixed in advance, as required for family-wise control. |
| `min_genes`, `max_genes` | 5--10 and 200 | Very small pathways are dominated by one gene; very large ones are unspecific. |
| `seed` | any fixed integer | Makes the null draws reproducible. |

Gene-length adjustment (`catfish_adjust_gene_p()`, `get_gene_lengths()`) is not
needed for MAGMA gene p-values, which are calibrated for gene size. It is for
gene p-values from methods that are not, such as an uncorrected minimum SNP p
per gene.

# 7. Output columns

## Westfall--Young (`catfish_wy_pathways()`, `res$pathways`)

| Column | Meaning |
|---|---|
| `pathway_id`, `pathway_name` | pathway |
| `z_std` | observed omnibus standardised against the pathway's own null; more negative is stronger. Results are sorted by it. |
| `p_fwer_std` | Westfall--Young family-wise p-value; the headline result |
| `p_calibrated` | per-pathway calibrated p-value |
| `n_genes` | genes of the pathway present in the gene results |
| `log10_p_omnibus`, `p_omnibus` | analytic (uncalibrated) omnibus |
| `p_fwer_wy` | family-wise p-value on the raw omnibus scale |
| `p_bonferroni` | Bonferroni on the analytic omnibus, for comparison |
| `gc_shift` | median gene z removed by genomic control (same in every row; 0 with `gc = FALSE`) |
| `acat_p`, `fisher_p`, `tfisher_p`, `minp_p`, `stouffer_p` | analytic component p-values (`catfish_run()` only) |

## `catfish_omni2_pathways()`

| Column | Meaning |
|---|---|
| `acat_p`, `fisher_p`, `tfisher_p_analytic`, `minp_p_analytic`, `stouffer_p_analytic` | analytic components |
| `omni_p_analytic` | analytic omnibus |
| `omni_p_cal`, `omni_p_fwer_wy_std`, `omni_p_fwer_wy`, `omni_z_std` | Westfall--Young results (as above) |
| `omni_p_global` | global-resampling p-value (`perm_mode = "global"` or `"both"`) |
| `omni_p_final` | best available: Westfall--Young calibrated, else global, else analytic |
| `*_BH`, `omni_p_final_q_storey`, `omni_p_final_q_bky10`, `omni_p_final_q_ihw10` | multiple-testing adjustments |
| `dominant_component` | component with the smallest p-value |
| `agreement_score`, `agreement_k`, `agreement_m` | fraction / number of components with p < 0.05 |
| `genes_used`, `gene_pvals_used` | genes of the pathway and their p-values |

## `gene_prioritize()`

| Column | Meaning |
|---|---|
| `rank`, `gene`, `symbol`, `score` | final rank and score |
| `n_layers` | layers (SNP, gene, pathway) that support the gene |
| `gwas_score`, `magma_score`, `pathway_score`, `hub_bonus` | the parts of the score |
| `gwas_neglog10p`, `gwas_rank` | strongest SNP of the gene |
| `magma_p`, `magma_rank` | MAGMA gene p-value and its rank |
| `best_pathway`, `pathway_stat`, `pathway_rank` | strongest tested pathway containing the gene |
| `n_pathways` | tested pathways containing the gene |

# 8. Function reference

Every exported function with all of its parameters, taken from the package
documentation (`?function_name` in R shows the same).

### `catfish_run()`

*Run the whole CATFISH analysis in one call*

```r
catfish_run(
  gwas_file = NULL,
  rename_columns = NULL,
  bfile = NULL,
  gene_loc = NULL,
  n_total = NULL,
  window = c(25, 25),
  gene_model = "snp-wise=mean",
  chroms = NULL,
  n_threads = 1L,
  nonhuman = TRUE,
  sep = NULL,
  chr_map_path = NULL,
  magma_dir = NULL,
  magma_pattern = NULL,
  gene_annot = NULL,
  pathways = NULL,
  species = NULL,
  pmn_gene_col = NULL,
  p_col = NULL,
  gene_regex = NULL,
  B = 10000L,
  use_tfisher = TRUE,
  tau_grid = c(0.01, 0.05, 0.5),
  min_genes = 10L,
  max_genes = 200L,
  seed = 20260706L,
  gc = TRUE,
  components = TRUE,
  prioritize = TRUE,
  targets = NULL,
  w_gwas = 1,
  w_magma = 1,
  w_pathway = 0.5,
  w_hub = 0.3,
  out_dir = "catfish_results",
  out_prefix = "catfish",
  verbose = TRUE
)
```

| Parameter | Description |
|---|---|
| `gwas_file` | GWAS summary statistics file (one SNP per row). |
| `rename_columns` | Named character vector mapping the standard names to the columns of `gwas_file`: `CHR`, `SNP`, `POS`, `PVALUE`, and optionally a per-SNP sample size (`N`, `NOBS`, or `NMISS` with `n_total`). Example: `c(CHR = "chr", SNP = "rs", POS = "ps", PVALUE = "p_wald", NMISS = "n_miss")`. |
| `bfile` | PLINK prefix of the LD reference genotypes (`.bed/.bim/.fam`). |
| `gene_loc` | MAGMA gene-location file. If `NULL`, the built-in file of `species` is used (maize, sorghum, arabidopsis). |
| `n_total` | Total sample size, when `rename_columns` has no per-SNP `N`/`NOBS`, or with `NMISS`. |
| `window` | SNP-to-gene window in kb, `c(upstream, downstream)`. |
| `gene_model` | MAGMA gene model (one). `"snp-wise=mean"` is recommended: see Details. |
| `chroms` | Chromosomes to analyse one MAGMA run at a time (e.g. `1:10`). `NULL` runs MAGMA once on the whole genome. |
| `n_threads` | Parallel MAGMA runs when `chroms` is given. |
| `nonhuman` | Passed to MAGMA annotation (`TRUE` for non-human data). |
| `sep` | Field separator of `gwas_file`; inferred from the extension if `NULL`. |
| `chr_map_path` | Optional chromosome-renaming map (see `gff3_to_geneloc`). |
| `magma_dir` | Folder with existing MAGMA output (`.genes.out` and `.genes.raw`). |
| `magma_pattern` | Regular expression selecting the `.genes.out` files in `magma_dir`; default all of them. |
| `gene_annot` | Optional MAGMA `.genes.annot` for the SNP layer of the candidate gene ranking when starting from `magma_dir` with a `gwas_file`. Without it, SNPs are assigned to genes by position. |
| `pathways` | Pathways: named list of gene vectors, or data.frame with `pathway_id`, `gene_id` (optional `pathway_name`). Give this or `species`. |
| `species` | Built-in pathways (and gene locations), e.g. `"sorghum"`; see `catfish_load_pathways`. |
| `pmn_gene_col` | Gene column of a built-in pathway file. |
| `p_col` | Gene p-value column. Default: `P_SNPWISE_MEAN` when present (a `multi=snp-wise` run), otherwise `P`. |
| `gene_regex` | Regular expression identifying gene IDs in the `.genes.raw` files; worked out from the gene IDs if `NULL`. |
| `B` | Westfall--Young null replicates. |
| `use_tfisher` | Include soft TFisher in the omnibus. |
| `tau_grid` | Fixed soft TFisher grid (also used for the component table). |
| `min_genes`, `max_genes` | Pathway size limits (genes present). |
| `seed` | Seed for the null draws. |
| `gc` | Gene-level genomic control (default `TRUE`); passed to `catfish_wy_pathways` and the component tests. |
| `components` | Add the five analytic component p-values to the pathway table. |
| `prioritize` | Rank candidate genes. |
| `targets` | Optional known genes; their ranks are reported. |
| `w_gwas`, `w_magma`, `w_pathway`, `w_hub` | Layer weights of `gene_prioritize`. |
| `out_dir` | Output folder. |
| `out_prefix` | Prefix of the output files. |
| `verbose` | Print progress. |

### `gff3_to_geneloc()`

*Convert a GFF3 file to a MAGMA gene location file (optionally recoding CHR to 1..K)*

```r
gff3_to_geneloc(
  gff,
  out,
  feature_type = "gene",
  id_fields = c("gene_id", "ID", "Name"),
  drop_scaffolds = TRUE,
  scaffold_regex = "^(scaf|scaffold|un|unplaced|chrUn|KI|GL|NW_|NT_)",
  chr_prefix = "chr",
  recode_chr = c("auto", "none", "order", "map"),
  chr_order = NULL,
  chr_map = NULL,
  strict_chr = FALSE,
  write_chr_map = TRUE,
  chr_map_out = NULL
)
```

| Parameter | Description |
|---|---|
| `gff` | Path to a GFF3 file. |
| `out` | Path to write the MAGMA gene location file (tab-delimited). |
| `feature_type` | GFF3 feature type(s) to use (default "gene"). |
| `id_fields` | Candidate attribute tags to use as gene IDs (searched in order), e.g. c("gene_id","ID","Name"). |
| `drop_scaffolds` | Logical; if TRUE, drop seqids that look like scaffolds. |
| `scaffold_regex` | Regex for scaffold/unplaced contigs to drop when drop_scaffolds=TRUE. |
| `chr_prefix` | Optional string to strip from start of seqid when reading (e.g. "chr"). |
| `recode_chr` | One of "none","auto","order","map". |
| `chr_order` | If recode_chr="order", the desired order of original chromosome labels. |
| `chr_map` | If recode_chr="map", a named vector mapping original->new numeric labels. |
| `strict_chr` | If TRUE, error when encountering CHR labels not in chr_order/chr_map. |
| `write_chr_map` | If TRUE, write a mapping file. |
| `chr_map_out` | Optional explicit path for the mapping file. Default: paste0(out, ".chr_map.tsv") |

### `magma_path()`

*Locate MAGMA binary*

```r
magma_path()
```

### `magma_annotate()`

*Run MAGMA SNP-to-gene annotation*

```r
magma_annotate(
  stats_file,
  rename_columns,
  gene_loc = NULL,
  out_prefix,
  out_dir = NULL,
  window = NULL,
  species = NULL,
  sep = NULL,
  nonhuman = NULL,
  filter_snps = NULL,
  chr_map_path = NULL,
  strict_chr = TRUE,
  verbose = TRUE
)
```

| Parameter | Description |
|---|---|
| `stats_file` | Path to the GWAS summary statistics file. |
| `rename_columns` | Named character vector mapping required MAGMA columns to the column names in `stats_file`; must contain at least `CHR`, `SNP`, and `POS` (e.g. `c(CHR = "chr", SNP = "rs", POS = "ps", PVALUE = "p")`). |
| `gene_loc` | Path to a MAGMA gene-location file (e.g. from `gff3_to_geneloc`). |
| `out_prefix` | Output file prefix. |
| `out_dir` | Optional output directory. Created if it does not exist. |
| `window` | Optional numeric vector `c(up, down)` giving the upstream/downstream window in kb around each gene. |
| `species` | Optional species tag (used to pick a built-in gene-location file when `gene_loc` is not supplied). |
| `sep` | Field separator for `stats_file`. If `NULL`, inferred from the file extension. |
| `nonhuman` | Logical or `NULL`; set `TRUE` for non-human data so MAGMA does not assume human chromosome coding. |
| `filter_snps` | Optional vector or file of SNPs to restrict the analysis to. |
| `chr_map_path` | Optional path to a two-column chromosome-renaming map. |
| `strict_chr` | Logical; if `TRUE`, error when the chromosome map is missing an entry rather than passing the value through unchanged. |
| `verbose` | Logical; print progress messages. |

### `magma_gene()`

*Run MAGMA gene analysis on SNP p-values*

```r
magma_gene(
  bfile,
  gene_annot,
  pval_file = NULL,
  stats_file = NULL,
  n_total,
  snp_col = "SNP",
  p_col = "P",
  rename_columns = NULL,
  out_prefix,
  out_dir = NULL,
  gene_model = "snp-wise=top",
  chroms = NULL,
  n_threads = 1,
  chr_keep = NULL,
  sep = NULL,
  chr_map_path = NULL,
  strict_chr = TRUE,
  verbose = TRUE
)
```

| Parameter | Description |
|---|---|
| `bfile` | PLINK prefix for LD reference. |
| `gene_annot` | MAGMA gene annotation file (.genes.annot). |
| `pval_file` | Optional; SNP p-value file for MAGMA's --pval. If NULL, `stats_file` must be provided. |
| `stats_file` | Optional; GWAS summary statistics file containing at least SNP ID and p-value columns. If provided, `pval_file` is ignored and a MAGMA p-value file is generated internally. |
| `n_total` | Total sample size (N=... in MAGMA) when using a constant N or computing per-SNP N from NMISS. Not required if `rename_columns` supplies an `N` or `NOBS` column. |
| `snp_col` | Column name containing SNP IDs (used only when `pval_file` is given and `rename_columns` is NULL). |
| `p_col` | Column name containing p-values (used only when `pval_file` is given and `rename_columns` is NULL). |
| `rename_columns` | Optional named character vector mapping standard names to the column names in `stats_file`, e.g. `c(CHR = "chr", SNP = "rs", POS = "ps", PVALUE = "p_wald", NMISS = "n_miss")` or `c(..., N = "N_eff")` or `c(..., NOBS = "Nobs")`. For `magma_gene()`, only `CHR`, `SNP`, `PVALUE`, `N`, `NOBS`, and `NMISS` are used. |
| `out_prefix` | Base output prefix (file name without directory). |
| `out_dir` | Optional output directory. If not NULL, the final prefix will be file.path(out_dir, out_prefix). The directory is created if it does not exist. For `stats_file` mode, the internal p-value file is also written into this directory. |
| `gene_model` | Character vector of one or more MAGMA gene models, e.g. "snp-wise=top", "snp-wise=mean", "multi". |
| `chroms` | Optional vector of chromosomes to analyse separately (e.g. 1:10), each written as `<out_prefix>_chr<k>`. Only used in `stats_file` mode. With `n_threads > 1` the chromosomes run in parallel, otherwise one after another. In a cluster array job give the job's single chromosome, e.g. `chroms = 5`. |
| `n_threads` | Number of parallel workers to use when `chroms` is provided. Defaults to 1 (no parallelism). Will be clamped to `min(length(chroms), detectCores() - 1)`. |
| `chr_keep` | Internal. Optional subset of chromosomes to retain. |
| `sep` | Field separator for `stats_file`: `","`, `" "`, or a tab. If `NULL`, inferred from the file extension. |
| `chr_map_path` | Optional path to a two-column chromosome-renaming map. |
| `strict_chr` | Logical; if `TRUE`, error when the chromosome map is missing an entry rather than passing the value through unchanged. |
| `verbose` | Logical; print progress messages. |

### `magma_read_genes()`

*Read and combine MAGMA gene results*

```r
magma_read_genes(
  files = NULL,
  dir = NULL,
  pattern = "\\.genes\\.out$",
  verbose = TRUE
)
```

| Parameter | Description |
|---|---|
| `files` | Character vector of `.genes.out` paths. Provide this or `dir`. |
| `dir` | Directory to search for files matching `pattern`. |
| `pattern` | Regular expression for the file names when `dir` is used. |
| `verbose` | Print the number of genes and chromosomes and the median of each gene p-value column. |

### `magma_genesraw_to_cor_pairs_banded()`

*Extract gene-gene correlations from MAGMA .genes.raw files*

```r
magma_genesraw_to_cor_pairs_banded(
  genes_raw_file,
  out_pairs_file,
  fixed_fields = NULL,
  gene_regex = "^FBgn",
  keep_abs_r_ge = 0,
  overwrite = TRUE,
  verbose = TRUE
)
```

| Parameter | Description |
|---|---|
| `genes_raw_file` | Path(s) to MAGMA `.genes.raw` files. Give one file per chromosome when MAGMA was run chromosome by chromosome; all are written to one output file. MAGMA computes gene--gene correlations only within a chromosome, so the files together hold the full correlation structure. |
| `out_pairs_file` | Output path for the correlation pairs file. |
| `fixed_fields` | Integer, or `NULL` (default) to read it from each file's `# COVAR` header line: the seven fixed MAGMA fields (GENE, CHR, START, STOP, NSNPS, NPARAM, N) plus one per covariate. For MAGMA 1.10 output (`# COVAR = NSAMP MAC`) this is 9 for every gene model. Used when a file has no `# COVAR` line. |
| `gene_regex` | Regular expression to identify gene IDs. Use organism-specific patterns: • Maize: `"^Zm"` • Fly: `"^FBgn"` • Arabidopsis: `"^AT"` • Sorghum: `"^SORBI"` |
| `keep_abs_r_ge` | Numeric. Only keep correlations with \|r\| >= this value. Default 0 keeps all correlations. |
| `overwrite` | Logical. Overwrite output file if it exists. |
| `verbose` | Logical. Print progress messages. |

### `catfish_load_pathways()`

*Load pathways as a data.frame*

```r
catfish_load_pathways(
  species = c("maize", "sorghum", "arabidopsis", "plant", "fly", "human",
    "human_reactome"),
  gene_col = NULL,
  drop_unknown = TRUE,
  file = NULL
)
```

| Parameter | Description |
|---|---|
| `species` | Character. One of "maize", "sorghum", "arabidopsis", "plant", "fly", "human" (HumanCyc) or "human_reactome" (Reactome, Entrez IDs, from MSigDB). |
| `gene_col` | Character or NULL. Column name to use for gene IDs. If NULL (default), the function will: • For PMN files: use "Gene-name" if present, else "Gene-id" • For fly: use "gene" |
| `drop_unknown` | Logical. If TRUE (default), drop rows where gene_id is NA, empty, or "unknown". |
| `file` | Optional path to your own pathway file in PMN / Pathway Tools tab-separated format (columns `Pathway-id`, `Pathway-name` and the gene column), e.g. a newer PMN release or another species from https://plantcyc.org. When given, `species` is ignored. |

### `catfish_pmn_file()`

*Path to built-in pathway file*

```r
catfish_pmn_file(
  species = c("maize", "sorghum", "arabidopsis", "plant", "fly", "human",
    "human_reactome")
)
```

| Parameter | Description |
|---|---|
| `species` | Which built-in database to use: • "maize": CornCyc pathways (corncyc_pathways.20230103) • "sorghum": SorghumBicolorCyc pathways with SORBI_ gene IDs, the IDs of the built-in sorghum.genes.loc (sorghumbicolorcyc_pathways.20230103.SORBI) • "arabidopsis": AraCyc pathways (aracyc_pathways.20230103) • "plant": PlantCyc general plant pathways (plantcyc_pathways.20230103) • "fly": FlyCyc Drosophila pathways (Fly_Cyc.tsv) • "human": HumanCyc pathways (humancyc_pathways.20261024) • "human_reactome": Reactome pathways for human, Entrez gene IDs (human_reactome.pw; 1,867 pathways from MSigDB v2026.2 C2:CP:REACTOME, c2.cp.reactome.v2026.2.Hs.entrez). MSigDB gene sets are distributed under the Creative Commons Attribution 4.0 International License; cite MSigDB and Reactome when you use them. One pathway per line: ID, then member genes. |

### `catfish_load_gmt()`

*Load pathways from a GMT file*

```r
catfish_load_gmt(file, min_size = 1L, max_size = Inf, drop_unknown = TRUE)
```

| Parameter | Description |
|---|---|
| `file` | Path to a .gmt file. |
| `min_size`, `max_size` | Integer. Keep only pathways whose gene count falls within these bounds after de-duplication. Defaults keep everything; 10 and 200 are common choices for GWAS pathway analysis. |
| `drop_unknown` | Logical. Drop blank and "unknown" gene identifiers. |

### `catfish_load_biocyc()`

*Load pathways from a BioCyc / Pathway Tools SmartTable export*

```r
catfish_load_biocyc(
  file,
  id_col = "Object ID",
  name_col = "Pathways",
  gene_col = "Genes of Pathways",
  min_size = 1L,
  max_size = Inf,
  drop_unknown = TRUE
)
```

| Parameter | Description |
|---|---|
| `file` | Path to the exported .csv. |
| `id_col`, `name_col` | Column names holding the pathway object ID and the pathway name. |
| `gene_col` | Column holding the gene identifier to use. For human analyses with MAGMA this should be an NCBI/Entrez gene ID column, since MAGMA's gene location files are Entrez-keyed. |
| `min_size`, `max_size` | Pathway size bounds, applied after de-duplication. |
| `drop_unknown` | Logical. Drop blank and "unknown" gene identifiers. |

### `catfish_write_set_annot()`

*Write a CATFISH pathway data.frame as a MAGMA --set-annot file*

```r
catfish_write_set_annot(pathways, file)
```

| Parameter | Description |
|---|---|
| `pathways` | data.frame from `catfish_load_gmt` or `catfish_load_biocyc`. |
| `file` | Output path. One row per pathway: set name followed by its gene identifiers, tab separated. |

### `catfish_acat_pathways()`

*ACAT pathway test on gene-level p-values*

```r
catfish_acat_pathways(
  gene_results,
  pathways = NULL,
  species = NULL,
  pmn_gene_col = NULL,
  gene_col = "GENE",
  p_col = "P",
  min_p = 1e-15,
  do_fix = TRUE,
  output = FALSE,
  out_dir = "catfish_acat"
)
```

| Parameter | Description |
|---|---|
| `gene_results` | data.frame with at least a gene ID column and a p-value column. For MAGMA `.genes.out`, these are commonly `GENE` and `P`. |
| `pathways` | Either: • A named list: each element is a character vector of gene IDs in that pathway, OR • A data.frame with columns `pathway_id`, `gene_id` (and optional `pathway_name`). |
| `species` | Character. One of "maize", "sorghum", "arabidopsis", "plant", or "fly". If provided, built-in pathways are loaded via `catfish_load_pathways`. Provide either `pathways` OR `species`, but not both. |
| `pmn_gene_col` | Character or NULL. Passed to `catfish_load_pathways(gene_col=...)` when `species` is used. If NULL, the loader uses default gene column for each species. |
| `gene_col` | Character. Column name in `gene_results` containing gene IDs. Default `"GENE"`. |
| `p_col` | Character. Column name in `gene_results` containing gene-level p-values. Default `"P"`. |
| `min_p` | Numeric. Lower cap for very small p-values before ACAT. Default 1e-15. |
| `do_fix` | Logical. If TRUE (default), clean/cap p-values with `fix_p_for_acat`. |
| `output` | Logical. If TRUE, write a CSV of results to `out_dir`. |
| `out_dir` | Character. Directory to write CSV when `output = TRUE`. Default `"catfish_acat"`. |

### `catfish_fisher_pathways()`

*Pathway-level Fisher test on gene-level p-values*

```r
catfish_fisher_pathways(
  gene_results,
  pathways = NULL,
  species = NULL,
  pmn_gene_col = NULL,
  gene_col = "GENE",
  p_col = "P",
  min_p = 1e-15,
  do_fix = TRUE,
  output = FALSE,
  out_dir = "catfish_fisher"
)
```

| Parameter | Description |
|---|---|
| `gene_results` | data.frame with at least a gene ID column and a p-value column. For MAGMA `.genes.out`, these are commonly `GENE` and `P`. |
| `pathways` | Either: • A named list: each element is a character vector of gene IDs in that pathway, OR • A data.frame with columns `pathway_id`, `gene_id` (and optional `pathway_name`). |
| `species` | Optional; one of "maize", "sorghum", "arabidopsis", "plant", "fly". If provided, built-in pathways are loaded via `catfish_load_pathways()`. Provide either `pathways` OR `species`, but not both. |
| `pmn_gene_col` | Optional; passed to `catfish_load_pathways(gene_col=...)` when `species` is used. If NULL, the PMN loader prefers "Gene-name" then "Gene-id". |
| `gene_col` | Column name in `gene_results` containing gene IDs (default "GENE"). |
| `p_col` | Column name in `gene_results` containing gene-level p-values (default "P"). |
| `min_p` | Lower cap for very small p-values (default 1e-15). |
| `do_fix` | If TRUE (default), clean/cap p-values using `fix_p_for_acat()` to ensure values lie strictly in (0,1). |
| `output` | If TRUE, write a CSV of results to `out_dir`. |
| `out_dir` | Directory to write CSV when `output = TRUE` (default "catfish_fisher"). |

### `catfish_soft_tfisher_adaptive_pathways()`

*Pathway-level *adaptive* soft TFisher (omnibus) on gene-level p-values*

```r
catfish_soft_tfisher_adaptive_pathways(
  gene_results,
  pathways = NULL,
  species = NULL,
  pmn_gene_col = NULL,
  gene_col = "GENE",
  p_col = "P",
  tau_grid = "auto",
  min_p = 1e-15,
  do_fix = TRUE,
  output = FALSE,
  out_dir = "catfish_tfisher_soft_adaptive"
)
```

| Parameter | Description |
|---|---|
| `gene_results` | data.frame with at least a gene ID column and a p-value column. For MAGMA `.genes.out`, these are commonly `GENE` and `P`. |
| `pathways` | Either: • A named list: each element is a character vector of gene IDs in that pathway, OR • A data.frame with columns `pathway_id`, `gene_id` (and optional `pathway_name`). |
| `species` | Optional; one of "maize", "sorghum", "arabidopsis", "plant", "fly". If provided, built-in PMN pathways are loaded via `catfish_load_pathways()`. Provide either `pathways` OR `species`, but not both. |
| `pmn_gene_col` | Optional; passed to `catfish_load_pathways(gene_col=...)` when `species` is used. If NULL, the PMN loader prefers "Gene-name" then "Gene-id". |
| `gene_col` | Column name in `gene_results` containing gene IDs (default "GENE"). |
| `p_col` | Column name in `gene_results` containing gene-level p-values (default "P"). |
| `tau_grid` | Adaptive soft-TFisher thresholds: either a numeric vector of `tau` values (internally sorted increasing, as required by `TFisher::stat.soft.omni()`), or `"auto"` (default) to choose the grid from the empirical quantiles of the pooled gene-level p-values. Values must satisfy 0 < tau < 1; the tau = 1 Fisher corner is excluded here because CATFISH already carries Fisher as a separate component test. |
| `min_p` | Lower cap for extremely small p-values before TFisher math (default 1e-15). |
| `do_fix` | If TRUE (default), clean/cap p-values using `fix_p_for_acat()` to ensure values lie strictly in (0,1). |
| `output` | If TRUE, write a CSV of results to `out_dir`. |
| `out_dir` | Directory to write CSV when `output = TRUE` (default "catfish_tfisher_soft_adaptive"). |

### `catfish_minp_pathways()`

*Pathway-level minP (Tippett/Wilkinson) on gene-level p-values*

```r
catfish_minp_pathways(
  gene_results,
  pathways = NULL,
  species = NULL,
  pmn_gene_col = NULL,
  gene_col = "GENE",
  p_col = "P",
  min_p = 1e-15,
  do_fix = TRUE,
  output = FALSE,
  out_dir = "catfish_minp"
)
```

| Parameter | Description |
|---|---|
| `gene_results` | data.frame with at least a gene ID column and a p-value column. For MAGMA `.genes.out`, these are commonly `GENE` and `P`. |
| `pathways` | Either: • A named list: each element is a character vector of gene IDs in that pathway, OR • A data.frame with columns `pathway_id`, `gene_id` (and optional `pathway_name`). |
| `species` | Optional; one of "maize", "sorghum", "arabidopsis", "plant", "fly". If provided, built-in PMN pathways are loaded via `catfish_load_pathways()`. Provide either `pathways` OR `species`, but not both. |
| `pmn_gene_col` | Optional; passed to `catfish_load_pathways(gene_col=...)` when `species` is used. If NULL, the PMN loader prefers "Gene-name" then "Gene-id". |
| `gene_col` | Column name in `gene_results` containing gene IDs (default "GENE"). |
| `p_col` | Column name in `gene_results` containing gene-level p-values (default "P"). |
| `min_p` | Lower cap for extremely small p-values before transformations (default 1e-15). |
| `do_fix` | If TRUE (default), clean/cap p-values using `fix_p_for_acat()` to ensure values lie strictly in (0,1). |
| `output` | If TRUE, write a CSV of results to `out_dir`. |
| `out_dir` | Directory to write CSV when `output = TRUE` (default "catfish_minp"). |

### `catfish_stoufferZ_pathways()`

*Pathway-level (weighted) Stouffer Z test on gene-level Z statistics*

```r
catfish_stoufferZ_pathways(
  gene_results,
  pathways = NULL,
  species = NULL,
  pmn_gene_col = NULL,
  gene_col = "GENE",
  z_col = "ZSTAT",
  weight_col = NULL,
  min_abs_w = 1e-08,
  min_p = 1e-15,
  alternative = c("greater", "two.sided", "less"),
  output = FALSE,
  out_dir = "catfish_stouffer_z"
)
```

| Parameter | Description |
|---|---|
| `gene_results` | data.frame with at least a gene ID column and a Z-statistic column. For MAGMA `.genes.out`, you may have `GENE` and `ZSTAT` (or similar). |
| `pathways` | Either: • A named list: each element is a character vector of gene IDs in that pathway, OR • A data.frame with columns `pathway_id`, `gene_id` (and optional `pathway_name`). |
| `species` | Optional; one of "maize", "sorghum", "arabidopsis", "plant", "fly". If provided, built-in PMN pathways are loaded via `catfish_load_pathways()`. Provide either `pathways` OR `species`, but not both. |
| `pmn_gene_col` | Optional; passed to `catfish_load_pathways(gene_col=...)` when `species` is used. If NULL, the PMN loader prefers "Gene-name" then "Gene-id". |
| `gene_col` | Column name in `gene_results` containing gene IDs (default "GENE"). |
| `z_col` | Column name in `gene_results` containing gene-level Z statistics (default "ZSTAT"). |
| `weight_col` | Optional column name in `gene_results` for per-gene weights (default NULL = equal weights). Examples: `sqrt(NSNPS)`, `1/SE`, etc. |
| `min_abs_w` | Non-finite / NA / <=0 weights are replaced by this small positive value (default 1e-8). |
| `min_p` | Lower cap for extremely small p-values to avoid underflow to 0 (default 1e-15). |
| `alternative` | One of `"greater"`, `"two.sided"`, `"less"`. `"greater"` tests for positive enrichment (default "greater"). |
| `output` | If TRUE, write a CSV of results to `out_dir`. |
| `out_dir` | Directory to write CSV when `output = TRUE` (default "catfish_stouffer_z"). |

### `catfish_wy_pathways()`

*Pathway testing with Westfall--Young family-wise calibration*

```r
catfish_wy_pathways(
  gene_results,
  pathways = NULL,
  species = NULL,
  pmn_gene_col = NULL,
  gene_col = "GENE",
  p_col = "P",
  chr_col = "CHR",
  z_col = NULL,
  cor_pairs = NULL,
  cor_file = NULL,
  B = 1000L,
  use_tfisher = TRUE,
  tau_grid = c(0.01, 0.05, 0.5),
  min_genes = 10L,
  max_genes = 200L,
  seed = 20260706L,
  gc = TRUE,
  verbose = TRUE,
  output = FALSE,
  out_dir = "catfish_wy"
)
```

| Parameter | Description |
|---|---|
| `gene_results` | data.frame of gene-level results (e.g. MAGMA `.genes.out`), with columns for gene ID, p-value and chromosome. |
| `pathways` | Named list of gene-ID vectors, or a data.frame with `pathway_id` and `gene_id`. Provide this or `species`. |
| `species` | Optional; built-in PMN pathways via `catfish_load_pathways()`. |
| `pmn_gene_col` | Passed to `catfish_load_pathways()` when `species` is used. |
| `gene_col`, `p_col`, `chr_col` | Columns of `gene_results` holding the gene ID, gene-level p-value and chromosome. |
| `z_col` | Optional column (e.g. MAGMA `ZSTAT`); if given, genes with a non-finite value are dropped, as MAGMA output is usually filtered. |
| `cor_pairs` | data.frame of gene pairs and correlations (gene1, gene2, r), e.g. from `magma_genesraw_to_cor_pairs_banded()`. |
| `cor_file` | Path to the same table on disk; used if `cor_pairs` is `NULL`. |
| `B` | Number of genome-wide null replicates. |
| `use_tfisher` | Include adaptive soft TFisher as a fifth component. |
| `tau_grid` | Fixed soft TFisher truncation grid. The default is the grid of Zhang et al. (2020) without its tau=1 corner, which is Fisher. |
| `min_genes`, `max_genes` | Pathway size limits, counted over genes present in `gene_results`. |
| `seed` | Random seed for the null draws (`NULL` leaves it unset). |
| `gc` | Gene-level genomic control: centre the gene z-scores at their median before testing (default `TRUE`). See Details. |
| `verbose` | Print progress. |
| `output` | If `TRUE`, write a CSV to `out_dir`. |
| `out_dir` | Output directory when `output = TRUE`. |

### `catfish_omni2_pathways()`

*CATFISH omnibus pathway p-values (analytic, global resampling, or Westfall--Young)*

```r
catfish_omni2_pathways(
  gene_results,
  pathways = NULL,
  species = NULL,
  pmn_gene_col = NULL,
  gene_col = "GENE",
  p_raw_col = "P",
  p_adj_col = NULL,
  z_col = NULL,
  weight_col = NULL,
  chr_col = NULL,
  tau_grid = "auto",
  tau_cap = 1,
  min_p = 1e-50,
  min_genes = 2L,
  max_genes = Inf,
  use_tfisher = TRUE,
  do_fix = TRUE,
  stouffer_min_abs_w = 1e-08,
  stouffer_alternative = c("greater", "two.sided", "less"),
  magma_out = NULL,
  include_magma_in_omni = TRUE,
  include_magma_in_perm = FALSE,
  omnibus = c("ACAT", "minP"),
  B_perm = NULL,
  perm_mode = c("none", "global", "wy", "both"),
  wy_tau_grid = c(0.01, 0.05, 0.5),
  tail_mode = c("empirical", "hybrid_gpd"),
  tail_switch_exceed = 10L,
  tail_gpd_k = 250L,
  tail_min_B = 10000L,
  tail_min_tail = 50L,
  B_global = 0L,
  B_wy = 0L,
  perm_pool = c("raw", "obs"),
  magma_cor_file = NULL,
  magma_cor_pairs = NULL,
  seed = NULL,
  gc = TRUE,
  output = FALSE,
  out_dir = "catfish_omni2"
)
```

| Parameter | Description |
|---|---|
| `gene_results` | data.frame of gene-level results (one row per gene recommended). Must contain `gene_col` and at least one p-value column (`p_raw_col` and/or `p_adj_col`). If you want Stouffer, it must contain `z_col`. If you want weighted Stouffer, it must contain `weight_col`. For Westfall--Young it must also contain `chr_col`. |
| `pathways` | Pathway definitions (provide `pathways` OR `species`). Either a data.frame with columns `pathway_id`, `gene_id` (optional `pathway_name`), or a named list mapping `pathway_id -> character vector of genes`. |
| `species` | Optional character (provide `species` OR `pathways`), one of `"maize"`, `"sorghum"`, `"arabidopsis"`, `"plant"`, `"fly"`, `"human"` (HumanCyc) or `"human_reactome"` (Reactome, Entrez IDs). When provided, pathways are loaded via `catfish_load_pathways()`. |
| `pmn_gene_col` | Optional character. When `species` is used, passed to `catfish_load_pathways(species=..., gene_col=pmn_gene_col)`. |
| `gene_col` | Character. Gene ID column in `gene_results`. Default `"GENE"`. |
| `p_raw_col` | Character. Raw p-value column in `gene_results`. Default `"P"`. |
| `p_adj_col` | Character or NULL. Preferred/adjusted p-value column. If supplied and present, used for ACAT/Fisher/TFisher/minP. Default NULL. |
| `z_col` | Character or NULL. Gene Z column for Stouffer. If NULL, Stouffer is skipped. Default NULL. |
| `weight_col` | Character or NULL. Optional weights for Stouffer. Default NULL (equal weights). |
| `chr_col` | Character or NULL. Chromosome column in `gene_results`. Required for `perm_mode = "wy"` / `"both"`; genes are drawn jointly within a chromosome and independently across chromosomes. Default NULL. |
| `tau_grid` | Adaptive soft-TFisher thresholds for the analytic component. Either a numeric vector of tau values, or `"auto"` (default) to select the grid from the empirical quantiles of the pooled gene-level p-values. |
| `tau_cap` | Numeric. Upper bound for `tau_grid` (sanity bound). Default 1. |
| `min_p` | Numeric. P-values are clipped to `[min_p, 1-min_p]` to avoid numeric issues. Default 1e-50. |
| `min_genes` | Integer. Minimum genes required per pathway. Default 2. |
| `max_genes` | Numeric. Maximum genes per pathway (counted over genes present). Default `Inf` (no upper limit). |
| `use_tfisher` | Logical. Include the adaptive soft-TFisher component. Default TRUE. |
| `do_fix` | Logical. If TRUE, uses your `fix_p_for_acat()` if available (else clips). Default TRUE. |
| `stouffer_min_abs_w` | Numeric. Non-finite / <=0 Stouffer weights replaced by this. Default 1e-8. |
| `stouffer_alternative` | One of `"greater"`, `"two.sided"`, `"less"` for Stouffer. Default "greater". |
| `magma_out` | Optional data.frame of MAGMA competitive results with columns `pathway_id`, `magma_pvalue`. |
| `include_magma_in_omni` | Logical. Include MAGMA p-values as a component in the omnibus (if provided). Default TRUE. |
| `include_magma_in_perm` | Logical. Include MAGMA p-values inside the resampling null. Generally not recommended. If FALSE and permutations are run, MAGMA is omitted from the permuted omnibus. Default FALSE. |
| `omnibus` | Character. `"ACAT"` or `"minP"` for combining method p-values. Default `"ACAT"`. The Westfall--Young omnibus always uses ACAT. |
| `B_perm` | Integer or NULL. Convenience parameter for number of resampling iterations. If provided, it overrides `B_global`/`B_wy` according to `perm_mode`. |
| `perm_mode` | Character. One of `"none"` (analytic only), `"global"` (global gene resampling), `"wy"` (Westfall--Young) or `"both"` (global and Westfall--Young; `omni_p_final` prefers the Westfall--Young calibrated p-value). |
| `wy_tau_grid` | Numeric. Fixed soft-TFisher grid used inside the Westfall--Young omnibus. Default `c(0.01, 0.05, 0.5)`: the grid of Zhang et al. (2020) without its tau = 1 corner, which is Fisher. |
| `tail_mode` | Character. `"empirical"` or `"hybrid_gpd"` (GPD tail extrapolation for global-resampling p-values below the Monte Carlo floor). Default `"empirical"`. |
| `tail_switch_exceed` | Integer. Switch to GPD extrapolation when the number of null exceedances falls below this. Default 10. |
| `tail_gpd_k` | Integer. Number of extreme null draws used to fit the GPD. Default 250. |
| `tail_min_B` | Integer. Minimum number of null draws required before GPD is attempted. Default 10000. |
| `tail_min_tail` | Integer. Minimum tail observations required for a stable GPD fit. Default 50. |
| `B_global` | Integer. Number of global-resampling iterations (used when `B_perm` is NULL). |
| `B_wy` | Integer. Number of genome-wide Westfall--Young null replicates (used when `B_perm` is NULL). |
| `perm_pool` | Character. For global resampling only: `"raw"` (use raw p-values) or `"obs"` (use p used by methods). |
| `magma_cor_file` | Character or NULL. Path to 3-col correlation pairs file: `gene1 gene2 r`. Used by Westfall--Young; without it (or `magma_cor_pairs`) genes are treated as independent. |
| `magma_cor_pairs` | data.frame or NULL. In-memory correlation pairs with columns `gene1`, `gene2`, `r`. If provided, it overrides `magma_cor_file`. |
| `seed` | Integer or NULL. Random seed. |
| `gc` | Gene-level genomic control (default `TRUE`): centre the gene p-values (and `z_col`) at their median before testing. See `catfish_wy_pathways` for when to set it to `FALSE`. |
| `output` | Logical. If TRUE, writes CSV to `out_dir` and attaches `"file"` attribute. Default FALSE. |
| `out_dir` | Character. Output directory. Default `"catfish_omni2"`. |

### `gene_prioritize()`

*Prioritise candidate genes from GWAS, MAGMA and pathway evidence*

```r
gene_prioritize(
  magma,
  pathway_results,
  pathways = NULL,
  species = NULL,
  pmn_gene_col = NULL,
  gwas = NULL,
  gene_annot = NULL,
  gene_loc = NULL,
  magma_gene_col = "GENE",
  magma_p_col = NULL,
  pathway_id_col = "pathway_id",
  pathway_stat_col = NULL,
  gwas_snp_col = "SNP",
  gwas_p_col = NULL,
  gwas_neglog10p_col = NULL,
  gwas_chr_col = "CHR",
  gwas_pos_col = "BP",
  symbol_col = NULL,
  window = 10000L,
  w_gwas = 1,
  w_magma = 1,
  w_pathway = 0.5,
  w_hub = 0.3,
  targets = NULL,
  chunk = 2000L,
  verbose = TRUE,
  output = FALSE,
  out_file = "gene_prioritize.csv"
)
```

| Parameter | Description |
|---|---|
| `magma` | MAGMA gene results: data.frame or file path. |
| `pathway_results` | Pathway results: data.frame or file path (CSV/TSV). |
| `pathways` | Pathway membership: list, data.frame or set-file path. Give this or `species`. |
| `species` | Built-in pathway set (see `catfish_load_pathways`); give this or `pathways`. |
| `pmn_gene_col` | Gene-ID column of the built-in set, passed to `catfish_load_pathways(gene_col = ...)` when `species` is used. |
| `gwas` | Optional SNP-level results: data.frame or file path. |
| `gene_annot` | Optional SNP-to-gene annotation: MAGMA `.genes.annot` path or data.frame with `gene` and `SNP`. |
| `gene_loc` | Optional gene locations: `.gene.loc` path or data.frame. |
| `magma_gene_col`, `magma_p_col` | Gene ID and p-value columns of `magma`. |
| `pathway_id_col`, `pathway_stat_col` | Pathway ID and statistic columns of `pathway_results`. |
| `gwas_snp_col`, `gwas_p_col`, `gwas_neglog10p_col`, `gwas_chr_col`, `gwas_pos_col` | Columns of `gwas`. |
| `symbol_col` | Column of `gene_loc` holding gene symbols: a name, or a number for a header-less `.gene.loc` file (MAGMA's NCBI files keep the symbol in column 6). `NULL` uses a `symbol` column if present, otherwise column 6 of a 6-column file. |
| `window` | Base pairs added on each side of a gene for the position route. Default 10000. |
| `w_gwas`, `w_magma`, `w_pathway`, `w_hub` | Layer weights. Defaults 1, 1, 0.5, 0.3. |
| `targets` | Optional character vector of known genes (IDs or symbols) to report as positive controls. |
| `chunk` | Genes read at a time from a `.genes.annot` file. |
| `verbose` | Print progress. |
| `output` | If `TRUE`, write the table to `out_file`. |
| `out_file` | CSV path used when `output = TRUE`. |

### `get_gene_lengths()`

*Extract gene lengths from a GFF3 file*

```r
get_gene_lengths(
  gff3_file,
  feature = "gene",
  id_key = "ID",
  keep_chr_prefix = NULL,
  output = FALSE,
  output_dir = ".",
  file_name = "gene_lengths.tsv"
)
```

| Parameter | Description |
|---|---|
| `gff3_file` | Path to a GFF3 file (.gff3 or .gff3.gz). |
| `feature` | Which feature type to treat as "gene" (default "gene"). |
| `id_key` | Attribute key used for gene ID (default "ID"). |
| `keep_chr_prefix` | If non-NULL, enforce a prefix in CHR (e.g. "chr"). If NULL, leaves seqid as-is. |
| `output` | Logical; if TRUE write TSV to output_dir/file_name. |
| `output_dir` | Directory for output if output=TRUE. |
| `file_name` | Output TSV file name if output=TRUE. |

### `catfish_adjust_gene_p()`

*Adjust MAGMA gene p-values using gene length and NSNPS*

```r
catfish_adjust_gene_p(
  gene_results,
  gene_lengths,
  gene_col = "GENE",
  nsnp_col = "NSNPS",
  p_col = "P",
  z_col = "ZSTAT",
  len_gene_col = "gene_id",
  len_col = "length",
  log1p_covars = TRUE
)
```

| Parameter | Description |
|---|---|
| `gene_results` | data.frame. Gene-level results from MAGMA (typically from reading a `.genes.out` file). Must contain gene ID, p-value, and NSNPS columns. |
| `gene_lengths` | data.frame. Gene lengths, typically from `get_gene_lengths`. Must contain gene ID and length columns. |
| `gene_col` | character. Column name in `gene_results` for gene IDs. Default "GENE". |
| `nsnp_col` | character. Column name in `gene_results` for the number of SNPs per gene. Default "NSNPS". |
| `p_col` | character. Column name in `gene_results` for raw p-values. Default "P". |
| `z_col` | character or NULL. Column name for Z-scores in `gene_results`. If NULL or not present, Z-scores are reconstructed from p-values (magnitude only). Default "ZSTAT". |
| `len_gene_col` | character. Column name in `gene_lengths` for gene IDs. Default "gene_id". |
| `len_col` | character. Column name in `gene_lengths` for gene length. Default "length". |
| `log1p_covars` | logical. If TRUE (default), apply log(1+x) transformation to NSNPS and gene length covariates before regression. |

### `fix_p_for_acat()`

*Clean p-values for ACAT*

```r
fix_p_for_acat(p, min_p = 1e-15)
```

| Parameter | Description |
|---|---|
| `p` | Numeric vector of p-values (NA allowed). Names are preserved. |
| `min_p` | Numeric. Lower cap for very small p-values (default 1e-15). |

