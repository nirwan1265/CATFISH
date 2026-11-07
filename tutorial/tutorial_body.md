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

