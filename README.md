# CATFISH

CATFISH is an R package for GWAS pathway analysis built on MAGMA gene-level outputs.  
It combines multiple pathway tests (ACAT, Fisher, adaptive soft TFisher, minP, Stouffer) into one final pathway-level result, with optional LD-aware MVN calibration using MAGMA gene-gene correlations.

## What CATFISH does

- Runs complementary pathway tests from gene-level GWAS statistics.
- Combines component tests into an omnibus p-value (`ACAT` or `minP`).
- Supports LD-aware calibration with MVN resampling using MAGMA correlation pairs.
- Includes built-in pathway loaders for `maize`, `sorghum`, `arabidopsis`, `plant`, and `fly`.

## Installation

### 1) Install MAGMA (external dependency)

Download MAGMA from the official site: https://ctg.cncr.nl/software/magma

### 2) Install CATFISH

```r
# install.packages("remotes")
remotes::install_github("nirwantandukar/CATFISH")
```

### 3) Set MAGMA path

```r
library(CATFISH)
options(magma.path = "/absolute/path/to/magma")
magma_path()
```

## End-to-end usage

### Step 1: Build gene location file from GFF3

```r
gene_loc <- gff3_to_geneloc(
  gff = "data/reference.gff3",
  out = "results/species.genes.loc"
)
```

### Step 2: MAGMA annotation and gene analysis

```r
ann <- magma_annotate(
  stats_file = "data/gwas.tsv",
  rename_columns = c(CHR = "CHR", SNP = "SNP", POS = "POS", PVALUE = "P"),
  gene_loc = "results/species.genes.loc",
  out_prefix = "trait",
  out_dir = "results/annot",
  window = c(10, 10),
  nonhuman = TRUE
)

magma_gene(
  bfile = "data/plink_prefix",
  gene_annot = ann$gene_annot,
  stats_file = "data/gwas.tsv",
  rename_columns = c(CHR = "CHR", SNP = "SNP", POS = "POS", PVALUE = "P"),
  n_total = 1000,
  out_prefix = "trait",
  out_dir = "results/magma_genes",
  gene_model = "multi=snp-wise",
  chroms = 1:10,
  n_threads = 10
)
```

Then combine per-chromosome `.genes.out` files into one `gene_results` data frame with at least:

- `GENE`
- `P`
- optional `ZSTAT` (needed for Stouffer)
- optional `NSNPS` (for p-adjustment step)

### Step 3 (optional): Adjust gene p-values for length/NSNPS

```r
gene_lengths <- get_gene_lengths("data/reference.gff3")
adj <- catfish_adjust_gene_p(
  gene_results = gene_results,
  gene_lengths = gene_lengths,
  gene_col = "GENE",
  nsnp_col = "NSNPS",
  p_col = "P",
  z_col = "ZSTAT"
)

gene_results$P_adj <- adj$p_adj[match(gene_results$GENE, adj$gene_id)]
```

### Step 4: Build gene-gene correlation pairs from MAGMA `.genes.raw`

```r
magma_genesraw_to_cor_pairs_banded(
  genes_raw_file = "results/magma_genes/trait_chr1.genes.raw",
  out_pairs_file = "results/magma_genes/gene_cor_pairs_chr1.txt",
  gene_regex = "^SORBI",   # set for your species
  keep_abs_r_ge = 0
)
```

For multi-chromosome analyses, run this per chromosome and concatenate (header only once), as shown in `usage2.R`.

### Step 5: Load pathways

Known species pathways:

```r
pathways <- catfish_load_pathways(species = "sorghum")
```

Custom pathways as long-format data frame:

```r
pathways <- data.frame(
  pathway_id = c("PWY1", "PWY1", "PWY2"),
  pathway_name = c("Pathway 1", "Pathway 1", "Pathway 2"),
  gene_id = c("GENE_A", "GENE_B", "GENE_C"),
  stringsAsFactors = FALSE
)
```

### Step 6: Run CATFISH omnibus pathway analysis

```r
omni <- catfish_omni2_pathways(
  gene_results = gene_results,
  pathways = pathways,                  # or species = "sorghum"
  gene_col = "GENE",
  p_raw_col = "P",
  p_adj_col = "P_adj",                  # optional
  z_col = "ZSTAT",                      # optional but recommended
  omnibus = "ACAT",                     # or "minP"
  B_perm = 10000,
  perm_mode = "mvn",                    # none/global/mvn/mvn_global/both
  magma_cor_file = "results/magma_genes/gene_cor_pairs_all.txt",
  mvn_marginal = "uniform",             # uniform or empirical
  mvn_calibrate_components = FALSE,
  min_genes = 2,
  n_threads = 8,
  seed = 123
)
```

For publication-grade analyses with GPD tail extrapolation:

```r
omni <- catfish_omni2_pathways(
  gene_results = gene_results,
  pathways = pathways,
  gene_col = "GENE",
  p_raw_col = "P",
  z_col = "ZSTAT",
  omnibus = "ACAT",
  B_perm = 1000000,                     # 1 million permutations
  perm_mode = "mvn",
  magma_cor_file = "results/magma_genes/gene_cor_pairs_all.txt",
  tail_mode = "hybrid_gpd",             # enable GPD extrapolation
  tail_switch_exceed = 10L,
  tail_gpd_k = 250L,
  min_p = 1e-50,
  n_threads = 12,                       # parallel processing
  seed = 123
)
```

## Core parameter guide

### `catfish_omni2_pathways()`

| Parameter | Meaning | Typical values |
|---|---|---|
| `gene_results` | Gene-level input table | data.frame with `GENE`, `P`, optional `ZSTAT` |
| `pathways` / `species` | Pathway source | custom df/list or built-in species |
| `p_raw_col`, `p_adj_col` | P-value columns for p-based methods | `P`, optional `P_adj` |
| `z_col` | Z-score column for Stouffer | `ZSTAT` |
| `omnibus` | Across-method combiner | `"ACAT"`, `"minP"` |
| `perm_mode` | Calibration mode | `"none"`, `"global"`, `"mvn"`, `"mvn_global"`, `"both"` |
| `B_perm` | Number of resampling draws | 5,000 to 1,000,000 |
| `magma_cor_file` / `magma_cor_pairs` | Gene-gene correlations for MVN | required for MVN modes |
| `mvn_marginal` | Marginal mapping in MVN mode | `"uniform"` (recommended), `"empirical"` |
| `mvn_calibrate_components` | Calibrate each component before omnibus | `FALSE` (faster), `TRUE` (stricter) |
| `tau_grid` | Adaptive TFisher threshold grid | numeric vector |
| `min_genes` | Minimum genes per pathway | usually `2` or more |
| `n_threads` | Parallel workers for MVN internals | integer (default: n_cores - 1) |
| `tail_mode` | P-value estimation mode | `"empirical"` or `"hybrid_gpd"` |
| `tail_switch_exceed` | Switch to GPD when exceedances below this | `10` (default) |
| `tail_gpd_k` | Number of extreme values for GPD fit | `250` (default) |
| `tail_min_B` | Minimum B required for GPD | `10000` (default) |
| `tail_min_tail` | Minimum tail observations for stable GPD | `50` (default) |
| `min_p` | Floor for reported p-values | `1e-50` (for GPD mode) |

### `magma_annotate()`

| Parameter | Meaning |
|---|---|
| `stats_file` | GWAS summary stats file |
| `rename_columns` | Map required columns (at least `CHR`, `SNP`, `POS`) |
| `gene_loc` | MAGMA gene location file (from `gff3_to_geneloc`) |
| `window` | Up/downstream kb window |
| `out_prefix`, `out_dir` | Output naming/location |

### `magma_gene()`

| Parameter | Meaning |
|---|---|
| `bfile` | PLINK prefix (`.bed/.bim/.fam`) |
| `gene_annot` | Output of `magma_annotate()` |
| `stats_file` + `rename_columns` | GWAS file and column mapping |
| `n_total` | Total sample size (when not provided per SNP in file) |
| `gene_model` | MAGMA gene model (e.g., `"multi=snp-wise"`) |
| `chroms`, `n_threads` | Per-chromosome parallel analysis |

### `catfish_load_pathways()`

| Parameter | Meaning |
|---|---|
| `species` | One of `maize`, `sorghum`, `arabidopsis`, `plant`, `fly` |
| `gene_col` | Optional gene column override for pathway source file |
| `drop_unknown` | Remove unknown/missing gene IDs |

## GPD tail extrapolation

Standard empirical p-values have a floor of approximately `1/(B+1)`. With `B = 10,000`, the minimum achievable p-value is ~10^-4, which limits detection of highly significant pathways.

CATFISH implements Generalized Pareto Distribution (GPD) tail extrapolation (Knijnenburg et al. 2009) to overcome this limitation. GPD fits an extreme value distribution to the tail of the null distribution, enabling p-value estimation far beyond the permutation floor.

### Enabling GPD mode

```r
omni <- catfish_omni2_pathways(
 gene_results = gene_results,
 pathways = pathways,
 B_perm = 1000000,
 perm_mode = "mvn",
 tail_mode = "hybrid_gpd",
 tail_switch_exceed = 10L,
 tail_gpd_k = 250L,
 tail_min_B = 10000L,
 tail_min_tail = 50L,
 min_p = 1e-50,
 n_threads = 12
)
```

### GPD parameters explained

- `tail_mode = "hybrid_gpd"`: Switches from empirical to GPD estimation when exceedances fall below `tail_switch_exceed`
- `tail_switch_exceed`: Threshold for switching to GPD (default: 10). When fewer than this many null values exceed the observed statistic, GPD extrapolation is used.
- `tail_gpd_k`: Number of extreme values used to fit the GPD (default: 250). Larger values give more stable fits.
- `tail_min_B`: Minimum permutations required before GPD is attempted (default: 10,000)
- `tail_min_tail`: Minimum tail observations needed for a stable GPD fit (default: 50)
- `min_p`: Floor for reported p-values (set to `1e-50` for GPD mode)

With `B = 1,000,000` and GPD, p-values as extreme as 10^-35 can be estimated reliably.

## Parallel processing

CATFISH supports multi-threaded parallelization for MVN resampling, significantly reducing runtime for large B values.

```r
omni <- catfish_omni2_pathways(
  ...,
  n_threads = 12   # Use 12 cores
)
```

- Uses `parallel::mclapply()` on Unix-based systems
- Parallelizes ACAT null, TFisher null, and omnibus null computations
- Only activates when `B >= 1000` to avoid overhead on small jobs
- Near-linear speedup with number of cores

**Example benchmark**: Sorghum analysis with 410 pathways and B = 1,000,000 completed in ~1 hour using 12 threads (vs. 50+ hours single-threaded).

## Main output columns

`catfish_omni2_pathways()` returns one row per pathway, including:

- component p-values: `acat_p`, `fisher_p`, `tfisher_p_analytic`, `minp_p_analytic`, `stouffer_p_analytic`
- optional `magma_pvalue` (if provided via `magma_out`)
- omnibus p-values: `omni_p_analytic`, `omni_p_global`, `omni_p_mvn`, `omni_p_final`
- diagnostics such as dominant component and calibration impact fields

## Reproducible examples in this repo

- Full workflow script: `usage2.R`
- Sorghum MVN + GPD run: `sorghum_catfish_results/run_catfish_GPD_1M.R`

## Notes

- CATFISH relies on MAGMA outputs; make sure MAGMA runs successfully before CATFISH steps.
- MVN calibration requires valid gene-gene correlation pairs (`gene1`, `gene2`, `r`).
- For publication-grade calibration, consider `perm_mode = "mvn"` and larger `B_perm`.
