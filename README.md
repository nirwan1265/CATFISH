# MAGCAT

**MAGCAT**: R interface to **MAGMA** with **ACAT / Fisher / TFisher**–based gene and pathway aggregation.

MAGCAT streamlines post-GWAS analysis by:

* running **MAGMA** (Multi-marker Analysis of Genomic Annotation) for SNP → gene,
* importing gene-level results into R,
* adjusting gene Z / p for **gene length** and **#SNPs per gene**,
* performing multiple **gene → pathway** combination tests:
  * ACAT (Aggregated Cauchy Association Test),
  * Fisher’s method (unweighted),
  * soft truncated Fisher (TFisher),
* combining them with an **ACAT-O–style omnibus ACAT**, and
* controlling FDR across pathways with BH and/or q-values.

It’s designed for reproducible, **tail-aware** pathway analysis on top of MAGMA.

---

## Installation

### 1. Install MAGMA (external dependency)

Download and install MAGMA from the official site:

* https://ctg.cncr.nl/software/magma

Make sure the `magma` executable is on your `PATH`, or note its full path.

### 2. Install MAGCAT in R

```r
# install.packages("devtools")  # if needed
devtools::install_github("YOUR_GITHUB_USERNAME/MAGCAT")
```

Replace `YOUR_GITHUB_USERNAME` with your GitHub username or org.

Then:

```r
library(MAGCAT)
```

### 3. Optional: set MAGMA path

If `magma` is not on your system `PATH`, configure it:

```r
MAGCAT::magma_set_path("/full/path/to/magma")
# see ?magma_path / ?magma_set_path for details
```

---

## Conceptual workflow

MAGCAT wraps MAGMA for SNP → gene analysis and then does flexible, R-based gene → pathway aggregation:

1. **SNP → gene (MAGMA)**  
   - Annotate SNPs to genes using GFF3-derived gene locations.  
   - Run MAGMA with `--gene-model multi=snp-wise` to get gene-level Z and p.

2. **Gene-level adjustment**  
   - Merge in **gene length** and **#SNPs per gene**.  
   - Regress raw MAGMA Z-scores on `log(gene_length)` and `log(NSNPS)`.  
   - Use residuals as **size/SNP-adjusted Z**, convert to adjusted p-values.

3. **Gene → pathway tests**  
   For each pathway, combine gene-level adjusted p-values (`P_adj`) using:
   - **ACAT** – sensitive to a few very small p’s,  
   - **Fisher’s method** – sensitive to many modest p’s,  
   - **Soft TFisher** – truncated/weighted Fisher, focusing on the tail.

4. **Omnibus combination (ACAT-O style)**  
   - Combine the three pathway p-values with ACAT again → **one omnibus p**.

5. **Multiple testing correction**  
   - Apply BH FDR and optional Storey q-values across pathways.

---

## Quick example

Minimal end-to-end example using bundled maize files in `inst/extdata/`:

```r
library(MAGCAT)

## 1. Load example MAGMA gene results -----------------------------

genes_file <- system.file(
  "magma_genes",
  "N_maize_MLM.snp_wise_top.genes.out",
  package = "MAGCAT"
)

genes_all <- read_magma_genes(genes_file)


## 2. Load gene lengths from GFF3-derived table -------------------

gene_len_file <- system.file(
  "extdata",
  "Zea_mays_gene_lengths.tsv",
  package = "MAGCAT"
)

maize_gene_len <- readr::read_tsv(gene_len_file, show_col_types = FALSE)


## 3. Adjust MAGMA gene p-values for gene length + NSNPS ----------

adj_out <- magcat_adjust_gene_p(
  gene_results = genes_all,      # MAGMA .genes.out
  gene_lengths = maize_gene_len, # from GFF3
  gene_col     = "GENE",
  nsnp_col     = "NSNPS",
  p_col        = "P",
  z_col        = "ZSTAT",
  len_gene_col = "gene_id",
  len_col      = "length"
)

genes_adj <- adj_out$genes   # Z_raw, Z_adj, P_adj, gene_length, log covariates, etc.


## 4. Load pathways (e.g. PMN maize pathways) ---------------------

pathways_maize <- magcat_load_pathways(species = "maize")
# or supply your own (data.frame with pathway_id, gene_id, pathway_name)


## 5. Pathway tests: ACAT, Fisher, soft TFisher -------------------

# (a) ACAT on adjusted gene p-values
res_acat <- magcat_acat_pathways(
  gene_results = genes_adj,
  pathways     = pathways_maize,
  gene_col     = "GENE",
  p_col        = "P_adj",
  B_perm       = 0L  # set >0 for permutation-calibrated ACAT
)

# (b) Fisher's method (unweighted) on adjusted gene p-values
res_fisher <- magcat_fisher_pathways(
  gene_results = genes_adj,
  pathways     = pathways_maize,
  gene_col     = "GENE",
  p_col        = "P_adj"
)

# (c) soft TFisher on adjusted gene p-values
res_tfsoft <- magcat_soft_tfisher_pathways(
  gene_results     = genes_adj,
  pathways         = pathways_maize,
  gene_col         = "GENE",
  p_col            = "P_adj",
  tau1             = 0.05,   # truncation / normalization parameter
  B_perm           = 1000L,  # gene-set permutations
  analytic_logical = TRUE
)


## 6. Omnibus ACAT across methods --------------------------------

res_omni <- magcat_omni_acat_pathways(
  res_acat    = res_acat,
  res_fisher  = res_fisher,
  res_tfisher = res_tfsoft,
  acat_col    = "acat_p",
  fisher_col  = "fisher_p",
  tfisher_col = "tfisher_p_analytic"
)

# Multiple testing (BH and qvalue)
res_omni$omni_p_BH <- p.adjust(res_omni$omni_p, method = "BH")

if (requireNamespace("qvalue", quietly = TRUE)) {
  res_omni$omni_p_q <- qvalue::qvalue(res_omni$omni_p)$qvalues
}

# Top pathways by omnibus p
head(res_omni[order(res_omni$omni_p), ], 10)
```

This gives, for each pathway:

* `acat_p` – ACAT p (gene → pathway),
* `fisher_p` – Fisher p,
* `tfisher_p_analytic` – soft TFisher p,
* `omni_p` – omnibus ACAT over the three,
* `omni_p_BH` – BH-adjusted FDR,
* `omni_p_q` – optional q-value.

---

## Detailed workflow and models

### 1. SNP → gene with MAGMA

MAGCAT doesn’t reimplement MAGMA; it wraps it and organizes the inputs/outputs.

#### 1.1 Gene locations from GFF3

Using `gff3_to_geneloc.R` and `gene_length_wrappers.R`, MAGCAT:

* Parses species GFF3 and creates a gene coordinate table:

  ```text
  gene_id   chr   start   end   length   strand
  ```

* Writes MAGMA-style gene location files `*.genes.loc`:

  ```text
  GENE   CHR   START   END
  ```

These live under `inst/extdata/` (e.g. `maize.genes.loc`, `sorghum.genes.loc`).

#### 1.2 MAGMA annotation and gene analysis

You supply:

* SNP location file `*.snp.loc`:

  ```text
  SNP_ID  CHR  BP
  ```

* GWAS p-value file `*.pval.txt`:

  ```text
  SNP_ID  P   [optional: BETA, SE, ...]
  ```

Then MAGCAT (via wrappers) runs, conceptually:

```bash
magma   --annotate   --snp-loc   <snp.loc>   --gene-loc  <genes.loc>   --out       <prefix>

magma   --bfile       <reference_panel>   --pval        <gwas.pval.txt> N=<N>   --gene-annot  <prefix>.genes.annot   --gene-model  multi=snp-wise   --out         <prefix>
```

Key points:

* `multi = snp-wise` combines:
  * snp-wise mean: sensitive to many small SNP effects,
  * snp-wise top: sensitive to a single strong SNP,
* while accounting for LD and #SNPs per gene.

Output: `<prefix>.genes.out` with columns like:

* `GENE`, `CHR`, `NSNPS`, `ZSTAT`, `P`, …

MAGCAT reads this via `read_magma_genes()` into a tidy `data.frame`.

---

### 2. Adjusting gene p-values for gene length and #SNPs

MAGMA already accounts for NSNPS/LD in its null, but gene-level Z and p can still show residual dependence on gene size and SNP density. MAGMA’s own gene-set test corrects for this using covariates; MAGCAT pushes that idea down one level by adjusting the gene Z-scores themselves.

#### 2.1 Gene length table

Using `gene_length_wrappers.R`, MAGCAT builds a gene length table (from GFF3):

```text
gene_id   chr   start   end   length   strand
```

For maize, this is stored e.g. as `Zea_mays_gene_lengths.tsv` in `inst/extdata/`.

#### 2.2 Regression model: `magcat_adjust_gene_p()`

Let:

* `Z_i` = MAGMA gene Z-stat for gene i,
* `L_i` = gene length (bp),
* `S_i` = number of SNPs in gene i (`NSNPS`).

MAGCAT fits the linear model:

```text
Z_i = β0 + β1 * log(L_i) + β2 * log(S_i) + ε_i
```

using all genes with valid data.

Then it defines adjusted Z and adjusted p as:

```text
Z_hat_i = β0 + β1 * log(L_i) + β2 * log(S_i)
Z_adj_i = Z_i - Z_hat_i
P_adj_i = 2 * Φ( -|Z_adj_i| )
```

where `Φ` is the CDF of N(0,1).

Interpretation:

* `P_adj_i` is the size- and SNP-adjusted gene p-value,
* removing linear dependence on log gene length and log #SNPs.

`magcat_adjust_gene_p()` returns:

* `genes`: original gene table augmented with:
  * `gene_length`, `log_gene_length`, `log_nsnp`,
  * `Z_raw`, `Z_adj`, `P_adj`,
* `fit`: the `lm` object (inspect coefficients, R², diagnostics, etc.).

---

### 3. Gene → pathway: combination tests

Given a pathway `g` containing genes `i ∈ g`, and adjusted gene p-values `p_i = P_adj_i`, MAGCAT provides three main pathway tests.

Let `k = |g|` be the number of genes in the pathway.

#### 3.1 ACAT: `magcat_acat_pathways()`

Aggregated Cauchy Association Test (ACAT) combines p-values using a Cauchy transform:

```text
T_ACAT = sum_{i=1}^k w_i * tan( π * (0.5 - p_i) )
```

where `w_i >= 0` and `sum w_i = 1`. In MAGCAT’s default, all weights are equal:

```text
w_i = 1 / k
```

Under the null (for independent or moderately dependent p-values), `T_ACAT` is approximately standard Cauchy:

```text
p_ACAT ≈ 0.5 - (1 / π) * atan(T_ACAT)
```

Properties:

* Very sensitive to a few very small p’s,
* Robust to dependence,
* Extremely fast (no permutations required).

MAGCAT also supports an optional permutation-calibrated ACAT:

* Permute gene labels across pathways (or within strata),
* Recompute ACAT for each permuted set,
* Empirical p:

```text
p_perm = (1 + # { b : p_ACAT^(b) <= p_ACAT_obs }) / (B + 1)
```

Output columns include `acat_p` and (if requested) `acat_p_perm`.

#### 3.2 Fisher’s method (unweighted): `magcat_fisher_pathways()`

Classic Fisher combination:

```text
T_F = -2 * sum_{i=1}^k log(p_i)
```

Under the null (independent p’s):

```text
T_F ~ chi-square with 2k degrees of freedom
p_Fisher = 1 - F_chisq_2k( T_F )
```

where `F_chisq_2k` is the chi-square CDF with `2k` degrees of freedom.

Properties:

* Sensitive to many modest p-values (polygenic signals within a pathway),
* Less dominated by a single tiny p than ACAT,
* Easy and deterministic.

`magcat_fisher_pathways()` computes `fisher_p` for each pathway using adjusted gene p-values.

#### 3.3 Soft TFisher (soft truncated Fisher): `magcat_soft_tfisher_pathways()`

TFisher generalizes Fisher’s method by truncating and/or weighting p-values, focusing on the lower tail (most significant genes). MAGCAT uses the soft-threshold version implemented in the `TFisher` R package.

Given a truncation/normalization parameter `tau1` (e.g. 0.05), soft TFisher emphasizes p-values below `tau1` and effectively downweights or nullifies very large p-values.

For pathway p-values `p_i`, MAGCAT computes:

1. The soft TFisher statistic:

   ```r
   stat <- TFisher::stat.soft(p = p_i, tau1 = tau1)
   ```

2. The null CDF (left-tail) via:

   ```r
   Fq <- TFisher::p.soft(q = stat, n = length(p_i), tau1 = tau1, M = NULL)
   ```

3. The right-tail (test) p-value:

   ```text
   p_TFisher = 1 - Fq
   ```

MAGCAT reports:

* `tfisher_stat` – the soft TFisher statistic,
* `tfisher_p_analytic` – analytic soft TFisher p-value (`1 - p.soft(...)`),
* optional `tfisher_p_perm` – permutation-calibrated p from gene-set permutations:

   ```text
   p_TFisher_perm = (1 + # { b : T_TFisher^(b) >= T_TFisher_obs }) / (B + 1)
   ```

Properties:

* Bridges Fisher and minP:
  * Good when there are a few strong genes among many nulls,
  * Retains power when there are multiple moderate signals.
* `tau1` tunes how tail-focused the test is.

---

### 4. Omnibus ACAT (ACAT-O style): `magcat_omni_acat_pathways()`

Different combination tests favor different signal architectures:

* ACAT: few very small p’s,
* Fisher: many modestly small p’s,
* TFisher: a few good genes while ignoring many nulls.

To avoid betting on a single model, MAGCAT constructs a pathway-level omnibus p by combining the three tests with ACAT again (ACAT-O style).

For each pathway, let:

* `p1 = p_ACAT`,
* `p2 = p_Fisher`,
* `p3 = p_TFisher`.

Define:

```text
T_omni = sum_{j=1}^3 w_j * tan( π * (0.5 - p_j) )
```

with equal weights `w_j = 1/3` by default, and

```text
p_omni = 0.5 - (1 / π) * atan(T_omni)
```

This `omni_p` is:

* small if any of the component tests is strongly significant,
* robust to unknown signal pattern (sparse vs polygenic vs tail-truncated).

`magcat_omni_acat_pathways()` returns a table with:

* `pathway_id`, `pathway_name`,
* input p’s (`acat_p`, `fisher_p`, `tfisher_p_analytic`),
* `omni_p` – omnibus ACAT p.

---

### 5. Multiple testing: BH and q-values

Across all pathways, MAGCAT allows standard multiple-testing correction:

1. Benjamini–Hochberg FDR:

   ```r
   res_omni$omni_p_BH <- p.adjust(res_omni$omni_p, method = "BH")
   ```

   Pathways with `omni_p_BH <= 0.05` (or 0.10) are often considered significant.

2. Storey q-values (if the `qvalue` package is installed):

   ```r
   res_omni$omni_p_q <- qvalue::qvalue(res_omni$omni_p)$qvalues
   ```

   - `omni_p_q` directly estimates FDR at each pathway.
   - You can distinguish:
     - `q <= 0.05` → high-confidence pathways,
     - `0.05 < q <= 0.10` → suggestive / discovery pathways.

---

## Function reference (high-level map)

### MAGMA integration

* `magma_set_path()`, `magma_path()` – configure/query the MAGMA binary path.
* `run_magma_annot()` – wrapper around `magma --annotate` to produce `.genes.annot` files.
* `run_magma_gene()` – wrapper around `magma --gene-analysis` with `--gene-model multi=snp-wise`.
* `read_magma_genes()` – read MAGMA `.genes.out` into a data.frame with `GENE`, `NSNPS`, `ZSTAT`, `P`, etc.

### Gene location & length from GFF3

* `get_gene_lengths_gff3()` – parse a GFF3 to a gene table: `gene_id`, `chr`, `start`, `end`, `length`, `strand`.
* `write_genes_loc()` – export GFF3-derived gene coordinates to MAGMA `*.genes.loc` format.

Example resources are stored under `inst/extdata/GFF3/` and `inst/extdata/Zea_mays_gene_lengths.tsv`.

### Gene-level p-value adjustment

* `magcat_adjust_gene_p()` (in `Adjust_pvalues.R`)
  * Fits: `Z_i ~ log(L_i) + log(S_i)`
  * Creates `Z_adj` and `P_adj`
  * Returns augmented gene table + `lm` fit.

### Pathway loading

* `magcat_load_pathways(species)` – load built-in PMN pathways for `"maize"`, `"sorghum"`, `"arabidopsis"`, `"plant"`.

Returns a data.frame or list mapping `pathway_id` → gene IDs. You can also pass your own pathways.

### Pathway-level tests

* `magcat_acat_pathways()` – ACAT-based gene → pathway test, plus optional permutation-calibrated ACAT.
* `magcat_fisher_pathways()` – unweighted Fisher’s method across gene p-values in each pathway.
* `magcat_soft_tfisher_pathways()` – soft TFisher (truncated/weighted Fisher) with analytic and/or permutation p-values.
* `magcat_ordmeta_pathways()` (optional) – ordmeta-based pathway test (minimum marginal p in joint order distribution).

### Omnibus

* `magcat_omni_acat_pathways()` – ACAT-O style omnibus combination of multiple pathway p’s into a single `omni_p`.

---

## Citation

If you use MAGCAT in a publication, please cite:

* MAGMA: de Leeuw et al., PLoS Comput Biol (2015).
* ACAT: Liu et al., Am J Hum Genet (2019).
* TFisher / truncated Fisher: Zaykin et al., Sheng & Yang, Zhang & Wu, etc.
* And this package (GitHub / Zenodo DOI once available).
