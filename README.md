# CATFISH : GWAS pathway analysis pipeline

This Markdown is structured into:

- **INTRODUCTION** — what the pipeline is and why multiple tests are needed  
- **METHODS** — paper‑ready subsections with explicit equations and assumptions  
- **USAGE** — installation + reproducible example commands and R snippets
- **RESULTS** — an example analysis for soil Nitrogen GWAS

---

## INTRODUCTION

### What CATFISH is

**CATFISH** (**C**ombining **A**CAT, **T**Fisher (soft), **F**isher, m**I**n-P, and **S**touffer for **H**olistic pathway analysis) is a model-averaged pathway framework that combines ACAT, soft TFisher, Fisher, Stouffer, and minP on top of LD-aware MAGMA gene-level GWAS statistics adjusted for gene length and SNP density. It then uses an omnibus test based on permutation-calibrated minP or ACAT to collapse these multiple pathway tests into a single, correlation-robust enrichment p-value that is sensitive to both sparse and polygenic pathway patterns.

CATFISH uses:

1. **MAGMA** for LD-aware **SNP → gene** inference (gene-level p-values).
2. **Multiple gene → pathway combination tests** (ACAT, Fisher, soft TFisher, Stouffer).
3. A **correlation-robust omnibus test** (permutation-calibrated minP and ACAT-O) that aggregates these correlated component tests into a single pathway-level p-value.

### Why multiple tests are needed

Pathways can be significant for statistically different reasons:

- one driver gene dominates,
- many genes show coordinated moderate enrichment,
- diffuse polygenic shift, and
- hybrids of these patterns.

No single gene-set statistic is uniformly most powerful across these various possibilities. Instead of relying on a single test, CATFISH runs several complementary pathway tests and combines them into one omnibus p-value. We summarize these patterns using a set of **pathway signal archetypes** (sparse driver, coordinated moderate enrichment, diffuse polygenic shift, hybrid driver–support, single-gene proxy), which describe different ways a pathway can appear significant.

# Pathway signal archetypes

CATFISH interprets each significant pathway by classifying the **rank profile** of its gene-level p-values into one of several archetypes and reporting which component test drove significance.

## Archetype I — Sparse Driver Architecture (SDA)

**Signature:** one or a few genes are extremely significant; most genes look null.

- The top gene p-value `p_(1)` is much smaller than alpha (e.g. 1e-6),
- The rest of the genes in the pathway have p-values that look roughly uniform.

**Best detectors in CATFISH:**

- **ACAT** (loves a few very small p’s)
- **minP / Tippett** (minimum p-value)
- **Hard / heavy truncation** (e.g. TFisher with small tau)

**Interpretation:**  
The pathway is significant because of **driver gene dominance**, not broad engagement. Biologically, this could be a core “bottleneck” gene whose annotation pulls in a whole pathway.

---

## Archetype II — Coordinated Moderate Enrichment (CME)

**Signature:** many genes show moderate association; no single gene is insanely extreme.

- A non-trivial fraction of genes have p in, say, `[1e-3, 0.05]`,
- The top p-value is not massively more extreme than the rest.

**Best detectors in CATFISH:**

- **Fisher’s method** (sum of log p’s)
- Optionally **wFisher / mean-Z** if you ever add weights

**Interpretation:**  
This reflects **collective functional engagement**: the pathway as a whole is involved, even if no single gene is a monster hit.

---

## Archetype III — Diffuse Polygenic Shift (DPS)

**Signature:** the pathway’s genes are, on average, slightly more associated than the genome-wide background, but almost none cross a conventional 0.05 threshold.

- Most p-values are > 0.05,
- But the **mean adjusted Z (`Z_adj`) is shifted** away from 0 in one direction.

**Best detectors in CATFISH:**

- **Stouffer / mean-Z on `Z_adj`** (unweighted; permutation-calibrated)
- Optionally **competitive regression models** (e.g. MAGMA competitive) if you include them

**Interpretation:**  
This is a **global pathway bias consistent with polygenicity** – lots of tiny pushes in the same direction, no obvious star gene. It’s the “many gnats, no dragon” scenario.

---

## Archetype IV — Hybrid Driver–Support (HDS)

**Signature:** a few very strong genes plus a supporting cast of moderately associated genes.

- The top one or few genes have very small p-values (e.g. < 1e-4),
- Several additional genes have p in `[1e-3, 0.05]`.

**Best detectors in CATFISH:**

- **Soft TFisher** (truncated / tail-focused, but not as extreme as minP)
- **Fisher** (picks up the moderate bulk)
- **Omnibus combo** (ACAT over ACAT/Fisher/TFisher/Stouffer)

**Interpretation:**  
This fits a **hierarchical pathway organization**: a few “driver” genes plus **supporting machinery**. It’s often what you expect for key biosynthetic or signaling pathways.

---

## Archetype V — Single-Gene Proxy Pathway (SGP)

**Signature:** the entire pathway signal is explained by one gene; remove that gene and the pathway is no longer significant.

Operationally you can flag it by:

- Removing the top gene `g*` from the pathway and recomputing the pathway p-value,
- Marking pathways where the p-value becomes non-significant after that removal.

**Best detectors / diagnostics in CATFISH:**

- **minP / Tippett** (will happily call these significant),
- **ACAT** (also loves a single extreme gene),
- **SGP check:** “recompute without top gene” as a *diagnostic*, not as another omnibus test.

**Interpretation:**  
The pathway is essentially a **proxy for a single gene-level association** (often due to dense annotation or overlapping pathway definitions). Biologically, the pathway may still be relevant, but you should interpret it as “this gene is driving everything.”

---

## Recommended reporting checklist (per significant pathway)

For each pathway that passes your FDR / q-value threshold, CATFISH should report:

- **Archetype call:**  
  SDA (Sparse driver) / CME (Coordinated moderate) / DPS (Diffuse polygenic) / HDS (Hybrid driver–support) / SGP (Single-gene proxy)

- **Driver test(s):**  
  Which component test(s) gave the strongest evidence?
  - ACAT vs Fisher vs soft TFisher vs Stouffer
  - Optionally, whether minP across methods was the omnibus winner.

- **SGP diagnostic:**  
  Does the pathway remain significant after removing the top gene (and maybe the top 2)?
  - If not → flag as **Single-Gene Proxy (SGP)**.

- **Bias controls checked:**
  - Gene size / gene length adjustment done (yes; via regression),
  - SNP density / NSNPs adjustment done (yes; via regression),
  - Pathway size (small vs huge),
  - LD reference panel used.

- **Calibration:**
  - Analytic vs permutation p-values used for each test,
  - Number of permutations `B` for each permutation-based p,
  - Any evidence of inflation / deflation in null simulations or permuted data.

---

### Bias warning

When calculating over‑representation or enrichment, pathway‑based analysis is subject to multiple biases; **gene size**, **pathway size**, **density of SNP coverage**, and **linkage disequilibrium (LD)** patterns are all factors that must be considered and appropriately addressed (White et al., 2020; PMC6391732), which we have tried to address in our pipeline.  
Link: https://pmc.ncbi.nlm.nih.gov/articles/PMC6391732/

---

## METHODS

### Notation

Let a pathway (gene set) be denoted by $S$, containing $G = |S|$ genes indexed by $g = 1,\dots,G$.

Let the (adjusted) gene-level p-values be:

$$
\mathbf{p}_S = (p_1, p_2, \dots, p_G).
$$

Let the ordered p-values be:

$$
p_{(1)} \le p_{(2)} \le \dots \le p_{(G)}.
$$

---

## 1) Gene-level association statistics (SNP → gene)

For each gene $g$, MAGMA produces a gene‑level association p‑value $p_g$ by aggregating SNP‑level signals within/near the gene while accounting for local LD using a reference panel.

Conceptually:

- SNPs are mapped to genes (gene boundaries with optional windows).
- A multi‑marker gene model accounts for LD among SNPs in the gene region.
- MAGMA outputs gene statistics (e.g., $Z_g$ and $p_g$).

> In CATFISH, these MAGMA gene p-values are treated as the primitive inputs to all pathway tests.

---

## 2) Competitive enrichment framing (pathway vs background)

Pathway interpretation is most stable under **competitive** enrichment logic: a pathway is significant if its member genes are *more associated* than genes outside the pathway (i.e., relative enrichment), rather than simply showing absolute polygenicity.

Practical implications:

- Prefer competitive gene‑set testing frameworks (e.g., MAGMA competitive gene‑set regression) when available.
- When using gene‑p combination tests (ACAT/Fisher/TFisher) as pathway detectors, interpret them as **enrichment summaries** and retain covariate controls / calibrations (below) to mitigate confounding.

---

## 3) Optional gene-level adjustment for gene size and SNP density

Even with LD-aware gene testing, gene‑level signals can exhibit residual dependence on gene size and SNP density. CATFISH optionally performs a post‑hoc adjustment at the gene level.

Let:

- $Z_g$ be the MAGMA gene Z‑statistic,
- $L_g$ be gene length (bp),
- $S_g$ be number of SNPs mapped to the gene (e.g., `NSNPS`).

Fit:

$$
Z_g = \beta_0 + \beta_1 \log(L_g) + \beta_2 \log(S_g) + \varepsilon_g.
$$

Define adjusted residual Z:

$$
Z^{\mathrm{adj}}_g = Z_g - \widehat{Z}_g,
\quad \widehat{Z}_g = \widehat{\beta}_0 + \widehat{\beta}_1 \log(L_g) + \widehat{\beta}_2 \log(S_g).
$$

Convert to a two‑sided adjusted p-value:

$$
p^{\mathrm{adj}}_g = 2\Phi\left(-|Z^{\mathrm{adj}}_g|\right),
$$

where $\Phi(\cdot)$ is the standard normal CDF.

---

## 4) Pathway-level test statistics (gene → pathway)

CATFISH computes multiple pathway statistics from $\{p_g\}_{g \in S}$.

### 4.1 ACAT (Aggregated Cauchy Association Test)

Define the Cauchy‑transformed score for each gene:

$$
t_g = \tan\left(\pi\left(\tfrac{1}{2} - p_g\right)\right).
$$

Define non‑negative weights $w_g \ge 0$ with $\sum_{g \in S} w_g = 1$ (default $w_g = 1/G$).

The ACAT statistic is:

$$
T_{\mathrm{ACAT}}(S) = \sum_{g \in S} w_g\, t_g
= \sum_{g \in S} w_g \tan\left(\pi\left(\tfrac{1}{2} - p_g\right)\right).
$$

The combined p-value is:

$$
p_{\mathrm{ACAT}}(S) = \tfrac{1}{2} - \frac{1}{\pi}\arctan\left(T_{\mathrm{ACAT}}(S)\right).
$$

**Key property (interpretation):** ACAT is asymptotically dominated by the smallest p-values, and is therefore sensitive to **sparse driver** architectures.

---

### 4.2 Fisher’s method (coordinated enrichment)

Fisher’s statistic is:

$$
T_{\mathrm{Fisher}}(S) = -2\sum_{g \in S} \log(p_g).
$$

Under independence,

$$
T_{\mathrm{Fisher}}(S) \sim \chi^2_{2G},
\quad
p_{\mathrm{Fisher}}(S) = 1 - F_{\chi^2_{2G}}\left(T_{\mathrm{Fisher}}(S)\right),
$$

where $F_{\chi^2_{2G}}(\cdot)$ is the $\chi^2$ CDF with $2G$ degrees of freedom.

**Key property (interpretation):** Fisher is sensitive to **coordinated moderate enrichment** (many moderately small p-values).

---

### 4.3 Weighted Fisher (optional)

A common weighted variant is:

$$
T_{\mathrm{wFisher}}(S) = -2\sum_{g \in S} w_g \log(p_g),
$$

with weights $w_g$ reflecting, for example, gene importance, expression specificity, or QC confidence (must be pre‑specified to avoid circularity).

---

### 4.4 Truncation / TFisher (tail-focused enrichment)

To emphasize only the lower tail (top genes), define a truncation threshold $\tau \in (0,1)$ and consider:

#### (A) Hard truncated Fisher (conceptual)

$$
T_{\mathrm{TF}}(S;\tau) = -2\sum_{g \in S: p_g \le \tau} \log(p_g).
$$

This interpolates between minP‑like behavior (very small $\tau$) and Fisher (large $\tau$).

#### (B) Soft TFisher (recommended)

Soft‑truncation uses a continuous tail‑weighting scheme that downweights large p-values without an abrupt cutoff. In practice, CATFISH/MAGCAT can compute **soft TFisher** using an analytic null CDF (and optionally permutation calibration).

**Key property (interpretation):** truncation/TFisher is powerful for **hybrid driver–support** architectures (a few strong genes + several moderate ones).

---

### 4.5 Stouffer’s method (mean-Z; diffuse polygenic shift)

Stouffer aggregates gene-level **Z** statistics (e.g., `Z_adj`) rather than p-values. For a pathway with genes `g = 1..G`, define:

- **Unweighted Stouffer:**
  - `Z_stouffer = (1/sqrt(G)) * sum_g Z_g`
  - Convert to p-value: `p_stouffer = 2 * Φ( -|Z_stouffer| )`

- **Weighted Stouffer (optional):**
  - Choose weights `w_g >= 0`
  - `Z_stouffer = (sum_g w_g * Z_g) / sqrt( sum_g w_g^2 )`
  - `p_stouffer = 2 * Φ( -|Z_stouffer| )`

Notes:
- Best for **Diffuse Polygenic Shift (DPS)** where the **average Z** is shifted but few genes pass 0.05.
- If you calibrate by permutation, compute `Z_stouffer` on permuted gene labels and use an empirical p-value.

---

### 4.6 minP / Tippett (single-gene proxy diagnostic)

Define the pathway’s minimum gene p-value:

- `p_min = min_g p_g`

Under independence, Tippett’s combined p-value is:

- `p_tippett = 1 - (1 - p_min)^G`

But because gene p-values are correlated (LD, shared biology), CATFISH typically treats **minP as a diagnostic** and/or calibrates it by permutation:

- `p_min_perm = (1 + #{b: p_min^(b) <= p_min_obs}) / (B + 1)`

Notes:
- Highly sensitive to **single extreme genes**, so it flags **Sparse Driver (SDA)** and especially **Single-Gene Proxy (SGP)** pathways.
- Use “remove top gene and recompute” as the practical SGP check.

---

## 5) Correlation among pathway tests

All pathway statistics above are functions of the **same** gene-level p-values $\{p_g\}$, and are therefore intrinsically correlated:

$$
\left(T_{\mathrm{ACAT}},\,T_{\mathrm{Fisher}},\,T_{\mathrm{TF}}(\tau)\right)
\text{ are dependent}.
$$

Examples of why:

- Fisher and TFisher share overlapping subsets of $\{p_g\}$.
- ACAT and TFisher are both strongly influenced by the extreme left tail (smallest $p_g$).
- LD and shared biology induce correlation among gene-level signals, further increasing dependence.

**Implication:** naïve combination rules assuming independence between component pathway tests would be invalid (often anti‑conservative).

---

## 6) Omnibus pathway testing (correlation-robust model averaging)

Because no single test is uniformly optimal across latent pathway architectures, for each pathway `S` we compute a small panel of component p-values:

- `p_ACAT(S)`       – ACAT on gene-level p-values
- `p_Fisher(S)`     – Fisher’s sum-of-logs combination
- `p_TF(S; tau)`    – soft truncated Fisher (TFisher) at threshold `tau`
- `p_Stouffer(S)`   – (optionally weighted) Stouffer / mean-Z combination
- `p_minP(S)`       – within-pathway minP (Tippett) across genes

We then summarize these into *two* omnibus statistics:
(i) an analytic **ACAT-O** combination across methods, and  
(ii) a permutation-calibrated **minimum-p (minP)** across methods.

---

### 6.1 Omnibus ACAT (ACAT-O across methods)

Let `p1, p2, p3, p4, p5` denote the five component p-values for pathway `S`, and let weights `v_j >= 0` satisfy `sum_j v_j = 1` (default `v_j = 1/5`).

Define the ACAT-O Cauchy statistic:

- `T_omni_ACAT(S) = sum_{j=1..5} v_j * tan( pi * (0.5 - p_j) )`

Under the global null, `T_omni_ACAT(S)` is approximately standard Cauchy, giving the analytic omnibus p-value:

- `p_omni_ACAT(S) = 0.5 - (1/pi) * atan( T_omni_ACAT(S) )`

This ACAT-O layer is most sensitive when **at least one** component test (e.g., ACAT for sparse drivers, Fisher / Stouffer for coordinated enrichment, TFisher for hybrid patterns, or minP for hard single-gene hits) is strongly significant, even if the others are not.

---

### 6.2 Omnibus minP across methods with permutation

To obtain a complementary, more conservative summary that explicitly accounts for correlation between component tests, we also compute a minimum-p omnibus statistic across methods:

- `T_omni_min(S) = min_{j in {1,2,3,4,5}} p_j`

equivalently:

- `T_omni_min(S) = min( p_ACAT(S), p_Fisher(S), p_TF(S; tau), p_Stouffer(S), p_minP(S) )`

Because the `p_j` are correlated (all are computed from the same set of gene-level p-values), we calibrate `T_omni_min(S)` by permutation.

For each of `B` permutations:

- randomize gene labels across pathways (preserving the empirical distribution of gene-level p-values),
- recompute all five component tests for each pathway,
- record, for permutation `b = 1..B`:

  - `T_omni_min^(b)(S) = min_j p_j^(b)(S)`

The permutation-based omnibus p-value is:

- `p_omni_min_hat(S) = (1 + #{ b : T_omni_min^(b)(S) <= T_omni_min(S) }) / (B + 1)`

We use `p_omni_min_hat(S)` as the **primary omnibus p-value** in the main analyses because it:

1. protects nominal type-I error under *arbitrary* dependence between component tests, and  
2. explicitly targets the “best” component test for each pathway while accounting for the fact that the best test is chosen post hoc (via the min across methods).

The analytic ACAT-O p-value `p_omni_ACAT(S)` is reported alongside as a **higher-power, model-based sensitivity analysis**, highlighting pathways that are consistently strong across methods or dominated by a single very informative test.


---

## 7) Optional permutation calibration (when you want empirical p-values)

Permutation can be used to calibrate pathway p-values (especially truncation-based statistics) under complex dependence and finite-sample quirks.

A generic permutation p-value is:

$$
p_{\mathrm{perm}} = \frac{1 + \sum_{b=1}^{B} \mathbf{1}\left(T^{(b)} \ge T^{\mathrm{obs}}\right)}{B + 1}.
$$

**Note:** the minimum achievable permutation p-value is $1/(B+1)$; increasing $B$ increases resolution but should not be interpreted as “making results more significant” unless the observed statistic truly lies in the extreme tail.

---

## 8) Multiple testing correction

Across all pathways, omnibus p-values $\{p_{\mathrm{omni}}(S)\}$ are adjusted using Benjamini–Hochberg FDR:

$$
q_{\mathrm{BH}}(S) = \mathrm{BH}\left(p_{\mathrm{omni}}(S)\right).
$$

Because each pathway yields a **single** omnibus p-value, no additional penalty is required for the number of component tests.

---

## USAGE

# MAGCAT (R package wrapper)

**MAGCAT** is the R interface implementation used to run CATFISH‑style workflows on top of MAGMA, and to compute ACAT/Fisher/TFisher + omnibus pathway statistics.

---

## Installation

### 1) Install MAGMA (external dependency)

Download MAGMA from the official site and make the `magma` executable available on your `PATH`:

- https://ctg.cncr.nl/software/magma

### 2) Install MAGCAT in R

```r
# install.packages("devtools")  # if needed
devtools::install_github("nirwan1265/MAGCAT")
library(MAGCAT)
```

### 3) Optional: set MAGMA path

```r
MAGCAT::magma_set_path("/full/path/to/magma")
```

---

## Conceptual workflow (end-to-end)

1. **SNP → gene (MAGMA)**
   - Prepare SNP locations (`*.snp.loc`) and gene locations (`*.genes.loc`).
   - Run MAGMA annotation and gene analysis to get gene Z and p.

2. **Gene-level adjustment (optional)**
   - Regress $Z_g$ on `log(gene_length)` and `log(NSNPS)`; derive $p^{adj}_g$.

3. **Gene → pathway tests**
   - Compute pathway p-values from adjusted gene p-values using:
     - ACAT,
     - Fisher,
     - soft TFisher (tail-focused).

4. **Omnibus**
   - Combine pathway p-values using ACAT to produce $p_{\mathrm{omni}}$.

5. **Multiple testing**
   - BH FDR (and optional Storey q-values).

---

## MAGMA commands (typical)

```bash
# 1) Annotate SNPs to genes
magma \
  --annotate \
  --snp-loc  <snp.loc> \
  --gene-loc <genes.loc> \
  --out      <prefix>

# 2) Gene analysis (LD-aware)
magma \
  --bfile      <LD_reference_panel_prefix> \
  --pval       <gwas.pval.txt> N=<N> \
  --gene-annot <prefix>.genes.annot \
  --gene-model multi=snp-wise \
  --out        <prefix>
```

---

## Quick R example (MAGCAT-style)

> **Note:** This is the exact end-to-end pipeline used (MAGMA → gene adjustment → pathway tests → omnibus). Paths, filenames, and column mappings should be edited to match your local files.

```r
############################################################
## MAGCAT PIPELINE: MAGMA → SIZE-ADJUSTED GENE P → PATHWAYS
## (with brief explanation of the key parameters)
############################################################

## MAGCAT gives you:
##  - R wrappers around the MAGMA binary (magma_annotate, magma_gene)
##  - Parallel chromosome-wise runs (chroms, n_threads)
##  - Plant-ready GFF3 → gene loc helpers
##  - Automatic gene-size/#SNP adjustment
##  - A suite of pathway tests + an omnibus layer
##
## Below is the end-to-end flow, with parameter explanations in comments.
############################################################


############################################################
## 1. Build MAGMA gene location file from GFF3
############################################################

gff_path     <- "/Users/.../Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3"
gene_loc_out <- "/Users/.../maize.genes.loc"

gff3_to_geneloc(
  gff        = gff_path,
  out        = "inst/extdata/maize.genes.loc",
  chr_prefix = "chr"  # strips "chr" prefix so MAGMA’s CHR matches PLINK CHR
)

## What this does / why:
## - Parses the GFF3, extracts gene features, and writes a MAGMA-ready
##   gene location file (GENE, CHR, START, STOP, STRAND…).
## - When you later call species = "maize" in MAGCAT, it automatically
##   uses inst/extdata/maize.genes.loc and you don’t have to think
##   about coordinates again.


############################################################
## 2. SNP → gene annotation with MAGMA (magma_annotate wrapper)
############################################################

stats_file1 <- "/Users/.../nitrogen_0-5cm_maize_LMM.txt"

magma_annotate(
  stats_file     = stats_file1,
  rename_columns = c(
    CHR    = "chr",    # your column "chr" → MAGMA expects "CHR"
    SNP    = "rs",     # your column "rs"  → "SNP"
    POS    = "ps",     # your column "ps"  → "POS"
    PVALUE = "p_wald"  # your p-value column → "PVALUE"
  ),
  species    = "maize",        # uses built-in maize.genes.loc from step 1
  out_prefix = "N_maize_MLM",  # prefix for MAGMA output files
  out_dir    = "annot",        # directory to write MAGMA .genes.annot
  window     = c(25, 25)       # +/- 25 kb around each gene for SNP mapping
)

## Under the hood:
## - MAGCAT builds the MAGMA command line and calls the MAGMA binary for you.
## - You just supply stats_file and column rename map; no manual shell scripting.


############################################################
## 3. Gene-level MAGMA (multi = snp-wise) with R wrapper
############################################################

bfile      <- "/Users/.../all_maize2"   # PLINK basename: .bed/.bim/.fam
stats_file2 <- "/Users/.../raw_GWAS_MLM_3PC_N.txt"

## 3A. Chromosome-wise run using NMISS (per-SNP sample size)
magma_gene(
  bfile      = bfile,
  gene_annot = "annot/N_maize_MLM.genes.annot",
  stats_file = stats_file2,
  n_total    = 3539,             # optional global N if NOBS/NMISS absent
  rename_columns = c(
    CHR    = "Chr",
    SNP    = "SNP",
    POS    = "Pos",
    PVALUE = "P.value",
    NMISS  = "n_miss"            # MAGMA interprets this as per-SNP sample size
  ),
  out_prefix = "N_maize_MLM",
  out_dir    = "magma_genes_by_chr",
  gene_model = c("multi=snp-wise"),  # multi-parameter SNP-wise model
  chroms     = 1:10,                 # run MAGMA separately for chr 1–10
  n_threads  = 10                    # run up to 10 MAGMA jobs in parallel
)

## 3B. Chromosome-wise run using NOBS (per-SNP N)
magma_gene(
  bfile      = bfile,
  gene_annot = "annot/N_maize_MLM.genes.annot",
  stats_file = stats_file2,
  rename_columns = c(
    CHR    = "Chr",
    SNP    = "SNP",
    POS    = "Pos",
    PVALUE = "P.value",
    NOBS   = "nobs"             # per-SNP N; preferred over n_total
  ),
  out_prefix = "N_maize_MLM",
  out_dir    = "magma_multi_snp_wise_genes_by_chr_N_maize",
  gene_model = c("multi=snp-wise"),
  chroms     = 1:10,
  n_threads  = 10
)

## Selling point:
## - magma_gene is a high-level R wrapper around the MAGMA binary:
##   * Handles all flags, temp files, and error checking.
##   * Accepts chroms + n_threads to run multiple chromosomes in parallel.
##   * Only requires a minimal rename_columns spec instead of reformatting
##     your GWAS files by hand.


############################################################
## 4. Combine per-chromosome MAGMA gene outputs
############################################################

files <- sprintf(
  "/Users/.../magma_multi_snp_wise_genes_by_chr_N_maize/N_maize_MLM_chr%d.multi_snp_wise.genes.out",
  1:10
)

gene_list <- lapply(files, function(f) {
  if (!file.exists(f)) stop("File not found: ", f)
  utils::read.table(f, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
})

genes_all_raw <- do.call(rbind, gene_list)

# Ensure the gene p column is named "P"
colnames(genes_all_raw)[9] <- "P"

# For genes appearing on multiple chromosomes or windows, keep smallest P
o         <- order(genes_all_raw$GENE, genes_all_raw$P)
genes_all <- genes_all_raw[o, ]
genes_all <- genes_all[!duplicated(genes_all$GENE), ]

# Optionally sort by CHR and START
if (all(c("CHR", "START") %in% names(genes_all))) {
  genes_all <- genes_all[order(genes_all$CHR, genes_all$START), ]
}

write.table(
  genes_all,
  file      = "/Users/.../magma_N_maize.txt",
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

## Why:
## - This consolidates 10 per-chromosome MAGMA outputs into a single
##   gene-level table with one row per gene, ready for downstream
##   adjustment and pathway analysis.


############################################################
## 5. Gene length extraction + bias adjustment
############################################################

## 5A. Extract gene lengths from GFF3
gff3_path <- system.file(
  "extdata", "GFF3",
  "Zea_mays.Zm-B73-REFERENCE-NAM-5.0.62.chr.gff3",
  package = "MAGCAT"
)

maize_gene_len <- get_gene_lengths(
  gff3_file  = gff3_path,
  output     = TRUE,
  output_dir = "inst/extdata",
  file_name  = "Zea_mays_gene_lengths.tsv"
)

## Output includes:
##   gene_id, chr, start, end, length
## and can be reused across traits.

## 5B. Adjust gene-level Z/P for gene length & #SNPs

genes_all <- read.table(
  "/Users/.../magma_N_maize.txt",
  header = TRUE, stringsAsFactors = FALSE
)

adj_out <- magcat_adjust_gene_p(
  gene_results = genes_all,
  gene_lengths = maize_gene_len,
  gene_col     = "GENE",
  p_col        = "P",
  z_col        = "ZSTAT",    # raw MAGMA Z
  len_gene_col = "gene_id",
  len_col      = "length"
  # nsnp_col   = "NSNPS"     # optional: adjusts for #SNPs per gene as well
)

genes_adj <- adj_out$genes    # includes Z_raw, Z_adj, P_adj, log_gene_length, log_nsnp
lm_fit    <- adj_out$fit      # lm(Z_raw ~ log_gene_length + log_nsnp)

write.csv(genes_adj, "genes_adj.csv", row.names = FALSE)

## Why:
## - MAGMA gene p-values are known to correlate weakly with gene size
##   and SNP density. magcat_adjust_gene_p:
##   * Fits a linear model: Z_raw ~ log(gene length) + log(#SNPs).
##   * Uses residuals (Z_adj) as “size/SNP-adjusted” Z-scores.
##   * Converts Z_adj back to P_adj = 2*pnorm(-|Z_adj|).
## - You then use P_adj in gene→pathway tests, reducing bias toward big genes.


############################################################
## 6. Load pathway definitions from PMN/CornCyc
############################################################

maize_pw <- magcat_load_pathways(
  species  = "maize",
  gene_col = "Gene-name"  # column in the PMN gene-set file that matches MAGMA gene IDs
)

## Features:
## - MAGCAT ships plant pathway collections (CornCyc, AraCyc, etc.)
## - You can also pass a custom list/data.frame of pathways if you want.


############################################################
## 7. Pathway-level tests (gene → pathway)
############################################################

## All tests use gene_results + species/pathways to:
##  - find genes per pathway,
##  - take their p-values (raw P or adjusted P_adj),
##  - compute a pathway-level p per method.

### 7A. ACAT per pathway

pw_res_acat_adj <- magcat_acat_pathways(
  gene_results = genes_adj,     # use size/SNP-adjusted P_adj
  species      = "maize",
  gene_col     = "GENE",
  p_col        = "P_adj",
  B            = 0L,            # B>0 → permutation-calibrated ACAT
  seed         = NULL,
  output       = TRUE,
  out_dir      = "acat_results"
)

### 7B. ordmeta (rank-based, minimum marginal p over gene rank)

ord_res <- magcat_ordmeta_pathways(
  gene_results = genes_adj,
  species      = "maize",
  gene_col     = "GENE",
  p_col        = "P_adj",
  effect_col   = "ZSTAT",   # needs signed effects for up/down
  is_onetail   = FALSE
)

### 7C. Weighted Fisher (wFisher) pathways

wf_res_raw <- magcat_wfisher_pathways(
  gene_results = genes_all,   # raw P
  species      = "maize",
  gene_col     = "GENE",
  p_col        = "P",
  effect_col   = "ZSTAT",
  weight_col   = "NSNPS",     # weight genes by #SNPs
  is_onetail   = FALSE
)

wf_res_adj <- magcat_wfisher_pathways(
  gene_results = genes_adj,   # adjusted P
  species      = "maize",
  gene_col     = "GENE",
  p_col        = "P_adj",
  effect_col   = "ZSTAT",
  is_onetail   = FALSE        # two-sided test with sign from ZSTAT
)

### 7D. Truncated Fisher / TFisher

tf_res_adj <- magcat_tfisher_pathways(
  gene_results     = genes_adj,
  species          = "maize",
  gene_col         = "GENE",
  p_col            = "P_adj",
  ptrunc           = 0.05,     # only genes with p <= 0.05 contribute
  B_perm           = 10000L,   # gene-set permutations for empirical p
  seed             = 123,
  analytic_logical = TRUE,     # also report chi-square approx p
  output           = TRUE,
  out_dir          = "magcat_tfisher"
)

### 7E. Soft TFisher (smooth truncation)

soft_tf_res_adj <- magcat_soft_tfisher_pathways(
  gene_results     = genes_adj,
  species          = "maize",
  gene_col         = "GENE",
  p_col            = "P_adj",
  tau1             = 0.05,     # soft truncation threshold
  B_perm           = 10000L,
  seed             = 123,
  analytic_logical = TRUE,
  output           = TRUE,
  out_dir          = "magcat_tfisher_soft"
)

### 7F. Stouffer (sum of Z across genes)

stouf_res <- magcat_stouffer_pathways(
  gene_results = genes_adj,
  species      = "maize",
  gene_col     = "GENE",
  p_col        = "P_adj",
  weight_col   = NULL,    # equal weights for all genes
  B_perm       = 100L,    # permutations for empirical Stouffer p
  seed         = 123,
  output       = TRUE,
  out_dir      = "magcat_stouffer"
)

### 7G. Gene-level minP per pathway

minp_res <- magcat_minp_pathways(
  gene_results = genes_adj,
  species      = "maize",
  gene_col     = "GENE",
  p_col        = "P_adj",
  B_perm       = 0L,      # analytic only (can turn on permutations)
  min_p        = 1e-15,
  do_fix       = TRUE,
  output       = TRUE,
  out_dir      = "magcat_minp_maize"
)


############################################################
## 8. Omnibus over methods (ACAT-O or method-level minP)
############################################################

omni_minp <- omni_pathways(
  gene_results      = genes_adj,
  species           = "maize",
  gene_col          = "GENE",
  p_col             = "P_adj",
  effect_col        = "Z_adj",     # use adjusted Z for direction
  #weight_col       = "NSNPS",    # optional, if you want wFisher to be weighted
  is_onetail        = FALSE,
  ptrunc            = 0.05,        # passed to internal TFisher component
  min_p             = 1e-15,       # floor for tiny p's (needed for ACAT stability)
  do_fix            = TRUE,
  omnibus           = "minP",      # "ACAT" = ACAT-O across methods, "minP" = minP over methods
  B_perm            = 10000L,      # permutations at the *omnibus* level
  seed              = 123,
  remove_singletons = TRUE,        # drop pathways with n_genes < 2
  output            = TRUE,
  out_dir           = "magcat_omni_full"
)

## In a single call omni_pathways gives you:
##   - acat_p       : ACAT over gene p’s per pathway
##   - wfisher_p    : wFisher per pathway
##   - tpm_p        : truncated Fisher per pathway
##   - stouffer_p   : Stouffer per pathway
##   - minp_gene_p  : gene-level minP per pathway
##   - omni_p       : combination of these methods (ACAT-O or method-level minP)
##   - omni_perm_p  : permutation-calibrated omnibus p
##   - BH FDR and q-values for omni_p and each component
##
## Selling point:
## - Users don’t have to orchestrate multiple scripts or re-run GWAS.
##   They:
##      (1) point MAGCAT to a GFF3 + PLINK + GWAS file,
##      (2) call magma_annotate + magma_gene once,
##      (3) call a single high-level pathway function (ACAT, TFisher, Stouffer,
##          minP, or omni_pathways),
##   and get a full panel of pathway p-values (plus permutation-calibrated
##   omnibus) with plant-aware defaults and parallel MAGMA integration.

```

---


## References

- White MJ et al. *Strategies for Pathway Analysis using GWAS and WGS Data*. Current Protocols in Human Genetics (2020)
  https://pmc.ncbi.nlm.nih.gov/articles/PMC6391732/
- de Leeuw CA et al. *MAGMA: Generalized Gene‑Set Analysis of GWAS Data*. PLoS Comput Biol (2015).
  https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004219
- Liu Y et al. *ACAT: A Fast and Powerful p Value Combination Method for Rare-Variant Analysis in Sequencing Studies*.  The American Journal of Human Genetics (2019).
  https://pubmed.ncbi.nlm.nih.gov/30849328/
- Zhang H et al. *TFisher: A powerful truncation and weighting procedure for combining p-values*. Annals of Applied Statistics (2020)
https://projecteuclid.org/journals/annals-of-applied-statistics/volume-14/issue-1/TFisher--A-powerful-truncation-and-weighting-procedure-for-combining/10.1214/19-AOAS1302.full
- Yoon S et al. *Powerful p-value combination methods to detect incomplete association*. Nature (2021)
  https://www.nature.com/articles/s41598-021-86465-y
- Tippett, L. H. C. *The Methods of Statistics*. London:Williams & Norgate (1931)
- Westfall, P. H., & Young, S. S. *Resampling-Based Multiple Testing: Examples and Methods for p-Value Adjustment*. New York: Wiley(1993)
  


