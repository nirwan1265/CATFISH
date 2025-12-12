# CATFISH : pathway analysis pipeline + signal archetypes

This Markdown is structured into:

- **INTRODUCTION** — what the pipeline is and why multiple tests are needed  
- **METHODS** — paper‑ready subsections with explicit equations and assumptions  
- **USAGE** — installation + reproducible example commands and R snippets  

---

## INTRODUCTION

### What CATFISH is

**CATFISH** is a gene‑set (pathway) analysis workflow that starts from SNP‑level association results and produces **pathway‑level significance** *plus* a **mechanistic interpretation** of *why* each pathway is significant.

CATFISH uses:

1. **MAGMA** for LD‑aware **SNP → gene** inference (gene-level p-values).
2. **Multiple gene → pathway combination tests** (ACAT, Fisher / wFisher, truncation/TFisher).
3. A **correlation‑robust omnibus test** that aggregates these correlated component tests.
4. A **signal archetype** framework to interpret pathway significance by its gene‑level p‑value rank profile.

### Why multiple tests are needed

Pathways can be significant for qualitatively different reasons:

- one **driver** gene dominates,
- many genes show **coordinated moderate** enrichment,
- or there is a **diffuse polygenic** shift.

No single gene‑set statistic is uniformly most powerful across these regimes. CATFISH therefore treats pathway detection as a **model‑averaging** problem over latent “signal architectures,” rather than a single one‑size‑fits‑all test.

# Pathway signal archetypes (interpretation layer)

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

Because no single test is uniformly optimal across latent pathway architectures, we compute multiple component p-values:

$$
\mathcal{P}(S) = \left\{ p_{\mathrm{ACAT}}(S),\, p_{\mathrm{Fisher}}(S),\, p_{\mathrm{TF}}(S;\tau) \right\}.
$$

We then combine them into a single omnibus p-value using ACAT again:

### 6.1 Omnibus ACAT statistic

Let $p_1, p_2, p_3$ denote the three component p-values and $v_j \ge 0$ be weights with $\sum_j v_j = 1$ (default $v_j = 1/3$).

Define:

$$
T_{\mathrm{omni}}(S) = \sum_{j=1}^{3} v_j \tan\left(\pi\left(\tfrac{1}{2} - p_j\right)\right).
$$

Then:

$$
p_{\mathrm{omni}}(S) =
\tfrac{1}{2} - \frac{1}{\pi}\arctan\left(T_{\mathrm{omni}}(S)\right).
$$

### 6.2 Why this is valid under dependence

ACAT is designed to be robust to **arbitrary dependence** among inputs in many practical settings because the Cauchy tail behavior yields stable combined p-values without needing an explicit covariance estimate.

**Interpretation:** omnibus significance indicates that *at least one plausible signal architecture* is supported by the data. It does **not** imply uniform enrichment across all genes.

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

> **Note:** function names below follow the style of your MAGCAT wrapper. If your package uses slightly different names/arguments, keep the structure and edit names to match your exported API.

```r
library(MAGCAT)

# 1) Read MAGMA gene results (.genes.out)
genes_file <- system.file(
  "magma_genes",
  "N_maize_MLM.snp_wise_top.genes.out",
  package = "MAGCAT"
)
genes_all <- read_magma_genes(genes_file)

# 2) Load gene lengths (from GFF3-derived table)
gene_len_file <- system.file("extdata", "Zea_mays_gene_lengths.tsv", package = "MAGCAT")
maize_gene_len <- readr::read_tsv(gene_len_file, show_col_types = FALSE)

# 3) Adjust gene p-values for gene length + NSNPS
adj_out <- magcat_adjust_gene_p(
  gene_results = genes_all,
  gene_lengths = maize_gene_len,
  gene_col     = "GENE",
  nsnp_col     = "NSNPS",
  p_col        = "P",
  z_col        = "ZSTAT",
  len_gene_col = "gene_id",
  len_col      = "length"
)
genes_adj <- adj_out$genes

# 4) Load pathways (built-in or user-supplied)
pathways_maize <- magcat_load_pathways(species = "maize")

# 5) Pathway tests
res_acat <- magcat_acat_pathways(
  gene_results = genes_adj,
  pathways     = pathways_maize,
  gene_col     = "GENE",
  p_col        = "P_adj",
  B_perm       = 0L
)

res_fisher <- magcat_fisher_pathways(
  gene_results = genes_adj,
  pathways     = pathways_maize,
  gene_col     = "GENE",
  p_col        = "P_adj"
)

res_tfsoft <- magcat_soft_tfisher_pathways(
  gene_results     = genes_adj,
  pathways         = pathways_maize,
  gene_col         = "GENE",
  p_col            = "P_adj",
  tau1             = 0.05,
  B_perm           = 1000L,
  analytic_logical = TRUE
)

# 6) Omnibus ACAT across component tests
res_omni <- magcat_omni_acat_pathways(
  res_acat    = res_acat,
  res_fisher  = res_fisher,
  res_tfisher = res_tfsoft,
  acat_col    = "acat_p",
  fisher_col  = "fisher_p",
  tfisher_col = "tfisher_p_analytic"
)

# 7) Multiple testing (BH)
res_omni$omni_p_BH <- p.adjust(res_omni$omni_p, method = "BH")

# Top pathways
head(res_omni[order(res_omni$omni_p), ], 10)
```

---

## Suggested output columns (pathway-level)

For each pathway, CATFISH/MAGCAT should report at minimum:

- `pathway_id`, `pathway_name`, `n_genes`
- `acat_p`
- `fisher_p` (and optionally `wfisher_p`)
- `tfisher_p_analytic` (and optionally `tfisher_p_perm`)
- `omni_p`
- `omni_p_BH` (and optional `omni_p_q`)
- `driver_component` (which test drove the omnibus)
- `archetype` (SDA/CME/DPS/HDS/SGP/HAA)

---

## References

- White MJ et al. *Strategies for Pathway Analysis using GWAS and WGS Data*. (PMC6391732)  
  https://pmc.ncbi.nlm.nih.gov/articles/PMC6391732/
- de Leeuw CA et al. *MAGMA: Generalized Gene‑Set Analysis of GWAS Data*. PLoS Comput Biol (2015).
- Liu Y et al. *ACAT / Cauchy combination test* (2019).
- TFisher / truncated Fisher family: Zaykin et al.; Sheng & Yang; Zhang & Wu (depending on implementation).


