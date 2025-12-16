# CATFISH : GWAS pathway analysis pipeline

This Markdown is structured into:

- **INTRODUCTION** — what the pipeline is and why multiple tests are needed  
- **METHODS** — paper‑ready subsections with explicit equations and assumptions  
- **USAGE** — installation + reproducible example commands and R snippets
- **RESULTS** — an example analysis for soil Nitrogen GWAS

---

## INTRODUCTION

### What CATFISH is

**CATFISH** (<ins>**C**</ins>ombining <ins>**A**</ins>CAT, <ins>**T**</ins>Fisher (soft), <ins>**F**</ins>isher, m<ins>**I**</ins>n-P, and <ins>**S**</ins>touffer for <ins>**H**</ins>olistic pathway analysis) is a model-averaged pathway framework that combines ACAT, soft TFisher, Fisher, Stouffer, and minP on top of LD-aware MAGMA gene-level GWAS statistics adjusted for gene length and SNP density. It then uses an omnibus test based on permutation-calibrated minP or ACAT to collapse these multiple pathway tests into a single, correlation-robust enrichment p-value that is sensitive to both sparse and polygenic pathway patterns. In short, CATFISH casts a wide net across complementary tests and reels in a single pathway p-value.

CATFISH uses:

1. **MAGMA** for LD-aware **SNP → gene** inference (gene-level p-values).
2. **Multiple gene → pathway combination tests** (ACAT, Fisher, soft TFisher, Stouffer, minP).
3. A **correlation-robust omnibus test** (permutation-calibrated minP and ACAT-O) that aggregates these correlated component tests into a single pathway-level p-value.

### Why multiple tests are needed

Pathways can be significant for statistically different reasons:

- one driver gene dominates,
- many genes show coordinated moderate enrichment,
- diffuse polygenic shift, and
- hybrids of these patterns.

No single gene-set statistic is uniformly most powerful across these various possibilities. Instead of relying on a single test, CATFISH runs several complementary pathway tests and combines them into one omnibus p-value. We summarize these patterns using a set of **pathway signal archetypes** (sparse driver, coordinated moderate enrichment, diffuse polygenic shift, hybrid driver–support, single-gene proxy), which describe different ways a pathway can appear significant.

# Pathway signal archetypes

We divide the pathway signals into a small set of archetypes that describe different ways a pathway can be enriched as explained above. We then use a combination of statistical tests chosen to be representative to each of these behaviors and a provide a biological example for each archetype.

## Archetype I — Sparse Driver Architecture (SDA)

**Signature:** one or a few genes are extremely significant.

**Gene-level p-value pattern:**

- The top gene(s) p-value is much smaller than the significance level α:
  $$p_{(1)} \ll \alpha \;\;(\text{e.g. } p_{(1)} \approx 10^{-6}\text{ or smaller}).$$

- The remaining gene p-values are approximately uniform:  
  $$(p_{(2)},\,p_{(3)},\,\dots,\,p_{(m)}) \sim \mathrm{Uniform}(0,1).$$

- A sharp “elbow” in the ranked p-values (one (or few) tiny hit, then a long flat tail).

**Interpretation:**  
An SDA pathway is significant because of **driver gene(s) dominance**, not broad engagement. One gene can carry the pathway signal when it sits at a mechanistic bottleneck (committed step, rate-limiting enzyme, key regulator, or essential transporter). Pathway annotations bundle a few true control points with many supporting enzymes and housekeeping genes; under SDA, trait-relevant variation concentrates in the control point(s) while the rest remain near-null. Statistically, this produces a tiny $p_{(1)}$ and a long tail of near-uniform values.

**Biological example:**  
**Aspartokinase in the aspartate-derived amino acid pathway**

The aspartate-derived amino-acid biosynthesis pathway converts **aspartate** into essential amino acids such as **lysine**, **threonine**, **methionine**, and **isoleucine**. In plants and bacteria, the first committed step is catalyzed by **aspartokinase (AK)**, which phosphorylates aspartate to **aspartyl-phosphate** and feeds a branched network that produces multiple end-products. Because AK sits at the top of this branchpoint, it acts as a **flux-controlling bottleneck**: changes in AK activity shift carbon and nitrogen flow through the entire network. Downstream enzymes (tailoring steps, dehydrogenases, transaminases) tend to have distributed or buffered roles. Hence the gene-level pattern matches SDA: one AK gene or a few genes downstream with very small $p$’s (e.g. $10^{-8}\text{-}10^{-10}$), while most others scatter over $(0,1)$.

**Best detectors in CATFISH:**
- **ACAT** — sensitive to a few extremely small p-values (RECOMMENDED).  
- **minP / Tippett** — targets the minimum p-value; optimal when one gene dominates.


---

## Archetype II — Coordinated Moderate Enrichment (CME)

**Signature:** many genes show moderate association; no single gene is insanely extreme.

**Gene-level p-value pattern:**
- A non-trivial fraction of genes have moderately small p-values, e.g.
  $$p_i \in [10^{-3},\,0.05]\quad\text{for many } i.$$
- The top signal is not orders-of-magnitude beyond the rest (no single-gene spike), e.g.
  $$p_{(1)} \not\ll p_{(k)}\quad\text{for small }k\ (\text{e.g., }k=5,10,20).$$
- The ranked p-values show a **broad shoulder** (many “pretty good” genes), not a sharp elbow.

**Interpretation:**  
CME reflects **collective functional engagement**: the pathway behaves like a coordinated module where many components contribute small-to-moderate effects. This is expected when phenotypes arise from distributed regulation, redundancy, and buffering/feedback (extreme single-gene effects are uncommon). Statistically, enrichment comes from many mildly informative genes rather than one dominating driver.

**Biological example:**  
**Cytokine / immune signaling cascades**

Cytokine and inflammatory signaling rarely operates through a single master switch. Instead, cascades such as **IL-6/JAK–STAT**, **TNF/NF-κB**, and **IL-1** function as multi-step circuits that sense damage or infection, propagate signals through adaptor/kinase layers, reprogram transcription, and then dampen the response through negative feedback. These modules contain many interacting parts that each tune inflammatory tone:

- **Cytokines** (e.g., *IL6*, *TNF*, *IL1B*) acting as soluble alarm signals.
- **Cell-surface receptors** (e.g., *IL6R*, TNFRSF family) that detect alarms and set sensitivity.
- **Intracellular adaptors and kinases** (e.g., JAKs, MAPKs, TRAFs, IKK complex) that relay signals through phosphorylation cascades.
- **Transcription factors** (e.g., STATs, RELA/NF-κB) that shift downstream gene expression programs.
- **Feedback and regulatory nodes** (inhibitors, suppressors, decoy receptors) that reshape, terminate, or localize responses.

Because pathway output emerges from **distributed tuning across many steps** (receptor abundance, kinase efficiency, adaptor binding, TF activation thresholds, and feedback strength), natural genetic variation tends to produce **many modest effects** rather than one overwhelming driver. This often appears as gene-level p-values clustered in a “pretty good” range across multiple pathway members, for example:
$$(p_{(1)},p_{(2)},\ldots)\  \text{contains many values around }10^{-3}\text{--}10^{-2}\ \text{without an extreme like }10^{-12}.$$
At the pathway level, the entire signaling circuit is slightly shifted, consistent with CME. In CATFISH, this pattern is best captured by tests that accumulate evidence across many genes, such as **Fisher** and **Stouffer/mean-Z** (and often **soft TFisher** with mild truncation)—because their power increases when many genes contribute moderate association rather than a single spike.

**Best detectors in CATFISH:**
- **Fisher’s method** (aggregates evidence across many moderately small p-values) (RECOMMENDED).
- **Stouffer / mean-Z** (gains power when many genes shift together).
- Optionally **wFisher / weighted Z** if you later add biologically informed weights.



---

## Archetype III — Diffuse Polygenic Shift (DPS)

**Signature:** the pathway’s genes are, on average, slightly more associated than the genome-wide background, but almost none cross a conventional significance threshold.

**Gene-level p-value / Z pattern:**
- Most gene p-values satisfy:
  $$p_i > 0.05\quad\text{for most } i.$$
Yet the pathway shows a small but consistent shift in adjusted gene-level Z-scores:

$$
\overline{Z}_{S,\mathrm{adj}} \;=\; \frac{1}{G}\sum_{i\in S} Z_{i,\mathrm{adj}} \;\neq\; 0,
$$

often with a coherent sign (bias in one direction).
- Equivalent distributional statement:
  $$\{p_i: i\in S\}\ \text{is subtly enriched toward smaller values relative to Uniform}(0,1),$$
  but without extreme outliers.
- Visual cue: no sharp elbow; instead, the ranked p-values show a gentle, global downward bend relative to null.

**Interpretation:**  
DPS reflects a **global pathway bias consistent with polygenicity**: many genes each contribute tiny effects in the same general direction, and the pathway is enriched because the entire module is subtly shifted rather than driven by a single “star” gene. This is the “many gnats, no dragon” regime. Spike-hunting tests (e.g., minP/Tippett or very aggressive truncation) are typically underpowered here because there is no single extreme p-value to exploit. In contrast, mean-/distribution-sensitive tests (Stouffer/mean-Z, Fisher, and competitive regression) are designed to detect exactly this low-amplitude, widespread deviation from the null.

**Biological example:**  

**Biological example – human height as a diffuse polygenic shift**

Human adult height is one of the clearest textbook examples of an extremely polygenic trait. Early GIANT meta-analyses in ~180k individuals identified 180 loci and “hundreds of variants” for height, yet these explained only about 10% of phenotypic variance, despite height’s ~80% heritability. Subsequent meta-analyses in ~700,000 Europeans expanded this to thousands of associated SNPs and hundreds of loci, reinforcing the view that height is influenced by very many common variants of tiny effect rather than a small set of large-effect genes. A recent “saturated map” analysis went further, reporting ~12,000 genome-wide significant SNPs in >7,000 genomic segments (covering ~21% of the genome), again consistent with Fisher’s original polygenic model for height proposed in 1918.

Biologically, these variants fall into multiple growth-related programs: growth plate chondrocyte proliferation and hypertrophy, extracellular matrix and cartilage organization, growth hormone and IGF-1 signalling, and morphogen pathways such as TGF-β and Hedgehog, along with broader developmental and endocrine regulators. Pathway and gene-set analyses of height GWAS loci have shown enrichment for signalling pathways like TGF-β and Hedgehog, and for genes involved in skeletal growth, growth plate regulation, and related Mendelian growth disorders. For example, Lango Allen et al. (2010) found that genes near height-associated variants cluster in biologically coherent pathways (TGF-β signalling, Hedgehog signalling, histone and growth/development gene sets) and that several SNPs near genes in these pathways “narrowly miss” genome-wide significance, implying many additional sub-threshold contributors within the same modules. Guo et al. (2018) similarly show that genes near height GWAS loci are enriched in growth-relevant pathways and tissues such as growth plate cartilage.

Viewed at the level of a single pathway (e.g. TGF-β signalling, Hedgehog signalling, or growth-plate extracellular matrix genes), this architecture naturally yields a **Diffuse Polygenic Shift (DPS)** pattern rather than a single screaming driver. Across the gene set, many genes carry one or more common variants with small effects on height; some genes achieve clear genome-wide significance, but many others show only modest or nominal association. The result is that, compared with random gene sets, the distribution of gene-level test statistics in these pathways is shifted toward stronger evidence overall—more genes with small or moderate p-values, fewer genes that look completely null—even though most individual genes would not, by themselves, justify a strong claim of association. This is exactly the situation where CATFISH’s DPS-oriented detectors (Stouffer/mean-Z on adjusted gene-level Z-scores, and optionally MAGMA-style competitive regression tests) are most appropriate: rather than hunting for a single driver gene, they test whether the **average** association signal across a biologically coherent pathway is subtly but consistently elevated relative to the genome-wide background.


**Best detectors in CATFISH:**
- **Stouffer / mean-Z on** `Z_adj` (unweighted; permutation-calibrated) (RECOMMENDED).
- Optionally **competitive regression-style gene-set models** (e.g., **MAGMA competitive**) if included as a component test (NOT INCLUDED IN CATFISH)


---

## Archetype IV — Hybrid Driver–Support (HDS)

**Signature:** a few very strong genes plus a supporting cast of moderately associated genes.

**Gene-level p-value pattern:**
- One or a few top genes are extremely significant, e.g.
  $$p_{(1)},\,p_{(2)} \ll 10^{-4}\quad(\text{often much smaller}).$$
- Beyond the top hits, several additional genes show moderate evidence:
  $$p_{(k)} \in [10^{-3},\,0.05]\quad\text{for multiple }k \text{ (support genes).}$$
- The remaining genes are near-null:
  $$p_{(j)} \sim \mathrm{Uniform}(0,1)\quad\text{for most other }j.$$
- Visual cue: a small “spike” at the top (drivers) plus a clear “shoulder” of moderately small p-values (support), then a flat tail.

**Interpretation:**  
HDS reflects a **hierarchical pathway organization**: a few “driver” genes dominate the strongest effects, while a set of **supporting machinery** genes contribute smaller but consistent association. This is common in pathways where flux or signal is gated by a small number of control points, yet successful pathway output also depends on coordinated activity across multiple downstream components. Statistically, this architecture sits between Sparse Driver (single spike) and Coordinated Moderate Enrichment (broad shoulder): there *is* a clear driver signal, but the pathway significance is strengthened by additional moderate signals rather than being purely carried by one gene.




**Biological example** 
**LDL cholesterol as a hybrid driver–support pathway**

LDL-cholesterol (LDL-C) regulation is a textbook case where a few genes act as major “drivers,” sitting on top of a broader polygenic background. Large GWAS and sequencing studies consistently highlight *LDLR*, *APOB* and *PCSK9* as core Mendelian hypercholesterolemia genes: rare coding or splice-altering variants in these genes can shift LDL-C by roughly half a standard deviation or more and underlie classical familial hypercholesterolemia, with greatly increased coronary artery disease risk. Recent whole-genome sequence analyses in >60,000 individuals extend this picture by showing that even rare *non-coding* variants near *LDLR* and *PCSK9* can have effects comparable to clinically reported exonic FH variants, reinforcing their role as dominant levers of LDL-C homeostasis.

Around these drivers sits a much larger supporting network of lipoprotein and cholesterol metabolism genes. GWAS of conventional lipids and NMR-based lipoprotein traits have identified dozens of additional loci influencing LDL particle size, concentration, and composition, including apolipoprotein clusters (*APOE/APOC*, *APOA1/A5*), hepatic lipase (*LIPC*), and transporters such as *ABCG5/ABCG8*. Individually, common variants at these loci tend to explain only a small fraction of LDL-C variance, but collectively they account for a substantial portion of the genome-wide polygenic signal, and large multi-ancestry meta-analyses now report hundreds of lipid-associated loci spread across this extended lipoprotein network.

If you translate this biology into gene-level association statistics for an LDL-related trait, the pathway does not look purely “sparse driver” (one screaming gene and nothing else), nor purely “coordinated moderate” (dozens of roughly equal signals). Instead, you typically see a handful of very strong gene-level signals at *LDLR*, *APOB*, *PCSK9*, and a few related loci, sitting on a broader base of moderately associated genes involved in lipoprotein assembly, remodeling, and cholesterol transport. That is precisely the **Hybrid Driver–Support (HDS)** pattern: a pathway whose enrichment reflects both a small set of dominant, large-effect genes and a genuine, but softer, contribution from the wider metabolic machinery. In CATFISH terms, this is the regime where **soft TFisher** (which emphasizes the lower tail while still counting moderate p-values), together with **Fisher** and the **omnibus combination**, is well matched to the underlying biology.


**Best detectors in CATFISH:**
- **Soft TFisher** (tail-focused; gains power from a few strong hits *plus* additional modest hits)
- **Fisher** (accumulates evidence across the moderate support set)
- **Omnibus combination** (e.g., ACAT across ACAT/Fisher/soft-TFisher/Stouffer; or permutation-calibrated minP across methods)


---

## Archetype V — Single-Gene Proxy Pathway (SGP)

**Signature:** the pathway looks significant only because it contains one very strong gene; the remaining members look essentially null.

**Gene-level pattern:**
- The top gene has an extremely small p-value (e.g. $10^{-6}\text{--}10^{-12}$).
- The rest of the genes have p-values that look noisy / $\mathrm{Uniform}(0,1)$, with no clear excess of small p’s.
- If you mentally “ignore” that top gene, there is no obvious enrichment left in the pathway.

**Interpretation:**  
The pathway is essentially acting as a proxy for a single gene-level association. This often happens when:
- an annotation term is very small (one gene plus a couple of weakly related neighbors), or
- multiple pathway definitions redundantly include the same driver gene, so several “different” pathways light up but all are tagging the same underlying gene.

Biologically, the driver gene can still be very important (and may also fit the Sparse Driver archetype), but the pathway-level claim carries no extra information beyond “this gene is strongly associated.” In CATFISH we therefore flag SGP patterns as a caution: interpret these as gene-centric hits with pathway labels attached, rather than as evidence that the entire pathway is collectively engaged.

**Biological example:**  
**PAH in phenylalanine metabolism**

Phenylalanine metabolism in humans is dominated by a single bottleneck enzyme, phenylalanine hydroxylase (PAH), a hepatic monooxygenase that hydroxylates phenylalanine to tyrosine in a tetrahydrobiopterin (BH₄)–dependent reaction (Scriver, 2007; Elhawary et al., 2022). Classical phenylketonuria (PKU) and related forms of hyperphenylalaninemia (HPA) arise when PAH activity is severely reduced: blood phenylalanine rises far above the normal range, tyrosine is relatively low, and neurotoxic metabolites accumulate, leading to the characteristic untreated PKU phenotype of very high phenylalanine, low tyrosine, and progressive neurological damage (Elhawary et al., 2022). Large clinical and molecular series show that the vast majority of HPA/PKU cases are caused by pathogenic variants in PAH itself, while only a minority are due to defects in BH₄ synthesis or recycling or in the PAH co-chaperone DNAJC12 (Blau et al., 2014; Himmelreich et al., 2021; Elhawary et al., 2022). In other words, within the broader “phenylalanine metabolism” pathway, PAH is the critical flux-controlling step: common and rare variation at PAH has large, direct effects on systemic phenylalanine levels, whereas most other pathway components (transporters, minor side-enzymes, cofactor recycling genes) exert much weaker or rarer effects at the population level. This is exactly the sparse driver architecture that CATFISH’s driver-oriented tests (ACAT, minP-like behavior, hard truncation) are designed to detect.

**Best detectors / diagnostics in CATFISH:**
- **minP / Tippett** (will happily call these significant)
- **ACAT** (also loves a single extreme gene)
- **SGP check:** “recompute without top gene” as a *diagnostic*, not as another omnibus test


---

### Bias warning
When calculating over-representation or enrichment, pathway-based analysis is subject to multiple biases; gene size, pathway size, density of SNP coverage, and linkage disequilibrium (LD) patterns all need to be handled explicitly (White et al., 2020; PMC6391732). In CATFISH, we address these in three steps: (i) SNP → gene is performed with MAGMA’s LD-aware multi–SNP model, so gene-level Z/P already account for local LD structure and SNP density; (ii) we then regress MAGMA gene Z-scores on log(gene length) and log(number of SNPs) and use the residual-based `P_adj` for all downstream gene → pathway tests, thereby removing remaining dependence on gene size and SNP density; and (iii) at the pathway level, we avoid naive over-representation tests and instead calibrate our omnibus statistics with gene-label permutations that preserve each gene’s adjusted p-value and the observed distribution of pathway sizes, providing enrichment p-values that are robust to these known biases.

Link: https://pmc.ncbi.nlm.nih.gov/articles/PMC6391732/

---

## METHODS

### Notation
## Gene-level covariate adjustment (gene length & SNP density)

Let a pathway (gene set) be denoted by $S$, containing $G = |S|$ genes indexed by $g = 1,\dots,G$.

For each gene $g \in S$, let:
- $Z_g$ denote the raw gene-level $Z$-statistic, and
- $p_g$ denote the corresponding two-sided $p$-value.

We collect these into vectors:

$$
\mathbf{Z}_S = (Z_1, Z_2, \dots, Z_G),
\qquad
\mathbf{p}_S = (p_1, p_2, \dots, p_G).
$$

To remove residual dependence on gene length and SNP density, we fit a linear model:

$$
Z_g = \beta_0 + \beta_1 \log(L_g) + \beta_2 \log(\mathrm{SNPs}_g) + \varepsilon_g.
$$

We define the adjusted $Z$-scores as the regression residuals:

$$
Z_{g,\mathrm{adj}} = \hat{\varepsilon}_g.
$$

The corresponding adjusted two-sided $p$-values are:

$$
p_{g,\mathrm{adj}} = 2\,\Phi\!\left(-\left|Z_{g,\mathrm{adj}}\right|\right),
$$

where $\Phi(\cdot)$ is the standard normal CDF.

We denote the adjusted vectors by:

$$
\mathbf{Z}_{S,\mathrm{adj}} = (Z_{1,\mathrm{adj}}, \dots, Z_{G,\mathrm{adj}}),
\qquad
\mathbf{p}_{S,\mathrm{adj}} = (p_{1,\mathrm{adj}}, \dots, p_{G,\mathrm{adj}}).
$$


---

## 1) Gene-level association statistics (SNP → gene)

For each gene $g$, MAGMA produces a gene‑level association p‑value $p_g$ by aggregating SNP‑level signals within/near the gene while accounting for local LD using a reference panel. We used MAGMA’s `multi=snp-wise` model, which combines a SNP-wise mean test (good for many small effects) and a SNP-wise top test (good for a single strong SNP) into a single LD-aware omnibus statistic per gene, making the gene p-values robust to different within-gene causal architectures.

The workflow includes:

- SNPs are mapped to genes (gene boundaries with optional windows).
- A multi‑marker gene model accounts for LD among SNPs in the gene region.
- MAGMA outputs gene statistics (e.g., $Z_g$ and $p_g$).

In CATFISH, these MAGMA gene p-values (or Z statistics) are treated as the inputs to all pathway tests.

---

## 2) Gene-level adjustment for gene size and SNP density

Even with LD-aware gene testing, gene‑level signals can exhibit residual dependence on gene size and SNP density. CATFISH performs a post‑hoc adjustment at the gene level.

Let:

- $Z_g$ be the MAGMA gene Z‑statistic,
- $L_g$ be gene length (bp),
- $S_g$ be number of SNPs mapped to the gene (e.g., $$NSNPS$$).

We fit a regression line as:

$$
Z_g = \beta_0 + \beta_1 \log(L_g) + \beta_2 \log(S_g) + \varepsilon_g.
$$

We define adjusted residual Z as:

$$
Z^{\mathrm{adj}}_g = Z_g - \widehat{Z}_g,
\quad \widehat{Z}_g = \widehat{\beta}_0 + \widehat{\beta}_1 \log(L_g) + \widehat{\beta}_2 \log(S_g).
$$

We then convert to a two‑sided adjusted p-value:

$$
p^{\mathrm{adj}}_g = 2\Phi\left(-|Z^{\mathrm{adj}}_g|\right),
$$

where $\Phi(\cdot)$ is the standard normal CDF.

---

## 3) Pathway-level test statistics (gene → pathway)

CATFISH computes multiple pathway statistics from $\{p_g\}_{g \in S}$ $\{Z_g\}_{Z \in S}$.

### 3.1 ACAT (Aggregated Cauchy Association Test)

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

### 3.2 Fisher’s method

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

### 3.3 Soft TFisher (tail-focused truncation and weighting)

To focus power on the lower tail of gene-level p-values while avoiding a hard cutoff, we use the **soft-thresholded TFisher** statistic of Zhang et al. (2020).

Fix a truncation / weighting parameter $\tau \in (0,1]$. For a pathway $S$ with adjusted gene-level p-values $\{p_g\}_{g \in S}$, the soft TFisher statistic is:

$$
T_{\mathrm{TF}}^{\mathrm{soft}}(S;\tau) = \sum_{g \in S} \left[ -2\log(p_g) + 2\log(\tau) \right]_{+}.
$$

Here $(x)_{+} = \max(x,0)$.

Equivalently, only genes with $p_g \le \tau$ contribute to the statistic, but their contribution is smoothly down-weighted as $p_g$ approaches $\tau$: very small $p_g$ behave like Fisher's $-2\log(p_g)$, $p_g$ near $\tau$ contribute little, and $p_g > \tau$ contribute exactly zero.

This corresponds to the $\tau_1 = \tau_2 = \tau$ special case of the general TFisher family (Zhang et al., 2020, Eq. 2.1–2.2).

Soft TFisher interpolates between:

- **Fisher's method** when $\tau = 1$, and  
- **hard truncated tests** (e.g. TPM / RTP) when $\tau < 1$, but with continuous, stable weighting instead of a sharp cutoff.

In CATFISH, we apply:

$$
T_{\mathrm{TF}}^{\mathrm{soft}}(S;\tau)
$$

to adjusted gene-level p-values $p_{g,\mathrm{adj}}$, using the analytic null distribution for TFisher with $\tau_1 = \tau_2 = \tau$ (Zhang et al., 2020, Theorem 1).

**Key property (interpretation):** Soft TFisher is especially sensitive to **Hybrid Driver–Support (HDS)** and **Coordinated Moderate Enrichment (CME)** architectures.


---

### 3.4 Stouffer's method (mean-Z; diffuse polygenic shift)

Stouffer aggregates gene-level **Z** statistics (e.g., $$Z_{g,\mathrm{adj}}$$) rather than p-values.  
For a pathway $$S$$ with $$G=|S|$$ genes, define:

**Unweighted Stouffer:**

$$
Z_{\mathrm{stouffer}}(S) = \frac{1}{\sqrt{G}} \sum_{g \in S} Z_g
$$

$$
p_{\mathrm{stouffer}}(S) = 2\,\Phi\!\left(-\left|Z_{\mathrm{stouffer}}(S)\right|\right)
$$

where $$\Phi(\cdot)$$ is the standard normal CDF.

**Weighted Stouffer (optional):** choose weights $$w_g \geq 0$$.

$$
Z_{\mathrm{stouffer}}^{(w)}(S) = \frac{\sum_{g \in S} w_g Z_g}{\sqrt{\sum_{g \in S} w_g^2}}
$$

$$
p_{\mathrm{stouffer}}^{(w)}(S) = 2\,\Phi\!\left(-\left|Z_{\mathrm{stouffer}}^{(w)}(S)\right|\right)
$$

**Notes:**
- Best for **Diffuse Polygenic Shift (DPS)**: many small positive (or negative) shifts in $$Z_g$$ that may not yield many $$p_g<0.05$$ individually.
- If calibrating by permutation, recompute $$Z_{\mathrm{stouffer}}$$ (or $$Z_{\mathrm{stouffer}}^{(w)}$$) under permuted gene labels and use an empirical p-value.

---
### 3.5 minP / Tippett (single-gene proxy diagnostic)

Define the pathway's minimum gene p-value:

$$
p_{\min} = \min_{g \in S} p_g
$$

Under independence, Tippett's combined p-value is:

$$
p_{\mathrm{tippett}} = 1 - (1 - p_{\min})^G
$$

But because gene p-values are correlated (LD, shared biology), CATFISH typically treats **minP as a diagnostic** and/or calibrates it by permutation:

Let:
- $$p_{\min}^{\mathrm{obs}}$$ be the observed minimum p-value
- $$p_{\min}^{(b)}$$ be the minimum p-value in permutation $$b$$
- $$B$$ be the total number of permutations

Then the permutation p-value is:

$$
p_{\min}^{\mathrm{perm}} = \frac{1 + \text{count}(p_{\min}^{(b)} \leq p_{\min}^{\mathrm{obs}})}{B + 1}
$$

where $$\text{count}(\cdot)$$ counts the number of permutations satisfying the condition.

**Notes:**
- Highly sensitive to **single extreme genes**, so it flags **Sparse Driver (SDA)** and especially **Single-Gene Proxy (SGP)** pathways.
- Use "remove top gene and recompute" as the practical SGP check.
---

## 4) Correlation among pathway tests

All pathway statistics above are functions of the same gene-level p-values $$\{p_g\}$$ and are therefore intrinsically correlated. For a given pathway $$S$$, define $$T_{\mathrm{ACAT}},\$$ $$T_{\mathrm{Fisher}},\$$ $$T_{\mathrm{TF}}^{\mathrm{soft}}(\tau),\$$ $$T_{\mathrm{Stouffer}},\$$ $$T_{\mathrm{minP}},\$$
where

- $$T_{\mathrm{ACAT}}$$ is the ACAT statistic on $$S$$ (Section 3.1),
- $$T_{\mathrm{Fisher}}$$ is the Fisher / sum–log–p statistic (Section 3.2),
- $$T_{\mathrm{TF}}^{\mathrm{soft}}(\tau)$$ is the soft–TFisher statistic with truncation parameter $$\tau$$ (Section 3.3),
- $$T_{\mathrm{Stouffer}}$$ is the (possibly weighted) Stouffer / mean-$$Z$$ statistic, and
- $$T_{\min} = \min_{g \in S} p_g$$ is the Tippett / minP statistic.

Under the null, this vector is **not** jointly independent, because each component is a deterministic function of the same set of p-values $$\{p_g\}$$ (or equivalently the same ordered p-values $$\{p_{(k)}\}$$).

Intuitively:

- **Fisher** and **soft TFisher** are both based on $$-\log(p_g)$$, with TFisher acting like a truncated / reweighted Fisher focused on $$p_g \le \tau$$. Whenever the sum of $$-\log(p_g)$$ is large (many small p-values), both statistics tend to be extreme in the same direction.

- **Stouffer** first maps p-values to $$Z$$-scores via $$Z_g = \Phi^{-1}(1 - p_g/2)$$ and then averages or sums them. This is a monotone transform of the same $$\{p_g\}$$, so pathways with many small p-values will simultaneously inflate Fisher, TFisher, and Stouffer.

- **ACAT** and **soft TFisher** both emphasize the **lower tail**: ACAT uses $$\tan\{\pi(0.5 - p_g)\}$$ (very heavy weight on the smallest p’s), while TFisher heavily weights small $$p_g$$ through $$-\log(p_g)$$ below $$\tau$$. A pathway with one or a few very small p-values will push both $$T_{\mathrm{ACAT}}$$ and $$T_{\mathrm{TF}}^{\mathrm{soft}}(\tau)$$ in the same direction.

- **minP** is literally the smallest p-value $$p_{(1)}$$. It is almost perfectly aligned with the extreme-tail behavior of ACAT and the most aggressive truncation versions of TFisher: whenever $$p_{(1)}$$ is tiny, ACAT, TFisher, and often Fisher / Stouffer will also look unusually extreme.

On top of this, linkage disequilibrium and shared biology induce correlation among the gene-level inputs themselves (the $$\{p_g\}$$ and corresponding $$\{Z_g\}$$ are not independent across genes), which further strengthens dependence among all of the pathway-level statistics.

**Implication.** Because $$T_{\mathrm{ACAT}},\$$ $$T_{\mathrm{Fisher}},\$$ $$T_{\mathrm{TF}}^{\mathrm{soft}}(\tau),\$$ $$T_{\mathrm{Stouffer}},\$$ $$T_{\mathrm{minP}},\$$
is strongly dependent, naïve combination rules that assume independence between component tests (e.g. simple Bonferroni or analytic minP under independence) would generally be mis-calibrated and often anti-conservative. This is why CATFISH uses **permutation-calibrated minP across methods** (and treats ACAT-O as a model-based sensitivity analysis), so that the joint null distribution of these correlated statistics is learned empirically rather than assumed.


---

## 5) Omnibus pathway p-value across methods

Each pathway in CATFISH is evaluated by a panel of complementary gene → pathway tests (ACAT, Fisher, soft TFisher, Stouffer, and minP), each tuned to a different pattern of gene-level signal. Because these five statistics are all functions of the same adjusted gene p-values, they are strongly correlated and can disagree in which pathways they call “most” enriched. Rather than choosing a single favorite test a priori, we combine their evidence into a **single omnibus pathway p-value** using two complementary strategies: (i) an analytic ACAT-O combination across method-level p-values, and (ii) a permutation-calibrated minP across methods that provides a dependence-robust, “best-of-tests” summary. The following subsections describe these two omnibus layers.


### 5.1 Omnibus ACAT (ACAT-O across methods)

Let $$p_1, \dots, p_5$$ denote the five component p-values for pathway $$S$$, and let weights $$v_j \ge 0$$ satisfy $$\sum_{j=1}^5 v_j = 1$$ (default $$v_j = 1/5$$).

Define the ACAT-O Cauchy statistic

$$T_{\mathrm{omni,ACAT}}(S) = \sum_{j=1}^5 v_j \,\tan\{\pi(0.5 - p_j)\},$$
which under the global null is approximately standard Cauchy. The analytic omnibus p-value is
$$p_{\mathrm{omni,ACAT}}(S)= 0.5 - \frac{1}{\pi}\arctan\{T_{\mathrm{omni,ACAT}}(S)\}.$$

This ACAT-O layer is most sensitive when **at least one** component test (e.g. ACAT for sparse drivers, Fisher/Stouffer for coordinated enrichment, soft TFisher for hybrid patterns, or minP for hard single-gene hits) is strongly significant, even if the others are only modest or null.

---

### 5.2 Omnibus minP across methods with permutation 

To obtain a complementary, more conservative summary that explicitly accounts for correlation between component tests, we also compute a minimum-p omnibus statistic across methods: 

$$T_{\mathrm{omni,min}}(S) = \min_{j \in \{1,\dots,5\}} p_j= \min\({ p_{\mathrm{ACAT}}(S), p_{\mathrm{Fisher}}(S),p_{\mathrm{TF}}(S;\tau), p_{\mathrm{Stouffer}}(S),p_{\mathrm{minP}}(S) \})$$ 

Because the $$p_j$$ are correlated (they are all computed from the same set of gene-level p-values), we calibrate $$T_{\mathrm{omni,min}}(S)$$ by gene-label permutation. 

For each of $$B$$ permutations: 

- randomize gene labels across pathways (preserving the empirical distribution of gene-level p-values and pathway sizes),
- recompute all five component tests for each pathway,
- record, for permutation $$b = 1,\dots,B$$, $$T_{\mathrm{omni,min}}^{(b)}(S) = \min_j p_j^{(b)}(S)$$

The permutation-based omnibus p-value is then 

$$\hat p_{\mathrm{omni,min}}(S)= \frac{1 + \text{count}\{\,b : T_{\mathrm{omni,min}}^{(b)}(S)\le T_{\mathrm{omni,min}}(S)\,\}}{B + 1}$$ 

We use $$\hat p_{\mathrm{omni,min}}(S)$$ as the **primary** omnibus p-value in the main analyses because it 

1. protects nominal type-I error under *arbitrary* dependence between component tests, and
2. explicitly targets the "best" component test for each pathway while accounting for the fact that this best test is effectively chosen post hoc via the min across methods.

The analytic ACAT-O p-value $$p_{\mathrm{omni,ACAT}}(S)$$ is reported alongside as a **higher-power, model-based sensitivity analysis**, highlighting pathways that are consistently strong across methods or dominated by a single very informative test.

---

## 6) Multiple testing correction

Across all pathways, omnibus p-values $\{p_{\mathrm{omni}}(S)\}$ are adjusted using Benjamini–Hochberg FDR:

$$
q_{\mathrm{BH}}(S) = \mathrm{BH}\left(p_{\mathrm{omni}}(S)\right).
$$

Because each pathway yields a **single** omnibus p-value, no additional penalty is required for the number of component tests.

---

## 7) Omnibus across methods with permutation (LD-aware upstream gene statistics)

To obtain a complementary, conservative omnibus summary that explicitly accounts for dependence among component pathway tests, we compute a minimum-p statistic across the five component p-values for each pathway $$S$$.
Let

$$\mathcal{P}(S)=\{p_{\mathrm{ACAT}}(S)\,p_{\mathrm{wFisher}}(S)\,p_{\mathrm{TPM}}(S;\tau)\,p_{\mathrm{Stouffer}}(S)\,p_{\mathrm{minP,gene}}(S)\}$$

The omnibus minimum statistic is

$$T_{\mathrm{omni,min}}(S)=\min \mathcal{P}(S)=\min\!\{p_{\mathrm{ACAT}}(S),\,p_{\mathrm{wFisher}}(S),\,p_{\mathrm{TPM}}(S;\tau),\,p_{\mathrm{Stouffer}}(S),\,p_{\mathrm{minP,gene}}(S)\}.$$

Because all $$p_j(S)$$ are computed from the *same* gene-level association evidence and therefore are correlated, $$T_{\mathrm{omni,min}}(S)$$ is not Uniform $$(0,1)$$ under the null. We therefore calibrate the minP omnibus by permutation.

#### 7.1 LD-aware gene-level inputs via MAGMA

All pathway tests operate on gene-level summary statistics computed using MAGMA's SNP-to-gene model with an LD reference panel. Specifically, GWAS summary statistics are mapped to genes (with a symmetric $$\pm 25$$ kb window), and MAGMA is run using the multi-model SNP-wise gene analysis. This step is LD-aware: within each gene, SNP associations are aggregated while accounting for linkage disequilibrium in the reference genotypes, producing gene-level test statistics (e.g., $$Z_g$$, $$P_g$$) that are appropriately calibrated under correlated SNP structure.

To reduce confounding by gene size and SNP density, we additionally adjust gene-level $$Z_g$$ (or $$-\log_{10} P_g$$) by regressing on log(gene length) and log(#SNPs) and using residual-based, size-adjusted gene p-values $$P_g^{\mathrm{adj}}$$. These adjusted gene p-values are the inputs to all pathway-level component tests.

This design ensures that LD is handled where it matters most and is most identifiable from data—at the SNP-to-gene aggregation stage—while keeping downstream pathway aggregation flexible across multiple alternative signal architectures.

#### 7.2 Component pathway tests

For each pathway $$S$$ with member genes $$g \in S$$, we compute five complementary pathway p-values from the (optionally adjusted) gene-level p-values $$\{P_g\}_{g\in S}$$:

1. **ACAT**: sensitive to sparse signals and a few very strong genes.
2. **Weighted Fisher (wFisher)**: aggregates evidence across genes with optional weights (e.g., SNP counts) and can incorporate effect direction when available.
3. **Truncated Product Method / truncated Fisher (TPM)** with truncation $$\tau$$: emphasizes moderate subsets of more significant genes while down-weighting null genes.
4. **Stouffer**: combines gene-level Z-scores (or Z reconstructed from p-values under the null) to target diffuse polygenic enrichment.
5. **Gene-level minP**: a diagnostic "single-gene proxy" detector, capturing whether the pathway is driven by one extreme gene.

These tests are deliberately correlated but have different power profiles across latent pathway architectures. The minP omnibus selects the best-performing component while controlling for this post hoc selection.

#### 7.3 Gene-label permutation calibration of the minP omnibus

We calibrate $$T_{\mathrm{omni,min}}(S)$$ using gene-label permutation, which preserves (i) the empirical distribution of gene-level p-values and (ii) the observed pathway sizes, while breaking pathway membership.

For each permutation $$b = 1,\dots,B$$ and each pathway $$S$$ with $$|S|$$ genes:

1. **Sample a null gene set** $$S^{(b)}$$ by selecting $$|S|$$ genes uniformly without replacement from the pool of all genes with non-missing gene-level p-values.
2. **Recompute all component tests** on the sampled gene set, yielding

   $$p_j^{(b)}(S), \quad j \in \{\mathrm{ACAT},\mathrm{wFisher},\mathrm{TPM},\mathrm{Stouffer},\mathrm{minP,gene}\}$$

3. Record the permutation min statistic

   $$T_{\mathrm{omni,min}}^{(b)}(S) = \min_j p_j^{(b)}(S)$$

The permutation-calibrated omnibus p-value is then

$$\hat p_{\mathrm{omni,min}}(S)=\frac{1 + \bigl|\{\,b : T_{\mathrm{omni,min}}^{(b)}(S)\le T_{\mathrm{omni,min}}(S)\,\}\bigr|}{B + 1}$$

In practice, this procedure is implemented in MAGCAT/CATFISH as follows:

- we construct the gene pool from MAGMA gene results (after optional size/SNP adjustment),
- for each pathway we generate $$B$$ random gene sets of identical size,
- for each permuted gene set we recompute the five component pathway tests using the same parameters (e.g., truncation $$\tau$$ for TPM),
- we compute $$\hat p_{\mathrm{omni,min}}(S)$$ using the standard "+1" correction to avoid zero p-values.

This permutation directly addresses correlation between the component tests because the full multivariate dependence induced by shared inputs and shared pathways is mirrored under permutation.

#### 7.4 Rationale and interpretation

We use $$\hat p_{\mathrm{omni,min}}(S)$$ as the **primary** omnibus p-value in the main analyses because it:

1. controls type-I error under *arbitrary dependence* among component tests (all derived from the same gene-level signals),
2. properly accounts for "winner's curse" induced by taking a minimum across multiple methods,
3. remains robust across diverse pathway signal architectures by adaptively selecting the most informative component per pathway.

The analytic ACAT-O omnibus $$p_{\mathrm{omni,ACAT}}(S)$$ is reported alongside as a **higher-power, model-based sensitivity analysis**, highlighting pathways with consistent support across methods or those dominated by a single exceptionally sensitive component.

#### 7.5 Note on LD-awareness and the permutation layer

A key design choice is that LD-awareness is enforced upstream at the SNP-to-gene stage via MAGMA's LD-informed gene model. The omnibus permutation described above is a gene-label permutation and therefore does not explicitly resimulate LD structure. Instead, it leverages the fact that the gene-level p-values already represent LD-adjusted gene evidence; the permutation step is used to calibrate dependence *between pathway-level test statistics* and to control for post hoc selection across multiple correlated pathway tests.

This tiered strategy allows computationally feasible, robust omnibus inference using only GWAS summary statistics and an LD reference panel, while retaining sensitivity to multiple plausible pathway architectures.


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

## Quick R example (CATFISH)

> **Note:** This is the exact end-to-end pipeline used (MAGMA → gene adjustment → pathway tests → omnibus). Paths, filenames, and column mappings should be edited to match your local files.

```r
############################################################
## CATFISH PIPELINE: MAGMA → SIZE-ADJUSTED GENE P → PATHWAYS
## (with brief explanation of the key parameters)
############################################################

## CATFISH gives you:
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

## - You can call species = "maize" in CATFISH , it automatically
##   uses inst/extdata/maize.genes.loc 

## - You can also use any gff3 file for your organism of choice
## - Parses the GFF3, extracts gene features,
## and writes a MAGMA-ready-to-use gene location file (GENE, CHR, START, STOP, STRAND…).


############################################################
## 2. SNP → gene annotation with MAGMA (magma_annotate wrapper)
############################################################

stats_file <- "/Users/.../raw_GWAS_MLM_3PC_N.txt"

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

## - magma_annotate() builds the MAGMA command line and calls the MAGMA binary for you.
## - You just supply stats_file and remember to rename your columns or do it here


############################################################
## 3. Gene-level MAGMA (multi = snp-wise) with R wrapper
############################################################

# Plink based bed/bim.fam files
bfile      <- "/Users/.../all_maize2"   # PLINK basename: .bed/.bim/.fam

## 3A. Chromosome-wise run using NMISS (per-SNP sample size)
## NMISS is the number of missing genotypes
magma_gene(
  bfile      = bfile,
  gene_annot = "annot/N_maize_MLM.genes.annot",
  stats_file = stats_file,
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
  gene_model = c("multi=snp-wise"),  # multi-parameter SNP-wise model. This works best for CATFISH. Combines top and mean SNPs. 
  chroms     = 1:10,                 # run MAGMA separately for chr 1–10
  n_threads  = 10                    # run up to 10 MAGMA jobs in parallel
)

## 3B. Chromosome-wise run using NOBS (per-SNP N)
## NOBS is the total number of genotypes used
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

## - magma_gene() is a R wrapper around the MAGMA binary.
##   * Handles all flags, temp files, and error checking.
##   * Additional can parallelize the analysis
##   * chroms + n_threads to run multiple chromosomes in parallel.
##   * Only requires a minimal rename_columns spec instead of reformatting


############################################################
## 4. Combine per-chromosome MAGMA gene outputs
############################################################

# Load the file
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

## - This combines 10 per-chromosome MAGMA outputs into a single file for ease of use


############################################################
## 5. Gene length extraction + SNP density bias adjustment
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


## 5B. Adjust gene-level Z/P for gene length & #SNPs

genes_all <- read.table(
  "/Users/.../magma_N_maize.txt", # from previous step
  header = TRUE, stringsAsFactors = FALSE
)

# Adjsut the pvalue based on gene length and number of snps
adj_out <- magcat_adjust_gene_p(
  gene_results = genes_all,
  gene_lengths = maize_gene_len,
  gene_col     = "GENE",
  p_col        = "P",
  z_col        = "ZSTAT",    # raw MAGMA Z
  len_gene_col = "gene_id",
  len_col      = "length"
  # nsnp_col   = "NSNPS"     # #SNPs per gene
)

genes_adj <- adj_out$genes    # includes Z_raw, Z_adj, P_adj, log_gene_length, log_nsnp
lm_fit    <- adj_out$fit      # lm(Z_raw ~ log_gene_length + log_nsnp)

write.csv(genes_adj, "genes_adj.csv", row.names = FALSE)

## - MAGMA gene p-values are known to correlate weakly with gene size and SNP density.
## - magcat_adjust_gene_p():
##   * Fits a linear model: Z_raw ~ log(gene length) + log(#SNPs).
##   * Uses residuals (Z_adj) as “size/SNP-adjusted” Z-scores.
##   * Converts Z_adj back to P_adj = 2*pnorm(-|Z_adj|).
## - You then use P_adj in gene→pathway tests, reducing bias toward big genes.


############################################################
## 6. Load pathway definitions from PMN/CornCyc or use saved files
############################################################

maize_pw <- magcat_load_pathways(
  species  = "maize",
  gene_col = "Gene-name"  # column in the PMN gene-set file that matches MAGMA gene IDs
)

############################################################
## 7. Omnibus combining methods (ACAT or minP)
############################################################

## All tests use gene_results + species/pathways to:
##  - find genes per pathway,
##  - take their p-values (raw P or adjusted P_adj),
##  - compute a pathway-level p per method.


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
  omnibus           = "minP",      # "minP" or "ACAT"
  B_perm            = 10000L,      # permutations at the *omnibus* level
  seed              = 123,
  remove_singletons = TRUE,        # drop pathways with n_genes < 2
  output            = TRUE,
  out_dir           = "magcat_omni_full"
)

## In a single call omni_pathways gives you:
##   - acat_p       : ACAT pvalue per pathway
##   - fisher_p    : Fisher pvalue per pathway
##   - tpm_p        : truncated Fisher (soft) pvalue per pathway
##   - stouffer_p   : Stouffer pvalue per pathway
##   - minp_gene_p  : minP pvalue per pathway
##   - omni_p       : combination of these methods (ACAT or minP) pvalue per pathway
##   - omni_perm_p  : permutation-calibrated omnibus pvalues
##   - BH FDR and q-values for omni_p and each component
##   - omni_perm_BH  : permutation-calibrated BH pvalues
##   - omni_perm_q  : permutation-calibrated q pvalues

############################################################
## 8. Pathway-level tests (gene → pathway)
############################################################

## You can run each pathway individually as well

### 8A. ACAT per pathway

pw_res_acat_adj <- magcat_acat_pathways(
  gene_results = genes_adj,     # adjusted Pvalue
  species      = "maize",
  gene_col     = "GENE",
  p_col        = "P_adj",
  B            = 10000L,            # 10,000 is good enough
  seed         = NULL,
  output       = TRUE,
  out_dir      = "acat_results"
)

### 8B. Fisher pathways

wf_res_raw <- magcat_wfisher_pathways(
  gene_results = genes_adj,   # raw P
  species      = "maize",
  gene_col     = "GENE",
  p_col        = "P_adj",
  effect_col   = "ZSTAT",
  #weight_col   = NULL, # If you have any other weights. 
  is_onetail   = FALSE
)

### 8C. Truncated Fisher (soft) / TFisher(soft)

soft_tf_res_adj <- magcat_soft_tfisher_pathways(
  gene_results     = genes_adj,
  species          = "maize",
  gene_col         = "GENE",
  p_col            = "P_adj",
  tau1             = 0.05,     # soft truncation threshold
  B_perm           = 10000L,  # permutations for empirical pvalue
  seed             = 123,
  analytic_logical = TRUE,
  output           = TRUE,
  out_dir          = "magcat_tfisher_soft"
)

### 8D. Stouffer (sum of Z across genes)

stouf_res <- magcat_stouffer_pathways(
  gene_results = genes_adj,
  species      = "maize",
  gene_col     = "GENE",
  p_col        = "P_adj",
  weight_col   = NULL,    # equal weights for all genes
  B_perm       = 10000L,    # permutations for empirical pvalue
  seed         = 123,
  output       = TRUE,
  out_dir      = "magcat_stouffer"
)

### 8E. Gene-level minP per pathway

minp_res <- magcat_minp_pathways(
  gene_results = genes_adj,
  species      = "maize",
  gene_col     = "GENE",
  p_col        = "P_adj",
  B_perm       = 10000L,      # permutations for empirical pvalue
  min_p        = 1e-15,
  do_fix       = TRUE,
  output       = TRUE,
  out_dir      = "magcat_minp_maize"
)

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
  


