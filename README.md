# CATFISH : GWAS pathway analysis pipeline

This Markdown is structured into:

- **INTRODUCTION** — what the pipeline is and why multiple tests are needed  
- **METHODS** — paper‑ready subsections with explicit equations and assumptions  
- **USAGE** — installation + reproducible example commands and R snippets
- **RESULTS** — an example analysis for soil Nitrogen GWAS

---

## INTRODUCTION

### What CATFISH is

**CATFISH** (**C**ombining **A**CAT, **T**Fisher (soft), **F**isher, m**I**n-P, and **S**touffer for **H**olistic pathway analysis) is a model-averaged pathway framework that combines ACAT, soft TFisher, Fisher, Stouffer, and minP on top of LD-aware MAGMA gene-level GWAS statistics adjusted for gene length and SNP density. It then uses an omnibus test based on permutation-calibrated minP or ACAT to collapse these multiple pathway tests into a single, correlation-robust enrichment p-value that is sensitive to both sparse and polygenic pathway patterns. In short, CATFISH casts a wide net across complementary tests and reels in a single pathway p-value.

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
  


