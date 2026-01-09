# CATFISH : GWAS pathway analysis pipeline

This Markdown is structured into:

- [**INTRODUCTION**](#introduction) — what the pipeline is and why multiple tests are needed  
- [**METHODS**](#methods) — paper‑ready subsections with explicit equations and assumptions  
- [**USAGE**](#usage) — installation + reproducible example commands and R snippets
- [**RESULTS**](#results) — an example analysis for soil Nitrogen GWAS

---

## INTRODUCTION

**CATFISH** (Combining <ins>**C**</ins>auchy combination (ACAT), <ins>**A**</ins>daptive TFisher (soft) test, <ins>**F**</ins>isher's test, m<ins>**I**</ins>n-P, and <ins>**S**</ins>touffer's method for <ins>**H**</ins>olistic pathway analysis) is a multi-test pathway framework built on LD-aware MAGMA gene-level GWAS statistics adjusted for gene length and SNP density that combines ACAT, soft TFisher, Fisher, Stouffer and minP. It then uses an omnibus test based on permutation-calibrated minP or ACAT to collapse these multiple pathway tests into a single, correlation-robust enrichment p-value that is sensitive to both sparse and polygenic pathway patterns. In short, CATFISH casts a wide net across complementary tests and reels in a single pathway p-value.

CATFISH uses:

1. **MAGMA** for LD-aware **SNP → gene** inference (gene-level p-values).
2. **Multiple gene → pathway combination tests** (ACAT, Fisher, soft TFisher, Stouffer, minP).
3. A **correlation-robust LD-aware omnibus test** (permutation-calibrated second minP and ACAT test) that aggregates these tests into a single pathway-level p-value.

### Why multiple tests are needed

Pathways can be significant for statistically different reasons:

- one driver gene dominates,
- many genes show coordinated moderate enrichment,
- diffuse polygenic shift, and
- hybrids of these patterns.

No single gene-set statistic is uniformly most powerful across these various possibilities. Instead of relying on a single test, CATFISH runs several complementary pathway tests and combines them into one pathway-level omnibus p-value. We define these patterns using a set of **pathway signal archetypes** (sparse driver, coordinated moderate enrichment, diffuse polygenic shift, hybrid driver–support, single-gene proxy), which describe different ways a pathway can appear significant.

# Pathway signal archetypes

We divide the pathway signals into a set of archetypes that describe different ways a pathway can be enriched in a biological sense explained in detail below. We then use a combination of statistical tests chosen to be representative to each of these behaviors and a provide a biological example for each archetype.

## Archetype I — Sparse Driver Architecture (SDA)

![Archetype I — Sparse Driver Architecture](Figures/Fig_archetypesI.png)
*Archetype I: Sparse Driver Architecture (SDA).*

**Signature:** A small number of genes are extremely significant; most genes look null.

**Gene-level p-value pattern:**

- There exists a small $K \ll G$ such that the top $K$ gene p-values are extremely small, $p_{(1)}, \dots, p_{(K)} \ll \alpha$ (e.g. $p_{(1)} \sim 10^{-7}$ or smaller),
- The remaining genes are approximately null, $p_{(K+1)}, \dots, p_{(G)} \sim \mathrm{Uniform}(0,1)$.
- This produces a sharp “elbow” in the ranked p-values (a few tiny hits followed by a long flat tail).

**Interpretation:**  
An SDA pathway is significant because **a small set of driver genes dominates the signal**, rather than broad involvement of most pathway members. This can occur when the trait-relevant biology passes through a **bottleneck** (committed step, rate-limiting enzyme, key regulator, or essential transporter) so that genetic variation concentrates its effect at a few control points. In contrast, pathway annotations typically include many additional enzymes, modifiers, and general “support” genes that may be necessary for pathway operation but do not carry strong association for the trait. Under SDA, association is therefore concentrated in the top $K$ genes, yielding very small ordered p-values $p_{(1)}, \ldots, p_{(K)}$ followed by a long tail $p_{(K+1)}, \ldots, p_{(G)}$ that is close to uniform.

**Biological example:**  
**Aspartokinase in the aspartate-derived amino acid pathway**

The aspartate-derived amino-acid biosynthesis pathway converts **aspartate** into essential amino acids, such as **lysine**, **threonine**, **methionine**, and **isoleucine**. In plants and bacteria, the first step is catalyzed by **aspartokinase (AK)**, which phosphorylates aspartate to **aspartyl-phosphate** that feeds multiple branched network to produces several end-products. Due to its position at the starting point of the pathway, AK acts as a **flux-controlling bottleneck**. Variations in AK activity alter carbon and nitrogen flow through the whole network. Downstream enzymes, including tailoring steps, dehydrogenases, and transaminases, typically exhibit dispersed or buffered roles. Thus, the gene-level pattern aligns with SDA. A single AK gene or a few genes downstream may exhibit very low $p$ values (e.g. $10^{-8}\text{-}10^{-10}$), while the majority are dispersed across the interval $(0,1)$.

**Best detectors in CATFISH:**
- **ACAT** — sensitive to a few extremely small p-values (RECOMMENDED).  
- **minP / Tippett** — targets the minimum p-value; optimal when one gene dominates.


---

## Archetype II — Coordinated Moderate Enrichment (CME)

![Archetype II — Coordinated Moderate Enrichment (CME)](Figures/Fig_archetypesII.png)
*Archetype II — Coordinated Moderate Enrichment (CME)*


**Signature:** many genes show moderate association; no single gene is extremely low.

**Gene-level p-value pattern:**
- A non-trivial fraction of genes have moderately small p-values, e.g.
  $$p_i \in [10^{-3}\,0.05]\quad\text{for many} i.$$
- The top signal is not orders-of-magnitude beyond the rest (no single-gene spike), e.g.
  $$p_{(1)} \not\ll p_{(k)}\quad\text{for small }k\ (\text{e.g., }k=5,10,20).$$
- The ranked p-values show a **broad shoulder** (many good genes) but not a sharp elbow.

**Interpretation:**  
CME signifies **collective functional engagement**. The pathway behaves like a coordinated module wherein numerous components exert small-to-moderate effects. This is expected when phenotypes emerge from distributed regulation, redundancy, and buffering/feedback as dramatic single-gene effects are rare. Statistically, enrichment arises from numerous mildly informative genes rather than single driver.

**Biological example:**  
**Cytokine / immune signaling cascades**

In numerous immunological pathways, the output is regulated not by a singular "master gene," but through distributed modulation across various tiers of a signaling circuit. A clear example is TNFα / IL-1β → NF-κB, where the activation of upstream receptors ultimately activates IKK complexes, which phosphorylate IκB inhibitors, facilitating the nuclear translocation of NF-κB family members. The temporal dynamics of NF-κB activation (rapid/transient versus slower/sustained) are significantly influenced by the stimulus class and receptor context  (Zhao et al., 2018). The significance of these dynamics lies in the variability of transcriptional outputs influenced by stimulus, NF-κB family composition, and cell type. Core feedback and marker targets encompass genes such as **NFKBIA (IκBα)** and **TNFAIP3 (A20)**, underscoring the notion that pathway behavior is governed by numerous regulatory nodes rather than a singular switch (Zhao et al., 2018).

This “many-knobs” architecture is also apparent one layer upstream in TNFRSF signaling. TNFRSF receptors bind trimeric TNFSF ligands, however, increasing evidence suggests that a single trimeric ligand–receptor complex fails to elicit complete signaling output. Rather, for certain TNFRSF pathways (notably the classical NF-κB pathway), successful activation may necessitate secondary interactions or clustering of multiple trimeric receptor complexes (Medler et al., 2019). Mechanistically, this indicates that pathway output relies on the coordinated effects of receptor assembly/avidity, adaptor recruitment, kinase activation thresholds, and the strength of negative feedback, precisely the type of system where typical genetic variation is anticipated to produce numerous modest perturbations rather than a singularly significant driver.

In CATFISH terminology, this yields a CME pattern. Within a cytokine/immune circuit, gene-level p-values may exhibit a surplus of moderate signals (e.g., numerous genes around 10⁻³–10⁻²) without a singular extreme outliers. The ordered gene-level p-values ($$p_{(1)}, p_{(2)}, \ldots$$) contains many values in the range of $$10^{-3} to 10^{-2},$$ without an extreme such as 10^{-12}. Consequently, CME pathways are optimally represented by evidence-accumulating tests (e.g., Fisher, Stouffer/mean-Z, and mild-truncation/soft-TFisher), whose efficacy is enhanced when numerous route members exhibit moderate associations, rather than depending on a singular peak.

**Best detectors in CATFISH:**
- **Fisher’s method** (aggregates evidence across many moderately small p-values) (RECOMMENDED).
- **Stouffer / mean-Z** (gains power when many genes shift together).
- Optionally **wFisher / weighted Z** if you later add biologically informed weights.


---

## Archetype III — Diffuse Polygenic Shift (DPS)

![Archetype III — Diffuse Polygenic Shift (DPS)](Figures/Fig_archetypesIII.png)
*Archetype III — Diffuse Polygenic Shift (DPS)*


**Signature:** the pathway’s genes are, on average, slightly more associated than the genome-wide background, but almost none cross a conventional significance threshold.

**Gene-level p-value / Z pattern:**
- Most gene p-values satisfy:
  $$p_i > 0.05\quad\text{for most } i.$$
Yet the pathway shows a small but consistent shift in adjusted gene-level Z-scores:

$$
\overline{Z}_{S,\mathrm{adj}} \;=\; \frac{1}{G}\sum_{i\in S} Z_{i,\mathrm{adj}} \;\neq\; 0,
$$

often with a coherent sign (bias in one direction).
- $$\{p_i: i\in S\}\ \text{is subtly enriched toward smaller values relative to Uniform}(0,1),$$
  but without extreme outliers.
- No sharp elbow, instead, the ranked p-values show a gentle, global downward bend relative to null.

**Interpretation:**  
DPS indicates a global pathway bias aligned with polygenicity. Numerous genes individually exert minimal impacts in a uniform direction, resulting in pathway enrichment due to the collective subtle shift of the entire module rather than the influence of a singular "star" gene. This is the regime characterized by numerous little issues and the absence of significant challenges. Spike-hunting tests, such as minP/Tippett or highly aggressive truncation, are generally underpowered in this context due to the absence of a singular extreme p-value to leverage. Conversely, mean-/distribution-sensitive tests (Stouffer/mean-Z, Fisher, and competitive regression) are specifically formulated to identify this subtle, pervasive divergence from the null hypothesis.

**Biological example:**  

**Biological example – human height as a diffuse polygenic shift**

The height of adult humans exemplifies a highly polygenic characteristic. Initial GIANT meta-analyses involving over 180,000 individuals identified 180 loci and "hundreds of variants" associated with height. Nevertheless, they accounted for merely 10% of the phenotypic variance, despite height exhibiting approximately 80% heritability. Subsequent meta-analyses involving about 700,000 Europeans identified thousands of linked SNPs and hundreds of locations, validating the perspective that height is regulated by several common variants of minimal effect rather than a limited number of high-impact genes. A recent investigation of a "saturated map" identified over 12,000 genome-wide significant SNPs across more than 7,000 genomic segments, encompassing around 21% of the genome, thereby reaffirming Fisher’s original polygenic hypothesis for height proposed in 1918.

These variants are classified into several growth-related processes, such as chondrocyte proliferation and hypertrophy in the growth plate, extracellular matrix and cartilage organization, growth hormone and IGF-1 signaling, and morphogen pathways including TGF-β and Hedgehog, as well as overarching developmental and endocrine regulators. Analyses of height GWAS loci through pathway and gene-set evaluations have demonstrated enrichment for signaling pathways such as TGF-β and Hedgehog, as well as for genes associated with skeletal growth, growth plate regulation, and pertinent Mendelian growth disorders. Lango Allen et al. (2010) discovered that genes adjacent to height-associated variants congregate in biologically coherent pathways, including TGF-β signaling, Hedgehog signaling, and histone and growth/development gene sets. Furthermore, several SNPs near these pathway genes "narrowly miss" genome-wide significance, suggesting numerous additional sub-threshold contributors within the same modules. Guo et al. (2018) demonstrate that genes adjacent to height GWAS loci are enriched in processes and tissues pertinent to growth, including growth plate cartilage.

When examined at the level of an individual pathway (e.g., TGF-β signaling, Hedgehog signaling, or growth-plate extracellular matrix genes), this structure inherently produces a DPS pattern. Throughout the gene set, numerous genes possess one or more prevalent variations that have minor impacts on height. Certain genes attain definitive genome-wide relevance, while many others exhibit relatively mild or nominal associations. The outcome indicates that, in contrast to random gene sets, the distribution of gene-level test statistics within these pathways exhibits a shift towards more robust evidence overall characterized by an increased number of genes with small or moderate p-values and a decreased number of genes appearing entirely null, despite the fact that most individual genes would not, on their own, substantiate a strong association claim. This scenario exemplifies the optimal application of CATFISH’s DPS-oriented detectors (Stouffer/mean-Z on adjusted gene-level Z-scores, and optionally MAGMA-style competitive regression tests). They assess whether the **average** association signal across a biologically coherent pathway is subtly yet consistently heightened in comparison to the genome-wide background.

**Best detectors in CATFISH:**
- **Stouffer / mean-Z on** `Z_adj` (unweighted; permutation-calibrated) (RECOMMENDED).
- Optionally **competitive regression-style gene-set models** (e.g., **MAGMA competitive**) if included as a component test (NOT INCLUDED IN CATFISH)


---

## Archetype IV — Hybrid Driver–Support (HDS)

![Archetype IV — Hybrid Driver–Support (HDS)](Figures/Fig_archetypesIV.png)
*Archetype IV — Hybrid Driver–Support (HDS)*


**Signature:** a few very strong genes plus a some moderately associated genes.

**Gene-level p-value pattern:**
- One or a few top genes are extremely significant, e.g.
  $$p_{(1)},\,p_{(2)} \ll 10^{-4}\quad(\text{often much smaller}).$$
- Beyond the top hits, several additional genes show moderate evidence:
  $$p_{(k)} \in [10^{-3},\,0.05]\quad\text{for multiple }k \text{ (support genes).}$$
- The remaining genes are near-null:
  $$p_{(j)} \sim \mathrm{Uniform}(0,1)\quad\text{for most other }j.$$
- A small “spike” at the top (drivers) plus a clear “shoulder” of moderately small p-values (support), then a flat tail.

**Interpretation:**  
Hybrid Driver–Support (HDS) exhibits a hierarchical pathway structure. A small number of “driver” genes exert the most significant effects, whilst a group of genes provides modest yet persistent associations. This is prevalent in pathways where flow or signal is regulated by a limited number of control points, whereas effective route output also relies on the synchronized activity of various downstream components. This architecture is statistically positioned between SDS and CME. A distinct driver signal exists, although the significance of the pathway is augmented by supplementary moderate signals.

**Biological example** – **LDL cholesterol as a hybrid driver–support pathway**

The regulation of LDL-cholesterol (LDL-C) exemplifies a scenario in which a limited number of genes serve as primary "drivers" within a more extensive polygenic framework. Extensive GWAS and sequencing investigations consistently identify *LDLR*, *APOB*, and *PCSK9* as fundamental genes associated with Mendelian hypercholesterolemia. Infrequent coding or splice-altering variants in these genes can alter LDL-C by approximately half a standard deviation or more and are responsible for classical familial hypercholesterolemia, significantly elevating the risk of coronary artery disease. Recent whole-genome sequencing investigations involving over 60,000 individuals reveal that even infrequent *non-coding* variations next to *LDLR* and *PCSK9* can exert effects comparable to clinically recognized exonic FH variants, hence supporting their significance as primary regulators of LDL-C homeostasis.

A supporting network of lipoprotein and cholesterol metabolism genes surrounds these drivers. GWAS of traditional lipids and nuclear magnetic resonance (NMR)-based lipoprotein characteristics have revealed numerous new loci affecting LDL particle size, concentration, and composition, including apolipoprotein clusters (*APOE/APOC*, *APOA1/A5*), hepatic lipase (*LIPC*), and transporters such as *ABCG5/ABCG8*. Individually, prevalent variants at these loci typically elucidate only a minor proportion of LDL-C variance, however, collectively, they constitute a significant segment of the genome-wide polygenic signal. Recent extensive multi-ancestry meta-analyses now identify hundreds of lipid-associated loci distributed throughout this extensive lipoprotein network.

Translating this biology into gene-level association statistics for an LDL-related trait reveals that the route is neither exclusively "sparse driver" nor entirely "coordinated moderate". Typically, one observes several robust gene-level signals at *LDLR*, *APOB*, *PCSK9*, and a few linked loci, supported by a wider array of modestly correlated genes implicated in lipoprotein assembly, remodeling, and cholesterol transport. The HDS pattern is characterized by a pathway whose enhancement is indicative of both a limited number of predominant, high-impact genes and a strong, albeit subtler, influence from the broader metabolic framework. In CATFISH terminology, this refers to the regime when **soft TFisher** (which prioritizes the lower tail while still considering moderate p-values), along with **Fisher** and the **omnibus combination**, aligns effectively with the underlying biology.

**Best detectors in CATFISH:**
- **Soft TFisher** (tail-focused; gains power from a few strong hits *plus* additional modest hits)
- **Fisher** (accumulates evidence across the moderate support set)
- **Omnibus combination** (e.g., ACAT across ACAT/Fisher/soft-TFisher/Stouffer; or permutation-calibrated minP across methods)

---

## Archetype V — Single-Gene Proxy Pathway (SGP)

![Archetype V — Single-Gene Proxy Pathway (SGP)](Figures/Fig_archetypesV.png)
*Archetype V — Single-Gene Proxy Pathway (SGP)*


**Signature:** the pathway looks significant only because it contains one very strong gene; the remaining members look essentially null.

**Gene-level pattern:**
- The top gene has an extremely small p-value (e.g. $10^{-7}\text{--}10^{-12}$).
- The rest of the genes have p-values that look noisy / $\mathrm{Uniform}(0,1)$, with no clear excess of small p’s.
- If you ignore the top gene, there is no obvious enrichment left in the pathway.

**Interpretation:**  
The pathway effectively serves as a proxy for a singular gene-level relationship. This frequently occurs when:
- an annotation term is very small (one gene plus a couple of weakly related neighbors), or
- numerous pathway definitions redundantly incorporate the same driver gene, resulting in several "distinct" pathways activating while all reference the same underlying gene.

Biologically, the driver gene remains significant (and may also align with the SDA), although the pathway-level assertion provides no additional information beyond indicating that "this gene is strongly associated to the pathway." In CATFISH, we consequently designate SGP patterns with a cautionary note, "view these as gene-centric findings with associated pathway designations, rather than as proof that the entire pathway is collectively activated."

**Biological example:**  
**PAH in phenylalanine metabolism**

In humans, phenylalanine metabolism is primarily regulated by a singular bottleneck enzyme, phenylalanine hydroxylase (PAH), a hepatic monooxygenase that catalyzes the conversion of phenylalanine to tyrosine in a process dependent on tetrahydrobiopterin (BH₄) (Scriver, 2007; Elhawary et al., 2022). Classical phenylketonuria (PKU) and associated forms of hyperphenylalaninemia (HPA) occur when phenylalanine hydroxylase (PAH) activity is significantly diminished, resulting in elevated blood phenylalanine levels, relatively decreased tyrosine levels, and the accumulation of neurotoxic metabolites. This culminates in the distinctive untreated PKU phenotype characterized by extremely high phenylalanine, low tyrosine, and progressive neurological impairment (Elhawary et al., 2022). Extensive clinical and molecular studies indicate that the predominant cause of HPA/PKU cases is pathogenic variants in PAH, whereas a minority results from deficiencies in BH₄ synthesis, recycling, or the PAH co-chaperone DNAJC12 (Blau et al., 2014; Himmelreich et al., 2021; Elhawary et al., 2022). In summary, within the overarching phenylalanine metabolism pathway, PAH represents the pivotal flux-controlling step. Both common and rare variations at PAH significantly influence systemic phenylalanine levels, while the majority of other pathway components (transporters, minor side-enzymes, cofactor recycling genes) have considerably weaker or less frequent effects at the population level.

**Best detectors in CATFISH:**
- **minP / Tippett** 
- **ACAT** 

---

## Archetype VI — Competitive Enrichment Above Background (CEAB) (OPTIONAL - Not used for CATFISH, just a check)

![Archetype VI — Competitive Enrichment Above Background (CEAB)](Figures/Fig_archetypesVI.png)
*Archetype VI — Competitive Enrichment Above Background (CEAB)*


**Signature:** the pathway is not merely “associated”, it is enriched above the genome-wide polygenic background. It passes a MAGMA competitive test ($\beta_s > 0$), indicating that genes inside the set exhibit, on average, greater associated than those outside the set.

**Gene-level pattern:**

Let $Z_g$ be the gene-level $Z$ used by MAGMA, derived from gene p-values via a probit transform (higher $Z$ = stronger association).  
In a CEAB pathway $s$:

- The distribution of $\{Z_g : g \in s\}$ is elevated in comparison to $\{Z_g : g \notin s\}$. Equivalently, $\mathrm{mean}(Z_{\mathrm{in\_set}}) > \mathrm{mean}(Z_{\mathrm{outside\_set}})$, rather than merely $\mathrm{mean}(Z_{\mathrm{in\_set}}) > 0$.
- The signal is generally not represented by a single-gene proxy, rather, one may see multiple modestly strong genes or a small top-tail alongside an elevated mean. The crucial aspect is that the *average* in-set relationship surpasses the background level
- The pathway retains its significance following MAGMA's default covariate adjustment for confounding gene characteristics (e.g., gene size and gene density, including log-transforms), indicating that it is not "large genes/dense genes" responsible for the observed effect.

**Interpretation:**  
CEAB signifies authentic enrichment rather than: (i) generic polygenicity (where multiple sets appear related under independent testing), (ii) artifacts from annotation overlap, or (iii) confounding effects due to gene size or density.  

Competitive tests possess a more generalized null hypothesis because they explicitly adjust for baseline associations inherent in polygenic characteristics. In MAGMA’s Crohn’s disease case study, numerous gene sets seemed related through self-contained testing. However, only one maintained significance under the competitive hypothesis demonstrating that CEAB identifies sets with associations that exceed expectations based on polygenicity.

**Biological example**
**From MAGMA’s Crohn’s Disease analysis**  

In the WTCCC Crohn’s disease dataset, MAGMA’s self-contained analysis identified 39 related gene sets. However the competitive analysis recognized only one of those 39 as enriched above background ("Regulation of AMPK activity via LKB1 (REACTOME))".

This exemplifies a standard CEAB pattern. The pathway meets the above-background enrichment criterion. Two supplementary sets (e.g., Cell adhesion molecules and ECM-receptor interaction) attained competitive significance solely when gene size/density correction was disabled, which MAGMA reads as (at least partially) confounding-induced inflating rather than genuine enrichment.

**Best detectors in CATFISH:**

- MAGMA competitive p-value (from `*.gsa.out`) as an external “above-background” anchor:  
  if $p_{\mathrm{comp}}$ is small and $\beta_s > 0$, interpret as CEAB-supported enrichment.

---

### Bias warning

Pathway-based analysis for assessing over-representation or enrichment is influenced by various biases, including gene size, pathway size, SNP coverage density, and linkage disequilibrium (LD) patterns, all of which must be addressed explicitly (White et al., 2020; PMC6391732). In CATFISH, we tackle these issues in three phases: (i) The SNP to gene analysis is conducted using MAGMA’s LD-aware multi-SNP model, ensuring that gene-level Z/P values account for local LD structure and SNP density; (ii) we subsequently regress MAGMA gene Z-scores against log(gene length) and log(number of SNPs), utilizing the residual-based $P_adj$ for all subsequent gene to pathway analyses, thereby eliminating any remaining dependence on gene size and SNP density; and (iii) at the pathway level, we avoid simplistic over-representation tests, opting instead to calibrate our omnibus statistics through gene-label LD-aware permutations that maintain each gene’s adjusted p-value and the observed distribution of pathway sizes, yielding enrichment p-values that are resilient to these established biases.

Link: https://pmc.ncbi.nlm.nih.gov/articles/PMC6391732/

---

## METHODS

### Notation

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

---

## 1) Gene-level association statistics (SNP → gene)

For each gene $g$, MAGMA generates a gene‑level association p‑value $p_g$ by aggregating SNP‑level signals within or within a window of the gene while accounting for local LD using a reference panel. We employed MAGMA’s `multi=snp-wise` model, which integrates a SNP-wise mean test (effective for numerous small effects) and a SNP-wise top test (effective for a single strong SNP) into a unified LD-aware omnibus statistic per gene, hence ensuring the robustness of gene p-values against varying within-gene causal structures.

The workflow emcompasses:

- SNPs are mapped to genes (gene boundaries with optional windows).
- A multi‑marker gene model accounts for LD among SNPs in the genic region.
- MAGMA generates gene statistics (e.g., $Z_g$ and $p_g$).

In CATFISH, the p-values (or Z statistics) of the MAGMA gene are utilized as inputs for all pathway analyses.

---

## 2) Gene-level adjustment for gene size and SNP density

Despite LD-aware gene testing, gene-level signals may still demonstrate residual dependency on gene size and SNP density. CATFISH executes a post-hoc correction at the gene level.

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

The corresponding adjusted two-sided $p$-values are:

$$
p_{g,\mathrm{adj}} = 2\,\Phi\!\left(-\left|Z_{g,\mathrm{adj}}\right|\right),
$$

where $\Phi(\cdot)$ is the standard normal CDF.

We denote the adjusted vectors by:

$$
\mathbf{Z}_{S\mathrm{adj}} = (Z_{1\mathrm{adj}}, \dots, Z_{G\mathrm{adj}}),
\qquad
\mathbf{p}_{S\mathrm{adj}} = (p_{1\mathrm{adj}}, \dots, p_{G\mathrm{adj}}).
$$


---

## 3) Pathway-level test statistics (gene → pathway)

CATFISH computes multiple pathway statistics from either unadjusted gene-level $p_g$ and $Z_g$, or the adjusted counterparts $p_{g,\mathrm{adj}}$ and $Z_{g,\mathrm{adj}}$. For convenience, we present definitions using the unadjusted inputs.

We define multiple pathway-level statistics as functionals of a common set of within-pathway, gene-level evidences $\{p_g\}_{g\in S}$ (and, when available, $\{Z_g\}_{g\in S}$).
The gene-level evidences for genes within the same pathway are typically dependent (e.g., due to linkage disequilibrium–induced correlation and shared genomic architecture). Consequently, the resulting pathway statistics are also mutually dependent, as they are derived from the same underlying inputs.

Therefore, closed-form reference calibrations that rely on independence of the gene-level tests (such as Fisher’s $\chi^2$ null or Tippett’s transformation for minP) are provided only as canonical or illustrative definitions. The final inferential $p$-values (both for each individual component statistic and for the omnibus statistic) are obtained via the unified null calibration procedure described in Section~4, which recomputes all statistics under a null-generating mechanism that preserves the dependence structure.


### 3.1 ACAT (Aggregated Cauchy Association Test)

We define the Cauchy‑transformed score for each gene:

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

ACAT is asymptotically dominated by the smallest p-values, and is therefore sensitive to SDAs.

---

### 3.2 Fisher’s method

The define the Fisher’s statistic as:

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

Fisher is sensitive to CMEs.

---

### 3.3 Adaptive Soft TFisher

A practical limitation of fixed- $\tau$ soft TFisher is that the appropriate tail-focus is contingent upon the *unknown* pathway topology (dense versus sparse, weak versus strong). TFisher proposes a **omnibus, data-adaptive** selection of truncation and weighting parameters, termed **oTFisher**, which autonomously identifies the most advantageous configuration for the observed p-value distribution (ref).

#### Soft TFisher family

For a pathway $S$ with adjusted gene-level p-values $\{p_g\}_{g\in S}$ and a soft-threshold parameter $\tau\in(0,1]$, the soft TFisher statistic is

$$W^{\mathrm{soft}}(S;\tau)=\sum_{g\in S}\left[-2\log(p_g)+2\log(\tau)\right]_{+},\qquad (x)_+=\max(x,0)$$

This is the $\tau_1=\tau_2=\tau$ special case of the TFisher family and implements a continuous down-weighting near the cutoff (soft vs hard truncation). 

#### Data-adaptive $\tau$ (oTFisher)

Let $\mathcal{T}=\{\tau_1,\dots,\tau_m\}$ be a small grid of candidate thresholds (e.g. a few small/medium/large values; TFisher shows that a sparse grid is usually sufficient). 

For each $\tau_j\in\mathcal{T}$, we compute:

1) the statistic $W^{\mathrm{soft}}(S;\tau_j)$, and  
2) its **null p-value**
$$p_{\tau_j}(S)=\Pr\!\left(W^{\mathrm{soft}}(\cdot;\tau_j)\ge W^{\mathrm{soft}}(S;\tau_j)\mid H_0\right),$$
using the TFisher null calculation for $W_n(\tau_1,\tau_2)$ with $\tau_1=\tau_2=\tau_j$. 

Then, we define the adaptive soft TFisher omnibus as the minimum across the grid as:

$$p_{\mathrm{aTF}}(S)=\min_{\tau\in\mathcal{T}} p_{\tau}(S)$$

This is the oTFisher concept, tailored to the soft-thresholding line $\tau_1=\tau_2$. Due to the dependence of the $p_{\tau_j}(S)$ (originating from the same ordered p-values), oTFisher offers an analytic calibration for the minimum across the grid by employing a multivariate normal approximation of the vector comprising component TFisher statistics (calculated through multivariate normal probabilities). In the soft scenario when $\tau_1=\tau_2=\tau$, the mean and covariance are simplified, allowing for efficient computation of the multivariate normal (MVN) probability, such as by Genz-style MVN cumulative distribution function (CDF). 

In CATFISH terminology, this provides a adaptive tail sensor without the necessity of hard-coding a single $\tau$ and without the need for extensive permutation for the within-method calibration (we still overlay our MVN/perm framework for the final omnibus across methods).

A minimal grid that usually works well in practice is something like:

$$\mathcal{T}=(\{0.01\,0.05\,0.5\,1\})$$

which oTFisher explicitly uses in its soft-thresholding omnibus examples, and which covers “rare-ish hits”, “moderate tail”, and “nearly Fisher”.

The oTFisher is more tailored towards HDS.

---

### 3.4 Stouffer's method

Stouffer's technique consolidates gene-level **Z** statistics instead of p-values and exhibits increased sensitivity to DPS. In CATFISH, the gene-level Z input is sourced directly from MAGMA’s gene output (`ZSTAT`). Significantly, MAGMA’s Z-scale is understood as a association-strength score, indicating that greater positive values signify stronger evidence of association, rather than indicating the direction of effect as trait-increasing or trait-decreasing. Consequently, the natural pathway-level Stouffer test in this context is one-sided (greater), assessing the enrichment of positive association strength inside the pathway.

**Default (unweighted) Stouffer:**

For a pathway $$S$$ with $$G=|S|$$ genes,

$$Z_{\mathrm{stouffer}}(S) = \frac{1}{\sqrt{G}} \sum_{g \in S} Z_g$$

$$p_{\mathrm{stouffer}}(S) = 1 - \Phi\!\left(Z_{\mathrm{stouffer}}(S)\right),$$

where $$\Phi(\cdot)$$ is the standard normal CDF.

**Optional (weighted) Stouffer:**

CATFISH optionally supports nonnegative per-gene weights $$w_g$$ (e.g., based on SNP count or other gene-level quantities). In that case,

$$Z_{\mathrm{stouffer}}^{(w)}(S)=\frac{\sum_{g\in S} w_g\, Z_g}{\sqrt{\sum_{g\in S} w_g^2}},\qquad p_{\mathrm{stouffer}}^{(w)}(S)=1-\Phi\!\left(Z_{\mathrm{stouffer}}^{(w)}(S)\right)$$

A **two-sided** option can be reported for completeness when a genuinely *signed* gene-level Z is available:

$$p_{\mathrm{stouffer}}^{(2\text{-sided})}(S)=2\,\Phi\!\left(-\left|Z_{\mathrm{stouffer}}(S)\right|\right),$$

This is not the default for MAGMA-style association-strength Z scores.

---

### 3.5 minP / Tippett test

For each pathway $S$ with $G$ genes, we define the minimum gene $p$-value as

$$
T_{\min}(S) = p_{\min}(S) = \min_{g \in S} p_g .
$$

An independence-based calibration is given by Tippett’s transform as: 

$$
p_{\mathrm{tippett}} = 1 - (1 - p_{\min})^{G}.
$$

However, CATFISH does not rely on this analytic mapping for inference because gene-level test statistics within a pathway are typically dependent (see Section~4). The minP statistic is emphasized not because it is uniquely sensitive to dependence (all constituent statistics are), but because it represents a qualitatively distinct mode of evidence that is driven almost entirely by the single most significant gene. Consequently, $(T_{\min}) serves primarily as a detector of sparse, single-gene–driven signals (SDA/SGP-type patterns), thereby complementing aggregate combination procedures (Fisher, Stouffer, softTFisher, ACAT) that are designed to capture more diffuse enrichment. We identify potential single-gene proxy pathways via a leave-one-gene-out diagnostic, in which the top-ranking gene is removed, and the test statistic is recomputed.

---

## 4) Unified null calibration (captures gene dependence and cross-method coupling)

For a pathway $S$, define the vector of component statistics, each computed from the same within-pathway
gene-level evidence:

$$
\mathbf{T}(S)=\big(T_{\mathrm{ACAT}},T_{\mathrm{Fisher}},T_{\mathrm{TF}}^{\mathrm{soft}},T_{\mathrm{Stouffer}},T_{\min}\big),
$$

where,

- $$T_{\mathrm{ACAT}}$$ is the ACAT statistic on $$S$$,
- $$T_{\mathrm{Fisher}}$$ is Fisher’s sum–log–p statistic,
- $$T_{\mathrm{TF}}^{\mathrm{soft}}(\tau)$$ is the soft–TFisher statistic at truncation $$\tau$$,
- $$T_{\mathrm{Stouffer}}$$ is the (optionally weighted) Stouffer Z-combination statistic, and
- $$T_{\min}=\min_{g\in S} p_g$$ is the Tippett/minP statistic.

Two sources of dependence are present: (i) gene-level inputs within $S$ are correlated (e.g., LD/shared genomic structure), and (ii) the components are mutually dependent because they are deterministic functions of the same $\{p_g\}$ (and, when used, $\{Z_g\}$).

To obtain valid inference without assuming independence at either level, CATFISH calibrates the omnibus under a single dependence-preserving null generator. We generate null replicates

$\{(\{p_g^{(b)}\},\{Z_g^{(b)}\})\}_{b=1}^{B}$ 

that preserve within-pathway gene dependence. For example, under an LD-aware MVN null,

$Z^{(b)} \sim \mathcal{N}(0, R_S)$,

followed by a Gaussian-copula mapping to $p$-values. Using the one-sided convention (larger $Z$ implies stronger association),

$p_g^{(b)} = 1 - \Phi(Z_g^{(b)})$.

(with a two-sided mapping used if the analysis adopts two-sided gene $p$-values). Alternatively, phenotype permutation may be used when feasible.

For each replicate $b$, we recompute all component pathway statistics to obtain $\mathbf{T}^{(b)}(S)$ and then compute the corresponding omnibus statistic using the same rule applied to the observed data:

$$
T_{\mathrm{omni}}^{(b)}(S)=\mathcal{O}\!\left(\mathbf{T}^{(b)}(S)\right),\qquad T_{\mathrm{omni}}^{\mathrm{obs}}(S)=\mathcal{O}\!\left(\mathbf{T}^{\mathrm{obs}}(S)\right).
$$

where $\mathcal{O}(\cdot)$ denotes the fixed omnibus decision rule (e.g., ACAT-across-methods, min-across-methods,
or another prespecified combination).

The omnibus $p$-value is estimated by the standard resampling tail probability:

$$
p_{\mathrm{omni}}(S)=\frac{1+\sum_{b=1}^{B}\mathbf{1}\!\left(T_{\mathrm{omni}}^{(b)}(S)\ge T_{\mathrm{omni}}^{\mathrm{obs}}(S)\right)}{B+1}
$$

Because each null replicate recomputes all components from an identically generated, dependence-preserving draw, this procedure jointly accommodates both intra-pathway gene dependence and inter-method coupling. Optionally, calibrated component-wise $p$-values can be derived in an analogous manner by applying the same tail-probability calculation to each component-specific test statistic $T_j$.



### 4.1 Deterministic coupling under the null (same $$\{p_g\}$$ reused)

Even under a pure null scenario where genes are independent, this vector is **not jointly independent** since each component is a deterministic function of the same multiset of p-values 
$p_g$ for $g \in S$ (equivalently the order statistics $P_{(k)}$). Specifically:

- Fisher and soft TFisher exhibit monotonicity in $$-\log(p_g)$$, with TFisher functioning as a shortened or reweighted Fisher that focuses weight on $$p_g \le \tau$$. Consequently, when numerous $$p_g$$ values are small (or the lower tail is relatively long), both $$T_{\mathrm{Fisher}}$$ and $$T_{\mathrm{TF}}^{\mathrm{soft}}(\tau)$$ typically exhibit extremity in the same direction.

- Stouffer integrates Z-scores, which in standard gene-set contexts, are monotonic transformations of p-values (or are directly supplied as gene-level association Z-scores). As a result, routes exhibiting consistent enrichment of tiny $$p_g$$ will also likely increase $$T_{\mathrm{Stouffer}}$$.

- ACAT and minP both are influenced by extreme tails. ACAT employs the heavy-tailed Cauchy transform $$\tan\{\pi(1/2 - p_g)\}$$, and minP is precisely $$P_{(1)}$$. Hence, when a pathway includes one (or more) unusually tiny p-values, both $$T_{\mathrm{ACAT}}$$ and $$T_{\min}$$ reach extreme values. Soft TFisher with a minimal $$\tau$$ exhibits analogous behavior as it assigns greater weight to the smallest p-values.

- Cross-method post hoc selection exacerbates dependence as any omnibus that utilizes the minimum across method p-values (or adaptively selects a $$\tau$$ in TFisher) give rise to “winner’s curse” behavior unless the null distribution is carefull adjusted.


### 4.2 Extra dependence from LD and shared gene-level correlation

In actual data, gene-level inputs exhibit correlation. The sets $$\{p_g\}$$ and $$\{Z_g\}$$ exhibit dependence across genes due to (i) the potential for adjacent genes to share SNPs or LD structures, (ii) the coupling of SNP-level evidence across genes through local LD, and (iii) the influence of a genome-wide polygenic background that may induce weak correlations in gene-level association metrics. This enhances reliance on pathway statistics beyond the aforementioned deterministic coupling.

In CATFISH we address this in two complementary ways:

1. **Upstream LD-aware SNP $$\rightarrow$$ gene aggregation (MAGMA)**

Gene-level $$p_g$$ and $$Z_g$$ are derived from MAGMA’s LD-informed SNP-to-gene model, which adjusts gene evidence considering the correlated SNP structure.

3. **Downstream empirical calibration of the omnibus across correlated tests**

Instead of supposing independence among $$\{p_j(S)\}$$, we calibrate the omnibus through two types of resampling:
   - Global gene-set resampling selects genes from the genome-wide repository and recalculates *all* component tests on the identical resampled gene sets, therefore capturing dependencies arising from shared inputs.
   - LD-aware MVN calibration simulates a singular correlated Gaussian draw $$Z \sim \mathcal{N}(0, R_S)$$ for each route (with $$R_S$$ constructed using MAGMA gene–gene correlations), and extracts p-based components from the identical draw by a Gaussian-copula mapping. This maintains both intra-pathway correlation and inter-test coherence (RECOMMENDED).

### 4.3 Implication for inference

Due to the strong dependence among $$\big(T_{\mathrm{ACAT}}, T_{\mathrm{Fisher}}, T_{\mathrm{TF}}^{\mathrm{soft}}(\tau), T_{\mathrm{Stouffer}}, T_{\min}\big)$$, naive combination rules that presume independence (such as analytic minP across methods under independence or straightforward Bonferroni adjustments across correlated method p-values) may be miscalibrated and frequently result in anti-conservatism. CATFISH regards the analytic across-method omnibus (ACAT-O or minP-O) as a useful summary and use resampling-calibrated omnibus p-values as the principal evidence. This calibration directly determines the null distribution of the omnibus using the identical recomputation methodology used to the observed pathway, thus regulating type-I error amidst dependency across component tests and LD among genes.

---

## 5) Omnibus pathway p-value across methods (analytic + resampling-calibrated; global + LD-aware MVN)

Each pathway in CATFISH is assessed by a panel of complementary gene-to-pathway tests (**ACAT**, **Fisher**, **adaptive soft TFisher**, **Stouffer**, and **gene-minP**), each calibrated to a distinct pattern of gene-level signal. Due to the fact that these five statistics derive from the same gene-level evidence ($$\{p_g\}$$ and optionally $$\{Z_g\}$$), they exhibit a **strong correlation** and may differ in their ranking of the "most" enriched pathways. Instead of selecting a singular test a priori, we consolidate the evidence for each pathway with a **single omnibus pathway p-value** employing two complementing methodologies:

1. **Analytic ACAT-O across methods** (a fast, sensitive summary across correlated method p-values), and  
2. **Resampling-calibrated omnibus** (global gene-set resampling and/or LD-aware MVN), which learns the null
   distribution of the omnibus under dependence and provides a conservative “best-of-tests” inference layer.

Throughout, let the five gene-derived component p-values for pathway $$S$$ be

$$
\mathcal{P}(S)=\{p_{\mathrm{ACAT}}(S),\,p_{\mathrm{Fisher}}(S),\,p_{\mathrm{TFisher}}(S),\,p_{\mathrm{minP}}(S),\,p_{\mathrm{Stouffer}}(S)\},
$$

with the Stouffer term included only when gene Z-scores are available.

> **MAGMA Competitive (Reported Separately):** If accessible, the MAGMA competitive gene-set p-value is disclosed.
> As an auxiliary pathway summary (`magma_pvalue`). It is **not included** in the resampling calibration by default.
> Utilize omnibus (`include_magma_in_perm=FALSE`) as the resampling layer beneath fails to produce a methodical null.
> for the MAGMA competitive regression without re-executing MAGMA (Section 5.6).


---

### 5.0 LD-aware SNP $$\rightarrow$$ gene inputs via MAGMA (upstream)

GWAS summary results were consolidated into gene-level association statistics utilizing MAGMA’s LD-aware SNP-to-gene model (multi-model SNP-wise gene analysis) with a reference LD panel. SNPs were assigned to genes utilizing a symmetric ± 25 kb frame. This produces per-gene association metrics such as $$p_g$$ and $$Z_g$$ that are adjusted to consider within-gene linkage disequilibrium (LD).

To mitigate confounding effects from gene size and SNP density, gene-level evidence may be adjusted by regressing on $$\log(L_g)$$ and $$\log(S_g)$$, utilizing residual-based size-adjusted gene p-values as inputs for pathways (where available as an alternative p-value column). In the implementation, p-based component tests utilize the preferred p-value column (e.g., an adjusted column if available, otherwise raw MAGMA gene p-values).

---

### 5.1 Component pathway tests (gene $$\rightarrow$$ pathway)

For each pathway $$S$$ with member genes $$g\in S$$, we compute:

1. **ACAT (gene-level)**
2. **Fisher (gene-level)**  
3. **Adaptive soft TFisher (gene-level)**  
4. **Gene-level minP (Sidák)**  
5. **Stouffer Z (optional; from gene Z-scores)**
6. 
**Numerical stability.** All p-values are clipped to $$[p_{\min},\,1-p_{\min}]$$ (e.g., $$p_{\min}=10^{-15}$$)
before applying $$\log(\cdot)$$ or $$\tan(\cdot)$$ transforms.

---

### 5.2 Omnibus ACAT across methods (ACAT-O)

Let $$p_1,\dots,p_K$$ denote the available component p-values for pathway $$S$$ (typically $$K=5$$), and let
weights $$v_j\ge 0$$ satisfy $$\sum_{j=1}^K v_j=1$$ (default $$v_j=1/K$$). Define the ACAT-O statistic

$$
T_{\mathrm{omni,ACAT}}(S)=\sum_{j=1}^{K} v_j \tan\{\pi(0.5-p_j)\}.
$$

The analytic omnibus p-value is

$$
p_{\mathrm{omni,ACAT}}(S)=0.5-\frac{1}{\pi}\arctan\{T_{\mathrm{omni,ACAT}}(S)\}.
$$

The ACAT-O layer has heightened sensitivity when at least one component test demonstrates great significance, regardless of the modest performance of other components (e.g., sparse drivers, coordinated enrichment, hybrid subsets).

---

### 5.3 “Best-of-tests” omnibus via minP across methods (minP-O)

To derive a supplementary, more conservative summary that just focuses on the *best-performing* component test for each pathway, we calculate a minimum-p omnibus across methodologies:

$$
T_{\mathrm{omni,min}}(S)=\min_{p\in\mathcal{P}(S)} p
=\min\{p_{\mathrm{ACAT}}(S),\,p_{\mathrm{Fisher}}(S),\,p_{\mathrm{TFisher}}(S),\,p_{\mathrm{Stouffer}}(S),\,p_{\mathrm{minP}}(S)\}.
$$

A simple analytic conversion under independence is

$$
p_{\mathrm{omni,min}}(S)=1-\bigl(1-T_{\mathrm{omni,min}}(S)\bigr)^K,
$$

but because the component tests are correlated, inference is based on **empirical calibration** (Section 5.4).

---

### 5.4 Resampling calibration of the omnibus (global and LD-aware MVN)

The five component pathway tests (ACAT, Fisher, soft TFisher, Stouffer, and gene-minP) are statistically dependent under the null hypothesis, as they all derive from the same gene-level inputs. To achieve proper significance levels that account for this reliance, CATFISH employs two complementary resampling frameworks:

1. **Global gene-set resampling**, which extracts null gene sets directly from the empirical MAGMA gene-level outcomes; and
2. **LD-aware multivariate normal (MVN) simulation**, which produces correlated null Z-scores utilizing the MAGMA gene–gene correlation framework for each pathway.

Both methodologies adhere to a fundamental principle: in every permutation duplicate, **the complete array of component tests is recalculated from a fresh sample of gene-level evidence**, so ensuring that the interdependence between methods is accurately maintained.

---

#### 5.4.1 Global gene-set resampling (`perm_mode="global"`)

**Concept.**  
The global resampling approach generates the null distribution of the omnibus by re-sampling genes instead of SNPs. It maintains the empirical distribution of gene-level p-values and Z-scores from MAGMA, while randomizing their allocation to pathways. This simulates a situation in which the genome-wide association landscape remains intact yet is not associated with any specific biological route designation.

**Step 1 – Define the global gene pool.**  
We define a gene pool $$\mathcal{G}$$ containing all genes with valid pathway inputs. Each gene contributes:
- its p-value $$p_g$$ from the preferred MAGMA p-value column (adjusted if provided, otherwise raw), and
- if Stouffer is used, its corresponding Z-score $$Z_g$$ from the `z_col` field in the same MAGMA gene results.

Both vectors $$\{p_g\}$$ and $$\{Z_g\}$$ are aligned so that each gene $$g$$ has a paired $$(p_g, Z_g)$$. This ensures that the p-value and Z-score for a given gene are always resampled together in every replicate.

**Step 2 – Sample null gene sets.**  
For each pathway $$S$$ of size $$|S| = d$$, and each permutation $$b = 1, \dots, B$$:

1. **Randomly draw gene indices**  
   Select $$d$$ unique gene indices $$I^{(b)} = (i_1, \dots, i_d)$$ uniformly from $$\{1, \dots, |\mathcal{G}|\}$$.  
   Sampling is *without replacement* when $$d \leq |\mathcal{G}|$$; if a pathway is extremely large, fallback to *with replacement* is allowed.

2. **Construct paired null evidence**  
   The resampled gene-level evidence is  

   $$P^{(b)} = (P_{i_1}, \dots, P_{i_d}), \quad Z^{(b)} = (Z_{i_1}, \dots, Z_{i_d}) \text{ (if Stouffer enabled)}$$

   This "paired resampling" ensures that each gene contributes its observed correlation between $$p_g$$ and $$Z_g$$ to the identical replicate, and that all component tests in replicate $$b$$ utilize the same foundational gene selection.
   
**Step 3 – Recompute all five component tests.**  
For every null draw $$(P^{(b)}, Z^{(b)})$$, CATFISH recomputes:
- $$p_{\mathrm{ACAT}}^{(b)}(S)$$ using the Cauchy combination on $$P^{(b)}$$;
- $$p_{\mathrm{Fisher}}^{(b)}(S)$$ via the sum–log–p statistic;
- $$p_{\mathrm{TFisher}}^{(b)}(S)$$ for each $$\tau$$ in the same grid as the observed test, recording the minimum;
- $$p_{\mathrm{minP}}^{(b)}(S)$$ from the smallest gene p-value;
- and $$p_{\mathrm{Stouffer}}^{(b)}(S)$$ using $$Z^{(b)}$$ under the same (one-sided) alternative.
  
The Stouffer null is typically regarded as unweighted for numerical stability; however, this does not influence dependence preservation, as the same genes are sampled collectively across all components.

**Step 4 – Combine resampled components into omnibus.**  
The set of replicate component p-values $$\{p_j^{(b)}(S)\}$$ are combined using the same omnibus rule (ACAT-O or minP-O) applied to the observed data:

$$p_{\mathrm{omni}}^{(b)}(S) = f_{\mathrm{omni}}\!\left(\{p_j^{(b)}(S)\}\right),$$

where $$f_{\mathrm{omni}}$$ denotes either the ACAT or minP operator.

**Step 5 – Empirical calibration.**  
The permutation-calibrated omnibus p-value is obtained as

$$\hat p_{\mathrm{omni,global}}(S) = \frac{1 + \left|\{\,b : p_{\mathrm{omni}}^{(b)}(S) \le p_{\mathrm{omni}}(S)\,\}\right|}{B + 1}$$

The "+1 correction" eliminates zero p-values and produces unbiased estimates, even with modest $$B$$ values.  
The cross-method correlation is inherently preserved because all five component statistics are recalculated on the identical resampled gene sets.

**Interpretation.**  
Global resampling offers a data-driven, LD-agnostic calibration of the omnibus, honoring the genome-wide distribution of gene-level evidence while randomizing route affiliation. It is computationally efficient and produces valid null hypotheses despite arbitrary reliance among component tests.

---

#### 5.4.2 LD-aware MVN calibration (`perm_mode="mvn"`)

**Concept.**  
Although global resampling maintains empirical p-value distributions, it fails to explicitly account for *within-pathway* linkage disequilibrium among genes. The LD-aware MVN calibration resolves this issue by generating correlated gene-level Z-scores from a multivariate normal distribution defined by MAGMA’s gene–gene correlation matrix $$R_S$$ for each pathway.


**Step 1 – Build the correlation matrix $$R_S$$.**  
For each pathway $$S=\{g_1,\dots,g_d\}$$, we extract pairwise gene correlations $$r_{ij}$$ from a MAGMA correlation file (three columns: `gene1`, `gene2`, `r`). We construct a symmetric $$d \times d$$ matrix $$R_S$$ with:
- diagonals set to 1 ($$R_{ii}=1$$),
- $$R_{ij}=r_{ij}$$ when available,
- missing pairs defaulted to 0,
- correlations clipped to $$|r| \le 0.999$$,
- and (optionally) enforced positive-definiteness using a *nearest positive definite* (nearPD) correction.

This ensures that $$R_S$$ is numerically stable and reflects LD and gene overlap within the pathway.

**Step 2 – Simulate correlated null Z-scores.**  
For each replicate $$b = 1, \dots, B$$:

$$Z^{(b)} \sim \mathcal{N}(0, R_S),$$

so that the null gene-level Z-scores have the same correlation structure as in the observed data.

**Step 3 – Derive null p-values from the same simulated $$Z^{(b)}$$.**  
To ensure that Stouffer and the p-based tests are coherent, both utilize the *same simulated draw* $$Z^{(b)}$$.  
We transform $$Z^{(b)}$$ into gene-level null p-values using a Gaussian copula:

- **Uniform marginals (default)**  

  $$U^{(b)} = \Phi(Z^{(b)}), \quad P^{(b)} = 2\min(U^{(b)},\,1-U^{(b)}),$$

  producing Uniform$$(0,1)$$ marginals while maintaining correlation via $$R_S$$.

- **Empirical marginals (optional)**  
  The same $$U^{(b)}$$ can be mapped through an empirical quantile function estimated from the observed genome-wide $$p_g$$ distribution. This preserves empirical tail inflation or deflation while keeping the correlation structure identical.

**Step 4 – Recompute component tests under MVN null.**  
Using the simulated p-values $$P^{(b)}$$ and the same Z-scores $$Z^{(b)}$$, CATFISH recomputes the full component panel:
- ACAT, Fisher, TFisher (utilizing the identical `tau_grid`), and gene-minP derived from $$P^{(b)}$$;
- Stouffer from $$Z^{(b)}$$, employing consistent weights and a one-sided alternative and integrates them via the chosen omnibus rule (ACAT-O or minP-O) to yield $$p_{\mathrm{omni}}^{(b)}(S)$$.

**Step 5 – Empirical calibration.**  
As in global resampling, the MVN-calibrated omnibus p-value is:

$$\hat p_{\mathrm{omni,mvn}}(S) = \frac{1 + \left|\{\,b : p_{\mathrm{omni}}^{(b)}(S) \le p_{\mathrm{omni}}(S)\,\}\right|}{B + 1}$$

This calibration explicitly accounts for the local LD structure encoded in $$R_S$$ and reproduces the multivariate dependence between genes under the null.

**Interpretation.**  
The MVN methodology offers a comprehensive LD-aware permutation layer: it transmits MAGMA’s gene–gene correlations into a multivariate normal framework and systematically reassesses all five component tests. By deriving both Z- and p-based methodologies from a singular correlated Gaussian sample, MVN resampling concurrently captures authentic gene-level dependence and cross-test correlation. It is more computationally intensive than global resampling but statistically more accurate to the true null in LD-rich areas.

---

### 5.5 Choice of final omnibus p-value (what we report)

Depending on the resampling mode used:

- $$p_{\mathrm{omni,analytic}}(S)$$: Analytic omnibus (ACAT-O or minP-O).  
- $$\hat p_{\mathrm{omni,global}}(S)$$: Global gene-set resampling calibrated omnibus.  
- $$\hat p_{\mathrm{omni,mvn}}(S)$$: LD-aware MVN calibrated omnibus.

The **final omnibus p-value** is chosen as:

$$
p_{\mathrm{omni,final}}(S) =
\begin{cases}
\hat p_{\mathrm{omni,mvn}}(S), & \text{if MVN calibration was performed;}\\
\hat p_{\mathrm{omni,global}}(S), & \text{else if global calibration was performed;}\\
p_{\mathrm{omni,analytic}}(S), & \text{otherwise.}
\end{cases}
$$

This final omnibus is then adjusted across pathways using Benjamini–Hochberg FDR to obtain
$$q_{\mathrm{omni,final}}(S)$$, reported as `omni_p_final_BH`.

---

### 5.6 Treatment of MAGMA competitive in the omnibus (explicit)

We also calculate and present the MAGMA competitive gene-set p-value (`magma_pvalue`) as an independent summary.
By default, it is **excluded** from the resampling-calibrated omnibus (`include_magma_in_perm=FALSE`) because the aforementioned resampling strategies provide null realizations just for **within-pathway** gene evidence ($$p_g$$ and $$Z_g$$ for genes $$g\in S$$). A principled null for the MAGMA competitive statistic necessitates rerunning a competitive regression (or MAGMA itself) for each duplicate on a suitable genome-wide null, which is not executed in this context. Thus, the resampling-calibrated omnibus is calculated exclusively for the five gene-derived component tests, and MAGMA competitive is analyzed in conjunction with the omnibus rather than being integrated into it.

---

## 6) Multiple testing correction

Across all pathways, the final omnibus p-values $$\{p_{\mathrm{omni,final}}(S)\}$$ are adjusted using the
Benjamini–Hochberg FDR procedure:

$$
q_{\mathrm{BH}}(S)=\mathrm{BH}\big(p_{\mathrm{omni,final}}(S)\big)
$$

Since each pathway produces a single final omnibus p-value, no supplementary penalty is necessary for the quantity of component tests. The post hoc "best-of-tests" selection is inherently addressed by the resampling calibration when activated.

---

## USAGE

# CATFISH (R package wrapper)

**CATFISH** is the R interface implementation used to run CATFISH‑style workflows on top of MAGMA, and to compute ACAT/Fisher/TFisher + omnibus pathway statistics.

---

## Installation

### 1) Install MAGMA (external dependency)

Download MAGMA from the official site and make the `magma` executable available on your `PATH`:

- https://ctg.cncr.nl/software/magma

### 2) Install CATFISH in R

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
     - soft TFisher (tail-focused),
     - Stouffer's test,
     - minP.

4. **Omnibus**
   - Combine pathway p-values using minP or ACAT to produce $p_{\mathrm{omni}}$.

5. **Multiple testing**
   - BH FDR (and optional Storey q-values).
  
6. Permutation
   - Use either random sampling or MVN (RECOMMENDED).

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
  effect_col        = "Z_adj",
  is_onetail        = FALSE,
  ptrunc            = 0.05,
  min_p             = 1e-15,
  do_fix            = TRUE,
  omnibus           = "ACAT",      # minP or "ACAT"
  B_perm            = 10000L,
  seed              = 123,
  perm_mode    = "mvn",       # mvn or 
  magma_genes_out = "/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_multi_snp_wise_genes_by_chr_N_maize/magma_N_maize.txt",
  remove_singletons = TRUE,
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
##   - BH FDR for omni_p and each component
##   - omni_perm_BH  : permutation-calibrated BH pvalues using perm pvalues

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

## RESULTS

### Drosophila melanogaster Starvation Response GWAS example 

### From sex-specific top genes to shared functional programs: pathway-level reproducibility and what the omnibus captures (Fig. A–F)

![Pathway analysis](Figures/Fig2/Fig2.png)
*Pathway analysis*


We first evaluated the highest-ranking genes from the sex-stratified GWAS to quantify SNP-level cross-sex reproducibility (Fig. A). Among the top 50 genes identified in each sex, the overlap was minimal, with only 4 genes shared between strata. Following LD adjusted SNP-to-gene aggregation using MAGMA (Fig. B), gene-level concordance also remained limited, with 5 of the top 50 genes shared across sexes. These findings indicate that the most strongly associated genes are largely sex-specific, consistent with a highly polygenic architecture in which small effect sizes, LD structure, and sampling variability can substantially reorder gene-level association rankings between sexes.

In this sex-stratified Drosophila GWAS, BH adjustment produced sparse or near-empty sets of discoveries at pathway level. This outcome is consistent with limited statistical power under sex stratification and mixed-model correction, together with the stringent multiple-testing burden imposed by thousands of genes and hundreds of pathways. Consequently, in addition to reporting BH-adjusted q-values, we characterized pathway prioritization using a rank-based Top-K framework (with a fixed K across sexes and methods) to enable direct comparisons when the number of FDR-significant findings is minimal. The smallest K at which a non-empty six-way intersection (across the five component tests and Omni_MVN) emerged was K = 20, which we adopt as the primary threshold for pathway-level summaries. Under this prioritization scheme, Omni_MVN increased the cross-sex concordance relative to several single-statistic summaries, consistent with the notion that weak, distributed association signals become more stable after aggregation from genes to pathways and integration across multiple statistics, rather than producing large numbers of FDR-significant pathways.

Substantial differences in cross-sex concordance were observed across the various component-level statistics (Fig. C). Stouffer’s method and Omni_MVN exhibited the greatest overlap in prioritized pathways, each sharing 7 of 20 pathways across sexes, whereas ACAT and Fisher’s method demonstrated limited overlap (1/20 and 2/20 shared pathways, respectively), and minP identified no shared pathways (0/20). This pattern is consistent with the distinct statistical properties of these methods. MinP is largely driven by the single most significant gene within a pathway, such that modest sex-specific variation in the leading gene can substantially alter pathway rankings even when the underlying pathway-level signal is comparable. By contrast, Stouffer’s method (and frequently TFisher) is more responsive to diffuse, pathway-wide enrichment, which tends to produce more stable cross-sex pathway prioritization under polygenic settings.

To evaluate whether cross-sex convergence observed for Omni_MVN primarily reflected effective replication of a single constituent statistic, as opposed to genuine multi-method integration and/or stabilization of pathway rankings, we decomposed the omnibus signal within each sex by intersecting the Omni_MVN Top-20 pathways with the Top-20 sets from each individual method (Fig. D). In females, the Omni_MVN Top-20 exhibited the greatest overlap with Fisher’s method (18/20 pathways), followed by TFisher and Stouffer (14/20 each), ACAT (12/20), and minP (10/20). In males, Omni_MVN overlapped most strongly with Stouffer (18/20) and TFisher (17/20), with more modest overlap with minP (11/20) and Fisher’s method (10/20), and the lowest overlap with ACAT (8/20). These patterns indicate that Omni_MVN preferentially tracks the constituent statistic whose signal is most concordant with the sex-specific enrichment architecture, while remaining non-identical to any individual method’s Top-20 list. This behavior is consistent with correlation-aware integration of multiple statistics rather than simple duplication of any single component method.

Within-sex UpSet visualizations further elucidate why integrative testing can yield greater stability than reliance on any single component method (Fig. E–F). In females, the Top-20 pathway set comprised a six-method consensus core of four pathways identified by all approaches, several intermediate-size partial overlaps (3–4 pathways), and numerous method-specific singletons (Fig. E). In males, the six-method consensus core was larger (six pathways), but the remaining overlap architecture was more fragmented (Fig. F). This pattern is consistent with the distinct statistical focus of the component tests. MinP preferentially detects pathways dominated by strong single-gene signals, Fisher’s method and ACAT emphasize pathways with a broad accumulation of small-to-moderate effects, and Stouffer’s method and TFisher are more sensitive to coordinated shifts in signal intensity (including truncated-tail enrichment). Overall, these results suggest that Omni_MVN improves cross-sex concordance mainly by integrating evidence across multiple enrichment statistics, which stabilizes pathway rankings under weak, distributed signal. Shared-gene overlap among pathways likely contributes to this stability, but the pattern is inconsistent with Omni_MVN simply reproducing any single component test.


--- 

### Pathway-level patterns of starvation tolerance are conserved across sexes, with sex-biased metabolic strategies


![Top 20 pathways enriched for female and male flies for starvation response](Figures/Fig_top20_bubbles_by_sex.png)
*Top 20 pathways enriched for female and male flies for starvation response*


We applied the CATFISH multi-statistic pathway framework to gene-level summary statistics from sex-stratified starvation resistance GWAS in *Drosophila melanogaster*, using pathway annotations derived from the BioCyc/MetaCyc databases. No pathway reached statistical significance after correction for multiple hypothesis testing at the prespecified false discovery rate (FDR) threshold. Accordingly, we treat these findings as exploratory pathway rankings rather than as definitive evidence of causal associations. Nonetheless, the highest-ranked pathways exhibit non-random biological clustering, implicating candidate metabolic and proteostasis-related processes that are broadly consistent with known aspects of starvation physiology. We must note that, the enriched pathways identified are most appropriately interpreted as an over-representation of genes annotated to the regulation and execution of gluconeogenic processes, rather than as direct evidence for altered metabolic flux through these pathways.

### Recurring cross-sex biological processes among the top-ranked pathways

A small set of pathway annotations appears in the top ranks for both sexes, suggesting a shared response to starvation:

**Glucose homeostasis and fuel partitioning.**  
Both sexes exhibited significant enrichment for *gluconeogenesis III*, consistent with the requirement to sustain circulating glucose concentrations as stored carbohydrate reserves become depleted. 

**NAD⁺ metabolism and redox balance.**  
*NAD biosynthesis* emerged as one of the most strongly enriched pathwy (*superpathway of NAD biosynthesis in eukaryotes* in both sexes), with females additionally showing enrichment for an NAD salvage pathway. This convergence suggests that genetic variation affecting pathways responsible for maintaining redox homeostasis and cofactor pools becomes particularly relevant under nutrient limitation.

**Membrane lipid remodeling.**  
Both sexes showed enrichment of sphingolipid-associated pathways (e.g., *ceramide degradation* and *glycosphingolipid biosynthesis*).

**Proteostasis and secretory pathway function.**  
*Protein N-glycosylation (initial phase)* was highly ranked in both sexes, suggesting that allelic variation in genes involved in endoplasmic reticulum (ER) and secretory pathway processing, including protein folding, glycoprotein maturation, and quality-control mechanisms, may modulate tolerance to starvation.

**Isoprenoid metabolism.**  
The *mevalonate pathway* was enriched in both sexes, indicating a shared annotation signal involving isoprenoid metabolism.


Beyond the shared core, the top-20 lists show sex-skewed emphases that may reflect differences in which biological modules are most genetically “visible” in each sex under starvation. **Females** exhibited an enhanced representation of lipid-associated annotations, including fatty acid– and PUFA-derived modules, together with pathways plausibly involved in stress-buffering processes (e.g., carbonyl detoxification and polyamine/methionine recycling). These patterns are consistent with a prominent lipid-handling component in the female-associated rankings. **Males** exhibited an enhanced representation of recycling- and salvage-associated annotations (e.g., nucleoside- and nitrogen-related modules), as well as additional pathways implicated in membrane lipid remodeling. Some pathway labels (e.g., those annotated as “mammals” or “yeast”) likely arise from cross-species naming conventions within MetaCyc.


Overall, despite the absence of FDR-significant pathways, the top-ranked results cluster around biologically plausible starvation-related categories—gluconeogenesis, NAD/redox metabolism, membrane lipid remodeling, proteostasis, and isoprenoid-related metabolism—with a small shared core and sex-skewed remainder.



---

### Arabidopsis thaliana Lowest Temperature of the Coldest Month GWAS example 

![Arabidopsis Lowest Temperature Environmental GWAS](Figures/Fig3/Fig3.png)
*[Arabidopsis Lowest Temperature Environmental GWAS*


### Multi-method agreement among top pathways for Arabidopsis cold-response GWAS and Omni_MVN overlaps broadly with each component test, while retaining a small distinct subset

For Arabidopsis, we analyzed a single phenotype representing cold-stress intensity, quantified as the coldest-month temperature (BIO6), and assessed pathway enrichment using the CATFISH framework (ACAT, Fisher, Adaptive TFisher, minP, and Stouffer), in conjunction with the omnibus correlation-aware synthesis method (Omni_MVN). To characterize the consistency with which pathways are prioritized across statistical approaches, we extracted the top 15 pathways from each method and visualized their overlap using an UpSet plot (Fig. 3A). A salient feature of this plot is a large, high-confidence “consensus core”: **six pathways are shared by all six approaches (the five individual tests plus Omni_MVN)**, indicating that a substantial portion of the strongest pathway-level signal is robust to the choice of enrichment statistic. Outside this consensus core, the remaining top-ranked pathways distribute into smaller intersections and method-specific sets, consistent with anticipated differences in sensitivity profiles among the tests (for example, methods that are more responsive to single extreme genes versus those that capture more diffuse shifts across many genes). Notably, Omni_MVN predominantly contributes to the major multi-method intersections, rather than yielding a large set of Omni-specific pathways, which is consistent with its role as an omnibus procedure that concentrates its highest-ranked pathways in regions of agreement among multiple methods.

We next quantified the relationship between the omnibus ranking and each constituent test by directly comparing the **top 15 Omni_MVN pathways** to the **top 15 pathways** identified by each individual method (Fig. 3B). Omni_MVN shared the majority of its leading pathways with every component statistic, with overlap sizes of **12/15 with ACAT**, **10/15 with Fisher**, **11/15 with Adaptive TFisher**, **11/15 with minP**, and **10/15 with Stouffer**. In each pairwise comparison, Omni_MVN also preserved a modest method-specific subset (the non-overlapping portion), indicating that it is not merely replicating the ranking of any single test. Rather, the consistent and substantial overlap across all five approaches suggests that Omni_MVN operates as an integrative summary statistic: it preferentially prioritizes pathways that receive convergent support across heterogeneous enrichment formulations, while still permitting a limited degree of omnibus-specific refinement. Taken together with the multi-set structure in Fig. 3A, these findings indicate that—in this Arabidopsis BIO6 application—the pathway-level signal is both (i) robust across multiple testing paradigms and (ii) efficiently consolidated by the omnibus framework into a stable, high-confidence top-ranked set.


### Pathway-level patterns of cold tolerance/adaptation

These top hits read like a **coherent cold-climate adaptation bundle** for Bio6 in *Arabidopsis*: **carbon buffering + photoprotection + wall mechanics + cuticle barrier**, with “support” metabolism that stabilizes **redox and nucleotide pools** when cold slows growth. Below I expand each module, and I’ll flag where the link is **directly cold/freezing-related** vs **plausible-but-not-specific**.

### Carbon storage and sugar-based cold protection

The most clearly Bio6-associated signal involves **carbohydrate partitioning**, specifically the pathways *starch biosynthesis* (PWY-622), *superpathway of sucrose and starch metabolism* (PWYQT-4466), and *UDP-α-D-glucose biosynthesis I* (PWY-7343). Cold acclimation typically shifts plants toward **soluble sugar accumulation and enhanced carbon-buffering capacity**, in part through modified starch turnover and sucrose metabolism. These soluble carbohydrates function as **cryoprotectants/osmoprotectants** and as readily mobilizable energy reserves under conditions where enzymatic activities and growth are constrained by low temperature (Fürtauer et al., 2019). Given that UDP-glucose is a central precursor for both **storage carbohydrates** and **cell-wall biosynthetic intermediates**, the observed enrichment of UDP-glucose biosynthesis is also consistent with the hypothesis that cold environments favor coordinated regulation of **carbon allocation** among storage, transport, and structural polysaccharides (Fürtauer et al., 2019). ([PubMed][1])

### Photosynthesis machinery and winter photoprotection

The pigment/cofactor cluster—*chlorophyll a* biosynthesis II (PWY-5064), carotenoid/neurosporene biosynthesis (CAROTENOID-PWY, PWY-6287), and phylloquinone (vitamin K1) biosynthesis (PWY-5863)—is consistent with a canonical “winter” constraint, in which **low temperature reduces metabolic capacity and repair rates while incident irradiance can remain high**, thereby generating an excitation pressure and reactive oxygen species (ROS) challenge. A key mechanism is that low temperature **exacerbates photoinhibition predominantly by decelerating repair processes** (e.g., D1 protein turnover in photosystem II), making the maintenance and optimization of photochemical efficiency and repair a critical dimension of cold tolerance (Murata et al., 2007; Bascuñán-Godoy et al., 2012). Carotenoids—particularly components of the xanthophyll cycle—play a central role in **dissipating excess excitation energy** via non-photochemical quenching and in mitigating oxidative damage under stress conditions (Demmig-Adams et al., 1996). Phylloquinone functions as an electron-transport cofactor in photosystem I (PSI); *Arabidopsis* mutants with defects in phylloquinone biosynthesis exhibit **reduced PSI activity**, providing a direct mechanistic rationale for why selection may act on this pathway under conditions of cold–light imbalance (Gross et al., 2006).

**Caveat (important for interpretation):** pigment/cofactor pathways can also track **radiation/photoperiod/latitude** correlated with Bio6, not only temperature per se (Huner et al., 1998). So I’d call these “highly plausible” but potentially **climate-covariate confounded** unless your top genes are clearly cold-response regulators. ([ScienceDirect][2])

### Cell wall remodeling and extracellular sugar polymers

The wall-associated pathway cluster—*xyloglucan biosynthesis* (PWY-5936), *UDP-galacturonate biosynthesis* (PWY-4), and *UDP-β-L-arabinose biosynthesis* (PWY-82)—is highly congruent with current knowledge of freezing and chilling biology. Cold and sub-zero acclimation are accompanied by extensive **apoplastic and cell-wall remodeling**, including alterations in wall composition and biophysical properties that affect tolerance to **freeze–thaw deformation and dehydration** (Takahashi et al., 2019). A growing body of evidence directly links discrete cell-wall constituents to acquired freezing tolerance, such as coordinated modifications in **pectic polymers** during acclimation and deacclimation (Kutsuno et al., 2022), and data indicating that the fine structure of pectic side chains contributes to the enhanced freezing tolerance conferred by cold acclimation (Takahashi et al., 2024). Consequently, the prominence of UDP-sugar precursor pathways in conjunction with hemicellulose- and pectin-related biosynthetic routes among the top-ranked pathways is not only mechanistically plausible but is increasingly substantiated by experimental evidence. ([Nature][3])

### Surface lipid barrier and cuticle-type traits

*Wax esters biosynthesis I* (PWY-5884) is plausibly linked to adaptation to cold climates via an “outer-shell” mechanism. Cold seasons frequently coincide with **low atmospheric humidity, elevated wind exposure, and freeze-induced dehydration risk**, conditions under which cuticular characteristics can exert a substantial influence on non-stomatal water loss and on the biophysical dynamics of tissue freezing. In *Arabidopsis thaliana*, genetic or physiological perturbation of cuticular wax deposition modulates sensitivity to **dehydration stress** and to **low-temperature/freezing outcomes** following cold acclimation. Comparisons between wax-deficient and wax-overproducing lines reveal quantifiable differences in these stress-response phenotypes (Rahman et al., 2021). Consequently, wax and cuticle biology represent a mechanistically credible axis of adaptation to Bio6-associated environments, particularly where the predominant selective regime reflects combined **cold and desiccation** stress rather than temperature effects in isolation. ([PubMed][4])

### Central energy, redox control, and nucleotide economy

The “support metabolism” block—*2-oxoglutarate → succinyl-CoA* (PWY-5084), *NAD de novo biosynthesis* (PYRIDNUCSYN-PWY), together with nucleotide salvage/dephosphorylation (*PWY-7193*, *PWY-7177*) and purine intermediate pathways (*PWY-6121/6122*)—can be interpreted as processes that contribute to maintaining metabolic homeostasis when low temperature constrains growth and perturbs cellular redox balance. Conceptually, exposure to cold can generate an “energy imbalance” and redox disequilibrium that necessitate coordinated adjustment of photosynthetic and respiratory activities, as well as of downstream metabolic networks (Huner et al., 1998). Among the pathways identified, *NAD de novo biosynthesis* is the most directly interpretable in a stress-physiological context: NAD/NADP homeostasis in plants is repeatedly associated with developmental regulation and stress tolerance, and de novo NAD synthesis from aspartate constitutes a well-established plant route explicitly discussed in relation to stress responses (Hashida et al., 2009).

In contrast, the linkage between nucleotide salvage/dephosphorylation pathways and cold exposure is likely to be more indirect. These routes may primarily reflect enhanced requirements for nucleic acid repair, transcriptional and translational capacity, and general metabolic “housekeeping” in cold-adapted genotypes, rather than representing cold-specific signaling or defense mechanisms. Unless the principal genes within these pathways are established components of canonical cold-response networks, it may be more appropriate to describe them as *enriched supportive metabolic processes* rather than as definitive “cold-response pathways” ([ScienceDirect][5]).


### Hormone and defense tuning

Finally, *IAA biosynthesis VII* (PWY-6066) and *salicylate glucosides biosynthesis II* (PWY-6623) are highly plausible candidates for regulatory adjustment points in cold environments, given that cold acclimation typically entails a combination of growth suppression and reprogramming of stress and defense responses. On the auxin side, low temperature can directly perturb auxin homeostasis and signaling via effects on transport: in *Arabidopsis* roots, cold rapidly inhibits basipetal auxin transport and disrupts the trafficking of auxin efflux carriers (Shibasaki et al., 2009). This makes enrichment of an “auxin module” consistent with chilling-associated shifts in growth strategies. On the salicylate side, prolonged chilling can promote accumulation of salicylic acid (SA), and SA has been shown to inhibit growth at low temperature in *Arabidopsis* (Scott et al., 2004). Thus, enrichment of a salicylate glucoside biosynthetic pathway—implicated in SA storage and turnover—is congruent with the hypothesis that cold climates select for genotypes that finely balance defense signaling against growth-related energetic and developmental costs. The direction of these effects can be context-dependent, but the underlying regulatory axis is well supported by experimental evidence ([PMC][6]).

---

## Papers behind the (Author, year) mentions

* Fürtauer et al., 2019 — cold acclimation drives broad metabolic rewiring including sugar/carbon partitioning. ([PubMed][1])
* Murata et al., 2007 — low temperature enhances photoinhibition largely by inhibiting PSII repair. ([ScienceDirect][2])
* Bascuñán-Godoy et al., 2012 — cold acclimation can improve resistance/recovery from low-temperature photoinhibition. ([PMC][7])
* Demmig-Adams et al., 1996 — carotenoids/xanthophyll cycle as core photoprotection under environmental stress. ([PubMed][8])
* Gross et al., 2006 — phylloquinone biosynthesis mutants show reduced PSI activity. ([PubMed][9])
* Huner et al., 1998 — framework for energy balance/acclimation to light + cold (helps explain pigment/cofactor enrichment). ([ScienceDirect][5])
* Takahashi et al., 2019 — cold and sub-zero acclimation induce cell-wall remodeling linked to freezing tolerance. ([Nature][3])
* Kutsuno et al., 2022 — temporal cell-wall polymer changes track freezing tolerance during acclimation/deacclimation. ([PMC][10])
* Takahashi et al., 2024 — pectic polymer structural changes contribute to acquired freezing tolerance. ([PubMed][11])
* Rahman et al., 2021 — cuticular wax levels affect dehydration + low-temperature/freezing responses. ([PubMed][4])
* Hashida et al., 2009 — NAD biosynthesis/homeostasis in plant development and stress responses. ([PMC][12])
* Shibasaki et al., 2009 — cold inhibits auxin transport via effects on trafficking of auxin carriers. ([PMC][6])
* Scott et al., 2004 — SA accumulation inhibits growth at chilling temperature in *Arabidopsis*. ([PMC][13])

If you paste the **top genes driving** each pathway (or the leading-edge genes from your MAGMA/CATFISH outputs), I can help you separate “**likely cold-specific**” vs “**likely climate-correlated**” signals, especially for the small gene sets (n=2–3) where rankings can be high-variance.

[1]: https://pubmed.ncbi.nlm.nih.gov/31671650/ "Dynamics of Plant Metabolism during Cold Acclimation - PubMed"
[2]: https://www.sciencedirect.com/science/article/pii/S0005272806003665 "Photoinhibition of photosystem II under environmental stress - ScienceDirect"
[3]: https://www.nature.com/articles/s41598-019-38688-3?utm_source=chatgpt.com "Both cold and sub-zero acclimation induce cell wall ..."
[4]: https://pubmed.ncbi.nlm.nih.gov/33557073/ "Dissecting the Roles of Cuticular Wax in Plant Resistance to Shoot Dehydration and Low-Temperature Stress in Arabidopsis - PubMed"
[5]: https://www.sciencedirect.com/science/article/abs/pii/S1360138598012485 "Energy balance and acclimation to light and cold - ScienceDirect"
[6]: https://pmc.ncbi.nlm.nih.gov/articles/PMC2814496/ "
            Auxin Response in Arabidopsis under Cold Stress: Underlying Molecular Mechanisms - PMC
        "
[7]: https://pmc.ncbi.nlm.nih.gov/articles/PMC3490872/ "
            Cold-acclimation limits low temperature induced photoinhibition by promoting a higher photochemical quantum yield and a more effective PSII restoration in darkness in the Antarctic rather than the Andean ecotype of Colobanthus quitensis Kunt Bartl (Cariophyllaceae) - PMC
        "
[8]: https://pubmed.ncbi.nlm.nih.gov/8647339/ "Carotenoids 3: in vivo function of carotenoids in higher plants - PubMed"
[9]: https://pubmed.ncbi.nlm.nih.gov/16617180/ "A plant locus essential for phylloquinone (vitamin K1) biosynthesis originated from a fusion of four eubacterial genes - PubMed"
[10]: https://pmc.ncbi.nlm.nih.gov/articles/PMC10107845/?utm_source=chatgpt.com "Temporal cell wall changes during cold acclimation and ..."
[11]: https://pubmed.ncbi.nlm.nih.gov/38335960/?utm_source=chatgpt.com "Structural changes in cell wall pectic polymers contribute to ..."
[12]: https://pmc.ncbi.nlm.nih.gov/articles/PMC2707885/ "
            The role of NAD biosynthesis in plant development and stress responses - PMC
        "
[13]: https://pmc.ncbi.nlm.nih.gov/articles/PMC514138/ "
            Salicylate Accumulation Inhibits Growth at Chilling Temperature in Arabidopsis - PMC
        "



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
  


