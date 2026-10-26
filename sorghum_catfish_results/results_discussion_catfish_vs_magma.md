# Comparative Pathway Analysis: CATFISH Omnibus versus MAGMA Gene-Set Analysis

## Results

### Pathway-Level Comparison Reveals Distinct Enrichment Profiles

To evaluate the performance of CATFISH relative to standard pathway analysis approaches, we compared pathway enrichment results from CATFISH omnibus testing with MAGMA gene-set analysis for sorghum stem volume. The two methods identified fundamentally different pathway signatures (Fig. X). Of the top 10 pathways ranked by each method, there was no overlap (0/10 pathways shared), indicating that CATFISH and MAGMA capture distinct aspects of genetic architecture underlying stem volume variation.

CATFISH identified 75 pathways exceeding the Bonferroni-corrected significance threshold (P < 1.22 × 10⁻⁴), whereas MAGMA failed to detect any pathways at this stringency (Table 1). The Spearman rank correlation between pathway rankings from the two methods was weak (ρ = 0.233), further demonstrating their divergent prioritization of biological processes.

**Table 1.** Summary of pathway analysis results comparing CATFISH omnibus and MAGMA gene-set analysis.

| Metric | CATFISH Omnibus | MAGMA |
|--------|-----------------|-------|
| Significant pathways (Bonferroni) | 75 | 0 |
| Top pathway P-value | 5.96 × 10⁻³⁵ | 9.13 × 10⁻⁴ |
| Top 10 overlap | 0/10 | — |
| Rank correlation (Spearman ρ) | 0.233 | — |

The most significant pathway identified by CATFISH was aerobic respiration I (cytochrome c; P = 5.96 × 10⁻³⁵), comprising 298 genes. Additional top-ranked CATFISH pathways included xylan biosynthesis (P = 1.27 × 10⁻²¹), geranylgeranyldiphosphate biosynthesis I via mevalonate (P = 1.23 × 10⁻³²), and phosphatidylcholine biosynthesis III (P = 1.22 × 10⁻²⁹). In contrast, MAGMA ranked D-galactose detoxification (P = 9.13 × 10⁻⁴), putrescine biosynthesis III (P = 1.12 × 10⁻³), and D-mannose degradation (P = 1.21 × 10⁻³) as the top pathways, none of which achieved genome-wide significance.

### Gene-Level Comparison Shows Partial Agreement

Despite the complete divergence at the pathway level, gene-level prioritization showed moderate concordance between methods. Comparison of the top 10 candidate genes ranked by CATFISH multi-layer scoring versus MAGMA gene P-values revealed 5 genes in common (50% overlap). The top-ranked gene by CATFISH scoring was SORBI_3009G238700 (score = 9.29), which achieved significance across all three evidence layers: GWAS association, MAGMA gene-level analysis, and membership in the top-ranked aerobic respiration pathway. Notably, the majority of top-ranked genes by both methods were located on chromosome 9, suggesting the presence of a major quantitative trait locus (QTL) influencing stem volume in this genomic region.

## Discussion

### CATFISH Detects Biologically Relevant Pathways Missed by Competitive Testing

The striking discordance between CATFISH and MAGMA pathway results reflects fundamental differences in their statistical frameworks. MAGMA employs a competitive test that asks whether genes in a pathway show stronger associations than genes outside the pathway. This approach is inherently conservative and may fail to detect pathways where the effect is driven by a subset of genes with strong signals—a pattern common in complex trait architecture. CATFISH, by contrast, combines multiple P-value aggregation methods (ACAT, Fisher, adaptive TFisher, minP, and Stouffer) through an omnibus framework, enabling detection of enrichment across diverse signal patterns including sparse driver effects, coordinated moderate signals, and diffuse polygenic contributions.

The pathways prioritized by CATFISH demonstrate clear biological relevance to stem volume and plant growth. Aerobic respiration provides the cellular energy required for biomass accumulation and stem elongation. Xylan biosynthesis directly contributes to secondary cell wall formation, a primary determinant of stem structural properties and lignocellulosic biomass composition in grasses (Rennie and Scheller, 2014). Geranylgeranyldiphosphate serves as a precursor for gibberellins, chlorophylls, and carotenoids—compounds essential for stem elongation and photosynthetic capacity (Ruiz-Sola and Rodriguez-Concepcion, 2012). Phosphatidylcholine and choline biosynthesis pathways support membrane biogenesis required for cell expansion during stem growth.

In contrast, the pathways ranked highest by MAGMA—D-galactose detoxification, putrescine biosynthesis, and D-mannose degradation—lack obvious mechanistic connections to stem volume determination. While these pathways may contribute to general metabolic homeostasis, their prioritization likely reflects statistical artifacts of the competitive testing framework rather than true biological enrichment.

### Multi-Layer Scoring Integrates Complementary Evidence

The partial overlap in top-ranked genes (50%) between CATFISH and MAGMA indicates that strong gene-level associations are detected by both approaches, but pathway context provides additional discriminatory power. The CATFISH multi-layer scoring framework, which integrates GWAS significance, MAGMA gene P-values, and pathway membership, identified candidate genes with convergent evidence across multiple analytical scales. The top candidate SORBI_3009G238700, a member of the aerobic respiration pathway, exemplifies how pathway-level analysis can elevate genes that might otherwise be overlooked based solely on single-gene association statistics.

The concentration of top-ranked genes on chromosome 9 by both methods suggests a major effect locus for stem volume in sorghum. This genomic region warrants further investigation through fine-mapping and functional validation approaches. The identification of pathway-supported candidates within this QTL region demonstrates the utility of integrating pathway analysis with positional information for candidate gene prioritization.

### Methodological Implications

These results highlight the complementary nature of different pathway analysis strategies and caution against relying solely on competitive gene-set tests for post-GWAS interpretation. The omnibus combination approach implemented in CATFISH, coupled with GPD tail extrapolation for precise P-value estimation, provides substantially greater power to detect pathway enrichment across diverse genetic architectures. For complex traits influenced by multiple genes with varying effect sizes, methods that accommodate heterogeneous signal patterns—rather than assuming uniform enrichment—may better capture the underlying biology.

The absence of significant MAGMA pathways despite 75 significant CATFISH pathways underscores the statistical power differential between these approaches. Researchers conducting pathway analysis of GWAS results should consider employing multiple complementary methods and evaluating biological plausibility of identified pathways rather than accepting results from any single approach uncritically.

## References

Rennie EA, Scheller HV (2014) Xylan biosynthesis. Curr Opin Biotechnol 26: 100-107

Ruiz-Sola MA, Rodriguez-Concepcion M (2012) Carotenoid biosynthesis in Arabidopsis: a colorful pathway. Arabidopsis Book 10: e0158
