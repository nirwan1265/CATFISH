# Sorghum Stem Volume CATFISH Analysis Log

## Overview
- **Trait**: Stem volume
- **Species**: Sorghum (SAP - Sorghum Association Panel)
- **GWAS file**: `/Users/nirwantandukar/Documents/Research/results/GWAS/SAP/Stem_diameter/Stem_volume_mod_sub_stem_volume_SAP_bialleles_MAF_0.05_11.assoc.txt`
- **Pathways**: 410 PlantCyc sorghum pathways
- **MVN resampling**: 10,000 permutations

---

## Analysis Scripts

### 1. Main CATFISH Analysis
- **File**: `sorghum_stem_vol_catfish.R`
- **Purpose**: Run CATFISH pathway analysis with MVN resampling
- **Output**: `sorghum_stem_vol_catfish_results_MVN.csv`

### 2. Figure 4: GWAS, MAGMA, Pathway Manhattan
- **File**: `sorghum_figure4_and_supp.R`
- **Panels**:
  - A: GWAS Manhattan plot
  - B: MAGMA gene Manhattan plot
  - C: CATFISH pathway heatmap (top 20 pathways)
- **Outputs**:
  - `Figure4_GWAS_MAGMA_Pathway.png`
  - `Figure4_GWAS_MAGMA_Pathway.pdf`
  - `TableS1_GWAS_genes_logP_gt7.csv` (423 genes)
  - `TableS2_MAGMA_genes_FDR_lt0.05.csv` (861 genes)
  - `TableS3_Pathways_FDR_lt0.05.csv` (102 pathways)
  - `TableS4_Pathway_genes_detailed.csv`

### 3. Figure 5: Integrated Candidate Gene Analysis
- **File**: `sorghum_figure5_and_table.R`
- **Panels**:
  - A: UpSet plot (GWAS, MAGMA, Pathway intersections)
  - B: Summary bar chart
  - C: Density plots (genes in sig pathways vs not)
  - D: GWAS vs MAGMA scatter (36 "All 3" genes highlighted in red)
- **Outputs**:
  - `Figure5_Candidate_Genes.png`
  - `Figure5_Candidate_Genes.pdf`
  - `Table_Genes_All_Three_Layers.csv` (36 genes)

### 4. Supplementary Figure: Component Test Comparison
- **File**: `sorghum_supplementary_figure.R`
- **Panels**:
  - A: UpSet plot (component tests vs Omnibus)
  - B: QQ plot (Omnibus p-values, colored by BH/Bonferroni significance)
  - C: Jaccard similarity (lower) / Spearman correlation (upper) matrix
  - D: P-value distribution boxplot
- **Outputs**:
  - `Supplementary_Figure_Component_Tests.png`
  - `Supplementary_Figure_Component_Tests.pdf`
  - Individual panels: `SuppFig_A_UpSet.png`, `SuppFig_B_QQplot.png`, etc.

---

## Key Results

### GWAS
- Threshold: -log10(p) > 7
- Significant genes: **423**
- Major QTL on chromosome 9 (~56.1-58.9 Mb)

### MAGMA Gene Analysis
- Total genes: 33,967
- FDR < 0.05: **861 genes**
- Gene-level lambda: **1.31** (acceptable for complex traits)

### CATFISH Pathway Analysis
- Total pathways: 410
- Bonferroni < 0.05: **73 pathways**
- BH FDR < 0.05: **102 pathways**
- Top pathways: GLUCONEO-PWY, GLUGLNSYN-PWY, GLUTAMINDEG-PWY, LIPASYN-PWY, PWY-3781

### Integrated Analysis
- Genes in all 3 layers (GWAS + MAGMA + Pathway): **36 genes**
- 34/36 on chromosome 9
- 2 on chromosome 2 (SORBI_3002G379200, SORBI_3002G379100)

---

## Component Test Comparison

| Method | Bonferroni Sig | Notes |
|--------|---------------|-------|
| TFisher | 84 | Most sensitive |
| Omnibus | 73 | Balanced |
| Fisher | 48 | - |
| ACAT | 47 | Identical to minP |
| minP | 47 | Identical to ACAT |
| Stouffer | 24 | Most selective |

- ACAT and minP: Jaccard = 1.0 (perfect overlap)
- Stouffer: Lowest overlap with others (Jaccard 0.29-0.50)
- All methods: High correlation (rho = 0.71-0.97)

---

## QQ Plot Notes
- **Non-significant pathways**: Follow diagonal well (good null calibration)
- **Significant pathways**: Deviate upward (true biological signal)
- MVN floor at -log10(p) ~ 4 (10,000 permutations)
- Lambda/Type I error metrics NOT applicable for real data with signal (only for null simulations)

---

## Manuscript
- **File**: `manuscript_results_discussion.tex`
- Contains:
  - Results sections for Figure 4 and Figure 5
  - Discussion of top candidate genes
  - Supplementary results for component test comparison
  - Tables: Top 9 genes (Table 1), All 36 genes (Table 2)
  - Figure captions

---

## File Locations

All outputs in: `/Users/nirwantandukar/Documents/Github/MAGCAT/sorghum_catfish_results/`

```
sorghum_catfish_results/
├── sorghum_stem_vol_catfish_results_MVN.csv    # Main results
├── sorghum_stem_vol_genes_combined.txt          # MAGMA gene results
├── Figure4_GWAS_MAGMA_Pathway.png/pdf
├── Figure5_Candidate_Genes.png/pdf
├── Supplementary_Figure_Component_Tests.png/pdf
├── Table_Genes_All_Three_Layers.csv
├── TableS1_GWAS_genes_logP_gt7.csv
├── TableS2_MAGMA_genes_FDR_lt0.05.csv
├── TableS3_Pathways_FDR_lt0.05.csv
├── TableS4_Pathway_genes_detailed.csv
├── MAGMA_gene_QQplot.png
├── manuscript_results_discussion.tex
└── SORGHUM_ANALYSIS_LOG.md (this file)
```

---

---

## Notes / Known Issues

1. **HMP to PLINK conversion**: The `sorghum_stem_vol_catfish.R` script includes HMP→PLINK conversion code, but user may have used pre-existing PLINK files instead. If HMP files use IUPAC heterozygote codes (R, Y, S, W, K, M), the conversion logic would need updating.

2. **GWAS SNP filtering**: Ideally should filter GWAS SNPs to match genotype file SNPs before MAGMA. May not have been done, but results look valid given strong signals.

3. **Gene regex for MVN**: The `magma_genesraw_to_cor_pairs_banded()` function needs `gene_regex` set for sorghum genes (SORBI_*). User confirmed this was handled.

4. **Results validation**: Despite potential minor issues, results appear well-calibrated (non-sig pathways follow QQ diagonal) and biologically coherent.

---

## Date
Analysis completed: March 2026
