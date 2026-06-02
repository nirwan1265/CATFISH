# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

CATFISH (Combining **C**auchy combination test (ACAT), **A**daptive TFisher, **F**isher, m**I**n-P, and **S**touffer for **H**olistic pathway analysis) is an R package for post-GWAS pathway enrichment analysis. It combines multiple p-value combination methods with LD-aware MAGMA gene-level statistics to detect pathway enrichment across different signal patterns.

## Build and Development Commands

```r
# Load package for development
devtools::load_all(".")

# Generate documentation (roxygen2)
devtools::document()

# Run tests
devtools::test()

# Run a single test file
testthat::test_file("tests/testthat/test-acat.R")

# Full package check
devtools::check()

# Install locally
devtools::install()

# Build vignettes
devtools::build_vignettes()
```

## External Dependencies

**MAGMA** must be installed separately (not an R package). Download from https://ctg.cncr.nl/software/magma

Set path via: `options(magma.path = "/path/to/magma")` or add to system PATH.

## Architecture

### Core Components

The package implements a multi-method pathway testing framework:

1. **Component Tests** (`R/*_wrappers.R`): Individual pathway tests
   - `catfish_acat_pathways()` - Cauchy combination (sparse signals)
   - `catfish_fisher_pathways()` - Fisher's method (moderate signals)
   - `catfish_minp_pathways()` - Minimum p-value (single driver)
   - `catfish_stoufferZ_pathways()` - Z-score aggregation (diffuse signals)
   - `catfish_soft_tfisher_adaptive_pathways()` - Truncated Fisher (hybrid patterns)

2. **Omnibus Integration** (`R/OMNIBUS_test.R`): Combines component tests
   - `catfish_omni2_pathways()` - Main entry point
   - Supports analytic, global resampling, and MVN resampling calibration

3. **MAGMA Interface** (`R/magma_wrappers.R`, `R/prepare_wrappers.R`):
   - `magma_gene()` - Run MAGMA gene analysis
   - `magma_genesraw_to_cor_pairs_banded()` - Extract gene correlations

4. **Data Preparation** (`R/gff3_to_geneloc.R`, `R/get_gene_lengths.R`):
   - `gff3_to_geneloc()` - Convert GFF3 to MAGMA format
   - `catfish_adjust_gene_p()` - Adjust for gene length/SNP density

### Pathway Signal Archetypes

The package is designed around different enrichment patterns:
- **SDA (Sparse Driver)**: Few extreme signals → ACAT, minP
- **CME (Coordinated Moderate)**: Many moderate signals → Fisher, Stouffer
- **DPS (Diffuse Polygenic)**: Subtle genome-wide shift → Stouffer
- **HDS (Hybrid Driver-Support)**: Mixed pattern → Soft TFisher, Omnibus

### Built-in Data

Located in `inst/extdata/`:
- `pathway/` - PMN pathways for maize, sorghum, arabidopsis, plant, drosophila
- `*.genes.loc` - Pre-computed gene location files

Load via: `catfish_load_pathways(species = "maize")`

### Gene Results Format

Gene-level results require columns: `GENE`, `P`, and optionally `ZSTAT`, `NSNPS`, `CHR`, `START`, `STOP`

## Key Files

- `R/OMNIBUS_test.R` - Core omnibus logic (~1,850 lines)
- `R/magma_wrappers.R` - MAGMA interface (~1,200 lines)
- `R/prepare_wrappers.R` - Data validation/preparation (~1,600 lines)
- `workflow_code.md` - Step-by-step usage guide
- `README.md` - Comprehensive documentation with methodology

## Testing

Tests use testthat edition 3. Test files in `tests/testthat/` cover each combination method. Many tests skip if optional packages (ACAT, TFisher, etc.) are not installed.

## Code Conventions

- Exported functions prefixed with `catfish_` or `magma_`
- All functions documented with roxygen2
- Gene ID matching is case-insensitive
- Pathways accepted as named list or data.frame with `pathway_id`, `gene_id`, `pathway_name` columns

## Recent Implementations (March 2026)

### GPD Tail Extrapolation

Implemented Knijnenburg et al. 2009 Generalized Pareto Distribution (GPD) tail extrapolation to overcome permutation floor limitations:

**Key Parameters in `catfish_omni2_pathways()`:**
- `tail_mode = "hybrid_gpd"` - Enable GPD tail extrapolation (vs `"empirical"` for pure Monte Carlo)
- `tail_switch_exceed = 10L` - Switch to GPD when exceedances fall below this threshold
- `tail_gpd_k = 250L` - Number of extreme values to fit GPD
- `tail_min_B = 10000L` - Minimum B required for GPD
- `tail_min_tail = 50L` - Minimum tail observations for stable GPD fit

**Why GPD?**
- Empirical p-values have floor of ~1/(B+1)
- With B=10,000, minimum p-value ≈ 1e-4
- GPD extrapolates beyond this floor by fitting extreme value distribution to tail
- Can achieve p-values of 10^-35 or smaller with B=1,000,000

### Parallel Processing

Added multi-threaded parallelization for inner null computations:

**Key Parameters:**
- `n_threads = NULL` - Number of threads (default: n_cores - 1)
- Uses `parallel::mclapply()` for Unix-based systems
- Parallelizes ACAT null, TFisher null, and omnibus null computations
- Only activates when B >= 1000 to avoid overhead

**Performance:**
- Near-linear speedup with number of cores
- Example: 50+ hours single-threaded → 2.4 hours with 12 threads for B=1M

### Progress Messages

Added progress tracking for long-running jobs:
- Shows `[MVN] Pathway X/Y: pathway_id` every 10 pathways
- Helps monitor job progress during large analyses

## Analysis Results

### Sorghum Stem Volume (B=1,000,000 GPD)

**Results Location:** `sorghum_catfish_results/`
- `sorghum_stem_vol_CATFISH_B1000000_GPD.csv` - Pathway results
- `candidate_genes_GPD_B1M_scored.csv` - Scored candidate genes
- `candidate_genes_top200_GPD_B1M.csv` - Top 200 candidates

**Key Findings:**
- Top pathway: aerobic respiration I (cytochrome c), p = 5.96e-35
- 75 pathways passed Bonferroni threshold
- Top gene: SORBI_3009G238700 (score 9.29, 3-layer support)
- 21 genes with 3-layer support (GWAS + MAGMA + pathway)

### Arabidopsis Cold Trait (B=1,000,000 GPD)

**Results Location:** `arabidopsis_catfish_results/`
- `arabidopsis_cold_CATFISH_B1000000_GPD.csv` - Pathway results
- `candidate_genes_GPD_B1M_scored_arabidopsis.csv` - Scored candidate genes
- `candidate_genes_top200_GPD_B1M_arabidopsis.csv` - Top 200 candidates

**Key Findings:**
- Top pathway: wax esters biosynthesis I, p = 7.70e-05
- Only 1 pathway passed Bonferroni threshold (weaker signals than sorghum)
- Top gene: AT1G74458 (score 8.90)
- 30 genes in significant pathway

## Example: Running GPD Analysis with Parallel Processing

```r
library(devtools)
load_all(".")

res <- catfish_omni2_pathways(
  gene_results = gene_results,
  pathways = pathways,
  gene_col = "GENE",
  p_raw_col = "P",
  z_col = "ZSTAT",
  B_perm = 1000000,           # 1 million permutations
  perm_mode = "mvn",
  magma_cor_pairs = cor_pairs,
  tail_mode = "hybrid_gpd",   # Enable GPD extrapolation
  tail_switch_exceed = 10L,
  tail_gpd_k = 250L,
  tail_min_B = 10000L,
  tail_min_tail = 50L,
  min_p = 1e-50,
  n_threads = 12              # Use 12 cores for parallelization
)
```

## Scoring Workflow

Multi-layer candidate gene scoring combines:
1. **GWAS layer**: Gene-level GWAS p-value rank
2. **MAGMA layer**: MAGMA gene p-value (Bonferroni threshold)
3. **Pathway layer**: Membership in significant pathways (GPD p-values)

Score formula:
```
score = -log10(gwas_rank/n_genes) - log10(magma_rank/n_genes) + pathway_bonus
pathway_bonus = -log10(best_pathway_p) / 10  # scaled down
```

See `sorghum_catfish_results/calculate_scores_GPD.R` for implementation.

## Computational Benchmarks

Block I benchmark added to `all_figs.R` for paper figure showing:
- Runtime vs B (permutation count) - linear scaling
- Runtime vs pathway size - linear scaling
- Memory usage remains stable
- 12 threads provides ~10x speedup over single-threaded

**Config:**
- Pathway sizes: 10, 25, 50, 100 genes
- B values: 1000, 5000, 10000, 50000, 100000
- 3 replicates per condition
- Results saved to `simulation_results/block_i_benchmark.rds`

## Manuscript Review (CATFISH.pdf) - Issues to Fix

### Critical Issues

1. **Missing/Placeholder References**: Multiple "(ref)" citations need actual references
   - Section 2.2: "The tail can be approximated (ref)"
   - Section 2.1.3: "P = 2×Φ(−|z|) (ref)"
   - Throughout Supplementary Materials

2. **Raw BibTeX in Text**: Section A2.9 has unrendered BibTeX entries like `\citep{gholami.2022}`, `\citet{zhu.2024}` that need fixing

3. **Inconsistent Archetype Naming**:
   - Main text uses "DPE" (Diffuse Polygenic Enrichment)
   - CLAUDE.md and some sections use "DPS" (Diffuse Polygenic Shift)
   - Need to standardize to one naming convention

4. **Outdated Results**:
   - Manuscript uses B=10,000 results (p-value floor ~10^-4)
   - Should update to B=1,000,000 GPD results (p-values to 10^-35)
   - This significantly strengthens the biological findings

### Sections Needing Completion

1. **Graphical Abstract**: Currently empty placeholder
2. **Data Availability**: Section exists but needs content
3. **Code Availability**: Section exists but needs content
4. **Author Contributions**: Section exists but needs content
5. **Declaration of Competing Interest**: Needs actual statement

### Suggestions for Improvement

1. Add comparison table: B=10k empirical vs B=1M GPD results
2. Include computational benchmark figure (Block I from all_figs.R)
3. Add runtime/scalability discussion for large-scale analyses
4. Consider adding Arabidopsis case study as validation

### File Locations
- Manuscript: `/Users/nirwantandukar/Downloads/CATFISH.pdf`
- Results to update from: `sorghum_catfish_results/sorghum_stem_vol_CATFISH_B1000000_GPD.csv`

## TFisher Tau Parameter Tuning

### Problem
With the default liberal `tau_grid = c(0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.001)`, TFisher was dominating the omnibus results:
- TFisher found 84 significant pathways while other methods found ~47-48
- TFisher dominated 100% of Bonferroni-significant pathways
- The omnibus was essentially just TFisher

### Understanding Tau
- **tau** is the threshold for including p-values in TFisher combination
- `tau = 0.2` is **liberal** (includes genes with p < 0.2)
- `tau = 1e-7` is **strict** (only includes genes with p < 0.0000001)
- For data with extreme signals (very small p-values), use strict tau

### Solution: Strict Tau Grid
Changed to `tau_grid = c(1e-5, 1e-6, 1e-7)` for sorghum (strong signals):
- TFisher significant: 84 → 29
- TFisher dominance: 100% → 56%
- Other methods (Fisher, ACAT) now contribute to significant pathways
- Omnibus: 75 → 52 (more conservative but balanced)

### Data-Dependent Behavior
- **Sorghum** (strong signals): TFisher dominates with liberal tau, use strict tau
- **Arabidopsis** (weak signals): TFisher gives p=0.73 with strict tau (no genes < 1e-5), Stouffer dominates

### Configuration
All scripts now have a `TAU_OPTION` at the top:
```r
TAU_OPTION <- "strict"  # "default" or "strict"
# default = c(0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.001)
# strict  = c(1e-5, 1e-6, 1e-7)
```

## Current State (March 2026)

### Completed Analyses
- ✅ Sorghum stem volume: B=1M GPD, strict tau, 12 threads (~1 hour, 410 pathways)
- ✅ Arabidopsis cold trait: B=1M GPD, 12 threads (1.86 hours, 321 pathways)
- ✅ Candidate gene scoring for both species
- ✅ Computational benchmark (Block I in all_figs.R)
- ✅ TFisher tau tuning investigation

### Sorghum Results (Strict Tau)
**File:** `sorghum_stem_vol_CATFISH_B1000000_GPD_strict_tau.csv`

| Method | Bonferroni Significant |
|--------|------------------------|
| ACAT | 47 |
| Fisher | 48 |
| TFisher | 29 |
| minP | 47 |
| Stouffer | 24 |
| Omnibus | 52 |

**Dominant component (Bonferroni significant):**
- TFisher: 29 (56%)
- Fisher: 13 (25%)
- ACAT: 7 (13%)
- Stouffer: 2 (4%)
- minP: 1 (2%)

**Candidate genes:**
- 21 genes with 3-layer support
- 186 genes with 2-layer support
- Top gene: SORBI_3009G231300 (score 9.06)

### Key Scripts
- `sorghum_catfish_results/run_catfish_GPD_1M.R` - Sorghum analysis (TAU_OPTION configurable)
- `sorghum_catfish_results/calculate_scores_GPD.R` - Sorghum scoring (TAU_OPTION configurable)
- `sorghum_figure4_and_supp.R` - Figure 4 + supplementary tables (TAU_OPTION configurable)
- `sorghum_figure5_and_table.R` - Figure 5 candidate genes (TAU_OPTION configurable)
- `sorghum_supplementary_figure.R` - Supplementary figure (TAU_OPTION configurable)
- `arabidopsis_catfish_results/run_catfish_GPD_1M_arabidopsis.R` - Arabidopsis analysis
- `all_figs.R` - Block I benchmark code

### Output Files (Strict Tau)
- `sorghum_stem_vol_CATFISH_B1000000_GPD_strict_tau.csv` - Pathway results
- `candidate_genes_GPD_B1M_scored_strict_tau.csv` - All scored genes
- `candidate_genes_top200_GPD_B1M_strict_tau.csv` - Top 200 candidates
- `Figure4_GWAS_MAGMA_Pathway_strict_tau.png/pdf`
- `Figure5_Candidate_Genes_strict_tau.png/pdf`
- `Supplementary_Figure_Component_Tests_strict_tau.png/pdf`
- `TableS3_Pathways_FDR_lt0.05_strict_tau.csv` - 96 pathways FDR < 0.05
- `TableS4_Pathway_genes_detailed_strict_tau.csv` - Gene-pathway associations
- `Table_Genes_All_Three_Layers_strict_tau.csv` - 35 genes with 3-layer support

### Notes
- B=5M was attempted for sorghum but killed - GPD with B=1M already achieves 10^-35 p-values
- Higher B provides diminishing returns when using GPD tail extrapolation
- Recommended: B=1M with GPD for production analyses
- Use strict tau for data with extreme gene-level signals

---

## DGRP Chili Pepper Project (April 2026)

### Overview

GWAS analysis of *Drosophila melanogaster* DGRP lines fed different chili peppers.

**Location:** `DGRP_chilipeppers/`

### Experimental Design

- **69 DGRP lines** (inbred, genotypes publicly available)
- **4 treatments:**
  - C = Control
  - B = Bell pepper (non-spicy)
  - S = Serrano (medium spicy)
  - H = Habanero (very spicy)
- **Phenotypes:**
  - Body Weight (BW) - measured separately for Males (M) and Females (F)
  - Triglycerides (TG) - sexes pooled
- **3 replicates per line/treatment/sex**

### Current Progress

**Step 1: BLUPs computed** ✅
- Script: `01_compute_BLUPs.R`
- Output: `BLUPs_all_traits.csv` (21 traits)
- Individual trait files in `blups_for_gwas/` folder (ready for DGRP web tool)

**Traits computed:**
- 8 BW traits: BW_C_F, BW_C_M, BW_B_F, BW_B_M, BW_S_F, BW_S_M, BW_H_F, BW_H_M
- 4 TG traits: TG_C, TG_B, TG_S, TG_H
- 9 derived traits:
  - BW_spicy_F/M, BW_nonspicy_F/M (averages)
  - BW_capsaicin_F/M = spicy - nonspicy (capsaicin effect contrast)
  - TG_spicy, TG_nonspicy, TG_capsaicin

**Heritability (Line variance %):**
| Trait | Heritability |
|-------|-------------|
| BW_S_F | 91.4% |
| BW_H_F | 87.6% |
| TG_H | 85.4% |
| TG_B | 84.5% |
| BW_H_M | 58.5% |

### Next Steps

**Step 2: Run GWAS**
- Option A: DGRP2 web tool (http://dgrp2.gnets.ncsu.edu/) - upload trait files
- Option B: Local GWAS with DGRP VCF + PLINK/GEMMA

**Step 3: Post-GWAS analysis**
- Manhattan plots
- Gene annotation
- Compare hits across treatments (spicy vs non-spicy)
- Potentially run CATFISH pathway analysis

### Key Questions

1. Do hot pepper treatments alter genetic architecture of TG/BW?
2. Are there capsaicin-responsive loci (GxE)?
3. Sex differences in response to spicy food?

### Files

```
DGRP_chilipeppers/
├── 01_compute_BLUPs.R          # BLUP computation script
├── phenotype.csv               # Raw phenotype data
├── BLUPs_all_traits.csv        # All 21 trait BLUPs
├── BLUP_summary_stats.csv      # Trait summaries
└── blups_for_gwas/             # Individual trait files for GWAS
    ├── BW_*.txt                # Body weight traits
    ├── TG_*.txt                # Triglyceride traits
    └── *_capsaicin.txt         # Capsaicin contrast traits
```

---

## Arabidopsis TuMV (Turnip Mosaic Virus) Project (April 2026)

### Overview

GWAS and CATFISH pathway analysis of *Arabidopsis thaliana* response to Turnip Mosaic Virus infection.

**Location:** `/Users/nirwantandukar/Documents/Research/results/GWAS/CATFISH/arabidopsis/`

### Phenotypes

- **TuMV_G_14** - Viral symptom severity at day 14 (G = symptom grade)
- **TuMV_G_21** - Viral symptom severity at day 21
- **TuMV_S_14** - (not analyzed yet)
- **TuMV_S_21** - (not analyzed yet)

### Analysis Completed

**MAGMA gene analysis** ✅
- 5 chromosomes per phenotype
- 28,435 genes per trait
- Gene correlation pairs extracted for MVN resampling

**CATFISH pathway analysis** ✅
- B = 100,000 permutations
- 264 AraCyc pathways tested
- 8 threads, GPD tail extrapolation

### Results Summary

| Trait | FDR < 0.05 | Bonferroni | MAGMA-sig genes |
|-------|------------|------------|-----------------|
| TuMV_G_14 | 13 pathways | 11 pathways | 8 genes |
| TuMV_G_21 | 2 pathways | 2 pathways | 4 genes |

**Top pathways (TuMV_G_14):**
1. cyanide detoxification I (p = 4.93e-17)
2. L-cysteine biosynthesis I (p = 1.09e-09)
3. volatile benzenoid biosynthesis I (p = 1.30e-08)
4. guanosine nucleotides degradation II (p = 2.74e-07)

**Top pathways (TuMV_G_21):**
1. choline biosynthesis III (p = 3.35e-06)
2. L-cysteine biosynthesis I (p = 2.51e-05)

### Candidate Genes

All MAGMA-significant genes are in a tight cluster on **chromosome 2**:

| Gene | TuMV_G_14 P-value | TuMV_G_21 P-value |
|------|-------------------|-------------------|
| AT2G14095 | 2.44e-14 | 1.10e-09 |
| AT2G14080 | 2.49e-14 | 4.39e-09 |
| AT2G14100 | 5.34e-14 | 3.27e-10 |
| AT2G14070 | 1.20e-12 | 1.69e-08 |
| AT2G14110 | 7.56e-09 | - |
| AT2G14115 | 1.27e-08 | - |
| AT2G14120 | 8.34e-08 | - |
| AT2G14060 | 5.25e-07 | - |

This cluster suggests a single causal locus with LD spreading across adjacent genes.

### Files

```
arabidopsis/
├── run_catfish_TuMV.R                    # Main analysis script
├── figure_TuMV_validation.R              # Figure generation script
├── TuMV_validation_summary.txt           # Validation story write-up
├── TuMV_G_14_validation_figure.png/pdf   # 4-panel validation figure
├── TuMV_G_14_CATFISH_B1e+05.csv          # Pathway results (264 pathways)
├── TuMV_G_14_candidate_genes_top200.csv  # Top 200 candidate genes
├── TuMV_G_14_candidate_genes_scored.csv  # All scored genes
├── TuMV_G_21_CATFISH_B1e+05.csv
├── TuMV_G_21_candidate_genes_top200.csv
├── TuMV_G_21_candidate_genes_scored.csv
├── TuMV_*.txt                            # GWAS summary stats
├── validate_arabidopsis/                 # Validation papers (PDFs)
└── magma_genes/
    ├── TuMV_G_14_chr*.genes.out          # MAGMA output per chromosome
    ├── TuMV_G_14_genes_combined.txt      # Combined gene results
    ├── TuMV_G_14_gene_cor_pairs.txt      # Gene correlations for MVN
    ├── TuMV_G_21_chr*.genes.out
    ├── TuMV_G_21_genes_combined.txt
    └── TuMV_G_21_gene_cor_pairs.txt
```

### Independent Validation

**Key finding:** Our candidate genes were independently validated by a separate paper!

- **Our GWAS:** Used symptom severity phenotypes (TuMV_G_14, TuMV_G_21)
- **Validation paper:** "Genetic basis of Arabidopsis thaliana responses to infection by naïve and adapted isolates of turnip mosaic virus"
- **Result:** Validated 2 genes from our AT2G14xxx cluster using different phenotypes

This is strong evidence because:
1. Different phenotypes → same causal region
2. Same virus system (TuMV)
3. Same 1001 Genomes resource
4. Consistent across timepoints (day 14 and day 21)

### Figures Generated

- `TuMV_G_14_validation_figure.png/pdf` - 4-panel figure:
  - A. GWAS Manhattan (3.3M SNPs)
  - B. MAGMA Manhattan (validated genes in red)
  - C. Top 15 pathways bubble plot
  - D. Top 20 candidate genes with scores
- `TuMV_validation_summary.txt` - Full validation write-up

### Notes

- Cyanide/cysteine pathways are classic plant defense responses - biologically relevant for viral infection
- Gene-pathway matching returned 0 matches due to ID mismatch (AraCyc uses gene names like "GAMT2", MAGMA uses "AT" locus IDs)
- Could re-run with `gene_col = "Gene-id"` in pathway loading to fix this
- Day 14 shows stronger signal than day 21 (more pathways, more genes)

---

## Maize Nitrogen CATFISH Analysis (May 2026)

### Overview

CATFISH pathway analysis of maize nitrogen use efficiency.

**Location:** `/Users/nirwantandukar/Documents/Research/results/GWAS/CATFISH/maize_nitrogen/`

**MAGMA results:** `/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_multi_snp_wise_genes_by_chr_N_maize/`

### Key Details

- **Sample size:** N = 3,107
- **Genes tested:** 34,873
- **Pathways tested:** 427 (CornCyc)
- **Permutations:** B = 100,000 with GPD tail extrapolation

### 3-Layer Scoring Pipeline

This analysis established the full 3-layer scoring methodology:

**Layer 1 - GWAS:** Uses `P_SNPWISE_TOP1` (best SNP p-value in gene)
**Layer 2 - MAGMA:** Uses `P` (gene-level multi-SNP model)
**Layer 3 - Pathway:** Membership in significant pathways (FDR < 0.10)

**Scoring Formula:**
```
score = -log10(gwas_rank/n) + -log10(magma_rank/n) + pathway_bonus
```

**Pathway Bonus (Multi-Pathway Boost):**
```
base_bonus = 0.1 × log2(n_all_pathways + 1)      # Any pathway
sig_bonus = (-log10(best_p)/10) × log2(n_sig_pathways + 1)  # Significant pathways
```

This rewards:
1. Genes in ANY annotated pathway (biological function known)
2. Genes in SIGNIFICANT pathways (trait-relevant)
3. Genes in MULTIPLE pathways (hub genes, stronger evidence)

### Results

**Top Pathways:**
| Pathway | P-value |
|---------|---------|
| Glutathione biosynthesis | 8.59e-08 |
| Glycine biosynthesis III | 1.14e-06 |
| L-homoserine biosynthesis | 2.02e-06 |

**Layer Support:**
- 3-layer: 1 gene (Zm00001eb220040)
- 2-layer: 4 genes
- 1-layer: 182 genes

**Top Candidate:** Zm00001eb220040 (Chr 5, score 8.84)
- 3-layer support (GWAS + MAGMA + Pathway)
- In 2 significant pathways

### Files

```
maize_nitrogen/
├── README.md                              # Full documentation
├── run_catfish_maize_N.R                  # CATFISH analysis script
├── calculate_scores_3layer.R             # 3-layer scoring script
├── maize_N_CATFISH_B1e+05_GPD.csv         # Pathway results
├── maize_N_candidate_genes_3layer_scored.csv  # All genes scored
├── maize_N_candidate_genes_3layer_top200.csv  # Top 200
├── maize_N_genes_in_sig_pathways.csv      # Genes in sig pathways
└── magma_gene_cor_pairs_maize_N.txt       # Gene correlations
```

### Notes

- Signals are weaker than sorghum/arabidopsis (no Bonferroni-significant genes)
- Used relaxed thresholds: GWAS < 1e-5, MAGMA < 1e-4
- Top candidates cluster on chromosome 5 (22-141 Mb region)
- Used default tau grid (appropriate for weaker signals)
