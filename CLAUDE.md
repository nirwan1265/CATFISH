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
