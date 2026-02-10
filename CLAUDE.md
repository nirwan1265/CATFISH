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
