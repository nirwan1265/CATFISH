# Spike-In Simulation Pipeline

This folder adds a semi-synthetic spike-in pipeline on top of the existing sorghum phenotype-permutation workflow.

The key idea is:

1. Start from the real BAP genotypes.
2. Build a null phenotype by permuting the observed dry-ttons phenotype among non-missing individuals.
3. Inject pathway-structured genetic signal into real genotype dosages for selected pathway genes.
4. Re-run the full `GWAS -> MAGMA -> CATFISH` pipeline.

This keeps real LD, MAF spectrum, SNP density, gene length, and pathway membership, while giving us known truth.

## Files

- `config.spikein.template.sh`
  Template config. Copy to `config.spikein.sh` on HPC and fill in the real paths.
- `01_spike_setup.R`
  Staged spike step 1. Chooses one real pathway and writes `causal_snps.txt` and `causal_map.tsv`.
- `01_spike_setup.sh`
  LSF wrapper for step 1.
- `02_spike_extract_geno.sh`
  Staged spike step 2. Extracts the chosen causal SNP dosages into `SPIKE_WORK/causal_geno.raw`.
- `03_make_spikein.R`
  Staged spike step 3. Builds synthetic phenotypes from the extracted causal dosage matrix.
- `03_make_spikein.sh`
  LSF wrapper for step 3.
- `04_gwas_array_LSF.sh`
  Runs `vcf2gwas` for each spiked phenotype x chromosome.
- `05_magma_array_LSF.sh`
  Runs MAGMA gene analysis via `05_magma_catfish.R`.
- `06_corpairs_once.sh`
  Reuses the existing gene-correlation extraction logic after MAGMA raw files exist.
- `07_catfish_array_LSF.sh`
  Runs CATFISH pathway analysis via `07_catfish_perm.R`.
- `08_aggregate_and_plot.sh`
  Pools the per-spike CATFISH outputs and creates null-style summary plots/tables.
- `09_spikein_eval.R`
  Staged spike evaluation script for power, un-spiked calibration, and dominant-component summaries.
- `09_spikein_eval.sh`
  LSF wrapper for staged spike evaluation.

## Recommended flow

Staged route:

1. Copy `config.spikein.template.sh` to `config.spikein.sh` and edit paths/settings.
2. Set a real `SPIKE_PATHWAY` in `config.spikein.sh`.
   If sorghum gene identifiers mismatch the MAGMA annotation, optionally set
   `SPIKE_GENE_COL` to `Gene-name` or `Gene-id`.
3. Build the staged spike input once:

```bash
bsub < 01_spike_setup.sh
bsub < 02_spike_extract_geno.sh
bsub < 03_make_spikein.sh
```

4. Submit the array jobs:

```bash
bsub < 04_gwas_array_LSF.sh
bsub < 05_magma_array_LSF.sh
bsub < 06_corpairs_once.sh
bsub < 07_catfish_array_LSF.sh
```

5. Once CATFISH outputs are done:

```bash
bsub < 08_aggregate_and_plot.sh
bsub < 09_spikein_eval.sh
```

## Outputs

- `PHENO_DIR/pheno_perm_###.csv`
  Spiked phenotype files in the same format used by `vcf2gwas -pf`.
- `PHENO_DIR/truth.tsv`
  Truth table for the staged spike run.
- `CATFISH_OUT/pathway_pvals_perm_###.csv`
  Per-replicate CATFISH pathway results.
- `SPIKE_SUMMARY_DIR/power_summary.csv`
  Recovery summary for the spiked pathway.
- `SPIKE_SUMMARY_DIR/calibration_unspiked.csv`
  Calibration summary for un-spiked pathways.

## Important note

These spike-ins are complementary to the existing idealized archetype simulations.
They do **not** replace the phenotype-permutation null calibration analysis.
