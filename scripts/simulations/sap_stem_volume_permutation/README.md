# SAP Stem-Volume Phenotype-Permutation Pipeline

This folder mirrors the BAP phenotype-permutation workflow, but for the sorghum SAP stem-volume trait.

The intended flow is:

1. Prepare genotype inputs once from the SAP chromosome VCF files and build per-chromosome MAGMA annotations.
2. Generate permuted stem-volume phenotype files.
3. Re-run the SAP GWAS once per permutation x chromosome.
4. Re-run MAGMA gene analysis.
5. Build the genotype-only MAGMA gene-correlation file once.
6. Run CATFISH pathway enrichment on each permutation.
7. Aggregate the null pathway p-values into QQ/histogram/lambda/type-I summaries.

## Chosen CATFISH settings

These scripts default to the newer settings we have been using:

- `mvn_calibrate_components = FALSE`
- adaptive/oTFisher tau grid = `c(0.1, 0.05, 0.01)`
- `tail_mode = "hybrid_gpd"`
- `B_perm = 10000`

## Phenotype file

The config now expects the phenotype CSV inside the pipeline run folder:

- `${PIPE_ROOT}/inputs/stem_volume.csv`

Source file on your laptop:

- `/Users/nirwantandukar/Documents/Research/data/SAP/GWAS_phenotype/Stem_volume/stem_volume.csv`

So on HPC the easiest flow is just:

```bash
mkdir -p /share/maize/ntanduk/CATFISH/simulation/sap_stem_volume_permutation_run/inputs
cp /path/to/stem_volume.csv /share/maize/ntanduk/CATFISH/simulation/sap_stem_volume_permutation_run/inputs/stem_volume.csv
```

## Important GWAS note

This SAP pipeline now assumes the GWAS reruns will use chromosome-specific VCFs through `vcf2gwas`, with a private per-job tmp copy of the VCF and its index. That avoids the `.csi` deletion collision you hit before when many jobs touched the same shared VCF.

## Files

- `config.sap_stem_volume.template.sh`
  Copy this to `config.sap_stem_volume.sh` on HPC and fill in the real paths.
- `00_prepare_genotype.R`
  Converts SAP chromosome VCFs to per-chromosome PLINK and builds MAGMA annotations from the real stem-volume GWAS.
- `00_prepare_genotype.sh`
  LSF wrapper for step 0.
- `01_make_permutations.R`
  Writes permuted stem-volume phenotypes aligned to the PLINK `.fam` sample order.
- `01_make_permutations.sh`
  LSF wrapper for step 1.
- `02_run_gwas_single.sh`
  Runs `vcf2gwas` for one permutation/chromosome using a private tmp VCF copy.
- `02_gwas_array_LSF.sh`
  Array wrapper for the SAP GWAS reruns.
- `03_magma_catfish.R`
  MAGMA gene analysis for one permutation/chromosome.
- `03_magma_array_LSF.sh`
  Array wrapper for MAGMA.
- `04_corpairs_once.R`
  Extracts the genotype-only MAGMA gene-correlation pairs once.
- `04_corpairs_once.sh`
  LSF wrapper for the correlation extraction.
- `05_catfish_perm.R`
  CATFISH pathway analysis for one permutation.
- `05_catfish_array_LSF.sh`
  Array wrapper for CATFISH.
- `06_aggregate_and_plot.R`
  Pools null pathway p-values and writes `lambda_table.csv`, `type1_error_table.csv`, QQ, and histogram.
- `06_aggregate_and_plot.sh`
  LSF wrapper for aggregation.

## Recommended run order

```bash
cp config.sap_stem_volume.template.sh config.sap_stem_volume.sh
# edit config.sap_stem_volume.sh

bsub < 00_prepare_genotype.sh
bsub < 01_make_permutations.sh
bsub < 02_gwas_array_LSF.sh
bsub < 03_magma_array_LSF.sh
bsub < 04_corpairs_once.sh
bsub < 05_catfish_array_LSF.sh
bsub < 06_aggregate_and_plot.sh
```

## Expected outputs

- `${PHENO_DIR}/pheno_perm_###.csv`
- `${GWAS_DIR}/perm_###/chrXX/assoc.txt`
- `${MAGMA_OUT}/perm_###/chrXX/*.genes.out`
- `${CORPAIRS}`
- `${CATFISH_OUT}/pathway_pvals_perm_###.csv`
- `${FINAL_DIR}/lambda_table.csv`
- `${FINAL_DIR}/type1_error_table.csv`
- `${FINAL_DIR}/QQ_omnibus_null.png`
- `${FINAL_DIR}/hist_omnibus_null.png`
