# Sorghum Reactome pathway data

- Bundled filename: `sorghum_reactome.csv`
- Original filename: `Reactome_pathways_sorghum_BTx623_maize_B73.csv`
- Added to CATFISH: 2026-11-27
- Species: *Sorghum bicolor*
- Assembly: BTx623 v3 (`GCA_000003195.3`)
- Gene identifiers: `SORBI_300...`
- Rows: 3,911 unique pathway–gene memberships
- Unique pathways: 237
- Unique genes: 2,864
- SHA-256: `f03dccb88d07d6a73e5e2f2aa224103575e2a2d11116b484e4e9059da5981c16`

The CSV is a long Reactome annotation table. CATFISH uses the `pathway_id`,
`pathway_name`, and `gene_id` columns; the remaining provenance, genomic
location, hierarchy, and reaction columns are retained in the bundled source.

Load it with:

```r
pw <- catfish_load_pathways(species = "sorghum_reactome")
```
