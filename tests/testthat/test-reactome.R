# Built-in human Reactome pathways (MSigDB C2:CP:REACTOME, Entrez IDs)

test_that("human_reactome loads as a long pathway table", {
  f <- catfish_pmn_file("human_reactome")
  expect_true(file.exists(f))
  expect_equal(basename(f), "human_reactome.pw")
  pw <- catfish_load_pathways(species = "human_reactome")
  expect_named(pw, c("pathway_id", "pathway_name", "gene_id"))
  expect_equal(length(unique(pw$pathway_id)), 1867L)
  expect_true(all(grepl("^REACTOME_", pw$pathway_id)))
  expect_true(all(grepl("^[0-9]+$", pw$gene_id)))      # Entrez IDs
  i <- pw$pathway_id == "REACTOME_PLASMA_LIPOPROTEIN_ASSEMBLY_REMODELING_AND_CLEARANCE"
  expect_true("3949" %in% pw$gene_id[i])                # LDLR
  expect_equal(pw$pathway_name[i][1],
               "plasma lipoprotein assembly remodeling and clearance")
})

test_that("the set file works directly as gene_prioritize pathways", {
  f <- catfish_pmn_file("human_reactome")
  pw <- catfish_load_pathways(species = "human_reactome")
  genes <- unique(pw$gene_id)[1:400]
  magma <- data.frame(GENE = genes, P = seq(1e-6, 1, length.out = 400))
  ids <- unique(pw$pathway_id[pw$gene_id %in% genes])[1:20]
  res <- data.frame(pathway_id = ids, log10_p_omnibus = -seq_along(ids))
  r <- gene_prioritize(magma, res, f, verbose = FALSE)
  expect_equal(attr(r, "layer_sizes")[["pathway"]],
               sum(genes %in% pw$gene_id[pw$pathway_id %in% ids]))
  r2 <- gene_prioritize(magma, res, species = "human_reactome", verbose = FALSE)
  expect_equal(r2, r)
  expect_error(gene_prioritize(magma, res, f, species = "human_reactome"),
               "not both")
})
