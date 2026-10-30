# Tests for gene_prioritize()

gp_toy <- function() {
  set.seed(1)
  genes <- paste0("g", 1:200)
  magma <- data.frame(GENE = genes, P = runif(200))
  magma$P[1:3] <- c(1e-12, 1e-6, 1e-4)
  gwas <- data.frame(SNP = paste0("rs", 1:2000), P = runif(2000))
  gwas$P[1:3] <- c(1e-20, 1e-8, 1e-5)
  gwas$NEGLOG10P <- -log10(gwas$P)
  gwas$CHR <- 1
  gwas$BP  <- seq_len(2000) * 1000
  annot <- data.frame(gene = rep(genes, each = 10), SNP = gwas$SNP)
  loc <- data.frame(gene = genes, chr = 1, start = (0:199) * 10000 + 1000,
                    end = (0:199) * 10000 + 9000, symbol = toupper(genes))
  pathways <- c(list(PWY_signal = genes[1:12]),
                setNames(lapply(1:14, function(i) sample(genes[-(1:12)], 12)),
                         paste0("PWY_", 1:14)))
  pw_res <- data.frame(pathway_id = names(pathways),
                       log10_p_omnibus = c(-15, log10(runif(14))))
  list(magma = magma, gwas = gwas, annot = annot, loc = loc,
       pathways = pathways, pw_res = pw_res)
}

test_that("the planted gene ranks first and all columns are present", {
  d <- gp_toy()
  r <- gene_prioritize(d$magma, d$pw_res, d$pathways, gwas = d$gwas,
                       gene_annot = d$annot, gene_loc = d$loc, verbose = FALSE)
  expect_equal(r$gene[1], "g1")
  expect_equal(r$symbol[1], "G1")
  expect_equal(r$best_pathway[1], "PWY_signal")
  expect_equal(r$rank, seq_len(nrow(r)))
  expect_true(all(c("rank", "gene", "symbol", "score", "n_layers", "gwas_score",
                    "magma_score", "pathway_score", "hub_bonus",
                    "gwas_neglog10p", "gwas_rank", "magma_p", "magma_rank",
                    "best_pathway", "pathway_stat", "pathway_rank",
                    "n_pathways") %in% names(r)))
  expect_false(is.unsorted(-r$score))
  w <- attr(r, "weights")
  expect_equal(r$score, w[["gwas"]] * r$gwas_score + w[["magma"]] * r$magma_score +
                 w[["pathway"]] * r$pathway_score + w[["hub"]] * r$hub_bonus)
  expect_equal(attr(r, "gwas_route"), "gene_annot")
})

test_that("-log10(p) is preferred to p and both give the same order here", {
  d <- gp_toy()
  a <- gene_prioritize(d$magma, d$pw_res, d$pathways, gwas = d$gwas,
                       gene_annot = d$annot, verbose = FALSE)
  b <- gene_prioritize(d$magma, d$pw_res, d$pathways,
                       gwas = d$gwas[, c("SNP", "P")],
                       gene_annot = d$annot, verbose = FALSE)
  expect_equal(a$gene, b$gene)
  expect_equal(a$gwas_neglog10p, b$gwas_neglog10p, tolerance = 1e-10)
})

test_that("position route assigns SNPs within the window", {
  d <- gp_toy()
  r <- gene_prioritize(d$magma, d$pw_res, d$pathways, gwas = d$gwas,
                       gene_loc = d$loc, window = 0, verbose = FALSE)
  expect_equal(attr(r, "gwas_route"), "position")
  ## gene g1 spans 1000-9000 bp, so SNPs rs1..rs9 (1000..9000 bp) fall inside
  i <- match("g1", r$gene)
  expect_equal(r$gwas_neglog10p[i], max(d$gwas$NEGLOG10P[1:9]))
})

test_that("files and R objects give the same result", {
  d <- gp_toy()
  td <- tempfile(); dir.create(td)
  write.table(d$magma, file.path(td, "x.genes.out"), quote = FALSE, row.names = FALSE)
  write.csv(d$pw_res, file.path(td, "pw.csv"), row.names = FALSE)
  writeLines(vapply(names(d$pathways), function(n)
    paste(c(n, d$pathways[[n]]), collapse = " "), ""), file.path(td, "sets.txt"))
  write.table(d$gwas[, c("SNP", "NEGLOG10P")], file.path(td, "gwas.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(vapply(unique(d$annot$gene), function(g)
    paste(c(g, "1:1:2", d$annot$SNP[d$annot$gene == g]), collapse = "\t"), ""),
    file.path(td, "x.genes.annot"))
  f <- gene_prioritize(file.path(td, "x.genes.out"), file.path(td, "pw.csv"),
                       file.path(td, "sets.txt"), gwas = file.path(td, "gwas.txt"),
                       gene_annot = file.path(td, "x.genes.annot"), verbose = FALSE)
  o <- gene_prioritize(d$magma, d$pw_res, d$pathways, gwas = d$gwas,
                       gene_annot = d$annot, verbose = FALSE)
  expect_equal(f$gene, o$gene)
  expect_equal(f$score, o$score)
})

test_that("missing columns give a clear error naming the argument", {
  d <- gp_toy()
  expect_error(gene_prioritize(d$magma, d$pw_res, d$pathways,
                               magma_gene_col = "SYMBOL", verbose = FALSE),
               "magma_gene_col")
  bad <- d$pw_res; names(bad)[2] <- "whatever"
  expect_error(gene_prioritize(d$magma, bad, d$pathways, verbose = FALSE),
               "pathway_stat_col")
})

test_that("without GWAS the score uses the MAGMA and pathway layers only", {
  d <- gp_toy()
  r <- gene_prioritize(d$magma, d$pw_res, d$pathways, verbose = FALSE)
  expect_true(all(r$gwas_score == 0))
  expect_equal(attr(r, "gwas_route"), "none")
})

test_that("targets are reported by ID or symbol", {
  d <- gp_toy()
  r <- gene_prioritize(d$magma, d$pw_res, d$pathways, gene_loc = d$loc,
                       targets = c("G1", "g2", "not_a_gene"), verbose = FALSE)
  t <- attr(r, "targets")
  expect_equal(nrow(t), 3)
  expect_equal(t$rank[1], r$rank[r$gene == "g1"])
  expect_true(is.na(t$rank[3]))
})

test_that("works on catfish_wy_pathways output directly", {
  d <- gp_toy()
  d$magma$CHR <- 1
  wy <- suppressWarnings(catfish_wy_pathways(d$magma, pathways = d$pathways,
                                             B = 30, min_genes = 5, seed = 1,
                                             verbose = FALSE))
  r <- gene_prioritize(d$magma, wy, d$pathways, verbose = FALSE)
  expect_equal(attr(r, "pathway_stat_col"), "z_std")
  r2 <- gene_prioritize(d$magma, wy, d$pathways, pathway_stat_col = "log10_p_omnibus",
                        verbose = FALSE)
  expect_equal(attr(r2, "pathway_stat_col"), "log10_p_omnibus")
  expect_equal(r$best_pathway[r$gene == "g1"], "PWY_signal")
})
