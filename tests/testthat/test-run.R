# catfish_run(): the whole analysis in one call, from MAGMA output or from
# GWAS summary statistics (MAGMA replaced by a stand-in that writes its files).

toy_genes <- function(chr, n = 30, pre = "Zm") {
  start <- seq(1000, by = 5000, length.out = n)
  set.seed(chr)
  p <- runif(n)
  if (chr == 1) p[1:5] <- c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4)
  data.frame(GENE = sprintf("%s%02d%04d", pre, chr, seq_len(n)), CHR = chr,
             START = start, STOP = start + 3000, NSNPS = 10, NPARAM = 3, N = 300,
             ZSTAT = qnorm(p, lower.tail = FALSE), P = p)
}
write_genes_out <- function(g, f) write.table(g, f, sep = " ", quote = FALSE, row.names = FALSE)
write_genes_raw <- function(g, f) {
  ln <- vapply(seq_len(nrow(g)), function(i) {
    r <- if (i > 1) "0.2" else character(0)           # correlation with the previous gene
    paste(c(g$GENE[i], g$CHR[i], g$START[i], g$STOP[i], 10, 3, 300, 300, 1, r), collapse = " ")
  }, character(1))
  writeLines(c("# VERSION = 110", "# COVAR = NSAMP MAC", ln), f)
}
toy_gwas <- function(genes, f) {
  snp <- do.call(rbind, lapply(seq_len(nrow(genes)), function(i)
    data.frame(chr = genes$CHR[i], rs = paste0(genes$GENE[i], "_s", 1:3),
               ps = genes$START[i] + c(100, 1500, 2900),
               p_wald = pmin(1, genes$P[i] * c(1, 3, 9)))))
  write.table(snp, f, sep = "\t", quote = FALSE, row.names = FALSE)
}
toy_pathways <- function(genes) {
  ids <- genes$GENE
  set.seed(9)
  c(list(PWY_signal = ids[1:8]),
    stats::setNames(lapply(1:9, function(i) sample(ids[-(1:8)], 8)), paste0("PWY_", 1:9)))
}

test_that("catfish_run from existing MAGMA output gives pathways, components and candidates", {
  td <- tempfile("run"); dir.create(file.path(td, "magma"), recursive = TRUE)
  g <- rbind(toy_genes(1), toy_genes(2))
  for (k in 1:2) {
    write_genes_out(g[g$CHR == k, ], file.path(td, "magma", sprintf("toy_chr%d.genes.out", k)))
    write_genes_raw(g[g$CHR == k, ], file.path(td, "magma", sprintf("toy_chr%d.genes.raw", k)))
  }
  gw <- file.path(td, "gwas.txt"); toy_gwas(g, gw)

  res <- catfish_run(magma_dir = file.path(td, "magma"), gwas_file = gw,
                     rename_columns = c(CHR = "chr", SNP = "rs", POS = "ps", PVALUE = "p_wald"),
                     pathways = toy_pathways(g), B = 200, min_genes = 5,
                     out_dir = file.path(td, "out"), out_prefix = "toy", verbose = FALSE)

  expect_named(res, c("pathways", "candidates", "genes", "pathway_sets", "settings", "files"))
  expect_equal(res$settings$p_col, "P")
  expect_equal(res$settings$gene_regex, "^Zm")
  expect_equal(res$pathways$pathway_id[1], "PWY_signal")
  expect_true(all(c("acat_p", "fisher_p", "tfisher_p", "minp_p", "stouffer_p") %in% names(res$pathways)))
  expect_false(anyNA(res$pathways$fisher_p))
  expect_true(res$pathways$p_fwer_std[1] <= 0.05)
  expect_equal(res$candidates$gene[1], g$GENE[1])
  for (f in res$files) expect_true(file.exists(f))
  ## correlations of both chromosomes in one file
  cp <- read.delim(res$files$cor_pairs)
  expect_equal(nrow(cp), 2 * 29)
})

test_that("catfish_run from GWAS statistics runs MAGMA per chromosome, then the rest", {
  skip_on_os("windows")
  td <- tempfile("run2"); dir.create(td)
  g <- rbind(toy_genes(1), toy_genes(2))
  gdir <- file.path(td, "truth"); dir.create(gdir)
  for (k in 1:2) {
    write_genes_out(g[g$CHR == k, ], file.path(gdir, sprintf("chr%d.genes.out", k)))
    write_genes_raw(g[g$CHR == k, ], file.path(gdir, sprintf("chr%d.genes.raw", k)))
  }
  ## stand-in for MAGMA: --annotate writes .genes.annot; a gene run copies the
  ## prepared output of the chromosome named in --out (..._chr<k>.<model>)
  fm <- file.path(td, "fake_magma")
  writeLines(c("#!/bin/sh",
               "out=''; ann=''",
               "while [ $# -gt 0 ]; do",
               "  case \"$1\" in --out) out=\"$2\"; shift;; --annotate) ann=1;; esac",
               "  shift",
               "done",
               "if [ -n \"$ann\" ]; then echo '# annot' > \"$out.genes.annot\"; exit 0; fi",
               "k=$(basename \"$out\" | sed 's/.*_chr\\([0-9]*\\)\\..*/\\1/')",
               sprintf("cp %s/chr$k.genes.out \"$out.genes.out\"", gdir),
               sprintf("cp %s/chr$k.genes.raw \"$out.genes.raw\"", gdir)), fm)
  Sys.chmod(fm, "755")
  old <- options(magma.path = fm); on.exit(options(old), add = TRUE)
  bfile <- file.path(td, "ref"); for (e in c(".bed", ".bim", ".fam")) file.create(paste0(bfile, e))
  loc <- file.path(td, "genes.loc")
  write.table(g[, c("GENE", "CHR", "START", "STOP")], loc, sep = "\t", quote = FALSE, row.names = FALSE)
  gw <- file.path(td, "gwas.txt"); toy_gwas(g, gw)

  res <- catfish_run(gwas_file = gw,
                     rename_columns = c(CHR = "chr", SNP = "rs", POS = "ps", PVALUE = "p_wald"),
                     n_total = 300, bfile = bfile, gene_loc = loc, chroms = 1:2,
                     pathways = toy_pathways(g), B = 200, min_genes = 5,
                     out_dir = file.path(td, "out"), out_prefix = "toy", verbose = FALSE)
  expect_equal(sort(basename(attr(res$genes, "files"))),
               c("toy_chr1.snp_wise_mean.genes.out", "toy_chr2.snp_wise_mean.genes.out"))
  expect_equal(nrow(res$genes), 60)
  expect_equal(res$pathways$pathway_id[1], "PWY_signal")
  expect_true(file.exists(res$files$gene_annot))
  expect_s3_class(res$candidates, "data.frame")
})

test_that("catfish_run checks its inputs", {
  expect_error(catfish_run(species = "sorghum"), "gwas_file")
  expect_error(catfish_run(magma_dir = tempdir()), "pathways")
  expect_error(catfish_run(gwas_file = "no_such_file", species = "sorghum"), "not found")
})

test_that("gene ID patterns are worked out from the IDs", {
  expect_equal(.catfish_gene_regex(c("SORBI_3001G000100", "SORBI_3010G1")), "^SORBI_")
  expect_equal(.catfish_gene_regex(c("Zm00001eb000010", "Zm00001eb1")), "^Zm")
  expect_equal(.catfish_gene_regex(c("1017", "348")), "^[0-9]+$")
  expect_equal(.catfish_gene_regex(c("AT1G01010", "ATCG00020", "ATMG1")), "^(AT|ATCG|ATMG)")
})
