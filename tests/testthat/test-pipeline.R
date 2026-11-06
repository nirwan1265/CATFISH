# Glue between the MAGMA steps and the pathway layer: per-chromosome MAGMA
# runs, reading them back, gene-gene correlations from several .genes.raw
# files, pathway loading, and the Westfall-Young output columns.

## ---- a stand-in for the MAGMA binary: records each call and writes the files
## MAGMA would write, so the R side can be tested without MAGMA installed
fake_magma <- function(dir) {
  f <- file.path(dir, "fake_magma")
  writeLines(c(
    "#!/bin/sh",
    "out=''; pval=''; ann=''",
    "for a in \"$@\"; do echo \"$a\"; done >> \"$(dirname \"$0\")/calls.log\"",
    "echo '---' >> \"$(dirname \"$0\")/calls.log\"",
    "while [ $# -gt 0 ]; do",
    "  case \"$1\" in --out) out=\"$2\"; shift;; --pval) pval=\"$2\"; shift;; --annotate) ann=1;; esac",
    "  shift",
    "done",
    "if [ -n \"$ann\" ]; then echo '# window' > \"$out.genes.annot\"; fi",
    "if [ -n \"$pval\" ]; then",
    "  n=$(($(wc -l < \"$pval\") - 1))",
    "  b=$(basename \"$out\")",
    "  printf 'GENE CHR START STOP NSNPS NPARAM N ZSTAT P\\n' > \"$out.genes.out\"",
    "  printf 'g_%s 1 1 100 %s 2 300 1.5 0.07\\n' \"$b\" \"$n\" >> \"$out.genes.out\"",
    "fi"
  ), f)
  Sys.chmod(f, "755")
  f
}

test_that("magma_gene with chroms and n_threads = 1 runs one MAGMA call per chromosome", {
  skip_on_os("windows")
  td <- tempfile("pipe"); dir.create(td)
  old <- options(magma.path = fake_magma(td)); on.exit(options(old), add = TRUE)

  bfile <- file.path(td, "ref")
  for (e in c(".bed", ".bim", ".fam")) file.create(paste0(bfile, e))
  annot <- file.path(td, "x.genes.annot"); writeLines("# annot", annot)
  gw <- data.frame(chr = rep(c(1, 2, 3), c(5, 3, 4)),
                   rs  = paste0("s", 1:12), ps = 1:12, p_wald = runif(12))
  stats <- file.path(td, "gwas.txt")
  write.table(gw, stats, sep = "\t", quote = FALSE, row.names = FALSE)

  res <- magma_gene(bfile = bfile, gene_annot = annot, stats_file = stats,
                    n_total = 300,
                    rename_columns = c(CHR = "chr", SNP = "rs", POS = "ps", PVALUE = "p_wald"),
                    out_prefix = "trait", out_dir = file.path(td, "magma_out"),
                    gene_model = "snp-wise=mean", chroms = c(1, 3),
                    n_threads = 1, verbose = FALSE)

  expect_named(res, c("1", "3"))
  outs <- sort(list.files(file.path(td, "magma_out"), "\\.genes\\.out$"))
  expect_equal(outs, c("trait_chr1.snp_wise_mean.genes.out",
                       "trait_chr3.snp_wise_mean.genes.out")[seq_along(outs)])
  expect_length(outs, 2)
  ## each run saw only its own chromosome's SNPs
  n1 <- nrow(read.table(file.path(td, "magma_out", "trait_chr1.pval.txt"), header = TRUE))
  n3 <- nrow(read.table(file.path(td, "magma_out", "trait_chr3.pval.txt"), header = TRUE))
  expect_equal(c(n1, n3), c(5, 4))
  calls <- readLines(file.path(td, "calls.log"))
  expect_equal(sum(calls == "--gene-model"), 2)

  ## and the per-chromosome outputs read back as one table
  g <- magma_read_genes(dir = file.path(td, "magma_out"), verbose = FALSE)
  expect_equal(nrow(g), 2)
  expect_equal(g$NSNPS, c(5, 4))
})

test_that("magma_read_genes combines files, keeps first duplicate, rejects mixed columns", {
  td <- tempfile("rg"); dir.create(td)
  mk <- function(f, d) write.table(d, file.path(td, f), sep = " ", quote = FALSE, row.names = FALSE)
  base <- function(g, chr, start, p) data.frame(GENE = g, CHR = chr, START = start,
                                               STOP = start + 10, NSNPS = 3, NPARAM = 2,
                                               N = 100, ZSTAT = qnorm(p, lower.tail = FALSE),
                                               P_MULTI = p, P_SNPWISE_MEAN = p)
  mk("t_chr10.genes.out", base(c("c", "d"), 10, c(1, 5), c(0.5, 0.6)))
  mk("t_chr2.genes.out",  base(c("a", "b"), 2,  c(1, 5), c(0.2, 0.9)))
  expect_message(g <- magma_read_genes(dir = td), "median P_SNPWISE_MEAN")
  expect_equal(g$GENE, c("a", "b", "c", "d"))         # numeric chromosome order
  expect_length(attr(g, "files"), 2)

  mk("t_dup.genes.out", base("a", 2, 1, 0.3))
  expect_message(g2 <- magma_read_genes(dir = td, verbose = FALSE), "more than one file")
  expect_equal(nrow(g2), 4)

  mk("t_bad.genes.out", data.frame(GENE = "z", CHR = 1, P = 0.1))
  expect_error(magma_read_genes(dir = td, verbose = FALSE), "Columns differ")
  expect_error(magma_read_genes(), "either 'files' or 'dir'")
})

test_that("cor pairs from several .genes.raw files equal the per-file results stacked", {
  td <- tempfile("raw"); dir.create(td)
  hdr <- c("# VERSION = 110", "# COVAR = NSAMP MAC")
  r1 <- file.path(td, "chr1.genes.raw"); r2 <- file.path(td, "chr2.genes.raw")
  writeLines(c(hdr,
               "SORBI_1 1 0 100 10 4 300 300 1",
               "SORBI_2 1 50 200 12 5 300 300 1 0.4",
               "SORBI_3 1 150 300 9 3 300 300 1 0.1 0.25"), r1)
  writeLines(c(hdr,
               "SORBI_4 2 0 100 10 4 300 300 1",
               "SORBI_5 2 50 200 12 5 300 300 1 -0.3"), r2)
  both <- file.path(td, "both.txt")
  magma_genesraw_to_cor_pairs_banded(c(r1, r2), both, gene_regex = "^SORBI", verbose = FALSE)
  x <- read.delim(both)
  expect_equal(nrow(x), 4)
  expect_equal(x$gene1, c("SORBI_2", "SORBI_3", "SORBI_3", "SORBI_5"))
  expect_equal(x$gene2, c("SORBI_1", "SORBI_1", "SORBI_2", "SORBI_4"))
  expect_equal(x$r, c(0.4, 0.1, 0.25, -0.3))

  one <- file.path(td, "one.txt"); two <- file.path(td, "two.txt")
  magma_genesraw_to_cor_pairs_banded(r1, one, gene_regex = "^SORBI", verbose = FALSE)
  magma_genesraw_to_cor_pairs_banded(r2, two, gene_regex = "^SORBI", verbose = FALSE)
  expect_equal(readLines(both), c(readLines(one), readLines(two)[-1]))

  ## the number of fixed fields follows the # COVAR header
  r3 <- file.path(td, "chr3.genes.raw")
  writeLines(c("# VERSION = 110", "# COVAR = NSAMP",
               "SORBI_6 3 0 100 10 4 300 300", "SORBI_7 3 50 200 12 5 300 300 0.5"), r3)
  three <- file.path(td, "three.txt")
  magma_genesraw_to_cor_pairs_banded(r3, three, gene_regex = "^SORBI", verbose = FALSE)
  expect_equal(read.delim(three)$r, 0.5)

  ## a pattern that matches no gene is an error, not an empty file
  expect_error(magma_genesraw_to_cor_pairs_banded(r1, one, gene_regex = "^FBgn", verbose = FALSE),
               "No gene record")
})

test_that("built-in sorghum pathways use the gene IDs of the built-in sorghum gene locations", {
  pw  <- catfish_load_pathways(species = "sorghum")
  loc <- read.table(system.file("extdata", "sorghum.genes.loc", package = "CATFISH"),
                    header = TRUE, stringsAsFactors = FALSE)
  expect_gt(mean(grepl("^SORBI_", pw$gene_id)), 0.99)
  expect_gt(mean(unique(pw$gene_id) %in% loc$GENE), 0.9)
})

test_that("catfish_load_pathways reads a user PMN-format file", {
  f <- tempfile(fileext = ".tsv")
  writeLines(c("Pathway-id\tPathway-name\tGene-name",
               "PWY-1\tfirst\tg1", "PWY-1\tfirst\tg2", "PWY-2\tsecond\tunknown",
               "PWY-2\tsecond\tg3"), f)
  pw <- catfish_load_pathways(file = f)
  expect_equal(pw$gene_id, c("g1", "g2", "g3"))
  expect_equal(unique(pw$pathway_name), c("first", "second"))
  expect_error(catfish_load_pathways(file = "no_such_file"), "not found")
})

test_that("catfish_wy_pathways reports pathway names", {
  set.seed(2)
  genes <- data.frame(GENE = paste0("g", 1:60), CHR = rep(1:2, each = 30),
                      P = runif(60))
  pw <- data.frame(pathway_id = rep(c("A", "B"), each = 10),
                   pathway_name = rep(c("alpha path", "beta path"), each = 10),
                   gene_id = paste0("g", c(1:10, 31:40)))
  res <- suppressWarnings(catfish_wy_pathways(genes, pathways = pw, B = 50,
                                              min_genes = 5, seed = 1, verbose = FALSE))
  expect_equal(names(res)[1:2], c("pathway_id", "pathway_name"))
  expect_equal(res$pathway_name[match(c("A", "B"), res$pathway_id)],
               c("alpha path", "beta path"))
  ## named-list input: the name is the ID
  res2 <- suppressWarnings(catfish_wy_pathways(genes,
            pathways = split(pw$gene_id, pw$pathway_id),
            B = 50, min_genes = 5, seed = 1, verbose = FALSE))
  expect_equal(res2$pathway_name, res2$pathway_id)
})
