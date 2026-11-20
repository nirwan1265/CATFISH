#!/usr/bin/env Rscript
## =============================================================================
## wy_logo.R -- leave-one-gene-out pathway evidence for the positive controls of
## the lipid gene score (Section 2.5). For each established gene g and each
## tested pathway P that contains it, P and P without g are scored together in
## one Westfall-Young run, so the change in the standardised statistic z_std
## from removing g is measured on the same null draws. 39_logo_recovery.R then
## gives each control the pathway evidence of its pathways without itself and
## recomputes its rank in the gene score: if the control still ranks high, the
## lift from the pathway layer comes from the other genes of its pathways.
##
## Inputs and settings are those of wy_gc_real.R (lipids): the same .genes.out,
## gene-level GC, correlation file, Reactome sets, B = 10000, 5 components,
## tau {0.01, 0.05, 0.5}, seed 20260706. Pathways: those in the paper's result
## (out/<A>_reactome_5comp_westfall_young.csv) that contain a control; the
## reduced sets may have one gene fewer than the 10-gene minimum.
##   Rscript wy_logo.R <analysis>      e.g. ldl_eur, hdl_afr
## Writes out_logo/<A>_logo.csv: pathway_id, base_pathway, removed_gene, z_std
## =============================================================================
suppressPackageStartupMessages(library(CATFISH))
A    <- commandArgs(trailingOnly = TRUE)[1]
ROOT <- Sys.getenv("GC_ROOT", "/share/maize/ntanduk/CATFISH/gc_real")
B    <- as.integer(Sys.getenv("WY_B", "10000"))
IN   <- file.path(ROOT, "inputs")
OUT  <- file.path(ROOT, "out_logo")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SEED <- 20260706L; TAU <- c(0.01, 0.05, 0.5); MIN <- 10L
if (!grepl("^(ldl|hdl)_(eur|afr|eas|sas|his)$", A)) stop("lipid analyses only, e.g. ldl_eur", call. = FALSE)

## positive controls of 06_candidate_genes.R, as Entrez IDs (from its output)
CONTROLS <- list(
  ldl = c(LDLR = "3949", PCSK9 = "255738", APOB = "338", APOE = "348", SORT1 = "6272",
          HMGCR = "3156", NPC1L1 = "29881", ABCG8 = "64241", LPA = "4018", CETP = "1071",
          LDLRAP1 = "26119", MYLIP = "29116", ANGPTL3 = "27329"),
  hdl = c(CETP = "1071", LIPC = "3990", LIPG = "9388", ABCA1 = "19", APOA1 = "335",
          APOA5 = "116519", LCAT = "3931", SCARB1 = "949", PLTP = "5360", LPL = "4023",
          ANGPTL4 = "51129", CUBN = "8029"))

## ---- gene statistics, GC, pathways: as wy_gc_real.R ----------------------------
g <- read.table(file.path(IN, "lipids", paste0(A, ".genes.out")), header = TRUE,
                comment.char = "#", stringsAsFactors = FALSE)
g$GENE <- as.character(g$GENE)
g <- g[is.finite(g$P) & is.finite(g$ZSTAT), ]
z <- stats::qnorm(log(pmax(g$P, 1e-320)), lower.tail = FALSE, log.p = TRUE)
g$P_GC <- stats::pnorm(z - stats::median(z), lower.tail = FALSE)
pw <- list()
for (ln in readLines(file.path(IN, "c2.cp.reactome.v2026.2.Hs.entrez.sets.annot"), warn = FALSE)) {
  f <- strsplit(ln, "\t", fixed = TRUE)[[1]]; if (length(f) < 2) next
  gs <- unique(f[-1]); pw[[f[1]]] <- gs[nzchar(gs)]
}
paper <- file.path(ROOT, "out", paste0(A, "_reactome_5comp_westfall_young.csv"))
if (!file.exists(paper)) stop("paper's result not found: ", paper, call. = FALSE)
tested <- utils::read.csv(paper, stringsAsFactors = FALSE)$pathway_id

## ---- each control's pathways, with and without it ----------------------------------
ctrl <- CONTROLS[[sub("_.*", "", A)]]
ctrl <- ctrl[ctrl %in% g$GENE]
sets <- list(); map <- list()
for (s in names(ctrl)) {
  id <- ctrl[[s]]
  for (P in tested[vapply(tested, function(p) id %in% pw[[p]], logical(1))]) {
    sets[[P]] <- pw[[P]]
    k <- paste0(P, "__minus__", id)
    sets[[k]] <- setdiff(pw[[P]], id)
    map[[k]] <- data.frame(pathway_id = c(P, k), base_pathway = P,
                           removed_gene = c(NA, id), stringsAsFactors = FALSE)
  }
}
map <- unique(do.call(rbind, map))
cat(sprintf("[%s] %d controls in %d tested pathways -> %d sets, B = %d\n",
            A, length(ctrl), length(unique(map$base_pathway)), length(sets), B))

t0 <- Sys.time()
res <- catfish_wy_pathways(g, pathways = sets, gene_col = "GENE", p_col = "P_GC", chr_col = "CHR",
                           z_col = "ZSTAT", cor_file = file.path(IN, "lipids", paste0(A, ".cor_pairs.txt")),
                           B = B, use_tfisher = TRUE, tau_grid = TAU, min_genes = MIN - 1L,
                           max_genes = 200L, seed = SEED, verbose = FALSE)
out <- merge(map, res[, c("pathway_id", "z_std", "n_genes")], by = "pathway_id", all.x = TRUE)
utils::write.csv(out, file.path(OUT, paste0(A, "_logo.csv")), row.names = FALSE)
cat(sprintf("[%s] done in %.1f min -> %s\n", A,
            as.numeric(difftime(Sys.time(), t0, units = "mins")), file.path(OUT, paste0(A, "_logo.csv"))))
