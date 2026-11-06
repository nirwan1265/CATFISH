#!/usr/bin/env Rscript
## =============================================================================
## 29_replication_tables.R -- candidate-gene tables for the cross-ancestry
## replication (Section 3.9.7): which evidence layer carries each gene in each
## ancestry group.
##
## Run from catfish2/scripts AFTER 06_candidate_genes.R for every trait and
## ancestry (TRAIT=ldl|hdl, ANC=eur|afr|eas|sas|his).
##
##   Rscript 29_replication_tables.R
##
## Writes
##   ../tex/auto/tab_replication_genes.tex          main text: genes in the top
##        TOP_G in all five groups, with combined, GWAS and MAGMA rank per group
##   ../tex/auto/tab_replication_candidates_ldl.tex supplementary: top ten
##   ../tex/auto/tab_replication_candidates_hdl.tex     per non-European group,
##        same columns as tab_lipid_candidates (21_lipid_results.R)
##   ../results/29_replication_gene_layers.csv      every gene in the top TOP_G
##        in at least one group, all layer ranks, long format
## =============================================================================
suppressPackageStartupMessages(library(data.table))

GP_DIR <- Sys.getenv("GP_DIR", "../results")
RES    <- Sys.getenv("RES_DIR", "../results")
AUT    <- Sys.getenv("AUT_DIR", "../tex/auto")
TOP_G  <- as.integer(Sys.getenv("TOP_G", "20"))      # same cut-off as Figure 17
N_SUPP <- as.integer(Sys.getenv("N_SUPP", "10"))     # rows per group, supp tables
for (d in c(RES, AUT)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

TRAITS <- c(ldl = "LDL-C", hdl = "HDL-C")
ANC    <- c("eur", "afr", "eas", "sas", "his")
ANC_L  <- c(eur = "EUR", afr = "AFR", eas = "EAS", sas = "SAS", his = "HIS")
## positive controls, named before prioritisation (identical to 06 and 21)
CONTROLS <- list(
  ldl = c("LDLR","PCSK9","APOB","APOE","SORT1","HMGCR","NPC1L1","ABCG8",
          "LPA","CETP","LDLRAP1","MYLIP","ANGPTL3"),
  hdl = c("CETP","LIPC","LIPG","ABCA1","APOA1","APOA5","LCAT","SCARB1",
          "PLTP","LPL","ANGPTL4","CUBN"))
LIPO_RX <- "^REACTOME_(PLASMA_LIPOPROTEIN|LDL_|HDL_|VLDL_|VLDLR_|CHYLOMICRON_)"

need <- function(f) if (!file.exists(f)) stop("missing input: ", f, call. = FALSE)
tex_esc <- function(x) gsub("([&%_#$])", "\\\\\\1", x)
fr <- function(x) ifelse(is.na(x), "--", formatC(round(x), format = "d", big.mark = "{,}"))
## readable Reactome name, as in 21_lipid_results.R
pretty_pw <- function(id, width = 40) {
  x <- tolower(gsub("_", " ", sub("^REACTOME_", "", id)))
  cut <- nchar(x) > width
  x[cut] <- substr(x[cut], 1, width - 1)
  fix <- c(rrna = "rRNA", mrna = "mRNA", trna = "tRNA", er = "ER", atp = "ATP",
           pparalpha = "PPAR$\\alpha$", ldl = "LDL", hdl = "HDL", vldl = "VLDL",
           vldlr = "VLDLR", mitf = "MITF", wdr5 = "WDR5", srebp = "SREBP",
           srebf = "SREBF", nr1h2 = "NR1H2", nr1h3 = "NR1H3", tp53 = "TP53",
           gtpase = "GTPase", gtpases = "GTPases", dna = "DNA", rna = "RNA",
           rac1 = "RAC1", egr2 = "EGR2", sox10 = "SOX10", upr = "UPR")
  x <- vapply(strsplit(x, " ", fixed = TRUE), function(w) {
    k <- match(w, names(fix)); w[!is.na(k)] <- fix[k[!is.na(k)]]
    paste(w, collapse = " ") }, character(1))
  first_fixed <- vapply(strsplit(x, " ", fixed = TRUE), function(w) w[1] %in% fix, logical(1))
  x[!first_fixed] <- sub("^([a-z])", "\\U\\1", x[!first_fixed], perl = TRUE)
  ifelse(cut, paste0(x, "\\ldots"), x)
}
gene_tex <- function(sym, t) {
  s <- sprintf("\\textit{%s}", tex_esc(sym))
  ifelse(sym %in% CONTROLS[[t]], sprintf("\\textbf{%s}", s), s)
}

## ---------------------------------------------------------------- inputs ---
G <- list()
for (t in names(TRAITS)) for (a in ANC) {
  f <- file.path(GP_DIR, if (a == "eur") sprintf("06_candidate_genes_%s.csv", t)
                         else sprintf("06_candidate_genes_%s_%s.csv", t, a))
  need(f)
  g <- fread(f, showProgress = FALSE)
  G[[paste(t, a)]] <- g[, .(trait = t, anc = a, symbol, rank, Score, rank_gwas,
                           rank_magma, rank_pathway, best_pathway)]
}

## ------------------------------------------------- long CSV of top genes ---
long <- rbindlist(lapply(names(TRAITS), function(t) {
  d <- rbindlist(G[paste(t, ANC)])
  top <- unique(d[rank <= TOP_G]$symbol)
  d <- d[symbol %in% top]
  d[, n_top := sum(rank <= TOP_G), by = symbol]
  d[, established := symbol %in% CONTROLS[[t]]]
  d[, lipo_pathway := any(grepl(LIPO_RX, best_pathway)), by = symbol]
  d
}))
fwrite(long[order(trait, -n_top, symbol, match(anc, ANC))],
       file.path(RES, "29_replication_gene_layers.csv"))
cat("29_replication_gene_layers.csv\n")

## ====================== main-text table: genes in the top TOP_G in all five ==
L <- c("\\begin{table*}[!t]", "\\centering",
       sprintf(paste0(
         "\\caption{Evidence layers behind the genes ranked in the top %d of the ",
         "combined candidate score in all five ancestry groups. For each group, ",
         "\\emph{Sc.} is the rank under the combined score, \\emph{GW} the rank of ",
         "the gene's best variant (GWAS layer) and \\emph{MA} the rank of its MAGMA ",
         "gene-level $p$-value, each among about $18{,}500$ genes (1 is the ",
         "strongest; tied ranks are averaged and rounded). Pre-specified ",
         "established genes in bold; $^\\dagger$~member of the plasma-lipoprotein ",
         "pathway subtree, the leading pathway of every replication set.}"), TOP_G),
       "\\label{tab:replication_genes}", "\\footnotesize",
       "\\setlength{\\tabcolsep}{3.5pt}",
       sprintf("\\begin{tabular}{l%s}", paste(rep("rrr", length(ANC)), collapse = "")),
       "\\toprule",
       paste0(" & ", paste(sprintf("\\multicolumn{3}{c}{\\textbf{%s}}", ANC_L[ANC]),
                           collapse = " & "), " \\\\"),
       paste(sprintf("\\cmidrule(lr){%d-%d}", 2 + 3 * (0:4), 4 + 3 * (0:4)), collapse = ""),
       paste0("\\textbf{Gene}", strrep(" & Sc. & GW & MA", length(ANC)), " \\\\"),
       "\\midrule")
for (t in names(TRAITS)) {
  d <- long[trait == t & n_top == length(ANC)]
  ## order: geometric-mean combined rank across groups
  ord <- d[, .(gm = exp(mean(log(rank)))), by = symbol][order(gm)]$symbol
  L <- c(L, sprintf("\\multicolumn{%d}{l}{\\textit{%s}} \\\\", 1 + 3 * length(ANC),
                    TRAITS[[t]]))
  for (s in ord) {
    r <- d[symbol == s][match(ANC, anc)]
    nm <- gene_tex(s, t)
    if (r$lipo_pathway[1]) nm <- paste0(nm, "$^\\dagger$")
    cells <- paste(sprintf("%s & %s & %s", fr(r$rank), fr(r$rank_gwas), fr(r$rank_magma)),
                   collapse = " & ")
    L <- c(L, sprintf("%s & %s \\\\", nm, cells))
  }
}
L <- c(L, "\\bottomrule", "\\end{tabular}", "\\end{table*}")
writeLines(L, file.path(AUT, "tab_replication_genes.tex"))
cat("tab_replication_genes.tex\n")

## ======================= supplementary: top ten per non-European group ======
for (t in names(TRAITS)) {
  L <- c("\\begin{table*}[!t]", "\\centering",
         sprintf(paste0(
           "\\caption{Ten highest-scoring candidate genes for %s in each ",
           "non-European ancestry group (GLGC 2021), in the format of ",
           "Table~\\ref{tab:lipid_candidates}. Rank columns give the gene's position ",
           "within each evidence layer (1 is the strongest): best variant within the ",
           "MAGMA gene annotation (GWAS), MAGMA gene-level $p$-value (MAGMA), and best ",
           "pathway by $z_{\\mathrm{std}}$ among pathway-annotated genes (Pathway). ",
           "Pre-specified established genes in bold.}"), TRAITS[[t]]),
         sprintf("\\label{tab:replication_candidates_%s}", t), "\\footnotesize",
         "\\begin{tabular}{rlrrrrl}", "\\hline",
         paste0("\\textbf{Rank} & \\textbf{Gene} & \\textbf{Score} & \\textbf{GWAS} & ",
                "\\textbf{MAGMA} & \\textbf{Pathway} & \\textbf{Best pathway} \\\\"),
         "\\hline")
  for (a in setdiff(ANC, "eur")) {
    g <- G[[paste(t, a)]][order(rank)][seq_len(N_SUPP)]
    L <- c(L, sprintf("\\multicolumn{7}{l}{\\textit{%s}} \\\\", ANC_L[[a]]))
    for (i in seq_len(nrow(g))) {
      bp <- if (is.na(g$best_pathway[i]) || g$best_pathway[i] == "") "--"
            else pretty_pw(g$best_pathway[i], 40)
      L <- c(L, sprintf("%d & %s & %.2f & %s & %s & %s & %s \\\\", g$rank[i],
                        gene_tex(g$symbol[i], t), g$Score[i], fr(g$rank_gwas[i]),
                        fr(g$rank_magma[i]), fr(g$rank_pathway[i]), bp))
    }
  }
  L <- c(L, "\\hline", "\\end{tabular}", "\\end{table*}")
  writeLines(L, file.path(AUT, sprintf("tab_replication_candidates_%s.tex", t)))
  cat(sprintf("tab_replication_candidates_%s.tex\n", t))
}
