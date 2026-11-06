#!/usr/bin/env Rscript
## =============================================================================
## make_tutorial.R -- builds tutorial/CATFISH_tutorial.md (and the .pdf, if
## pandoc and xelatex are installed) from tutorial_body.md plus a function
## reference generated from the package's own Rd files, so every parameter is
## listed exactly as documented. Rerun after devtools::document(), from the
## package root:
##   Rscript tutorial/make_tutorial.R
## =============================================================================
pkg <- "."
ns  <- readLines(file.path(pkg, "NAMESPACE"))
exp <- sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", ns, value = TRUE))
tag <- function(x) { t <- attr(x, "Rd_tag"); if (is.null(t)) "" else t }
render <- function(x) {
  if (is.character(x)) return(paste(x, collapse = ""))
  t <- tag(x)
  inner <- function() paste(vapply(x, render, ""), collapse = "")
  switch(t,
    "\\code" = paste0("`", gsub("`", "", inner()), "`"),
    "\\link" = inner(), "\\linkS4class" = inner(),
    "\\eqn" = { s <- render(x[[1]]); s <- gsub("\\\\tau", "tau", s); s <- gsub("[{}\\\\]", "", s); s },
    "\\emph" = paste0("*", inner(), "*"),
    "\\strong" = paste0("**", inner(), "**"),
    "\\url" = inner(), "\\href" = render(x[[2]]),
    "\\item" = { if (length(x) >= 2) paste0(" • ", render(x[[1]]), ": ", render(x[[2]])) else " • " },
    "\\itemize" = , "\\enumerate" = , "\\describe" = paste0(" ", inner()),
    "\\dots" = "...", "\\ldots" = "...", "\\R" = "R",
    "\\cr" = " ", "\\sQuote" = paste0("'", inner(), "'"), "\\dQuote" = paste0('"', inner(), '"'),
    "\\pkg" = inner(), "\\var" = inner(), "\\file" = paste0("`", inner(), "`"),
    "\\samp" = paste0("`", inner(), "`"), "\\verb" = paste0("`", inner(), "`"),
    "\\deqn" = render(x[[1]]),
    "COMMENT" = "",
    inner())
}
clean <- function(s) { s <- gsub("\\s+", " ", s); s <- gsub("\\|", "\\\\|", s); trimws(s) }
out <- character()
files <- list.files(file.path(pkg, "man"), "\\.Rd$", full.names = TRUE)
db <- lapply(files, tools::parse_Rd)
names(db) <- files
alias_of <- function(rd) unlist(lapply(rd[vapply(rd, tag, "") == "\\alias"], render))
order_fun <- c("catfish_run",
  "gff3_to_geneloc", "magma_path", "magma_annotate", "magma_gene", "magma_read_genes",
  "magma_genesraw_to_cor_pairs_banded",
  "catfish_load_pathways", "catfish_pmn_file", "catfish_load_gmt", "catfish_load_biocyc",
  "catfish_write_set_annot",
  "catfish_acat_pathways", "catfish_fisher_pathways", "catfish_soft_tfisher_adaptive_pathways",
  "catfish_minp_pathways", "catfish_stoufferZ_pathways",
  "catfish_wy_pathways", "catfish_omni2_pathways", "gene_prioritize",
  "get_gene_lengths", "catfish_adjust_gene_p", "fix_p_for_acat")
if (!setequal(order_fun, exp))
  stop("Exported functions and the tutorial list differ: ",
       paste(c(setdiff(exp, order_fun), setdiff(order_fun, exp)), collapse = ", "),
       " -- add new functions to order_fun.", call. = FALSE)
for (fn in order_fun) {
  rd <- Filter(function(r) fn %in% alias_of(r), db)[[1]]
  sec <- function(t) rd[vapply(rd, tag, "") == t]
  title <- clean(render(sec("\\title")[[1]]))
  usage <- paste(vapply(sec("\\usage")[[1]], render, ""), collapse = "")
  usage <- trimws(gsub("\n{2,}", "\n", usage))
  out <- c(out, paste0("### `", fn, "()`"), "", paste0("*", title, "*"), "",
           "```r", usage, "```", "")
  args <- sec("\\arguments")
  if (length(args)) {
    items <- args[[1]][vapply(args[[1]], tag, "") == "\\item"]
    out <- c(out, "| Parameter | Description |", "|---|---|")
    for (it in items) {
      nm <- clean(render(it[[1]])); ds <- clean(render(it[[2]]))
      nm <- paste(paste0("`", trimws(strsplit(nm, ",")[[1]]), "`"), collapse = ", ")
      out <- c(out, paste0("| ", nm, " | ", ds, " |"))
    }
    out <- c(out, "")
  }
}
body <- readLines(file.path(pkg, "tutorial", "tutorial_body.md"))
md   <- file.path(pkg, "tutorial", "CATFISH_tutorial.md")
writeLines(c(body, out), md)
message("wrote ", md)
if (nzchar(Sys.which("pandoc")) && nzchar(Sys.which("xelatex"))) {
  pdf <- sub("\\.md$", ".pdf", md)
  st <- system2("pandoc", c(shQuote(md), "-o", shQuote(pdf), "--pdf-engine=xelatex",
                            "--toc", "--toc-depth=2", "-V", "geometry:margin=2cm",
                            "-V", "fontsize=10pt", "-V", "colorlinks=true",
                            "--highlight-style=tango"))
  if (st == 0) message("wrote ", pdf) else warning("pandoc failed; only the .md was written")
}
