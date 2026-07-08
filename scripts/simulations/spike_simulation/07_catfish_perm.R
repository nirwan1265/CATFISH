#!/usr/bin/env Rscript

suppressMessages({ library(data.table) })
data.table::setDTthreads(1L)

getcfg <- function(k){ v<-Sys.getenv(k); if(v=="") stop(paste("config missing:",k)); v }
parse_num_grid <- function(x, default) {
  if (!nzchar(x)) return(default)
  parts <- unlist(strsplit(gsub("[,;]", " ", x), "\\s+"))
  vals <- suppressWarnings(as.numeric(parts))
  vals <- vals[is.finite(vals) & !is.na(vals)]
  if (!length(vals)) default else vals
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("need permutation index")
B <- as.integer(args[1]); BPAD <- sprintf("%03d", B)

MAGMA_OUT    <- getcfg("MAGMA_OUT")
CORPAIRS     <- getcfg("CORPAIRS")
CATFISH_OUT  <- getcfg("CATFISH_OUT")
CATFISH_REPO <- getcfg("CATFISH_REPO")
GENE_MODEL   <- getcfg("GENE_MODEL")
CHROMS       <- strsplit(trimws(getcfg("CHROMS")), "\\s+")[[1]]
SEED         <- as.integer(getcfg("MASTER_SEED")) + B

MVN_CALIBRATE_COMPONENTS <- tolower(Sys.getenv("MVN_CALIBRATE_COMPONENTS", "false")) %in% c("1","true","yes")
PATCH_EMP_NULL_LOWER <- tolower(Sys.getenv("PATCH_EMP_NULL_LOWER", "false")) %in% c("1","true","yes")
TAU_GRID <- parse_num_grid(Sys.getenv("TAU_GRID", ""), c(0.1, 0.05, 0.01))
B_PERM <- suppressWarnings(as.integer(Sys.getenv("CATFISH_B_PERM", "10000")))
if (!is.finite(B_PERM) || is.na(B_PERM) || B_PERM < 1L) B_PERM <- 10000L
TAIL_MODE <- trimws(Sys.getenv("TAIL_MODE", "empirical"))
if (!nzchar(TAIL_MODE)) TAIL_MODE <- "empirical"
TAIL_SWITCH_EXCEED <- suppressWarnings(as.integer(Sys.getenv("TAIL_SWITCH_EXCEED", "10")))
if (!is.finite(TAIL_SWITCH_EXCEED) || is.na(TAIL_SWITCH_EXCEED) || TAIL_SWITCH_EXCEED < 1L) TAIL_SWITCH_EXCEED <- 10L
TAIL_GPD_K <- suppressWarnings(as.integer(Sys.getenv("TAIL_GPD_K", "250")))
if (!is.finite(TAIL_GPD_K) || is.na(TAIL_GPD_K) || TAIL_GPD_K < 10L) TAIL_GPD_K <- 250L
TAIL_MIN_B <- suppressWarnings(as.integer(Sys.getenv("TAIL_MIN_B", "10000")))
if (!is.finite(TAIL_MIN_B) || is.na(TAIL_MIN_B) || TAIL_MIN_B < 1L) TAIL_MIN_B <- 10000L
TAIL_MIN_TAIL <- suppressWarnings(as.integer(Sys.getenv("TAIL_MIN_TAIL", "50")))
if (!is.finite(TAIL_MIN_TAIL) || is.na(TAIL_MIN_TAIL) || TAIL_MIN_TAIL < 5L) TAIL_MIN_TAIL <- 50L

PATHWAY_FILE <- Sys.getenv("PATHWAY_FILE",
  file.path(CATFISH_REPO, "inst/extdata/pathway/sorghumbicolorcyc_pathways.20230103.SORBI"))

suppressMessages(devtools::load_all(CATFISH_REPO, quiet = TRUE))
model_tag <- gsub("[^A-Za-z0-9]+", "_", GENE_MODEL)

cat(sprintf(
  "[catfish] perm=%d repo=%s mvn_calibrate_components=%s patch_emp_null_lower=%s threads=%s tau_grid=%s B_perm=%d tail_mode=%s tail_min_B=%d\n",
  B, CATFISH_REPO, MVN_CALIBRATE_COMPONENTS, PATCH_EMP_NULL_LOWER,
  Sys.getenv("CATFISH_THREADS", "1"),
  paste(TAU_GRID, collapse = ","), B_PERM, TAIL_MODE, TAIL_MIN_B
))

if (PATCH_EMP_NULL_LOWER) {
  ns <- asNamespace("CATFISH")
  unlockBinding(".emp_p_null_lower", ns)
  assign(".emp_p_null_lower", function(p_null) {
    p2 <- as.numeric(p_null)
    B <- length(p2)
    if (B < 1L) return(numeric(0))
    bad <- !is.finite(p2) | is.na(p2)
    if (any(bad)) p2[bad] <- Inf
    rmax <- rank(p2, ties.method = "max")
    as.numeric(rmax) / B
  }, envir = ns)
  lockBinding(".emp_p_null_lower", ns)
  cat("[catfish] patched .emp_p_null_lower for simulation run\n")
}

gene_files <- file.path(MAGMA_OUT, sprintf("perm_%s", BPAD), sprintf("chr%s", CHROMS),
                        sprintf("perm_%s_chr%s.%s.genes.out", BPAD, CHROMS, model_tag))
missing <- CHROMS[!file.exists(gene_files)]
if (length(missing)) stop(sprintf("perm %s INCOMPLETE: missing chr %s", BPAD, paste(missing, collapse=", ")))

read_one <- function(path) {
  x <- read.table(path, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
  if (!"GENE"  %in% names(x)) stop("no GENE col in ", path)
  if (!"ZSTAT" %in% names(x)) stop("no ZSTAT col in ", path)
  pcol <- if ("P_MULTI" %in% names(x)) "P_MULTI" else if ("P" %in% names(x)) "P" else stop("no P col")
  x$P <- x[[pcol]]; x
}
genes <- do.call(rbind, lapply(gene_files, read_one))
genes <- genes[order(genes$GENE, genes$P), , drop=FALSE]
genes <- genes[!duplicated(genes$GENE), , drop=FALSE]; rownames(genes) <- NULL

if (!file.exists(PATHWAY_FILE)) stop("PATHWAY_FILE not found: ", PATHWAY_FILE)
pr <- utils::read.delim(PATHWAY_FILE, header=TRUE, stringsAsFactors=FALSE, check.names = FALSE)
pathways <- unique(data.frame(
  pathway_id   = pr[["Pathway-id"]],
  pathway_name = pr[["Pathway-name"]],
  gene_id      = pr[["Gene-name"]],
  stringsAsFactors = FALSE))

omni <- catfish_omni2_pathways(
  gene_results              = genes,
  species                   = NULL,
  pathways                  = pathways,
  gene_col                  = "GENE",
  p_raw_col                 = "P",
  p_adj_col                 = NULL,
  z_col                     = "ZSTAT",
  tau_grid                  = TAU_GRID,
  min_p                     = 1e-15,
  do_fix                    = TRUE,
  stouffer_alternative      = "greater",
  omnibus                   = "ACAT",
  perm_mode                 = "mvn",
  B_perm                    = B_PERM,
  magma_cor_file            = CORPAIRS,
  make_PD                   = TRUE,
  mvn_marginal              = "uniform",
  mvn_calibrate_components  = MVN_CALIBRATE_COMPONENTS,
  tail_mode                 = TAIL_MODE,
  tail_switch_exceed        = TAIL_SWITCH_EXCEED,
  tail_gpd_k                = TAIL_GPD_K,
  tail_min_B                = TAIL_MIN_B,
  tail_min_tail             = TAIL_MIN_TAIL,
  n_threads                 = as.integer(Sys.getenv("CATFISH_THREADS", "1")),
  min_genes                 = 2L,
  seed                      = SEED,
  output                    = FALSE
)

dir.create(CATFISH_OUT, recursive=TRUE, showWarnings=FALSE)
keep <- intersect(c("pathway_id","pathway_name","n_genes",
                    "omni_p_final","omni_p_mvn","omni_p_analytic",
                    "omni_p_mvn_compcal",
                    "acat_p","fisher_p","minp_p_analytic",
                    "stouffer_p_analytic","tfisher_p_analytic","tfisher_p_min","tau_hat",
                    "acat_p_mvn_cal","fisher_p_mvn_cal","tfisher_p_mvn_cal",
                    "minp_p_mvn_cal","stouffer_p_mvn_cal"), names(omni))
out <- as.data.table(omni)[, ..keep]; out[, perm := B]
fwrite(out, file.path(CATFISH_OUT, sprintf("pathway_pvals_perm_%s.csv", BPAD)), nThread = 1L)
cat(sprintf("[catfish] perm %d done: %d pathways\n", B, nrow(out)))
