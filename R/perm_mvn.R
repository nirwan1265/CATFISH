# ============================================================
# 1) Read MAGMA gene-level table (genes.out) for gene order + chr
# ============================================================
magma_read_genes_out <- function(genes_out_file,
                                 gene_col = "GENE",
                                 chr_col  = "CHR") {
  tab <- utils::read.table(genes_out_file, header = TRUE, stringsAsFactors = FALSE)
  if (!all(c(gene_col, chr_col) %in% names(tab))) {
    stop("genes.out is missing expected columns: ", gene_col, ", ", chr_col, call. = FALSE)
  }
  tab[[gene_col]] <- tolower(as.character(tab[[gene_col]]))
  tab
}

# ============================================================
# 2) Read MAGMA gene-gene correlations (BEST CASE)
#    - If you have a file that explicitly lists gene pairs + r
#      (e.g. a .corr/.genes.corr/.genes.r2-like), parse it here.
#    - If you don't have such a file, see fallback below.
# ============================================================
magma_read_gene_cor_pairs <- function(cor_file,
                                      gene1_col = 1,
                                      gene2_col = 2,
                                      r_col     = 3,
                                      sep       = "",
                                      header    = FALSE) {
  tab <- utils::read.table(cor_file, header = header, sep = sep, stringsAsFactors = FALSE)
  if (ncol(tab) < max(gene1_col, gene2_col, r_col)) {
    stop("cor_file doesn't have enough columns.", call. = FALSE)
  }
  out <- data.frame(
    g1 = tolower(as.character(tab[[gene1_col]])),
    g2 = tolower(as.character(tab[[gene2_col]])),
    r  = as.numeric(tab[[r_col]]),
    stringsAsFactors = FALSE
  )
  out <- out[is.finite(out$r) & !is.na(out$g1) & !is.na(out$g2), , drop = FALSE]
  out
}

# ============================================================
# 3) Build R_S for a pathway from gene order + correlation pairs
#    - Missing pairs => 0 correlation (MAGMA also assumes distant=0) :contentReference[oaicite:2]{index=2}
# ============================================================
magma_build_R_for_pathway <- function(genes_S,
                                      genes_out_tab,
                                      cor_pairs = NULL,
                                      gene_col  = "GENE",
                                      chr_col   = "CHR",
                                      make_PD   = TRUE,
                                      pd_eps    = 1e-6) {
  genes_S <- tolower(as.character(genes_S))
  genes_S <- intersect(genes_S, genes_out_tab[[gene_col]])
  G <- length(genes_S)
  if (G < 2L) return(NULL)

  # Keep in genes.out order (often positional-ish by chr)
  ord <- match(genes_S, genes_out_tab[[gene_col]])
  genes_S <- genes_S[order(ord)]

  R <- diag(1, G)
  rownames(R) <- colnames(R) <- genes_S

  if (!is.null(cor_pairs) && nrow(cor_pairs)) {
    # fill symmetric correlations
    idx1 <- match(cor_pairs$g1, genes_S)
    idx2 <- match(cor_pairs$g2, genes_S)
    ok <- which(!is.na(idx1) & !is.na(idx2))
    if (length(ok)) {
      i <- idx1[ok]; j <- idx2[ok]; r <- cor_pairs$r[ok]
      # clip just in case
      r <- pmax(pmin(r, 0.999999), -0.999999)
      R[cbind(i, j)] <- r
      R[cbind(j, i)] <- r
    }
  }

  if (make_PD) {
    # Ensure positive definite for MVN
    # If Matrix is available, use nearPD; otherwise do ridge fix.
    if (requireNamespace("Matrix", quietly = TRUE)) {
      npd <- Matrix::nearPD(R, corr = TRUE)
      R <- as.matrix(npd$mat)
    } else {
      diag(R) <- diag(R) + pd_eps
    }
  }

  R
}

# ============================================================
# 4) Simulate null Z and p for a pathway (block diagonal by chr)
# ============================================================
magma_simulate_null_Zp <- function(genes_S,
                                  genes_out_tab,
                                  R_S,
                                  B,
                                  gene_col = "GENE",
                                  chr_col  = "CHR") {
  genes_S <- rownames(R_S)
  G <- length(genes_S)

  # rmvnorm needs PD covariance; here covariance = correlation matrix
  if (requireNamespace("mvtnorm", quietly = TRUE)) {
    Zmat <- mvtnorm::rmvnorm(n = B, mean = rep(0, G), sigma = R_S)
  } else if (requireNamespace("MASS", quietly = TRUE)) {
    Zmat <- MASS::mvrnorm(n = B, mu = rep(0, G), Sigma = R_S)
  } else {
    stop("Need mvtnorm or MASS for MVN simulation.", call. = FALSE)
  }

  # two-sided p from Z
  Pmat <- 2 * stats::pnorm(-abs(Zmat))

  list(Z = Zmat, P = Pmat, genes = genes_S)
}
