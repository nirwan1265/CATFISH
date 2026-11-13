#!/bin/bash
# =============================================================================
# Step 30: install the ARTP2 R package (sARTP) into ARTP_RLIB under usrapps.
# github.com/zhangh12/ARTP3 serves the ARTP2 repository (version 0.9.49);
# needs a C++11 compiler and the R packages Formula and data.table.
#   bsub < 30_artp_install_LSF.sh
# If the compute node has no internet, this job stops at the download and
# says so; nothing is installed in the home directory.
# =============================================================================
#BSUB -J cmp_artp_install
#BSUB -q sara
#BSUB -n 1
#BSUB -W 120
#BSUB -o logs/artp_install.%J.out
#BSUB -e logs/artp_install.%J.err
set -eo pipefail
source /share/maize/ntanduk/CATFISH/compare/config.sh
mkdir -p "$ARTP_RLIB"
export R_LIBS_USER="$ARTP_RLIB"
SRC="${ARTP_RLIB}/src_ARTP2"; mkdir -p "$SRC"
TGZ="${SRC}/ARTP2-master.tar.gz"
if [ ! -s "$TGZ" ]; then
  curl -fsSL -o "$TGZ" https://github.com/zhangh12/ARTP2/archive/refs/heads/master.tar.gz \
    || { echo "download failed (no internet on the compute node?)" >&2; exit 1; }
fi
TGZ="$TGZ" Rscript -e '
lib <- Sys.getenv("ARTP_RLIB"); .libPaths(c(lib, .libPaths()))
need <- c("Formula", "data.table")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) install.packages(miss, lib = lib, repos = "https://cloud.r-project.org")
install.packages(Sys.getenv("TGZ"), lib = lib, repos = NULL, type = "source")
suppressPackageStartupMessages(library(ARTP2, lib.loc = lib))
cat("ARTP2", as.character(packageVersion("ARTP2", lib.loc = lib)), "installed in", lib, "\n")
print(args(ARTP2::sARTP))
'
