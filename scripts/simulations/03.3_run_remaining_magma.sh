#BSUB -J "redomagma[1-350]%80"
#BSUB -n 1
#BSUB -W 5000
##BSUB -q sara
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=32GB]"
#BSUB -o stdout.%J.%I
#BSUB -e stderr.%J.%I
set -euo pipefail
PIPE_DIR=/share/maize/ntanduk/CATFISH/simulation/permutation_pipeline
source "${PIPE_CONFIG:-$PIPE_DIR/config.sh}"

LINE=$(sed -n "${LSB_JOBINDEX}p" "$PIPE_DIR/missing_jobs.txt")
B=$(echo "$LINE"  | awk '{print $1}')
CHR=$(echo "$LINE" | awk '{print $2}')
echo "[redomagma] perm=$B chr=$CHR"

# safety: confirm the GWAS input exists before MAGMA
WORKDIR="${GWAS_DIR}/perm_${B}/chr${CHR}"
[ -s "${WORKDIR}/assoc.txt" ] || { echo "ERROR: no assoc.txt for perm $B chr $CHR -- run GWAS first" >&2; exit 1; }

Rscript "$PIPE_DIR/03_magma_catfish.R" "$((10#$B))" "$CHR"
