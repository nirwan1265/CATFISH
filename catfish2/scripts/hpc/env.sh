# =============================================================================
# env.sh -- shared location config, sourced by every script in this directory.
#
# These scripts run in two places with different layouts:
#   cluster : /share/maize/ntanduk/CATFISH        (panel, data and MAGMA here)
#   laptop  : ~/Documents/Github/MAGCAT_data
# Detect rather than assume, and let DATA override either way.
#
#   DATA=/some/other/path bash 13_munge_trait.sh
# =============================================================================

CLUSTER_ROOT="${CLUSTER_ROOT:-/share/maize/ntanduk/CATFISH}"
LOCAL_ROOT="${LOCAL_ROOT:-$HOME/Documents/Github/MAGCAT_data}"

if [ -z "${DATA:-}" ]; then
  if [ -d "$CLUSTER_ROOT" ]; then DATA="$CLUSTER_ROOT"; WHERE="cluster"
  elif [ -d "$LOCAL_ROOT" ]; then DATA="$LOCAL_ROOT";   WHERE="laptop"
  else
    echo "cannot find the data root." >&2
    echo "  looked for: $CLUSTER_ROOT" >&2
    echo "              $LOCAL_ROOT" >&2
    echo "  set it explicitly:  DATA=/path/to/root $0" >&2
    exit 1
  fi
else
  WHERE="DATA override"
fi

BIM="${BIM:-$DATA/g1000_eur.bim}"
LIPIDS="${LIPIDS:-$DATA/lipids}"
MAGMA_BIN="${MAGMA_BIN:-/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/magma_v1.10/magma}"

echo "data root : $DATA  ($WHERE)"
