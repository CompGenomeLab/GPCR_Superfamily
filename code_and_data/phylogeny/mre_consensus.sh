#!/usr/bin/env bash
set -euo pipefail

PHYLOGENY_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_NAME="All_woC-Ste3_clust0.7_cov1_Oct24_minsize2"
BASE_NAME="${RUN_NAME}_linsi_ep0.123_trimends_treein_ginsi_ep0.123"

ML_TREES="${1:-$PHYLOGENY_DIR/data/trees/${BASE_NAME}.raxml.mlTrees}"
BOOTSTRAPS="${2:-$PHYLOGENY_DIR/data/trees/${BASE_NAME}.raxml.bootstraps}"
CONSENSUS_PREFIX="${3:-$PHYLOGENY_DIR/data/trees/${BASE_NAME}_consensus}"
SUPPORT_PREFIX="${4:-$PHYLOGENY_DIR/data/trees/${BASE_NAME}_consensus_TBE}"
RAXML_BIN="${RAXML_BIN:-raxml-ng}"
THREADS="${THREADS:-1}"

mkdir -p "$(dirname "$CONSENSUS_PREFIX")" "$(dirname "$SUPPORT_PREFIX")"

"$RAXML_BIN" \
  --consense MRE \
  --tree "$ML_TREES" \
  --prefix "$CONSENSUS_PREFIX"

"$RAXML_BIN" \
  --redo \
  --support \
  --tree "${CONSENSUS_PREFIX}.raxml.consensusTreeMRE" \
  --bs-trees "$BOOTSTRAPS" \
  --bs-metric TBE \
  --threads "$THREADS" \
  --prefix "$SUPPORT_PREFIX"

echo "MRE consensus: ${CONSENSUS_PREFIX}.raxml.consensusTreeMRE"
echo "TBE-supported MRE consensus: ${SUPPORT_PREFIX}.raxml.support"
