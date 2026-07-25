#!/usr/bin/env bash
set -euo pipefail

PHYLOGENY_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_NAME="All_woC-Ste3_clust0.7_cov1_Oct24_minsize2"
GUIDE_PREFIX="$PHYLOGENY_DIR/data/alignments/${RUN_NAME}_linsi_ep0.123_trimends_treein"

GUIDE_STEPS="${1:-${GUIDE_PREFIX}.mafft_steps.txt}"
REORDERED_FASTA="${2:-${GUIDE_PREFIX}.reordered.fasta}"
OUTPUT_FASTA="${3:-$PHYLOGENY_DIR/data/alignments/${RUN_NAME}_linsi_ep0.123_trimends_treein_ginsi_ep0.123.fasta}"
MAFFT_BIN="${MAFFT_BIN:-mafft}"
THREADS="${THREADS:-40}"

mkdir -p "$(dirname "$OUTPUT_FASTA")"

"$MAFFT_BIN" \
  --thread "$THREADS" \
  --treein "$GUIDE_STEPS" \
  --globalpair \
  --maxiterate 1000 \
  --ep 0.123 \
  "$REORDERED_FASTA" > "$OUTPUT_FASTA"

echo "Tree-guided G-INS-i alignment: $OUTPUT_FASTA"
