#!/usr/bin/env bash
set -euo pipefail

PHYLOGENY_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_NAME="All_woC-Ste3_clust0.7_cov1_Oct24_minsize2"

INPUT_FASTA="${1:-$PHYLOGENY_DIR/data/representatives/${RUN_NAME}.fasta}"
OUTPUT_FASTA="${2:-$PHYLOGENY_DIR/data/alignments/${RUN_NAME}_linsi_ep0.123.fasta}"
MAFFT_BIN="${MAFFT_BIN:-mafft}"
THREADS="${THREADS:-40}"

mkdir -p "$(dirname "$OUTPUT_FASTA")"

"$MAFFT_BIN" \
  --thread "$THREADS" \
  --localpair \
  --maxiterate 1000 \
  --ep 0.123 \
  "$INPUT_FASTA" > "$OUTPUT_FASTA"

echo "L-INS-i alignment: $OUTPUT_FASTA"
