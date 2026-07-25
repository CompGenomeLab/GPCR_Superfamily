#!/usr/bin/env bash
set -euo pipefail

STAGE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_NAME="All_clust0.7_cov1_Oct25_minsize2"

ACCESSION_LIST="${1:-$STAGE_DIR/data/representatives/${RUN_NAME}.txt}"
SOURCE_DB="${2:-${UNIPROT_MMSEQS_DB:?Set UNIPROT_MMSEQS_DB or provide the source MMseqs2 database as argument 2.}}"
OUTPUT_DB="${3:-$STAGE_DIR/data/representatives/${RUN_NAME}}"
OUTPUT_FASTA="${4:-$STAGE_DIR/data/representatives/${RUN_NAME}.fasta}"
MMSEQS_BIN="${MMSEQS_BIN:-mmseqs}"

mkdir -p "$(dirname "$OUTPUT_DB")" "$(dirname "$OUTPUT_FASTA")"

"$MMSEQS_BIN" createsubdb "$ACCESSION_LIST" "$SOURCE_DB" "$OUTPUT_DB"
"$MMSEQS_BIN" convert2fasta "$OUTPUT_DB" "$OUTPUT_FASTA"
sed -i 's/\*//g' "$OUTPUT_FASTA"

echo "Representative FASTA: $OUTPUT_FASTA"
