#!/usr/bin/env bash
set -euo pipefail

PHYLOGENY_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_NAME="All_woC-Ste3_clust0.7_cov1_Oct24_minsize2"

INPUT_ALIGNMENT="${1:-$PHYLOGENY_DIR/data/alignments/${RUN_NAME}_linsi_ep0.123_trimends.fasta}"
OUTPUT_PREFIX="${2:-$PHYLOGENY_DIR/data/alignments/${RUN_NAME}_linsi_ep0.123_trimends_kpicgappy}"
MODEL_PREFIX="${3:-$PHYLOGENY_DIR/data/model_selection/${RUN_NAME}_linsi_ep0.123_trimends_kpicgappy_bestmodel}"
TREE_PREFIX="${4:-$PHYLOGENY_DIR/data/trees/${RUN_NAME}_linsi_ep0.123_trimends_kpicgappy}"

CLIPKIT_BIN="${CLIPKIT_BIN:-clipkit}"
IQTREE_BIN="${IQTREE_BIN:-iqtree2}"
RAXML_BIN="${RAXML_BIN:-raxml-ng}"
THREADS="${THREADS:-40}"

mkdir -p "$(dirname "$OUTPUT_PREFIX")" "$(dirname "$MODEL_PREFIX")" "$(dirname "$TREE_PREFIX")"

TRIMMED_ALIGNMENT="${OUTPUT_PREFIX}.fasta"
"$CLIPKIT_BIN" "$INPUT_ALIGNMENT" -o "$TRIMMED_ALIGNMENT" -m kpic-gappy

"$IQTREE_BIN" \
  -s "$TRIMMED_ALIGNMENT" \
  -m MF -mfreq FO -cmin 5 -cmax 15 -T 10 \
  --mset JTT --mrate R \
  -pre "$MODEL_PREFIX"

MODEL="$(awk '/Best-fit model:/ {print $3; exit}' "${MODEL_PREFIX}.log")"
test -n "$MODEL"

"$RAXML_BIN" \
  --search \
  --tree 'pars{20}' \
  --msa "$TRIMMED_ALIGNMENT" \
  --model "$MODEL" \
  --threads "$THREADS" \
  --workers AUTO \
  --prefix "$TREE_PREFIX"

echo "Preliminary tree: ${TREE_PREFIX}.raxml.bestTree"
