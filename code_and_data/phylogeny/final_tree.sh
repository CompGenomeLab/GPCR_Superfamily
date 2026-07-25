#!/usr/bin/env bash
set -euo pipefail

PHYLOGENY_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_NAME="All_woC-Ste3_clust0.7_cov1_Oct24_minsize2"
BASE_NAME="${RUN_NAME}_linsi_ep0.123_trimends_treein_ginsi_ep0.123"

INPUT_ALIGNMENT="${1:-$PHYLOGENY_DIR/data/alignments/${BASE_NAME}.fasta}"
TRIMMED_ALIGNMENT="${2:-$PHYLOGENY_DIR/data/alignments/${BASE_NAME}_kpicgappy.fasta}"
MODEL_PREFIX="${3:-$PHYLOGENY_DIR/data/model_selection/${BASE_NAME}_bestmodel}"
TREE_PREFIX="${4:-$PHYLOGENY_DIR/data/trees/${BASE_NAME}}"

CLIPKIT_BIN="${CLIPKIT_BIN:-clipkit}"
IQTREE_BIN="${IQTREE_BIN:-iqtree2}"
RAXML_BIN="${RAXML_BIN:-raxml-ng}"
THREADS="${THREADS:-40}"

mkdir -p "$(dirname "$TRIMMED_ALIGNMENT")" "$(dirname "$MODEL_PREFIX")" "$(dirname "$TREE_PREFIX")"

"$CLIPKIT_BIN" "$INPUT_ALIGNMENT" -o "$TRIMMED_ALIGNMENT" -m kpic-gappy

"$IQTREE_BIN" \
  -s "$TRIMMED_ALIGNMENT" \
  -m MF -mfreq FO -cmin 5 -cmax 15 -T 10 \
  --mset JTT --mrate R \
  -pre "$MODEL_PREFIX"

MODEL="$(awk '/Best-fit model:/ {print $3; exit}' "${MODEL_PREFIX}.log")"
test -n "$MODEL"

"$RAXML_BIN" \
  --all \
  --tree 'pars{100}' \
  --msa "$TRIMMED_ALIGNMENT" \
  --model "$MODEL" \
  --bs-trees 'autoMRE{500}' \
  --bs-metric TBE \
  --threads "$THREADS" \
  --workers AUTO \
  --prefix "$TREE_PREFIX"

echo "Best ML tree: ${TREE_PREFIX}.raxml.bestTree"
echo "TBE-supported tree: ${TREE_PREFIX}.raxml.support"
