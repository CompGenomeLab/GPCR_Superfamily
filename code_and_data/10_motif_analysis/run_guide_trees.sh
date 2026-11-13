#!/usr/bin/env bash
# Guide trees (STE2, (Glutamate, rest)) and MAFFT merge steps from the two Oct08 all-family best trees.
set -euo pipefail
cd "$(dirname "$0")"
PREP="../../../GPCR_Superfamily/code_and_data/superfamily_alignment/prepare_mafft_treein.py"
FASTA="../../step8_representative_selection/fasta/representatives_all_Oct08/representatives_all_Oct08_TMregion.fasta"
for aln in ginsi_ep0.123 ginsi; do
  name="representatives_all_Oct08_${aln}"
  mkdir -p "$name"
  python make_guide_tree.py "../${name}_kpicgappy.raxml.bestTree" "$name/${name}_guide.nwk"
  python "$PREP" --fasta "$FASTA" --tree "$name/${name}_guide.nwk" --out-prefix "$PWD/$name/${name}_treein"
done
