#!/bin/bash
#SBATCH --job-name=treein_ginsi
#SBATCH --cpus-per-task=40
#SBATCH --mem=32GB
#SBATCH --time=2-0:0:0
#SBATCH --output=%j_treein_ginsi.out

# Tree-guided MAFFT G-INS-i of the 590 Oct08 representatives, guide tree (STE2, (Glutamate, rest)),
# run twice: --ep 0.123 with the guide tree from the --ep 0.123 best tree, default --ep with the one
# from the default best tree. Inputs (guide_trees/ in the archive) are the MAFFT
# merge steps and the sequences already in the guide tree's leaf order with their original headers
# (--treein numbers sequences by input order). Writes treein_alignments/representatives_all_<setting>_treein.fasta.
# Submit from this folder (sbatch run_treein_ginsi.sh).
set -euo pipefail

source "$(pwd)/../../config/paths.sh"
DIR="$GPCR_ZENODO_DIR/10_motif_analysis/guide_trees"
OUT="$GPCR_ZENODO_DIR/10_motif_analysis/treein_alignments"
mkdir -p "$OUT"
cd "$DIR"
THREADS=${SLURM_CPUS_PER_TASK:-40}

for setting in ginsi_ep0.123 ginsi; do
  name=representatives_all_Oct08_${setting}
  steps=$(find . -name "${name}_treein.mafft_steps.txt" | head -1)
  seqs=$(find . -name "${name}_treein_input.fasta" | head -1)
  [ -n "$steps" ] && [ -n "$seqs" ] || { echo "input files for $name not found under $DIR"; exit 1; }
  ep=""
  [ "$setting" = ginsi_ep0.123 ] && ep="--ep 0.123"
  echo "$name: $steps, $seqs, $ep"
  mafft --thread "$THREADS" --treein "$steps" --globalpair --maxiterate 1000 $ep "$seqs" > "$OUT/representatives_all_${setting}_treein.fasta"
done
