#!/usr/bin/env bash
#SBATCH --job-name=gpcr_superfamily_mafft
#SBATCH --mem=160GB
#SBATCH --cpus-per-task=40
#SBATCH --nodes=1
#SBATCH --time=0-3:0:0
#SBATCH --array=0-3
#SBATCH --output=%x_%A_%a.out

set -euo pipefail

# MAFFT G-INS-i alignment of the 7TM cores of each sequence set (one set per array task),
# once with --ep 0.123 (main) and once with the default gap extension penalty, and
# ClipKIT kpic-gappy trimming of both. Run after build_representative_fastas.py.

# USER CONFIGURATION ----------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config/paths.sh"
STAGE_DIR="${GPCR_ZENODO_DIR}/08_superfamily_alignment"
THREADS="${SLURM_CPUS_PER_TASK:-40}"

sets=(representatives_all representatives_no_ste2 representatives_no_ste2ste3classCV2 representatives_and_literature)
# END USER CONFIGURATION ------------------------------------------------------

set_name="${sets[${SLURM_ARRAY_TASK_ID:-0}]}"
seqs="${STAGE_DIR}/fasta/${set_name}/${set_name}_TMregion.fasta"
out="${STAGE_DIR}/alignments/${set_name}"
mkdir -p "${STAGE_DIR}/alignments"
echo "== ${set_name}: $(grep -c '>' "${seqs}") sequences"

mafft --thread "${THREADS}" --maxiterate 1000 --globalpair --ep 0.123 "${seqs}" > "${out}_ginsi_ep0.123.fasta"
mafft --thread "${THREADS}" --maxiterate 1000 --globalpair "${seqs}" > "${out}_ginsi.fasta"

for aln in ginsi_ep0.123 ginsi; do
    clipkit "${out}_${aln}.fasta" -o "${out}_${aln}_kpicgappy.fasta" -m kpic-gappy
done
