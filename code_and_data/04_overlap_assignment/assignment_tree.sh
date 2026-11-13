#!/usr/bin/env bash
#SBATCH --job-name=gpcr_assignment_tree
#SBATCH --mem=160GB
#SBATCH --cpus-per-task=40
#SBATCH --nodes=1
#SBATCH --time=0-1:0:0
#SBATCH --output=%x_%j.out

set -euo pipefail

# Assignment tree: the contested sequences and the ten highest-scoring members of each
# family involved (overlaps_and_reps_accession.txt from build_overlap_rep_set.py) are
# fetched, aligned with FAMSA, trimmed at 90% gaps and a tree is inferred with FastTree.

# USER CONFIGURATION ----------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config/paths.sh"
TREE_DIR="${GPCR_ZENODO_DIR}/04_overlap_assignment/assignment_tree"
THREADS="${SLURM_CPUS_PER_TASK:-40}"
# END USER CONFIGURATION ------------------------------------------------------

overlaps="${TREE_DIR}/overlaps_and_reps_accession.txt"
overlap_db="${TREE_DIR}/overlaps_and_reps_accession_db"

mmseqs createsubdb --id-mode 1 "${overlaps}" "${UNIPROT_MMSEQS_DB}" "${overlap_db}"
mmseqs convert2fasta "${overlap_db}" "${overlap_db}.fasta"
sed -i 's/\*//g' "${overlap_db}.fasta"

famsa "${overlap_db}.fasta" "${overlap_db}_FAMSA.fasta" -t "${THREADS}"
trimal -in "${overlap_db}_FAMSA.fasta" -out "${overlap_db}_FAMSA_gt0.1.fasta" -gt 0.1

export OMP_NUM_THREADS="${THREADS}"
FastTreeMP "${overlap_db}_FAMSA_gt0.1.fasta" > "${overlap_db}_FAMSA_gt0.1.tree"
