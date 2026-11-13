#!/usr/bin/env bash
#SBATCH --job-name=gpcr_topology_alignment
#SBATCH --mem=160GB
#SBATCH --cpus-per-task=40
#SBATCH --nodes=1
#SBATCH --time=0-3:0:0
#SBATCH --array=0-18
#SBATCH --output=%x_%A_%a.out

set -euo pipefail

# Topology alignment of each family's single-core sequences: the 7TM core (+5 residues) of
# every sequence is aligned with FAMSA, every aligned residue is replaced by its topology
# state and columns with 90% or more gaps are removed. Run after write_core_membership.py.

# USER CONFIGURATION ----------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config/paths.sh"
STAGE_DIR="${GPCR_ZENODO_DIR}/03_topology_filter"
WORK_DIR="${WORK_DIR:-${STAGE_DIR}/work}"
THREADS="${SLURM_CPUS_PER_TASK:-40}"

families=(
    classA Olfactory Vomeronasal1 classT Nematode GPR1 cAMP classB1 classB2 Mth
    classF FSLB GP143 TM116 GP157 STE3 classC Vomeronasal2 STE2
)
# END USER CONFIGURATION ------------------------------------------------------

family="${families[${SLURM_ARRAY_TASK_ID:-0}]}"
mkdir -p "${WORK_DIR}" "${STAGE_DIR}/topology_alignments"

threeline="${WORK_DIR}/deeptmhmm_predictions_oriented.3line"
[[ -s "${threeline}" ]] || gzip -dc "${STAGE_DIR}/deeptmhmm_predictions_oriented.3line.gz" > "${threeline}"

tm_db="${WORK_DIR}/${family}_genes_7TMfilter_db"
gzip -dc "${STAGE_DIR}/core_sequences/${family}.fasta.gz" > "${tm_db}.fasta"

python3 "${SCRIPT_DIR}/../08_superfamily_alignment/extract_tm_region.py" \
    --fasta "${tm_db}.fasta" \
    --topology "${threeline}" \
    --out "${tm_db}_TMonly.fasta" \
    --boundaries "${tm_db}_TMregion_boundaries.tsv" \
    --flank 5

export OMP_NUM_THREADS="${THREADS}"
famsa "${tm_db}_TMonly.fasta" "${tm_db}_TMonly_FAMSA.fasta" -t "${THREADS}"

python3 "${SCRIPT_DIR}/project_topology_onto_alignment.py" \
    --alignment "${tm_db}_TMonly_FAMSA.fasta" \
    --topology "${threeline}" \
    --boundaries "${tm_db}_TMregion_boundaries.tsv" \
    --out "${tm_db}_FAMSA_topology_converted.fasta"

trimal -in "${tm_db}_FAMSA_topology_converted.fasta" \
    -out "${STAGE_DIR}/topology_alignments/${family}_genes_7TMfilter_db_FAMSA_topology_converted_gt0.1.fasta" \
    -gt 0.1
