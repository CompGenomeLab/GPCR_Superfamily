#!/usr/bin/env bash
#SBATCH --job-name=gpcr_deeptmhmm
#SBATCH --mem=36GB
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --nodes=1
#SBATCH --time=0-6:0:0
#SBATCH --array=0-18
#SBATCH --output=%x_%A_%a.out

set -euo pipefail

# DeepTMHMM topology prediction of every clade member, one family per array task.
# Families too large to finish within the wall limit are split with split_fasta_shards.py
# and run with deeptmhmm_predict_shards.sh instead.
# DeepTMHMM is used under its academic license and is not part of the repository; set
# DEEPTMHMM_DIR to the folder holding its predict.py.

# USER CONFIGURATION ----------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config/paths.sh"
DEEPTMHMM_DIR="${DEEPTMHMM_DIR:-/path/to/DeepTMHMM-Academic-License-v1.0}"
DEEPTMHMM_OUT="${DEEPTMHMM_OUT:-${GPCR_ZENODO_DIR}/03_topology_filter/deeptmhmm_output}"

families=(
    classA Olfactory Vomeronasal1 classT Nematode GPR1 cAMP classB1 classB2 Mth
    classF FSLB GP143 TM116 GP157 STE3 classC Vomeronasal2 STE2
)
# END USER CONFIGURATION ------------------------------------------------------

family="${families[${SLURM_ARRAY_TASK_ID:-0}]}"
fasta="${DEEPTMHMM_OUT}/inputs/${family}.fasta"
outdir="${DEEPTMHMM_OUT}/${family}"

mkdir -p "$(dirname "${fasta}")"
gzip -dc "${GPCR_ZENODO_DIR}/02_clade_extraction/clade_sequences/${family}.fasta.gz" > "${fasta}"
rm -rf "${outdir}"

cd "${DEEPTMHMM_DIR}"
python3.9 predict.py --fasta "${fasta}" --output-dir "${outdir}"
