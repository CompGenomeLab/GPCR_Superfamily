#!/usr/bin/env bash
#SBATCH --job-name=gpcr_deeptmhmm_shard
#SBATCH --mem=50GB
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --nodes=1
#SBATCH --time=0-15:0:0
#SBATCH --output=%x_%A_%a.out

set -euo pipefail

# For large families only, where the whole FASTA does not finish within the wall limit.
# Split the family's FASTA first (shards are <family>_shard<N>.fasta with N starting at 1):
#   python split_fasta_shards.py $DEEPTMHMM_OUT/inputs/<family>.fasta --shard-size 25000
# then submit with the family as an argument and the array range covering its shards:
#   sbatch --array=1-6 deeptmhmm_predict_shards.sh classA
# Each shard writes to <family>/shard<N>, because predict.py always writes the same
# filenames into --output-dir and parallel shards would overwrite one another.

# USER CONFIGURATION ----------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config/paths.sh"
DEEPTMHMM_DIR="${DEEPTMHMM_DIR:-/path/to/DeepTMHMM-Academic-License-v1.0}"
DEEPTMHMM_OUT="${DEEPTMHMM_OUT:-${GPCR_ZENODO_DIR}/03_topology_filter/deeptmhmm_output}"
# END USER CONFIGURATION ------------------------------------------------------

family="${1:-}"
shard="${SLURM_ARRAY_TASK_ID:-1}"
if [[ -z "${family}" ]]; then
    echo "ERROR: no family given. Usage: sbatch --array=1-N deeptmhmm_predict_shards.sh <family>" >&2
    exit 1
fi

fasta="${DEEPTMHMM_OUT}/inputs/shards/${family}_shard${shard}.fasta"
outdir="${DEEPTMHMM_OUT}/${family}/shard${shard}"
if [[ ! -f "${fasta}" ]]; then
    echo "ERROR: shard not found: ${fasta}" >&2
    exit 1
fi
rm -rf "${outdir}"

echo "family: ${family}"
echo "shard:  ${shard}"
echo "fasta:  ${fasta} ($(grep -c '^>' "${fasta}") sequences)"
echo "outdir: ${outdir}"

cd "${DEEPTMHMM_DIR}"
python3.9 predict.py --fasta "${fasta}" --output-dir "${outdir}"
