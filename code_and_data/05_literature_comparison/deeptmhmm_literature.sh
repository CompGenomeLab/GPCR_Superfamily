#!/usr/bin/env bash
#SBATCH --job-name=gpcr_deeptmhmm_literature
#SBATCH --mem=36GB
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --nodes=1
#SBATCH --time=0-1:0:0
#SBATCH --output=%x_%j.out

set -euo pipefail

# DeepTMHMM topology prediction of the 365 literature sequences tested for deeper origins;
# predicted_topologies.3line is the input of literature_7tm_gate.py.
# DeepTMHMM is used under its academic license and is not part of the repository; set
# DEEPTMHMM_DIR to the folder holding its predict.py.

# USER CONFIGURATION ----------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config/paths.sh"
DEEPTMHMM_DIR="${DEEPTMHMM_DIR:-/path/to/DeepTMHMM-Academic-License-v1.0}"
OUT_DIR="${GPCR_ZENODO_DIR}/05_literature_comparison"
# END USER CONFIGURATION ------------------------------------------------------

outdir="${OUT_DIR}/deeper_origin_tested_sequences_deeptmhmm"

cd "${DEEPTMHMM_DIR}"
python3.9 predict.py --fasta "${OUT_DIR}/deeper_origin_tested_sequences.fasta" --output-dir "${outdir}"
cp "${outdir}/predicted_topologies.3line" "${OUT_DIR}/predicted_topologies.3line"
