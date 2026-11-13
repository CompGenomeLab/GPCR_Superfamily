#!/usr/bin/env bash
#SBATCH --job-name=gpcr_literature_mmseqs
#SBATCH --mem=160GB
#SBATCH --cpus-per-task=40
#SBATCH --nodes=1
#SBATCH --time=0:60:0
#SBATCH --array=0-3
#SBATCH --output=%x_%A_%a.out

set -euo pipefail

# MMseqs2 search of the three literature query sets (written by write_unmatched_fasta.py,
# build_remaining_fasta.py and write_unsearched_fasta.py) and of the sequences tested for
# deeper origins against the searched proteomes, one set per array task.

# USER CONFIGURATION ----------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config/paths.sh"
OUT_DIR="${GPCR_ZENODO_DIR}/05_literature_comparison"
THREADS="${SLURM_CPUS_PER_TASK:-40}"

query_sets=(unmatched_literature_sequences remaining_literature_sequences unsearched_literature_sequences
            deeper_origin_tested_sequences)
# END USER CONFIGURATION ------------------------------------------------------

name="${query_sets[${SLURM_ARRAY_TASK_ID:-0}]}"
query="${OUT_DIR}/${name}.fasta"
querydb="${OUT_DIR}/${name}_db"
outdb="${OUT_DIR}/${name}_search_db"

mmseqs createdb "${query}" "${querydb}" --compressed 0 -v 3

mmseqs search "${querydb}" "${UNIPROT_MMSEQS_DB}" "${outdb}" "${OUT_DIR}/tmp_${name}" --threads "${THREADS}" \
    -e 10 --sub-mat blosum62.out --seed-sub-mat VTML80.out \
    -s 7.5 --max-seqs 100 --remove-tmp-files 1

mmseqs convertalis "${querydb}" "${UNIPROT_MMSEQS_DB}" "${outdb}" "${OUT_DIR}/${name}.m8" --threads "${THREADS}" \
    --format-output query,target,theader,pident,fident,nident,alnlen,qlen,tlen,mismatch,raw,bits,qcov,tcov,evalue
