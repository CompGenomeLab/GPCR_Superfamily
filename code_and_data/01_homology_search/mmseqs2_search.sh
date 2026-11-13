#!/usr/bin/env bash
#SBATCH --job-name=gpcr_mmseqs
#SBATCH --mem=160GB
#SBATCH --cpus-per-task=40
#SBATCH --nodes=1
#SBATCH --time=5-00:00:00
#SBATCH --array=0-1
#SBATCH --output=%x_%A_%a.out

set -euo pipefail

# USER CONFIGURATION ----------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config/paths.sh"
SEARCH_RESULTS_DIR="${GPCR_ZENODO_DIR}/01_homology_search"
THREADS="${SLURM_CPUS_PER_TASK:-40}"

# Families searched from organism-specific FASTA seeds rather than HMM profiles.
families=(GPR1_YEAST FSLB_DICDI)
# END USER CONFIGURATION ------------------------------------------------------

task_id="${SLURM_ARRAY_TASK_ID:-0}"
if (( task_id >= ${#families[@]} )); then
    echo "Array index ${task_id} is outside the families list" >&2
    exit 2
fi

family="${families[$task_id]}"
family_dir="${SEARCH_RESULTS_DIR}/${family}"
querydb="${family_dir}/${family}_query_db"
outdb="${family_dir}/${family}_UniProt_s7.5_filter0_eprofile0.1_db"
outm8="${outdb}.m8"
headers="${outdb}_headers10000.txt"
retrieved_db="${family_dir}/${family}_db"

mkdir -p "${family_dir}" "${family_dir}/tmp"
mmseqs createdb \
    "${SCRIPT_DIR}/input_sequences/${family}/${family}.fasta" \
    "${querydb}" --compressed 0
mmseqs search "${querydb}" "${UNIPROT_MMSEQS_DB}" "${outdb}" "${family_dir}/tmp" \
    --threads "${THREADS}" -e 1000 --e-profile 0.1 \
    --sub-mat blosum62.out --seed-sub-mat VTML80.out -s 7.5 \
    --max-seqs 1000000 --remove-tmp-files 1 \
    --prefilter-mode 0 --num-iterations 5

mmseqs convertalis "${querydb}" "${UNIPROT_MMSEQS_DB}" "${outdb}" "${outm8}" \
    --threads "${THREADS}" \
    --format-output query,target,theader,pident,fident,nident,alnlen,qlen,tlen,mismatch,raw,bits,qcov,tcov,evalue

head -n 10000 "${outm8}" | cut -f2 > "${headers}"
mmseqs createsubdb --id-mode 1 "${headers}" "${UNIPROT_MMSEQS_DB}" "${retrieved_db}"
mmseqs convert2fasta "${retrieved_db}" "${retrieved_db}.fasta"
sed -i 's/\*//g' "${retrieved_db}.fasta"

famsa "${retrieved_db}.fasta" "${retrieved_db}_FAMSA.fasta" -t "${THREADS}"
trimal -in "${retrieved_db}_FAMSA.fasta" \
    -out "${retrieved_db}_FAMSA_gt0.1.fasta" -gt 0.1
export OMP_NUM_THREADS="${THREADS}"
FastTreeMP "${retrieved_db}_FAMSA_gt0.1.fasta" > "${retrieved_db}_FAMSA_gt0.1.tree"
