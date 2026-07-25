#!/usr/bin/env bash
#SBATCH --job-name=familydb
#SBATCH --mem=160GB
#SBATCH --cpus-per-task=40
#SBATCH --nodes=1
#SBATCH --output=%j_familydb.out
#SBATCH --time=0-2:0:0
#SBATCH --account=PAS1794
#SBATCH --array=0-17

set -euo pipefail

PROJECT_DIR="${GPCR_ANALYSIS_DIR:-/path/to/GPCR_Superfamily}"
UNIPROT_MMSEQS_DB="${UNIPROT_MMSEQS_DB:-/path/to/UniprotRefEuk}"
THREADS="${SLURM_CPUS_PER_TASK:-40}"

gene_names=(Olfactory Vomeronasal1 Vomeronasal2 Nematode STE2 STE3 Mth GP143 classF classC classT classB1 classB2 GPR1 FSLB classA cAMP GP157)
task_id="${SLURM_ARRAY_TASK_ID:-0}"
if (( task_id >= ${#gene_names[@]} )); then
    echo "Array index ${task_id} is outside the family list" >&2
    exit 2
fi
gene_name="${gene_names[$task_id]}"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
family_dir="${GPCR_FAMILY_DIR:-${PROJECT_DIR}/families}"
base_dir="${family_dir}/${gene_name}"
querydb="${base_dir}/${gene_name}_genes_db"
filtered_db="${base_dir}/${gene_name}_genes_filtered_db"

mkdir -p "${base_dir}"

accession_file="${GPCR_ACCESSION_TXT_DIR:-${family_dir}/accessions}/${gene_name}_accessions.txt"
test -s "${accession_file}" || {
    echo "Missing accession list: ${accession_file}" >&2
    exit 3
}
test "${UNIPROT_MMSEQS_DB}" != "/path/to/UniprotRefEuk" || {
    echo "Set UNIPROT_MMSEQS_DB before running family processing" >&2
    exit 3
}

# Retrieve the filtered family accession set.
mmseqs createsubdb --id-mode 1 \
    "${accession_file}" \
    "${UNIPROT_MMSEQS_DB}" "${querydb}"
mmseqs convert2fasta "${querydb}" "${querydb}.fasta"

# Align, trim gap-rich columns, and remove sequences shorter than 200 residues.
export OMP_NUM_THREADS="${THREADS}"
famsa "${querydb}.fasta" "${querydb}_FAMSA.fasta" -t "${THREADS}"
trimal -in "${querydb}_FAMSA.fasta" \
    -out "${querydb}_FAMSA_gt0.1.fasta" -gt 0.1
python3 "${script_dir}/fasta_filter_min_length.py" \
    --input "${querydb}_FAMSA_gt0.1.fasta" --min-ungapped-length 200

# Build and align the length-filtered family database.
mmseqs createsubdb --id-mode 1 \
    "${querydb}_FAMSA_gt0.1_filtered_ids.txt" \
    "${UNIPROT_MMSEQS_DB}" "${filtered_db}"
mmseqs convert2fasta "${filtered_db}" "${filtered_db}.fasta"
famsa "${filtered_db}.fasta" "${filtered_db}_FAMSA.fasta" -t "${THREADS}"

# Cluster at 70% identity and export cluster representatives.
mkdir -p "${base_dir}/tmp_cluster"
mmseqs cluster "${filtered_db}" "${filtered_db}_clust0.7_clust0_cov1" \
    "${base_dir}/tmp_cluster" --min-seq-id 0.7 -s 7.5 \
    --threads "${THREADS}" --cluster-reassign --remove-tmp-files 1 \
    --cluster-mode 0 --cov-mode 1
mmseqs createtsv "${filtered_db}" "${filtered_db}" \
    "${filtered_db}_clust0.7_clust0_cov1" \
    "${filtered_db}_clust0.7_clust0_cov1.tsv"

cut -f1 "${filtered_db}_clust0.7_clust0_cov1.tsv" | sort -u \
    > "${filtered_db}_clust0.7_clust0_cov1_reps.tsv"
mmseqs createsubdb --id-mode 1 \
    "${filtered_db}_clust0.7_clust0_cov1_reps.tsv" "${filtered_db}" \
    "${base_dir}/${gene_name}_genes_filtered_reps_db"
mmseqs convert2fasta "${base_dir}/${gene_name}_genes_filtered_reps_db" \
    "${base_dir}/${gene_name}_genes_filtered_reps_db.fasta"
sed -i 's/\*//g' "${base_dir}/${gene_name}_genes_filtered_reps_db.fasta"

echo "Completed family processing for ${gene_name}"
