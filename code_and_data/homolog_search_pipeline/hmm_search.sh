#!/usr/bin/env bash
#SBATCH --job-name=gpcr_hmmsearch
#SBATCH --mem=160GB
#SBATCH --cpus-per-task=40
#SBATCH --nodes=1
#SBATCH --time=05:00:00
#SBATCH --account=PAS1794
#SBATCH --array=0-19
#SBATCH --output=%x_%A_%a.out

set -euo pipefail

# USER CONFIGURATION ----------------------------------------------------------
PROJECT_DIR="${GPCR_ANALYSIS_DIR:-/path/to/GPCR_Superfamily}"
SEARCH_RESULTS_DIR="${GPCR_SEARCH_RESULTS_DIR:-${PROJECT_DIR}/search_results}"
UNIPROT_FASTA="${UNIPROT_FASTA:-/path/to/Uniprot_Ref_Eukaryota_Proteome.fasta}"
UNIPROT_MMSEQS_DB="${UNIPROT_MMSEQS_DB:-/path/to/UniprotRefEuk}"
THREADS="${SLURM_CPUS_PER_TASK:-40}"

# HMM-profile searches. Seed FASTA searches containing DICDI/YEAST are handled
# separately by mmseqs2_search.sh.
families=(
    STE3 STE2 Vomeronasal1 Vomeronasal2 Nematode Mth
    classB classC classF classT classA cAMP GP143
    Galpha Ggamma Gbeta Olfactory GP157 TBP RBP1
)

# Maximum number of top-ranked hits retained for each family.
declare -A max_hits=(
    [STE3]=10000
    [STE2]=10000
    [Vomeronasal1]=10000
    [Vomeronasal2]=10000
    [Nematode]=10000
    [Mth]=10000
    [classB]=150000
    [classC]=50000
    [classF]=50000
    [classT]=20000
    [classA]=500000
    [cAMP]=50000
    [GP143]=10000
    [Galpha]=500000
    [Ggamma]=30000
    [Gbeta]=30000
    [Olfactory]=200000
    [GP157]=10000
    [TBP]=50000
    [RBP1]=50000
)
# END USER CONFIGURATION ------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
task_id="${SLURM_ARRAY_TASK_ID:-0}"
if (( task_id >= ${#families[@]} )); then
    echo "Array index ${task_id} is outside the families list" >&2
    exit 2
fi

family="${families[$task_id]}"
if [[ ! -v "max_hits[${family}]" ]]; then
    echo "No retrieval limit configured for ${family}" >&2
    exit 3
fi
family_dir="${SEARCH_RESULTS_DIR}/${family}"
querydb="${family_dir}/${family}_hmm_E1000_db"
hmm_output="${family_dir}/${family}_hmmsearch_E1000.txt"

mkdir -p "${family_dir}"
shopt -s nullglob
hmm_models=("${SCRIPT_DIR}/input_sequences/${family}"/*.hmm)
if (( ${#hmm_models[@]} != 1 )); then
    echo "Expected exactly one HMM for ${family}; found ${#hmm_models[@]}" >&2
    exit 4
fi

hmmsearch --textw 50000 --noali --max -E 1000 \
    --cpu "${THREADS}" -o "${hmm_output}" "${hmm_models[0]}" "${UNIPROT_FASTA}"

python3 "${SCRIPT_DIR}/hmm_header_get.py" \
    "${hmm_output}" "${max_hits[$family]}"
mmseqs createsubdb --id-mode 1 \
    "${family_dir}/${family}_hmmsearch_E1000_headers.txt" \
    "${UNIPROT_MMSEQS_DB}" "${querydb}"

if [[ "${family}" == "classA" ]]; then
    # The Class A search is too large to align directly. Cluster the search
    # hits and use representatives for curation; the TSV is later used to
    # expand selected representatives back to the complete protein set.
    classa_cluster="${querydb}_clust0.8"
    classa_reps="${classa_cluster}_reps"

    mmseqs cluster "${querydb}" "${classa_cluster}" \
        "${family_dir}/tmp_classA_cluster" --min-seq-id 0.8 -s 7.5 \
        --threads "${THREADS}" --cluster-reassign --remove-tmp-files 1
    mmseqs createtsv "${querydb}" "${querydb}" \
        "${classa_cluster}" "${classa_cluster}.tsv"
    mmseqs createsubdb "${classa_cluster}" "${querydb}" "${classa_reps}"
    mmseqs convert2fasta "${classa_reps}" "${querydb}.fasta"
else
    mmseqs convert2fasta "${querydb}" "${querydb}.fasta"
fi
sed -i 's/\*//g' "${querydb}.fasta"

python3 "${SCRIPT_DIR}/trim_sequences_using_hmm.py" \
    "${hmm_output}" "${querydb}.fasta"

export OMP_NUM_THREADS="${THREADS}"
for variant in "${querydb}_hmmregions" "${querydb}"; do
    famsa "${variant}.fasta" "${variant}_FAMSA.fasta" -t "${THREADS}"
    trimal -in "${variant}_FAMSA.fasta" \
        -out "${variant}_FAMSA_gt0.1.fasta" -gt 0.1
    FastTreeMP "${variant}_FAMSA_gt0.1.fasta" > "${variant}_FAMSA_gt0.1.tree"
done
