#!/usr/bin/env bash
#SBATCH --job-name=gpcr_order_scores
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=3
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --array=0-21
#SBATCH --output=%x_%A_%a.out

set -euo pipefail

# USER CONFIGURATION ----------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config/paths.sh"
SEARCH_RESULTS_DIR="${GPCR_ZENODO_DIR}/01_homology_search"

# One search type (hmm or m8) is required for every entry. DICDI/YEAST seed
# searches use MMseqs2; all other entries use HMMER.
families=(
    STE3 STE2 Vomeronasal2 Nematode Mth
    GPR1_YEAST classB classC classF classT classA cAMP GP143 TM116
    Galpha Ggamma Gbeta Olfactory GP157 FSLB_DICDI TBP RBP1
)
search_types=(
    hmm hmm hmm hmm hmm
    m8 hmm hmm hmm hmm hmm hmm hmm hmm
    hmm hmm hmm hmm hmm m8 hmm hmm
)
# END USER CONFIGURATION ------------------------------------------------------

task_id="${SLURM_ARRAY_TASK_ID:-0}"
if (( ${#families[@]} != ${#search_types[@]} )); then
    echo "families and search_types must have equal lengths" >&2
    exit 2
fi
if (( task_id >= ${#families[@]} )); then
    echo "Array index ${task_id} is outside the families list" >&2
    exit 2
fi

family="${families[$task_id]}"
search_type="${search_types[$task_id]}"
base_dir="${SEARCH_RESULTS_DIR}/${family}"

order_variant() {
    local score_flag="$1" score_file="$2" variant="$3" suffix="$4"
    local ordered_tree="${variant}_FAMSA_gt0.1_${suffix}ordered.tree"

    python3 "${SCRIPT_DIR}/tree_order_by_hmm.py" \
        "${score_flag}" "${score_file}" --tree "${variant}_FAMSA_gt0.1.tree"
    python3 "${SCRIPT_DIR}/fasta_order_by_tree.py" \
        --tree "${ordered_tree}" --fasta "${variant}_FAMSA_gt0.1.fasta"
    python3 "${SCRIPT_DIR}/line_plot_hmm_and_index_treeorder.py" \
        "${score_flag}" "${score_file}" --tree "${ordered_tree}" \
        --weighting-type none --window-size 50 \
        --html-out "${variant}_FAMSA_gt0.1_${suffix}ordered_score_lineplot.html"
}

case "${search_type}" in
    hmm)
        querydb="${base_dir}/${family}_hmm_E1000_db"
        score_file="${base_dir}/${family}_hmmsearch_E1000.txt"
        order_variant --hmm "${score_file}" "${querydb}_hmmregions" hmm
        order_variant --hmm "${score_file}" "${querydb}" hmm
        ;;
    m8)
        order_variant --m8 \
            "${base_dir}/${family}_UniProt_s7.5_filter0_eprofile0.1_db.m8" \
            "${base_dir}/${family}_db" m8
        ;;
    *)
        echo "Unknown search type: ${search_type}" >&2
        exit 3
        ;;
esac
