#!/bin/bash

#SBATCH --job-name=scf_summary
#SBATCH --mem=40GB
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --output=%j_scf_summary.out
#SBATCH --time=0-0:30:0


# ===========================================================================
# run_scf_gather.sh
#
# Collect sCF across the 100 ML trees of each of the 8 datasets (4 rep sets x 2
# alignments), compare the rep sets within each alignment, and plot. Run AFTER the sCF_calculate.sh array has finished.
#
#     sbatch run_scf_gather.sh      (or: bash run_scf_gather.sh)
#
# Submit it FROM this folder: gather_scf.py, compare_scf_trees.py and
# discover_groupings.py must sit next to it.
# ===========================================================================
set -euo pipefail

# --------------------------- EDIT THESE IF NEEDED --------------------------
source "$(pwd)/../../config/paths.sh"
BASE="$GPCR_ZENODO_DIR/09_superfamily_phylogeny/sCF_analysis"
RESULTS_ROOT="$GPCR_ZENODO_DIR/09_superfamily_phylogeny/sCF"

sets=(representatives_all representatives_no_ste2 representatives_no_ste2ste3classCV2 representatives_and_literature)
labels=("all" "no STE2" "core (no STE2/STE3/C/V2)" "core + literature")
# literature sequences are kept out of split names like the divergent families; the
# divergent_placement column of the _long.tsv shows where they sit
DIVERGENT=(STE2 STE3 classC LIT_GP157 LIT_classT LIT_Vomeronasal1 LIT_GPR1like)
# ---------------------------------------------------------------------------

PY=$(command -v python3 || command -v python)
if [[ -z "$PY" ]]; then echo "no python found on PATH" >&2; exit 1; fi
echo "python : $PY"
echo "base   : $BASE"
echo

# One pass per alignment, each comparing the 4 rep sets (the plots take at most 5 datasets).
for aln in ginsi_ep0.123 ginsi; do
    DIRS=()
    for s in "${sets[@]}"; do DIRS+=("$BASE/scf_${s}_${aln}"); done
    RESULTS="$RESULTS_ROOT/$aln"
    mkdir -p "$RESULTS"
    echo "===================== $aln ====================="

    # --------------------- sanity: did every array task finish? ----------------
    for k in "${!DIRS[@]}"; do
        d=${DIRS[$k]}
        if [[ ! -d "$d" ]]; then echo "MISSING directory: $d" >&2; continue; fi
        n=$(ls "$d"/*.cf.tree 2>/dev/null | wc -l)
        echo "$(basename "$d") : $n / 100 cf.tree files"
        if [[ "$n" -ne 100 ]]; then
            echo "  !! expected 100 - rerun the missing trees of $(basename "$d") (sCF_calculate.sh array)" >&2
        fi
    done
    echo

    # --------------------------- 1. per-tree summaries -------------------------
    "$PY" gather_scf.py "${DIRS[@]}" -o "$RESULTS/scf_per_tree_summary.tsv"
    echo

    # ------------------- 2. comparison, p-values, figure -----------------------
    "$PY" compare_scf_trees.py "$RESULTS/scf_per_tree_summary.tsv" \
          --out-stats "$RESULTS/scf_tree_stats.tsv" \
          --out-pairs "$RESULTS/scf_tree_pairwise.tsv" \
          --out-fig   "$RESULTS/scf_tree_comparison" \
          --labels "${labels[@]}"
    echo

    # ------- 3. every multi-family split, its recovery and sCF --------
    #  (a) divergent families (STE2, STE3, classC incl. Vomeronasal2) left out of split names
    "$PY" discover_groupings.py "${DIRS[@]}" -o "$RESULTS/grouping_discovery_divergent_excluded" \
          --labels "${labels[@]}" --mode excluded --divergent "${DIVERGENT[@]}" \
          --min-families 2 --min-trees 10 --threshold 33 --plot-min-share 0.5
    echo

    #  (b) every family counted as it is, divergent ones included
    "$PY" discover_groupings.py "${DIRS[@]}" -o "$RESULTS/grouping_discovery_divergent_included" \
          --labels "${labels[@]}" --mode included --divergent "${DIVERGENT[@]}" \
          --min-families 2 --min-trees 10 --threshold 33 --plot-min-share 0.5
    echo
done

echo "done. results in $RESULTS_ROOT/"
ls -1R "$RESULTS_ROOT"
