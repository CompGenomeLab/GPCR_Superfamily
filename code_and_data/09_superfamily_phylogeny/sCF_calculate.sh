#!/bin/bash
#
#
#SBATCH --job-name=scf_array
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --output=%A_%a_scf.out
#SBATCH --time=0-1:0:0
#SBATCH --array=0-799

# One task = one ML search tree of one dataset (8 datasets x 100 trees).
#   dataset d = task / 100 (0-7), same order as the tree-making array:
#   set = d / 2, alignment = d % 2
#   tree i = task % 100 + 1


file_names=(representatives_all representatives_no_ste2 representatives_no_ste2ste3classCV2 representatives_and_literature)
alignments=(ginsi_ep0.123 ginsi)

d=$((SLURM_ARRAY_TASK_ID / 100))
suffix=${file_names[$((d / 2))]}
aln=${alignments[$((d % 2))]}
i=$((SLURM_ARRAY_TASK_ID % 100 + 1))

# Submit from this folder (sbatch sCF_calculate.sh); IQ-TREE 2.3.6 is the version used.
source "$(pwd)/../../config/paths.sh"
iqtree2="${IQTREE2:-iqtree2}"
scfdir="$GPCR_ZENODO_DIR/09_superfamily_phylogeny/sCF_analysis"

pre="$GPCR_ZENODO_DIR/09_superfamily_phylogeny/runs/${suffix}/${suffix}_${aln}"
trimmed_msa="$GPCR_ZENODO_DIR/08_superfamily_alignment/alignments/${suffix}_${aln}_kpicgappy.fasta"
ml_trees=$pre"_kpicgappy.raxml.mlTrees"
# run still in its bootstrap phase -> only the temporary file exists
if [[ ! -s "$ml_trees" && -s "$ml_trees.TMP" ]]; then
    ml_trees="$ml_trees.TMP"
fi
bmodel=$(grep "Best-fit model:" $pre"_kpicgappy_bestmodel.log" | awk '{print $3}')

if [[ -z "$bmodel" ]]; then
    echo "no Best-fit model in ${pre}_kpicgappy_bestmodel.log" >&2; exit 1
fi

treedir="$scfdir/trees_${suffix}_${aln}"
outdir="$scfdir/scf_${suffix}_${aln}"
mkdir -p "$treedir" "$outdir"

# tree i = line i of the .raxml.mlTrees file (ML search #i)
tree="$treedir/tree_${i}.treefile"
sed -n "${i}p" "$ml_trees" > "$tree"
if [[ ! -s "$tree" ]]; then
    echo "no tree $i in $ml_trees" >&2; exit 1
fi

pre_out="$outdir/scf_${i}"

echo "task $SLURM_ARRAY_TASK_ID : $suffix $aln tree $i from $(basename "$ml_trees"), model $bmodel"

$iqtree2 -te "$tree" \
         -s "$trimmed_msa" \
         -m "$bmodel" \
         --scfl 100 \
         --prefix "$pre_out" \
         -T "$SLURM_CPUS_PER_TASK" \
         -redo
