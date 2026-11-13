#!/usr/bin/env bash
#SBATCH --job-name=gpcr_superfamily_tree
#SBATCH --mem=160GB
#SBATCH --cpus-per-task=40
#SBATCH --nodes=1
#SBATCH --time=2-0:0:0
#SBATCH --array=0-7
#SBATCH --output=%x_%A_%a.out

set -euo pipefail

# One task per trimmed alignment (4 sets x 2 alignments; task i -> set i/2, alignment i%2):
# ModelFinder (IQ-TREE 2.3.6) selects the model, RAxML-NG 1.2.2 infers the tree from 100
# parsimony starting trees with autoMRE bootstrapping (up to 500 replicates) and TBE
# support, and booster gives the transfer index of every leaf in the best tree.

# USER CONFIGURATION ----------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config/paths.sh"
iqtree2="${IQTREE2:-iqtree2}"
raxmlng="${RAXMLNG:-raxml-ng}"
booster="${BOOSTER:-booster}"
STAGE_DIR="${GPCR_ZENODO_DIR}/09_superfamily_phylogeny"
STATS_DIR="${STAGE_DIR}/booster_stats"
CUTOFF=0.3                # 0.3 = branches with TBE >= 0.7

sets=(representatives_all representatives_no_ste2 representatives_no_ste2ste3classCV2 representatives_and_literature)
alignments=(ginsi_ep0.123 ginsi)
# END USER CONFIGURATION ------------------------------------------------------

task="${SLURM_ARRAY_TASK_ID:-0}"
set_name="${sets[$((task / 2))]}"
aln="${alignments[$((task % 2))]}"
echo "== task ${task}: ${set_name}, ${aln}"

trimmed_msa="${GPCR_ZENODO_DIR}/08_superfamily_alignment/alignments/${set_name}_${aln}_kpicgappy.fasta"
pre="${STAGE_DIR}/runs/${set_name}/${set_name}_${aln}"
mkdir -p "${STAGE_DIR}/runs/${set_name}" "${STATS_DIR}"

$iqtree2 -s "${trimmed_msa}" -m MF -mfreq FO -cmin 5 -cmax 15 -T 10 --mset JTT --mrate R -pre "${pre}_kpicgappy_bestmodel"
bmodel=$(grep "Best-fit model:" "${pre}_kpicgappy_bestmodel.log" | awk '{print $3}')

$raxmlng --all --tree pars{100} --msa "${trimmed_msa}" --model "${bmodel}" --bs-trees autoMRE{500} \
    --bs-metric TBE --threads 40 --workers AUTO --prefix "${pre}_kpicgappy"

# booster: per-taxon transfer index for this tree
PREFIX="${pre}_kpicgappy"
OUT="${PREFIX}.booster"
echo "Bootstraps: $(wc -l < "${PREFIX}.raxml.bootstraps") trees"
$booster -a tbe \
    -i "${PREFIX}.raxml.bestTree" \
    -b "${PREFIX}.raxml.bootstraps" \
    -@ "${SLURM_CPUS_PER_TASK:-1}" \
    -d "${CUTOFF}" \
    -o "${OUT}.tbe.nwk" \
    -r "${OUT}.tbe_raw.nwk" \
    -S "${OUT}.tbe_stats.txt" \
    -c

cp "${OUT}.tbe_stats.txt" "${STATS_DIR}/$(basename "${PREFIX}").tbe_stats.txt"
