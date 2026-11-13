#!/usr/bin/env bash
#SBATCH --job-name=gpcr_family_alignment
#SBATCH --mem=160GB
#SBATCH --cpus-per-task=40
#SBATCH --nodes=1
#SBATCH --time=0-1:0:0
#SBATCH --array=0-18
#SBATCH --output=%x_%A_%a.out

set -euo pipefail

# Family-wide alignment of the 7TM cores of each family's final members (after the Dikarya
# reassignment, final_membership_reassigned.tsv), one family per array task: the members are
# fetched, cut to their 7TM core with five flanking residues, aligned with FAMSA, and the
# alignment's topology projection is trimmed at 90% gaps. The TMonly FAMSA alignments are
# the family alignments read by extract_family_representatives.py, family_column_counts.py,
# trim_family_alignments_consensus.py and map_motif_residues.py.

# USER CONFIGURATION ----------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config/paths.sh"
STAGE_DIR="${GPCR_ZENODO_DIR}/10_motif_analysis"
FAMILIES_DIR="${FAMILIES_DIR:-${STAGE_DIR}/family_alignments/families}"
THREADS="${SLURM_CPUS_PER_TASK:-40}"

families=(
    classA Olfactory Vomeronasal1 classT Nematode GPR1 cAMP classB1 classB2 Mth
    classF FSLB GP143 TM116 GP157 STE3 classC Vomeronasal2 STE2
)
# END USER CONFIGURATION ------------------------------------------------------

family="${families[${SLURM_ARRAY_TASK_ID:-0}]}"
family_dir="${FAMILIES_DIR}/${family}"
mkdir -p "${family_dir}" "${STAGE_DIR}/family_alignments/topology_alignments"

threeline="${FAMILIES_DIR}/deeptmhmm_predictions_oriented.3line"
[[ -s "${threeline}" ]] || gzip -dc "${GPCR_ZENODO_DIR}/03_topology_filter/deeptmhmm_predictions_oriented.3line.gz" > "${threeline}"

# the family's column of the wide membership table
accessions="${family_dir}/${family}_filtered_accessions.txt"
awk -F'\t' -v fam="${family}" 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == fam) c = i; next }
    c && $c != "" { print $c }' \
    "${STAGE_DIR}/phylogeny_reassignment/final_membership_reassigned.tsv" > "${accessions}"

querydb="${family_dir}/${family}_filtered_genes_db"
mmseqs createsubdb --id-mode 1 "${accessions}" "${UNIPROT_MMSEQS_DB}" "${querydb}"
mmseqs convert2fasta "${querydb}" "${querydb}.fasta"
sed -i 's/\*//g' "${querydb}.fasta"

python3 "${SCRIPT_DIR}/../08_superfamily_alignment/extract_tm_region.py" \
    --fasta "${querydb}.fasta" \
    --topology "${threeline}" \
    --out "${querydb}_TMonly.fasta" \
    --boundaries "${querydb}_TMregion_boundaries.tsv" \
    --flank 5

export OMP_NUM_THREADS="${THREADS}"
famsa "${querydb}_TMonly.fasta" "${querydb}_TMonly_FAMSA.fasta" -t "${THREADS}"

python3 "${SCRIPT_DIR}/../03_topology_filter/project_topology_onto_alignment.py" \
    --alignment "${querydb}_TMonly_FAMSA.fasta" \
    --topology "${threeline}" \
    --boundaries "${querydb}_TMregion_boundaries.tsv" \
    --out "${querydb}_FAMSA_topology_converted.fasta"

trimal -in "${querydb}_FAMSA_topology_converted.fasta" \
    -out "${STAGE_DIR}/family_alignments/topology_alignments/${family}_filtered_genes_db_FAMSA_topology_converted_gt0.1.fasta" \
    -gt 0.1
