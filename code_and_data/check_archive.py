"""Check that the Zenodo archive is in place and holds the files each stage reads.

Usage:
    python code_and_data/check_archive.py
"""
import sys

from gpcr_paths import DATA, FAMILIES, exists, stage

SEARCHES = ["classA", "Olfactory", "classT", "Nematode", "cAMP", "classB", "Mth", "classF", "GP143",
            "TM116", "GP157", "STE3", "classC", "Vomeronasal2", "STE2", "Galpha", "Gbeta", "Ggamma",
            "TBP", "RBP1", "FSLB_DICDI", "GPR1_YEAST"]
REQUIRED = {
    "01_homology_search": [f"{s}" for s in SEARCHES] + ["classA/classA_hmm_E1000_db_clust0.8.tsv"],
    "02_clade_extraction": ["extracted_clades/alignments", "extracted_clades/trees", "clade_membership.tsv",
                            "clade_membership_controls.tsv"] + [f"clade_sequences/{f}.fasta.gz" for f in FAMILIES],
    "03_topology_filter": ["deeptmhmm_predictions_oriented.3line", "core_membership.tsv", "topology_alignments",
                           "topology_helix_scores.tsv", "topology_membership.tsv"],
    "04_overlap_assignment": ["overlap_diagnostics.tsv", "assignment_tree", "final_membership.tsv",
                              "acc_taxid_map.tsv", "tax_counts_per_family.tsv"],
    "05_literature_comparison": ["literature_sequences.tsv", "literature_comparison.tsv",
                                 "unmatched_literature_sequences.m8", "remaining_literature_sequences.m8",
                                 "unsearched_literature_sequences.m8", "literature_all_outcomes.tsv",
                                 "Supplementary_Table_literature_receptors.tsv", "deeper_origin_candidates.tsv",
                                 "deeper_origin_7TM_gate.tsv", "literature_reps_7TM_flank5.fasta",
                                 "deeper_origin_tested_sequences.fasta", "predicted_topologies.3line"],
    "06_family_origins": ["taxonomy/taxonomy_ncbi_searched.json", "taxonomy/ncbi_taxdump/names.dmp",
                          "taxonomy/ncbi_taxdump/nodes.dmp", "taxonomy_ncbi_all.nwk",
                          "family_presence_absence_ncbi_all.csv", "dollo_origins.tsv",
                          "lineage_origin_figure.html"],
    "07_representative_selection": ["representatives_all_families.tsv", "representatives_all.txt",
                                    "representatives_no_ste2.txt", "representatives_no_ste2ste3classCV2.txt",
                                    "taxonomy_eukaryotes_filtered.json"],
    "08_superfamily_alignment": ["fasta", "alignments/representatives_all_ginsi_ep0.123_kpicgappy.fasta"],
    "09_superfamily_phylogeny": ["runs", "trees", "booster_stats/tindex_mean_across_final_trees.tsv", "sCF",
                                 "representatives_family_map.tsv", "tree_support_tables/S_tree_runs.tsv"],
    "10_motif_analysis": ["treein_alignments", "guide_trees", "phylogeny_reassignment/final_membership_reassigned.tsv",
                          "phylogeny_reassignment/inputs_Sep22/best_trees", "phylogeny_reassignment/inputs_Sep22/ml_trees",
                          "phylogeny_reassignment/inputs_Sep22/representatives_all_Sep22_topo0.8_ginsi.fasta",
                          "family_alignments/family_reps_TMonly_FAMSA", "family_alignments/consensus_trimmed_msas",
                          "anchors/bridge_representatives_ranksum.tsv", "website_inputs", "consensus_mapping"],
    "11_structure_superposition": ["structures", "af_models", "motif_residue_numbers.tsv", "member_sequences.tsv",
                                   "superposition_core_fit_gp157.tsv", "local_backbone_check_gp157.tsv"],
}


def main():
    if not DATA.is_dir():
        sys.exit(f"No archive at {DATA}; unpack it beside the repository or set GPCR_ZENODO_DIR.")
    missing = [f"{name}/{item}" for name, items in REQUIRED.items() for item in items
               if not exists(stage(name, item))]
    print(f"archive: {DATA}")
    if missing:
        print(f"{len(missing)} missing:")
        for m in missing:
            print("  " + m)
        sys.exit(1)
    print(f"all {sum(len(v) for v in REQUIRED.values())} checked files and folders present")


if __name__ == "__main__":
    main()
