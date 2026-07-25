# -*- coding: utf-8 -*-
"""
Created on Fri Sep  5 15:02:38 2025

@author: selcuk.1
"""

import os

from Bio import SeqIO


PROJECT_DIR = os.environ.get("GPCR_ANALYSIS_DIR", os.environ.get("GPCR_ZENODO_DIR", "/path/to/GPCR_Superfamily"))
SEARCH_RESULTS_DIR = os.environ.get(
    "GPCR_SEARCH_RESULTS_DIR", os.path.join(PROJECT_DIR, "search_results")
)
EXTRACTED_CLADE_DIR = os.environ.get(
    "GPCR_EXTRACTED_CLADE_DIR", os.path.join(PROJECT_DIR, "extracted_clades")
)


def filter_fasta_by_ids(fasta_file, id_list, output_file):
    """
    Read a FASTA file, filter sequences based on provided IDs, and write the filtered sequences to an output file.

    Parameters:
    fasta_file (str): Path to the input FASTA file.
    id_list (list): List of IDs to filter.
    output_file (str): Path to the output FASTA file.
    """
    # Read the FASTA file and store sequences in a dictionary
    sequence_dict = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}

    # Write the filtered sequences to the output FASTA file
    with open(output_file, 'w') as fasta_output:
        for seq_id in id_list:
            if seq_id in sequence_dict:
                fasta_output.write(f">{seq_id}\n{sequence_dict[seq_id]}\n")
            else:
                print(f"ID {seq_id} not found in the FASTA file.")

from ete3 import Tree

def find_leaf_nodes_within_boundaries(boundaries, newick_file):
    """
    boundaries: list of items, each either
      [left_id, right_id]            -> include leaves in this MRCA clade
      [left_id, right_id, -1]        -> EXCLUDE leaves in this MRCA clade
    """
    # Read the phylogenetic tree
    tree = Tree(newick_file, format=1)
    # print(tree)

    included_order = []   # maintain order of first inclusion encounter
    included_seen = set() # for de-duping while preserving order
    excluded = set()      # collect all leaves to subtract at the end

    for b in boundaries:
        if len(b) < 2:
            continue  # ignore malformed entries silently

        left_id, right_id = b[0], b[1]
        node = tree.get_common_ancestor(left_id, right_id)
        leaves = node.get_leaf_names()

        # If a third entry exists and equals -1 -> exclusion
        if len(b) == 3 and b[2] == -1:
            excluded.update(leaves)
        else:
            for name in leaves:
                if name not in included_seen:
                    included_order.append(name)
                    included_seen.add(name)

    # subtract excluded while preserving the original inclusion order
    final_list = [name for name in included_order if name not in excluded]

    print("Total number of identified orthologs:", len(final_list))
    return final_list

    
def prune_tree_by_ids(tree_file, keep_ids, output_tree_file):
    """
    Read a Newick tree, prune away all leaves except those in keep_ids,
    and write the resulting subtree to output_tree_file.
    """
    # Load the tree (assuming format=1 for standard unrooted Newick)
    tree = Tree(tree_file, format=1)

    # Prune the tree in place, keeping only the specified leaf names.
    # Any leaf not in keep_ids will be removed; internal nodes will be collapsed as needed.
    tree.prune(keep_ids, preserve_branch_length=True)

    # Write out the pruned tree
    tree.write(outfile=output_tree_file)

def _resolve_paths(gene_name: str, gene_data: dict):
    hmmregion = gene_data["hmmregion"]
    if hmmregion == 1:
        hmmregion_txt = "hmmregions_"
    elif hmmregion == 0:
        hmmregion_txt = ""
    else:
        hmmregion_txt = "NA"

    if "tree" in gene_data:
        tree_path  = gene_data["tree"]
        fasta_path = gene_data["fasta"]
    else:
        # Zenodo preserves one folder per search, using the exact family name.
        is_m8 = str(hmmregion).upper() == "NA"
        search_name = gene_data.get("search_name", gene_name)
        folder_names = gene_data.get("search_folders", [search_name])
        if is_m8:
            stem = f"{search_name}_db_FAMSA_gt0.1"
            ordered_suffix = "m8ordered"
        else:
            stem = f"{search_name}_hmm_E1000_db_{hmmregion_txt}FAMSA_gt0.1"
            ordered_suffix = "hmmordered"

        direct_candidates = [
            (
                os.path.join(SEARCH_RESULTS_DIR, folder, f"{stem}_{ordered_suffix}.tree"),
                os.path.join(SEARCH_RESULTS_DIR, folder, f"{stem}_treeordered.fasta"),
            )
            for folder in folder_names
        ]
        legacy_candidate = (
            os.path.join(PROJECT_DIR, "hmmorder_trees", f"{stem}_{ordered_suffix}.tree"),
            os.path.join(PROJECT_DIR, "treeorder_msa", f"{stem}_treeordered.fasta"),
        )
        candidates = direct_candidates + [legacy_candidate]
        tree_path, fasta_path = next(
            (
                (candidate_tree, candidate_fasta)
                for candidate_tree, candidate_fasta in candidates
                if os.path.exists(candidate_tree) and os.path.exists(candidate_fasta)
            ),
            direct_candidates[0],
        )

    return hmmregion, hmmregion_txt, tree_path, fasta_path

def ortholog_ids_for(gene_name: str, gene_boundaries: dict):
    """
    Returns (ortholog_ids, hmmregion, hmmregion_txt, tree_path, fasta_path)
    using the same logic you had in the loop.
    """
    gene_data = gene_boundaries[gene_name]
    boundary  = gene_data["boundary"]
    hmmregion, hmmregion_txt, tree_path, fasta_path = _resolve_paths(gene_name, gene_data)
    if not os.path.exists(tree_path) or not os.path.exists(fasta_path):
        raise FileNotFoundError(
            f"Missing ordered inputs for {gene_name}. Expected tree '{tree_path}' "
            f"and FASTA '{fasta_path}'. Set GPCR_SEARCH_RESULTS_DIR to the "
            "Zenodo search_results directory or provide explicit tree/fasta paths "
            "in gene_boundaries."
        )
    ortholog_ids = find_leaf_nodes_within_boundaries(boundary, tree_path)
    return ortholog_ids, hmmregion, hmmregion_txt, tree_path, fasta_path



#%%
gene_boundaries = {
    "STE3": {"boundary": [["tr|A0A9P4MTP4|A0A9P4MTP4_9PLEO", "tr|F8Q0V9|F8Q0V9_SERL3"],
                          ["tr|G4TTQ9|G4TTQ9_SERID","tr|A0A1U7QRA7|A0A1U7QRA7_MESAU",-1]], "hmmregion": 1},
    "STE2": {"boundary": [["tr|B6H2Y5|B6H2Y5_PENRW","tr|A0A0P1BRD9|A0A0P1BRD9_9BASI"]], "hmmregion": 1},
    "Vomeronasal1": {"boundary": [["tr|H9GQ48|H9GQ48_ANOCA", "tr|A0A6P8SKI3|A0A6P8SKI3_GEOSA"]],
                      "hmmregion": 0, "search_name": "classT"},
    "Vomeronasal2": {"boundary": [["tr|A0A8C4SHY2|A0A8C4SHY2_ERPCA", "tr|A0A8J1JL86|A0A8J1JL86_XENTR"]], "hmmregion": 1},
    "Nematode": {"boundary": [["tr|A0A9P1EYY0|A0A9P1EYY0_9PELO","tr|A0A9Q1CMJ0|A0A9Q1CMJ0_HOLLE"]],"hmmregion":0},
    "Mth": {"boundary": [["tr|A0A3Q0J1M1|A0A3Q0J1M1_DIACI","tr|A0A8J2KMS7|A0A8J2KMS7_9HEXA"]],"hmmregion":0},
    "GPR1_YEAST": {"boundary": [["tr|Q6C142|Q6C142_YARLI","tr|L8GQK5|L8GQK5_ACACF"]],"hmmregion":"NA",
                   "search_folders": ["GPR1_YEAST"]},
    "classB": {"boundary": [["tr|A0A401SL16|A0A401SL16_CHIPU","tr|A0A074Y717|A0A074Y717_AURSE"]],"hmmregion":1,
               "output_name": "classB1"},
    "classB2": {"boundary": [
        ["tr|A0A914AT92|A0A914AT92_PATMI", "tr|A0A179FJJ4|A0A179FJJ4_METCM"],
        ["tr|A0A401SL16|A0A401SL16_CHIPU", "tr|A0A821QRT1|A0A821QRT1_9BILA", -1],
    ], "hmmregion": 1, "search_name": "classB"},
    "classC": {"boundary": [["tr|A0A553N6M9|A0A553N6M9_TIGCA","tr|A0A168Q443|A0A168Q443_ABSGL"],
                            ["tr|A0A8C5RW22|A0A8C5RW22_LATLA","tr|A0A9W7WV79|A0A9W7WV79_TRIRA"]],"hmmregion":1},
    "classF": {"boundary": [["tr|A0A9J7KIY8|A0A9J7KIY8_BRAFL","tr|A0A1D1VV00|A0A1D1VV00_RAMVA"]],"hmmregion":0},
    "classT": {"boundary": [["tr|A0A1U7TIF8|A0A1U7TIF8_CARSF","tr|H3A2X9|H3A2X9_LATCH"]],"hmmregion":0},
    "classA": {"boundary": [["tr|A0A6P4ZN87|A0A6P4ZN87_BRABE","tr|A0A8J6GE46|A0A8J6GE46_MICOH"],["tr|H3AS51|H3AS51_LATCH","tr|A0DLB5|A0DLB5_PARTE"],
                            ["tr|A0A8C3XM65|A0A8C3XM65_CHESE","tr|A0A482WL87|A0A482WL87_LAOST",-1]],"hmmregion":1},
    "cAMP": {"boundary": [["sp|P35352|CAR3_DICDI","tr|A0A152A2K3|A0A152A2K3_TIELA"],["tr|A0A9Q0LWC9|A0A9Q0LWC9_ANAIG","tr|L8H1T1|L8H1T1_ACACF"],
                          ["tr|A0AA89BWC4|A0AA89BWC4_PINIB","tr|A0A9Q0LWC9|A0A9Q0LWC9_ANAIG",-1]], "hmmregion":1},
    "GP143": {"boundary": [["tr|A0A8D0GDM2|A0A8D0GDM2_SPHPU","tr|A0A0D2U8X7|A0A0D2U8X7_CAPO3"]],"hmmregion":0},
    "Galpha": {"boundary":[["tr|A0A0M3JM16|A0A0M3JM16_ANISI","tr|A0AAD1XBA8|A0AAD1XBA8_EUPCR"]],"hmmregion":0},
    "Ggamma": {"boundary":[["tr|A0AAD5UE64|A0AAD5UE64_9FUNG","tr|A0A6J8BV00|A0A6J8BV00_MYTCO"],["tr|A0A8H7DX24|A0A8H7DX24_PLEOS","tr|A0A5N5PJW3|A0A5N5PJW3_PANHP"]],"hmmregion":1},
    "Gbeta": {"boundary":[["tr|A0A674NTF0|A0A674NTF0_TAKRU","tr|A0A5E8BDK1|A0A5E8BDK1_9ASCO"]],"hmmregion":0},
    "Olfactory": {"boundary":[["tr|A0A6P8R8A3|A0A6P8R8A3_GEOSA","tr|A0A1S3HJI5|A0A1S3HJI5_LINAN"]],"hmmregion":0},
    "GP157": {"boundary":[["tr|A0A8C3RV18|A0A8C3RV18_CHESE","tr|A0AAD3DB90|A0AAD3DB90_9STRA"]],"hmmregion":1},
    "FSLB_DICDI": {"boundary": [["sp|Q55DT7|FSLD_DICDI","tr|A0A1U7LQE6|A0A1U7LQE6_NEOID"]],"hmmregion":"NA",
                   "search_folders": ["FSLB_DICDI"]},
    "LRLA_DICDI": {"boundary": [["sp|Q54MC6|LRLA_DICDI","tr|A0A152A1B2|A0A152A1B2_TIELA"]],"hmmregion":"NA",
                   "search_folders": ["LRLA_DICDI"]},
    "TBP": {"boundary": [["tr|A0A3P8S939|A0A3P8S939_AMPPE","tr|A0A9Q1RRU0|A0A9Q1RRU0_9SOLA"],["tr|A0A8K0LVE6|A0A8K0LVE6_9PEZI","tr|A0A6A5ZXD8|A0A6A5ZXD8_9PLEO"]],"hmmregion":1},
    "RBP1": {"boundary": [["tr|H2Z9E6|H2Z9E6_CIOSA","tr|A0A6G0W177|A0A6G0W177_APHCR"]],"hmmregion":0},
    
}#%%


def main():
    os.makedirs(os.path.join(EXTRACTED_CLADE_DIR, "alignments"), exist_ok=True)
    os.makedirs(os.path.join(EXTRACTED_CLADE_DIR, "trees"), exist_ok=True)

    for gene_name, gene_data in gene_boundaries.items():
        print("Extracting clade", gene_name)
        ortholog_ids, hmmregion, hmmregion_txt, tree, fasta_file = ortholog_ids_for(gene_name, gene_boundaries)

        is_m8 = str(hmmregion).upper() == "NA"
        output_name = gene_data.get("output_name", gene_name)
        if is_m8:
            out_file = os.path.join(
                os.path.join(EXTRACTED_CLADE_DIR, "alignments"),
                f"{output_name}_db_FAMSA_gt0.1_treeordered_extracted_clade.fasta",
            )
        else:
            out_file = os.path.join(
                EXTRACTED_CLADE_DIR, "alignments",
                f"{output_name}_hmm_E1000_db_FAMSA_gt0.1_treeordered_extracted_clade.fasta",
            )

        filter_fasta_by_ids(fasta_file, ortholog_ids, out_file)
        if is_m8:
            out_tree = os.path.join(
                os.path.join(EXTRACTED_CLADE_DIR, "trees"),
                f"{output_name}_db_FAMSA_gt0.1_m8ordered_extracted_clade.tree",
            )
        else:
            out_tree = os.path.join(
                os.path.join(EXTRACTED_CLADE_DIR, "trees"),
                f"{output_name}_hmm_E1000_db_{hmmregion_txt}FAMSA_gt0.1_hmmordered_extracted_clade.tree",
            )
        prune_tree_by_ids(tree, ortholog_ids, out_tree)


if __name__ == "__main__":
    main()
