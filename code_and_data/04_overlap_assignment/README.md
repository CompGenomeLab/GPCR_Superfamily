# 04 · Assigning sequences claimed by two families, and per-genome counts

Methods: *Controls and Filtering* (last two paragraphs); Table S2.

After the topology filter, 382 sequences are still in the clades of two families.
`analyze_family_overlaps.py` lists them and scores each in every family claiming it by its
median-normalized search score: its bit score divided by the median bit score of that family's members,
so scores from different searches can be compared (`overlap_diagnostics.tsv`). `build_overlap_rep_set.py`
adds the ten highest-scoring members of each family involved in an overlap, the set is aligned with
FAMSA, columns with 90% or more gaps are removed and a tree is inferred with FastTree (`assignment_tree.sh`).
`annotate_tree_with_families.py` labels the leaves with their families and `build_picker.py` builds the
page the assignments were made in (`serve_picker.py` serves it and records every decision in
`mrca_selections.tsv`). Each contested sequence was assigned to the family whose top-scoring members it
clusters with, either by a clade of the tree (`sweep`) or individually (`manual`, using the
median-normalized scores where the tree gives no clear signal). `write_final_members.py` applies
`mrca_selections.tsv` to the topology-filtered membership.

`build_count_tables.py` then maps every final member to its genome (taxID, the `OX=` field of its
UniProt header) and counts members per genome and family; the controls are counted from their clade
membership, since the topology filter does not apply to them.

```bash
source ../../config/paths.sh
python analyze_family_overlaps.py                       # overlap_diagnostics.tsv
python build_overlap_rep_set.py --diagnostics $GPCR_ZENODO_DIR/04_overlap_assignment/overlap_diagnostics.tsv \
    --after-tsv $GPCR_ZENODO_DIR/03_topology_filter/topology_membership.tsv \
    --search-dir $GPCR_ZENODO_DIR/01_homology_search \
    --output $GPCR_ZENODO_DIR/04_overlap_assignment/assignment_tree/overlaps_and_reps_accession.txt
sbatch assignment_tree.sh                               # FAMSA, trimAl -gt 0.1, FastTree -> assignment_tree/
python annotate_tree_with_families.py ...               # see the script's usage
python build_picker.py ... && python serve_picker.py    # the decisions are in mrca_selections.tsv
python write_final_members.py --after-tsv $GPCR_ZENODO_DIR/03_topology_filter/topology_membership.tsv \
    --selections mrca_selections.tsv --output-dir $GPCR_ZENODO_DIR/04_overlap_assignment/final_families \
    --merged-tsv $GPCR_ZENODO_DIR/04_overlap_assignment/final_membership.tsv
python build_count_tables.py
```

`mrca_selections.tsv` is the curated record of the 382 decisions (and of the 150 family representatives
of the tree). Verified: applied to the topology-filtered membership it gives the final membership of the
paper exactly (391,856 sequences), and `build_count_tables.py` reproduces the per-genome counts of every
family and control without a single differing cell.

**Archive** (`04_overlap_assignment/`):

| File | Contents |
|---|---|
| `overlap_diagnostics.tsv` | the 382 contested sequences with their median-normalized scores |
| `assignment_tree/` | the sequences, FAMSA alignment, trimmed alignment and FastTree tree of the assignment set, and the family-labelled tree |
| `final_membership.tsv`, `final_families/` | final membership, one column (one file) per family |
| `acc_taxid_map.tsv` | every final member with its taxID |
| `tax_counts_per_family.tsv` | members per genome and family, controls included |
