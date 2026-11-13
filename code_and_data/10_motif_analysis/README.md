# 10 · Tree-guided superfamily alignment and shared motifs

Methods: *Tree-Guided Superfamily Alignment* (connecting family-wide and superfamily-wide MSAs, logo
generation); Figure 4; Tables S8 and S9; Fig. S2.

**Tree-guided alignment.** `make_guide_tree.py` turns the best ML tree of all 590 representatives into a
guide tree with STE2 joined at the root and the rest rooted between the Glutamate-related families and
the others; `prepare_mafft_treein.py` writes the MAFFT merge steps (`run_guide_trees.sh`), and
`run_treein_ginsi.sh` realigns the 7TM cores with G-INS-i on that tree (`--ep 0.123`, and with the
default penalty for the second alignment).

**Dikarya Class A to GPR1** (`phylogeny_reassignment/`). `tree_side_per_sequence.py` counts, for every
non-animal Class A representative, the trees and ML searches in which it falls on the GPR1 side of the
GPR1 / Class A-like split; `map_fungal_classA_to_clades.py` assigns every fungal Class A sequence to the
clade of its best local-alignment match; `apply_phylogeny_reassignments.py` moves the sequences listed in
`phylogeny_reassignments.tsv` (136 Dikarya sequences) to GPR1 and writes the membership the family-wide
alignments are built from; `fungal_motif_states.py` checks their conserved residues.

**Family columns.** Each family's final members are aligned (7TM cores, FAMSA; `family_alignment.sh`) and the rows of
its representatives are taken out (`extract_family_representatives.py`). `rank_bridge_representatives.py`
ranks the representatives by the sum of their topology-score ranks within the family and across the
superfamily alignments, and `write_logo_anchors.py` writes the anchors and the GPCRdb numbering references
(`gpcrdb_references.tsv`). `consensus_column_mapping.py` with `consensus_mapping_top_reps.py` maps every
superfamily column to the family column the top five representatives vote for, when at least a third
agree, one-to-one. `trim_family_alignments_consensus.py` keeps the mapped family columns (cluster),
`consensus_logo_mappings.py` builds the logo data of every family (and of the animal, fungal and plant
members of cAMP and Class C), and `label_by_human_receptors.py` with `apply_voted_labels.py` labels
every column with the GPCRdb number voted by the family's human receptors. `run_treein_logos.sh` runs
the mapping, logo and label steps.

**Conservation.** `superfamily_motif_positions.py` holds the rule: a family conserves a position when its
most frequent residue plus residues of the same similarity group (D/E, F/W/Y/H, K/R, N/Q, I/V, L/M)
reach 60% of its sequences, gaps counted; positions conserved in at least 12 of the 18 GPCR families
(STE2 excluded) are selected (28 positions). `make_fig4.py` composes Figure 4 (scatter, structure view
from `render_fig4b_structure.py`, and the logo exported from gpcrevolution.org/superfamily-logo) and names
the motifs over their logo columns (MOTIFS, each defined in its class's GPCRdb numbering as in the cited
papers; corrected 21 Oct: HETx at E3x50 only, NPGQ at Q7x49 only, the class F extended switch at W3x43,
F6x36 and Y6x40).
`ste_consensus_similarity.py` counts the conserved residues STE2 and STE3 share with each family
(Table S9).

```bash
sbatch family_alignment.sh              # family-wide alignments -> family_alignments/families/<family>/
#   then, where those alignments are (run from $GPCR_ZENODO_DIR/10_motif_analysis/family_alignments/):
#   export FAMILIES_DIR=$GPCR_ZENODO_DIR/10_motif_analysis/family_alignments/families
#   python extract_family_representatives.py $GPCR_ZENODO_DIR/08_superfamily_alignment/fasta/representatives_all/representatives_all.fasta $FAMILIES_DIR family_reps_TMonly_FAMSA
#   python family_column_counts.py $FAMILIES_DIR <out_dir> --groups family_column_counts_kingdom_groups.tsv
#   python trim_family_alignments_consensus.py --columns consensus_trim_columns_allreps.tsv --out-dir consensus_trimmed_msas
bash run_guide_trees.sh                 # then on the cluster: sbatch run_treein_ginsi.sh
bash run_treein_logos.sh
python superfamily_motif_positions.py
python ste_consensus_similarity.py
python ste_consensus_similarity.py --mappings superfamily_logo_mappings_lower_gap_treein --out ste_consensus_similarity_lower_gap
python render_fig4b_structure.py   # needs PyMOL; figure4/fig3b_motif_positions_2rh1.png
python make_fig4.py                # figure4/fig4_motifs.svg
```

Verified: rerun from the archive, all 25 logo data files (19 families, the six kingdom rows) are
identical to those of gpcrevolution.org, and the STE2/STE3 tables are identical to the published ones.

**Archive** (`10_motif_analysis/`): `guide_trees/`, `treein_alignments/`, `phylogeny_reassignment/`
(decisions, per-sequence readouts, the reassigned membership and counts), `family_alignments/`
(the family-wide alignments, `families/<family>/<family>_filtered_genes_db_TMonly_FAMSA.fasta.gz`, which
the cluster steps above read uncompressed; representative rows, mapped family alignments, topology alignments, kingdom groups), `anchors/`,
`website_inputs/` (the files of gpcrevolution.org the mapping reads: superfamily alignments with the
logo anchors, GPCRdb numbering files of the human receptors, and the logo files used as templates),
`consensus_mapping/` (column mappings and logo data), `gpcrdb_labels/` and `figure4/` (the structure
session and render, the exported logo, and Figure 4 as SVG and PNG).

The Dikarya reassignment was made on an earlier set of trees, from the September 22 representatives
(the eight best trees and the 100 ML search trees of five runs); these trees, the representative table and
their alignment are in `phylogeny_reassignment/inputs_Sep22/`, so `tree_side_per_sequence.py` and
`map_fungal_classA_to_clades.py` rerun from the archive. The trees only define the two reference groups;
each sequence is assigned by its best pairwise match.
