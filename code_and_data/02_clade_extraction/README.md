# 02 · Extracting the homologs of each family

Methods: *Extracting Inferred Homologs for Each Family* (MSA generation, ordering the gene tree and
MSA, extraction criteria); Table S1 ("Clade members").

The leaves of every search tree are ordered by search score: the tree is rooted on its lowest-scoring
leaf and, recursively, the sister clade with the higher mean score is placed on the left
(`tree_order_by_hmm.py`). The alignment is put in the same order (`fasta_order_by_tree.py`) and the
scores are drawn along the leaf order (`line_plot_hmm_and_index_treeorder.py`), so that high-scoring
clades, their alignment patterns and the drops in score where other families begin can be read side
by side. The clades chosen this way are recorded in `identify_gene_clade.py` as the leftmost and
rightmost leaf of each clade (`gene_boundaries`; a third element `-1` marks a clade subtracted from
the family).

`write_clade_membership.py` turns the extracted clades into family membership:

1. Class A representatives are expanded to all members of their 80% identity cluster (stage 01).
2. A sequence retrieved both by a broader search and by the search of a more specific family is given
   to the specific family: GPR1 and Vomeronasal 1 are removed from Class A, Vomeronasal 2 from Class C.
3. Five Olfactory accessions are added from `curated_accession_additions.tsv`.

```bash
source ../../config/paths.sh
sbatch tree_fasta_order_and_line_plot.sh     # ordered trees, alignments and line plots (stage 01 folders)
python identify_gene_clade.py                # extracted_clades/
python write_clade_membership.py             # clade_membership.tsv, clade_membership_controls.tsv
```

Verified: the clade membership written from the archive equals the table used in the paper for all 19
families (492,350 sequences).

**Archive** (`02_clade_extraction/`):

| File | Contents |
|---|---|
| `extracted_clades/alignments/`, `extracted_clades/trees/` | alignment rows and pruned tree of every extracted clade |
| `clade_membership.tsv` | one column per family, the accessions of its clades (input of stage 03) |
| `clade_membership_controls.tsv` | the same for the five controls |
| `clade_sequences/<family>.fasta.gz` | full-length sequences of every clade member, with their UniProt headers |

**Website data** (`website/build_family_assignment_data.py`): writes the per-family data of the
gpcrevolution.org family-assignment page (search tree, binned MSA, search score and the clade after each
filter step) from the archive: `python website/build_family_assignment_data.py --out <folder>`.
