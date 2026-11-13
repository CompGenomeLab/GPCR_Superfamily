# Analysis code

One folder per stage, in the order of the Methods. Each stage reads the outputs of the stages before
it from the Zenodo archive and writes its own outputs to the archive folder of the same name
(`gpcr_paths.py` resolves both; see the [top-level README](../README.md)). Each folder's README gives
the method, the scripts in the order they run, and the archive files.

| Stage | Methods section | Main result |
|---|---|---|
| [`01_homology_search/`](01_homology_search/README.md) | Homology Searches and Sequence Retrieval | search hits, alignments and gene trees of every search |
| [`02_clade_extraction/`](02_clade_extraction/README.md) | Extracting Inferred Homologs: MSA, tree ordering, extraction criteria | clade membership of the 19 families (492,350 sequences) |
| [`03_topology_filter/`](03_topology_filter/README.md) | Controls and Filtering: 7TM pattern and topology score | 408,489 single-core and 391,856 topology-filtered sequences |
| [`04_overlap_assignment/`](04_overlap_assignment/README.md) | Controls and Filtering: overlap assignment and controls | final membership and per-genome counts |
| [`05_literature_comparison/`](05_literature_comparison/README.md) | Comparison with Published Classifications | Tables S3a-d, the 18 literature sequences |
| [`06_family_origins/`](06_family_origins/README.md) | Taxonomic Distribution Visualization; Inference of Family Origins | Figure 2, Dollo origins, Table S4 |
| [`07_representative_selection/`](07_representative_selection/README.md) | Representative Protein Sequence Selection | the 590 representatives and three sets |
| [`08_superfamily_alignment/`](08_superfamily_alignment/README.md) | Superfamily Alignment | eight superfamily alignments |
| [`09_superfamily_phylogeny/`](09_superfamily_phylogeny/README.md) | Phylogenetic Inference of the GPCR Superfamily Hierarchy | eight ML trees, support, Tables S6-S7 |
| [`10_motif_analysis/`](10_motif_analysis/README.md) | Tree-Guided Superfamily Alignment | family-column mapping, logos, motif positions, Figure 4, Tables S8-S9 |
| [`11_structure_superposition/`](11_structure_superposition/README.md) | Superimposition of Identified Motifs | Figures 5 and 6, Fig. S3, Table S10 |
| [`supplementary_tables/`](supplementary_tables/README.md) | | Tables S1-S10 and the supplementary-figure document |

Shell scripts marked as cluster jobs are SLURM arrays and source `config/paths.sh`; everything else
runs locally with the packages in `requirements.txt` (stage 11 also needs PyMOL).
