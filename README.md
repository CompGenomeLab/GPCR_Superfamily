# GPCR Superfamily Is Defined by Shared Motifs

Code for the identification of GPCR families across eukaryotes, their filtering and assignment, the
comparison with published classifications, the inference of family origins, the superfamily trees,
the shared sequence motifs and their structural test.

The data every stage reads and writes are in the accompanying Zenodo archive,
[doi:10.5281/zenodo.21573037](https://doi.org/10.5281/zenodo.21573037). Interactive versions of the
figures are on [gpcrevolution.org](https://gpcrevolution.org): the phyletic distribution, the
superfamily trees and the superfamily logo.

## Layout

```text
project/
|-- GPCR_Superfamily/          this repository
|   |-- code_and_data/         one folder per stage, in the order of the Methods
|   |-- config/paths.sh        archive and database locations for the shell (SLURM) scripts
|   `-- requirements.txt
`-- zenodo_data/               the unpacked archive, one folder per stage with the same names
```

| Stage | What it does |
|---|---|
| `01_homology_search` | HMMER and MMseqs2 searches of the UniProt 2025_02 eukaryotic reference proteomes, FAMSA alignments and FastTree gene trees |
| `02_clade_extraction` | score-ordered trees and alignments, the extracted clades of every family, clade membership |
| `03_topology_filter` | DeepTMHMM, the single GPCR-like 7TM core, the topology score and its 0.8 cutoff |
| `04_overlap_assignment` | sequences claimed by two families, final membership, counts per genome |
| `05_literature_comparison` | receptors of 15 published surveys matched to our families |
| `06_family_origins` | taxonomy of the searched genomes, Dollo origins, Figure 2 |
| `07_representative_selection` | 590 representatives by taxonomic diversity and topology score |
| `08_superfamily_alignment` | 7TM cores of the four sequence sets, G-INS-i alignments, ClipKIT |
| `09_superfamily_phylogeny` | ModelFinder, RAxML-NG, TBE, transfer index, stems, recovery, sCF |
| `10_motif_analysis` | tree-guided alignment, family-column mapping, logos, conserved positions, Figure 4 |
| `11_structure_superposition` | structure superimposition, motif distances, Figures 5-6 |
| `supplementary_tables` | Tables S1-S10 and the supplementary-figure document |

## Getting started

```bash
cd project/GPCR_Superfamily
python -m pip install -r requirements.txt
python code_and_data/check_archive.py        # confirms the archive is complete and in place
```

The archive is found beside the repository; set `GPCR_ZENODO_DIR` when it is elsewhere. Scripts read
gzip-compressed archive files directly, so nothing has to be unpacked except the archive itself.
Python scripts run locally; the shell scripts are SLURM jobs for the steps run on a cluster (searches,
DeepTMHMM, alignments, tree inference, site concordance factors) and expect HMMER, MMseqs2, FAMSA,
trimAl, FastTree, MAFFT, ClipKIT, IQ-TREE 2 and RAxML-NG on the path. Rerunning the searches also needs
the UniProt 2025_02 eukaryotic reference proteome FASTA and its MMseqs2 database (`config/paths.sh`).

Every stage can be run on its own from the archive. Each stage README says which of its results were
checked by rerunning the stage from the archive against the files used in the paper.

## Family names

Files use short internal names: `classA`, `classB1`, `classB2`, `classC`, `classF`, `classT` (Taste-2),
`FSLB` (FSL), `Vomeronasal1`, `Vomeronasal2`, `GP143`, `GP157`, `GPR1`, `TM116`, `Mth`, `Nematode`,
`Olfactory`, `cAMP`, `STE2`, `STE3`; the controls are `Galpha`, `Gbeta`, `Ggamma`, `TBP` and `RBP1`
(RPB1). The search `classT` is shared by Taste-2 and Vomeronasal 1, and `classB` by Class B1 and B2.
