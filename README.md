# GitHub repository for the manuscript "A Unified Evolutionary Framework for the GPCR Superfamily"

Code and small input files supporting the identification of GPCR homologs across eukaryotes, curation of family-level homolog sets, taxonomic mapping, representative selection, superfamily phylogeny, motif analysis, and structural comparisons.

## Quick start: use GitHub code with Zenodo data

The GitHub repository contains the analysis code, documentation, and small
inputs. The accompanying Zenodo archive contains the larger retained search,
curation, and family-processing datasets. The code reads the archive through
environment variables; do not copy Zenodo files into the GitHub repository.

Place the two downloads beside one another:

```text
project/
|-- GPCR_Superfamily/   # GitHub repository
`-- zenodo_data/        # unpacked Zenodo archive
    |-- search_results/
    |-- extracted_clades/
    `-- families/
```

From the GitHub repository on a Linux or HPC system, install the Python
requirements, connect the archive, and validate the downloaded archive:

```bash
cd project/GPCR_Superfamily
python -m pip install -r requirements.txt
source config/paths.sh
python code_and_data/validate_zenodo_archive.py
```

FASTA alignments and search-result text files in the Zenodo archive are
gzip-compressed. Before running a pipeline stage, decompress them in the local
working copy so every input has its standard filename:

```bash
find "$GPCR_ZENODO_DIR" -type f -name '*.gz' -exec gunzip {} +
python code_and_data/validate_zenodo_archive.py
```

This replaces the compressed files in the local working copy with their
uncompressed equivalents. Retain the original downloaded archive separately
if you need to preserve the deposited compressed files.

If the unpacked Zenodo archive is stored elsewhere, set its location before
sourcing the resolver:

```bash
export GPCR_ZENODO_DIR=/path/to/zenodo_data
source config/paths.sh
python code_and_data/validate_zenodo_archive.py
```

The validation command accepts either compressed or uncompressed archive files
and confirms that the archive is complete and matches the expected pipeline
handoffs. Most users can work directly from the retained Zenodo results.
Re-running the earliest search stages also requires the external UniProt
FASTA/MMseqs2 database and the Linux/HPC software listed in the stage-specific
READMEs.

## Repository layout

```text
GPCR_Superfamily/
|-- code_and_data/
|   |-- homolog_search_pipeline/   # Step 1 - HMMER/MMseqs2 searches, alignments, exploratory trees
|   |-- post_processing/           # Step 2 - score ordering and homolog-clade extraction
|   |-- family_processing/         # Step 3 - accessions, E-values, taxonomy, clustering
|   |-- phylogeny/                 # Step 4 - representative selection and superfamily phylogeny
|   |   `-- data/                  # alignments, trees, and MSA visualization
|   |-- superfamily_alignment/     # Step 5 - extended alignment and motif/logo mappings
|   |-- structure_visualization/   # Step 6 - PyMOL structural comparisons
|-- config/
|-- requirements.txt
`-- .gitignore
```

All analysis scripts and folder-level documentation live under `code_and_data/`. Each major stage has a README describing its role, inputs, outputs, analysis settings, and relationship to the next stage.

## File checksums

A checksum is a short fingerprint calculated from a file's exact bytes. It
does not change the analysis and is not another data format. It lets a user
verify that a downloaded file is byte-for-byte identical to the deposited
copy and was not truncated or corrupted.

After the Zenodo ZIP is final, create `SHA256SUMS.txt` from the repository
root with:

```bash
python code_and_data/create_checksums.py ../zenodo_data.zip \
  code_and_data/structure_visualization/superfamily_logo_positions_overall_look.pse \
  code_and_data/structure_visualization/motifs_across_families.pse
```

Generate the manifest only after the release files stop changing; editing or
recompressing a file correctly produces a different checksum.

## Workflow (in order)

1. **Homolog discovery** - Search a UniProtKB-derived eukaryotic protein database with HMMER or iterative MMseqs2 and retrieve top-scoring sequences. HMMER searches produce alignments and trees from both complete sequences and HMM-recognized regions. MMseqs2 searches produce alignments and trees from complete sequences. Alignments are generated with FAMSA, trimmed with trimAl, and used to infer FastTree trees.
2. **Homolog-clade curation** - Order trees and alignments by HMMER/MMseqs2 score, inspect score transitions, and extract manually selected homologous clades.
3. **Family processing and taxonomy** - Expand Class A cluster representatives, remove predefined overlaps, use the original search results to apply E-value filtering and resolve remaining overlaps, map retained accessions to taxIDs, construct family-count tables, apply the minimum-length filter, cluster family sequences, and build family-wide alignments.
4. **Representative selection and phylogeny** - Rank candidates using cluster size and taxonomic diversity, retain the study representatives, construct and refine the superfamily alignment, infer the primary phylogeny, calculate support and consensus trees, and render the final alignment used in the manuscript figure.
5. **Extended superfamily alignment and motif analysis** - Add the divergent family representatives omitted from the primary tree, construct the 550-sequence alignment, connect it to family-wide MSAs, and generate motif/logo mappings.
6. **Structural analysis** - Compare conserved positions in experimental structures or AlphaFold models.

Taxonomy scripts are stored directly in `family_processing/` so the complete accession-to-family-dataset workflow remains in one unit. The `superfamily_alignment/` stage supports the manuscript superfamily logo and the interactive [gpcrevolution.org/superfamily-logo](https://gpcrevolution.org/superfamily-logo) page.

## Search strategy and retained controls

- `GPR1_YEAST`, `FSLB_DICDI`, and `LRLA_DICDI` use iterative MMseqs2 searches initiated from seed FASTA sequences.
- All other configured searches use HMMER profile searches.
- STE2 remains included in the search, clade-curation, and family-processing stages.
- The heterotrimeric G-protein subunits (`Galpha`, `Gbeta`, and `Ggamma`) and the `TBP`/`RBP1` housekeeping controls remain included.
- RAMP, RGS, Arrestin, and Dvl are not part of the final configured analysis.

## Where to read more

| Folder | What it covers |
|---|---|
| [homolog_search_pipeline](code_and_data/homolog_search_pipeline/README.md) | Search inputs, HMMER/MMseqs2 SLURM scripts, retrieval limits, alignments, and exploratory trees |
| [post_processing](code_and_data/post_processing/README.md) | Score-ordered trees/FASTAs, line plots, manual clade boundaries, and extracted-clade outputs |
| [family_processing](code_and_data/family_processing/README.md) | Family accession tables, E-value/length filtering, taxonomic mapping, clustering, and family-wide MSAs |
| [phylogeny](code_and_data/phylogeny/README.md) | Representative selection, MAFFT/tree scripts, verified alignments and trees, and the MSA figure generator |
| [superfamily_alignment](code_and_data/superfamily_alignment/README.md) | Construction of the 550-sequence extended alignment and family/superfamily motif mappings |
| [structure_visualization](code_and_data/structure_visualization/README.md) | PyMOL sessions used for Figure 3b and Figure 4, conserved-position selections, and structural provenance |

Start with `homolog_search_pipeline`, continue through `post_processing` and `family_processing`, and then use `phylogeny`. The `superfamily_alignment` stage extends the representative set for motif/logo analysis before structural interpretation.

## Configuration

The Zenodo archive is a separate directory beside the GitHub repository:

```text
superfamily_paper/
|-- GPCR_Superfamily/
`-- zenodo_data/
    |-- search_results/
    |-- extracted_clades/
    `-- families/
```

With this layout, source the shared resolver before running operational scripts:

```bash
source config/paths.sh
```

If the archive is elsewhere, copy the example, replace only the archive and
source-database paths, and source it:

```bash
cp config/paths.env.example config/paths.env
source config/paths.env
```

`GPCR_ZENODO_DIR` points to the archive root. The resolver derives
`GPCR_ANALYSIS_DIR`, `GPCR_EXTRACTED_CLADE_DIR`, `GPCR_FAMILY_TABLE_DIR`,
`GPCR_ACCESSION_TXT_DIR`,
`GPCR_CLUSTER_DIR`, and the other stage paths from that one location.
`GPCR_SEARCH_RESULTS_DIR` points to the downloaded Zenodo `search_results`
tree, where every family retains its own search folder, and `GPCR_FAMILY_DIR`
points to the downloaded or generated `families` tree. The Class A expansion
table is selected with `GPCR_CLASSA_CLUSTER_TSV`. Source database paths are set
with `UNIPROT_FASTA` and `UNIPROT_MMSEQS_DB`; executable/module settings remain
overridable in the corresponding shell scripts.

## Data on Zenodo

The large pipeline archive is available on Zenodo: [doi:10.5281/zenodo.21573037](https://doi.org/10.5281/zenodo.21573037).

GitHub contains scripts, documentation, small search inputs, curation tables, the verified primary phylogeny, and the extended superfamily-alignment materials. Larger search, family-processing, and supporting outputs are archived on Zenodo.

The archive has three top-level data directories:

| Directory | Contents |
|---|---|
| `search_results/` | Per-family HMMER/MMseqs2 output, retrieved sequences, exploratory alignments/trees, ordered files, and the Class A search-stage cluster |
| `extracted_clades/` | Locally curated clade alignments and matching pruned trees produced from the ordered search results |
| `families/` | Family accession/taxonomy tables, `accessions/<family>_accessions.txt` inputs, filtered datasets, 70% clusters, family MSAs, representative sets, manually trimmed `trimends` alignments, and final superfamily analysis outputs |

Preserve this layout after downloading the archive. If it is not beside the
repository, set `GPCR_ZENODO_DIR` before sourcing `config/paths.sh`.

## Notes

- Only the Zenodo root and external source-database paths need to be supplied;
  derived analysis paths should not be edited individually.
- Dataset-specific filenames and manually curated clade boundaries are intentionally retained because they record the analysis used for the manuscript.
- If a required large file is not present on GitHub, obtain it from the Zenodo archive or regenerate it using the preceding stage.
