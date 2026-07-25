# Homolog search pipeline

Full HMMER/MMseqs2 output, retrieved sequence databases, alignments, and exploratory trees are distributed through the Zenodo archive described in the [top-level README](../../README.md#data-on-zenodo). This folder contains the search scripts, helper programs, and the exact HMM/FASTA query inputs used to identify GPCR-family homologs in a UniProtKB-derived eukaryotic protein database.

Two search scripts are retained because the analyses begin from different input types:

- `hmm_search.sh` starts from a profile HMM and runs HMMER.
- `mmseqs2_search.sh` starts from a protein seed FASTA and runs iterative MMseqs2.

Both workflows retrieve target sequences, align them with FAMSA, trim gap-rich columns with trimAl, and infer exploratory trees with FastTree.

- **HMMER searches** produce one alignment and tree from the complete retrieved sequences and another alignment and tree from the HMM-recognized regions.
- **MMseqs2 searches** produce one alignment and tree from the complete retrieved sequences.

## Search inputs

Exact inputs are stored under [`input_sequences/`](input_sequences/README.md):

- HMMER profile searches: 20 HMM files, one per configured entry.
- MMseqs2 seed searches: `GPR1_YEAST`, `FSLB_DICDI`, and `LRLA_DICDI` FASTAs.

All other configured entries use HMMER, including STE2, the GPCR families, the three heterotrimeric G-protein subunits, and the TBP/RBP1 controls.

## Required databases

Both shell scripts require a UniProtKB-derived eukaryotic protein database prepared for MMseqs2:

```bash
UNIPROT_MMSEQS_DB=/path/to/UniprotRefEuk
```

The HMMER workflow also requires the same protein collection as FASTA:

```bash
UNIPROT_FASTA=/path/to/Uniprot_Ref_Eukaryota_Proteome.fasta
```

The study MMseqs2 database was created on the HPC system with:

```bash
mmseqs createdb "$UNIPROT_FASTA" "$UNIPROT_MMSEQS_DB"
```

The manuscript records a retrieval date of May 9, 2025. UniProt release
2025_02 was released on April 23, 2025, and was the active release on that
date; release 2025_03 followed on June 18, 2025. The study database is
**UniProt 2025_02**.

Source `config/paths.sh` before submission. Search outputs then go to
`$GPCR_ZENODO_DIR/search_results` by default, while the HMM and seed inputs are
read from this GitHub repository. An external archive location can be selected
by setting `GPCR_ZENODO_DIR` before sourcing the resolver.

## Step 1A - HMMER profile searches

Script: `hmm_search.sh`

For the family selected by the SLURM array index, the script:

1. Runs `hmmsearch --max -E 1000` against the eukaryotic protein FASTA.
2. Uses `hmm_header_get.py` to retain the configured number of top-ranked hits.
3. Retrieves those sequences from the MMseqs2 database.
4. Clusters the 500,000-hit Class A result at 80% identity and retains cluster representatives for exploratory alignment/tree construction.
5. Aligns the complete retrieved sequences with FAMSA, trims the alignment with trimAl, and infers a FastTree tree.
6. Uses `trim_sequences_using_hmm.py` to extract the HMMER domain envelope from
   each sequence. The envelope spans `envfrom` through `envto`: it contains the
   directly aligned domain and any short flanking residues HMMER assigns to the
   same domain. This preserves the coordinate choice used in the study.
7. Aligns the HMM-region sequences with FAMSA, trims the alignment with trimAl, and infers a second FastTree tree.

Family-specific retrieval limits are encoded in the `max_hits` array. Examples include 500,000 for Class A, 200,000 for Olfactory, 20,000 for Class T/Taste-2, and 10,000 for several smaller GPCR families.

Submit after setting database and project paths:

```bash
sbatch code_and_data/homolog_search_pipeline/hmm_search.sh
```

## Step 1B - Iterative MMseqs2 seed searches

Script: `mmseqs2_search.sh`

For `GPR1_YEAST`, `FSLB_DICDI`, and `LRLA_DICDI`, the script:

1. Creates an MMseqs2 query database from the family seed FASTA.
2. Runs five search iterations with `-e 1000 --e-profile 0.1 -s 7.5`, prefilter mode 0, and up to 1,000,000 candidates.
3. Exports a 15-column `.m8` table containing scores, coverage, and E-values.
4. Retrieves the configured top target sequences.
5. Generates a FAMSA alignment, trimAl-trimmed alignment, and FastTree tree.

Submit with:

```bash
sbatch code_and_data/homolog_search_pipeline/mmseqs2_search.sh
```

## Output layout

Each array task writes to:

```text
search_results/<family>/
```

The Zenodo archive uses this same layout for both HMMER and MMseqs2 searches,
so downstream scripts can read it directly through `GPCR_SEARCH_RESULTS_DIR`.

### HMMER outputs

| Pattern | Contents |
|---|---|
| `<family>_hmmsearch_E1000.txt` | Standard HMMER output used for ranking hits and extracting HMM coordinates |
| `<family>_hmmsearch_E1000_headers.txt` | Top-ranked accessions used for MMseqs2 target retrieval |
| `<family>_hmm_E1000_db` | Retrieved MMseqs2 target database |
| `classA_hmm_E1000_db_clust0.8.tsv` | Class A 80%-identity cluster membership used for later expansion |
| `<family>_hmm_E1000_db.fasta` | Complete retrieved protein sequences |
| `<family>_hmm_E1000_db_FAMSA.fasta` | FAMSA alignment of complete sequences |
| `<family>_hmm_E1000_db_FAMSA_gt0.1.fasta` | TrimAl alignment of complete sequences |
| `<family>_hmm_E1000_db_FAMSA_gt0.1.tree` | FastTree tree from complete sequences |
| `<family>_hmm_E1000_db_hmmregions.fasta` | HMM-recognized region from each sequence |
| `<family>_hmm_E1000_db_hmmregions_FAMSA.fasta` | HMM-region FAMSA alignment |
| `<family>_hmm_E1000_db_hmmregions_FAMSA_gt0.1.fasta` | TrimAl alignment of HMM regions |
| `<family>_hmm_E1000_db_hmmregions_FAMSA_gt0.1.tree` | FastTree tree from HMM regions |

### MMseqs2 outputs

MMseqs2 searches do not create HMM-region files. Their main outputs are:

| Pattern | Contents |
|---|---|
| `*.m8` | MMseqs2 tabular search output |
| `<family>_db.fasta` | Complete retrieved protein sequences |
| `<family>_db_FAMSA.fasta` | FAMSA alignment of complete sequences |
| `<family>_db_FAMSA_gt0.1.fasta` | TrimAl alignment of complete sequences |
| `<family>_db_FAMSA_gt0.1.tree` | FastTree tree from complete sequences |

## Helper scripts

### `hmm_header_get.py`

Reads standard HMMER text output, keeps the first requested hits, normalizes UniProt-style identifiers to accessions, and writes `<search>_headers.txt` for `mmseqs createsubdb`.

```bash
python hmm_header_get.py family_hmmsearch_E1000.txt 10000
```

### `trim_sequences_using_hmm.py`

Parses HMM domain coordinates and writes one continuous HMM-aligned span for each retrieved sequence. The span begins at the minimum `ali-from` and ends at the maximum `ali-to` reported for that target.

```bash
python trim_sequences_using_hmm.py family_hmmsearch_E1000.txt family_db.fasta
```

## Required software

- SLURM
- HMMER 3.3.2
- MMseqs2
- FAMSA 2.2.2
- trimAl 1.4
- FastTree 2.1.11 (`FastTreeMP`)
- Python 3 with Biopython

The next stage is [`post_processing/`](../post_processing/README.md), which orders the exploratory trees and alignments by search score and extracts the curated homologous clades.
