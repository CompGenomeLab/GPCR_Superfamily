# Family processing and motif/logo preparation

Large family FASTAs, MMseqs2 databases, cluster tables, and family-wide alignments are distributed through the Zenodo archive described in the [top-level README](../../README.md#data-on-zenodo). This folder contains the scripts that convert manually extracted homolog sets into the filtered and clustered family datasets used for taxonomic summaries, representative selection, and motif analysis.

## Upstream inputs

The principal inputs come from `post_processing/identify_gene_clade.py`:

```text
extracted_clades/alignments/<family>_*_treeordered_extracted_clade.fasta
```

Class A additionally uses the 80%-identity search-stage cluster table so selected Class A representatives can be expanded back to the complete hit set:

```text
search_results/classA/classA_hmm_E1000_db_clust0.8.tsv
```

Set its exact downloaded location explicitly. This avoids copying the table
out of the Class A search folder:

```bash
export GPCR_CLASSA_CLUSTER_TSV=/path/to/zenodo/search_results/classA/classA_hmm_E1000_db_clust0.8.tsv
```

## Step 1 - Build the family accession table

Script: `prepare_family_accession_table.py`

The script reads every extracted family FASTA, converts headers to UniProt accessions, and writes a table with one column per family. It also:

1. Expands Class A cluster representatives to the full cluster membership.
2. Removes accessions from broad families when they were assigned to a more specific family.
3. Writes one-accession-per-line files for later MMseqs2 sequence retrieval.
4. Optionally writes an interactive family-overlap plot.

The study configuration removes `GPR1` and `Vomeronasal1` accessions from Class A and removes `Vomeronasal2` accessions from Class C.

### Representative table structure

The accession table is a wide TSV: each column is a family and each non-empty cell is one UniProt accession. Families contain different numbers of accessions, so blank cells are expected near the bottom of shorter columns.

| classA | classB1 | classC | STE2 | GPR1 | FSLB |
|---|---|---|---|---|---|
| A0A123 | A0A456 | A0A789 | Q12345 | A0A111 | Q55DT7 |
| A0A124 | A0A457 | A0A790 | Q12346 | A0A112 | A0A222 |
| A0A125 |  | A0A791 |  | A0A113 |  |

The accessions above are illustrative. The actual table contains the accessions extracted for the study.

`LRLA_DICDI` is retained in `search_results/` as a control search but is not
propagated into the 23-family accession table. Five Olfactory accessions that
were present in the archived study assignment but fall outside the retained
tree-boundary FASTA are recorded explicitly in
`curated_accession_additions.tsv`. The script applies this small curation table
automatically; reproducing the study does not require editing Python code.

Important outputs:

| Output | Contents |
|---|---|
| `tables/accessions_by_family_*.tsv` | One accession column per extracted family |
| `accessions/<family>_accessions.txt` | Family-specific accessions for sequence retrieval |
| `tables/overlap_plot_*.html` | Pairwise overlap diagnostic |

Olfactory is not expanded with this table. It was not clustered during the
search stage; its 70% cluster table was generated later together with the other
17 final families.

Run after sourcing `config/paths.sh`:

```bash
python code_and_data/family_processing/prepare_family_accession_table.py
```

## Step 2 - Apply E-value and taxonomic mapping

The accession table does not contain E-values. It records which accession is assigned to which family and acts as a whitelist for the next step.

This step requires the original search-output files in addition to the accession table. Set the root of the per-family Zenodo search folders first:

```bash
export GPCR_SEARCH_RESULTS_DIR=/path/to/zenodo/search_results
```

- HMMER text output (`*_hmmsearch_E1000.txt`) for every HMM-based family.
- MMseqs2 tabular output (`*.m8`) for `GPR1_YEAST`, `FSLB_DICDI`, and `LRLA_DICDI`.

These are the original outputs produced by `homolog_search_pipeline/`; the alignments and trees are not used for the E-value lookup.

The taxonomy scripts stored directly in this folder scan the immediate family
folders under that root, so HMMER and MMseqs2 files remain in place. They read
these search results to obtain the E-value and taxID for every family-accession
assignment:

- HMMER output is used for HMM-based families.
- MMseqs2 `.m8` output is used for seed-search families.

For each accession in a family column, the taxonomy script finds the matching result from that family's search. If more than one acceptable result is available within that family search, it keeps the lowest E-value within that source search.

Class B1 and Class B2 are the exception to the one-folder-per-family rule.
They are separate extracted clades from the shared Class B search, so both
columns obtain their E-values and taxIDs from
`search_results/classB/classB_hmmsearch_E1000.txt`. The output `prefix` column
continues to record the assigned family (`classb1` or `classb2`), matching the
archived mapping tables.

The archived family databases themselves were constructed from
`acc_taxid_map_prefilter_Oct24_e1000_all.tsv`, before this downstream E-value
filter. `prepare_family_accession_files.py` converts that table into the exact
`accessions/<family>_accessions.txt` inputs whose sequence counts match the
archived `<family>_genes_db.fasta` files. The E-value-0.1 mapping is retained
separately for filtered taxonomic summaries and representative selection.

### Resolving accessions assigned to multiple families

Overlap is handled at two points:

1. **Predefined broad-family overlaps are removed in Step 1.** For example, accessions assigned to GPR1 or Vomeronasal 1 are removed from Class A, and Vomeronasal 2 accessions are removed from Class C.
2. **Remaining cross-family overlaps are resolved by E-value in Step 2.** The study configuration uses `dedupe_by_best_evalue=True`. If the same accession still has acceptable assignments to more than one family, the mapping rows are sorted by E-value and the accession is retained only for the family with the lowest E-value.

The resulting filtered mapping therefore contains at most one family assignment per accession in the configured study run.

The wide accession table is converted into a filtered, long-format mapping table such as:

| family | accession | taxID | evalue | prefix |
|---|---|---:|---:|---|
| classA | A0A123 | 9606 | 2.4e-35 | classa |
| classB1 | A0A456 | 7227 | 7.1e-18 | classb1 |
| STE2 | Q12345 | 559292 | 0.004 | ste2 |

Assignments with E-value greater than `0.1`, no search hit, or no recoverable taxID are not included in this mapping. They are reported separately in `missing_acc_taxid_*.tsv` for auditing.

Each retained mapping is also summarized into a taxID-by-family count table. In simplified form:

| taxID | classA_count | classB1_count | STE2_count |
|---:|---:|---:|---:|
| 9606 | 25 | 8 | 0 |
| 7227 | 12 | 3 | 0 |
| 559292 | 4 | 0 | 6 |

This is how the Step 1 accession table is used for E-value filtering: it defines the assignments to test, while the original HMMER/MMseqs2 files provide the values used to keep or remove them.

The current pipeline writes three canonical count views in one run:

| Output | Filtering |
|---|---|
| `tax_counts_per_family_unfiltered.tsv` | No downstream E-value cutoff and no global best-family deduplication. It is still limited to the extracted-clade accession assignments and the retained E≤1000 search output. |
| `tax_counts_per_family_e0.1_BestEval.tsv` | E-value ≤0.1; when an accession remains in multiple families, retain only its lowest-E-value family. |
| `tax_counts_per_family_e0.1_BestEval_allControls.tsv` | The preceding table restricted to taxa with `Galpha`, `Gbeta`, `Ggamma`, `RBP1`, and `TBP` counts all ≥1. |

The first two views also have matching long-format `acc_taxid_map_*.tsv`
files. The control-complete view is a row subset of the BestEval count table,
so a duplicate long-format mapping is not written.

Each view also has a matching pairwise accession-overlap plot:

- `overlap_plot_unfiltered.html`
- `overlap_plot_e0.1_BestEval.html`
- `overlap_plot_e0.1_BestEval_allControls.html`

These preserve the original overlap definition: an overlap is the same protein
accession assigned to two families. The unfiltered view can contain overlaps;
both BestEval views should report none because BestEval retains only one family
per accession.

Scripts:

- `create_taxID_family_count_table.py` creates the unfiltered, E-value/BestEval,
  and all-controls-present canonical views.
- `combine_tables_and_filter.py` is retained for reproducing the older October
  table merge, where GPCR and control counts came from separate historical files.

`combine_tables_and_filter.py` reads and writes tax-count tables under
`GPCR_FAMILY_TABLE_DIR` (by default, `families/tables`). Its output can be
overridden with `GPCR_CONTROL_FILTERED_TABLE`. The `families/accessions`
directory is not used here; it contains inputs for sequence retrieval rather
than tax-count tables.

```bash
python code_and_data/family_processing/create_taxID_family_count_table.py
python code_and_data/family_processing/combine_tables_and_filter.py
```

After taxonomy/E-value filtering, run `prepare_family_accession_files.py` to regenerate the accession lists used by the family database script:

```bash
python code_and_data/family_processing/prepare_family_accession_files.py
```

The preferred input is the mapping table written by `family_processing/create_taxID_family_count_table.py`, with at least `family` and `accession` columns.

## Step 3 - Retrieve, length-filter, cluster, and align each family

Script: `familydb_cluster_align.sh`

For each family listed in `gene_names`, the complete workflow is intended to:

1. Retrieve sequences from the UniProtKB MMseqs2 database using the filtered accession list.
2. Generate a FAMSA family alignment.
3. Trim columns containing at least 90% gaps with trimAl (`-gt 0.1`).
4. Remove sequences shorter than 200 ungapped residues in the trimmed alignment.
5. Rebuild the filtered MMseqs2 database.
6. Cluster filtered sequences at 70% identity using `-s 7.5 --cluster-reassign`.
7. Write cluster membership and representative-accession tables.
8. Generate the filtered family-wide MSA used for representative selection and downstream superfamily motif mapping.

Zenodo may store very gap-rich family-wide MSAs as `.fasta.gz`. The downstream
`msa_trim_ref.py` reader accepts either the uncompressed `.fasta` produced by
the family-processing job or its gzip-compressed archive form.

The shell script is a SLURM array job. Source `config/paths.sh`, set the UniProt
MMseqs2 database path, and use an array range
matching `gene_names`. It runs the complete workflow in order; no code blocks
need to be enabled or disabled:

```bash
sbatch code_and_data/family_processing/familydb_cluster_align.sh
```

For each array family it retrieves the filtered accessions, constructs and
trims the initial family alignment, applies the 200-residue filter, creates the
filtered MMseqs2 database and family-wide alignment, performs 70%-identity
clustering, and exports the cluster representatives. Existing outputs should be
moved to a separate retained-results directory before a clean reproduction run.

## Downstream superfamily alignment and motif mapping

The filtered family-wide MSAs are passed to [`superfamily_alignment/`](../superfamily_alignment/README.md). That stage constructs the extended representative alignment and connects its columns back to the family-wide MSAs for the manuscript and website motif/logo analysis.

## Required software

- Python 3 with pandas and Biopython
- MMseqs2
- FAMSA 2.2.2
- trimAl 1.4
- SLURM for the supplied array script

## Main outputs passed downstream

| Output category | Used by |
|---|---|
| Filtered accession/taxID tables | taxonomy summaries and website phyletic distributions |
| 70%-identity cluster tables | representative selection for phylogenetic inference |
| Filtered family-wide MSAs | extended superfamily alignment and motif/logo mapping |
| Cluster representative FASTAs | superfamily alignment and phylogenetic inference |
