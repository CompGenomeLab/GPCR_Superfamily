# GPCR superfamily representative selection and phylogeny

This folder contains the representative-selection code, the verified primary phylogeny inputs and outputs, and the alignment visualization used for the manuscript figure.

The primary phylogeny contains 455 representatives. Class C, Vomeronasal 2, and STE3 were omitted from tree inference because of their high divergence. The [`superfamily_alignment/`](../superfamily_alignment/README.md) stage adds these families back and constructs the extended alignment used for motif/logo analysis.

The scripts are stored directly in this folder, following the same organization as the other analysis stages. The verified outputs are grouped under [`data/`](data/README.md).

## Pipeline

The scripts use the included primary phylogeny files by default. Paths can be supplied as positional arguments when running another dataset. Executable names can be overridden with `MAFFT_BIN`, `CLIPKIT_BIN`, `IQTREE_BIN`, or `RAXML_BIN`.

### Step 1 - Select representatives

Script: `select_superfamily_representatives.py`

The script uses family cluster tables, the filtered accession/taxID table, and the taxonomy JSON produced upstream. It ranks candidate representatives using cluster size and taxonomic diversity. Clusters containing only one sequence are excluded. The complete selection supplies the extended 550-sequence alignment described in [`superfamily_alignment/`](../superfamily_alignment/README.md).

With `GPCR_ZENODO_DIR` set to the downloaded archive root, the default inputs are:

- `families/clusters/*_genes_filtered_db_clust0.7_clust0_cov1.tsv`
- `families/tables/acc_taxid_map_BestEval_May22_e0.1.tsv`
- `families/reference_data/taxonomy_eukaryotes_filtered.json`

The staged cluster tables under `families/clusters/` must be byte-identical to their per-family copies under `families/<family>/`. Run `python code_and_data/validate_zenodo_archive.py` from the repository root before representative selection; the validator reports any disagreement. The retained 550-accession result is `superfamily_alignment/data/representatives/All_clust0.7_cov1_Oct25_minsize2.txt`. The separately written `manual_superfamily_representatives.tsv` contains only the 17 manually curated family-linking references used for the logo analysis and is not the 550-sequence selection output.

Class C, Vomeronasal 2, and STE3 are removed from that complete set for primary tree inference, leaving 455 representatives. The filtered accession list and corresponding unaligned FASTA are stored in `data/representatives/`.

### Step 2 - Initial MAFFT L-INS-i alignment

Input: the 455-sequence representative FASTA.

Output: the initial full-length L-INS-i alignment.

```bash
bash initial_linsi_alignment.sh
```

The script uses `--localpair --maxiterate 1000 --ep 0.123`.

### Step 3 - Manual 7TM trimming

The initial alignment was inspected in Jalview. Variable N- and C-terminal regions and poorly aligned regions outside the shared seven-transmembrane core were removed manually.

In the retained analysis, filenames containing `trimends` identify these
manually trimmed alignments. The Zenodo manifest labels them as manual outputs
and preserves their original filenames.

Input: the initial L-INS-i alignment in `data/alignments/`.

Output: the manually trimmed 7TM alignment in the same folder.

### Step 4 - Preliminary tree

Input: the manually trimmed 7TM alignment.

Outputs: the initial ClipKIT `kpic-gappy` alignment, ModelFinder results, and preliminary RAxML-NG tree.

```bash
bash preliminary_tree.sh
```

ModelFinder uses:

```text
-m MF -mfreq FO -cmin 5 -cmax 15 -T 10 --mset JTT --mrate R
```

ModelFinder selected `JTT+FO+R8`. The preliminary RAxML-NG search used 20 parsimony starting trees.

### Step 5 - Prepare the MAFFT guide tree

Input: the manually trimmed 7TM FASTA and preliminary RAxML-NG best tree.

Outputs: a reordered FASTA, numbered Newick tree, MAFFT merge steps, and mismatch reports.

```bash
python prepare_mafft_treein.py
```

The saved preliminary tree and FASTA contain the same 455 sequence names.

### Step 6 - Tree-guided MAFFT G-INS-i alignment

Input: the reordered 7TM FASTA and MAFFT guide-tree steps.

Output: the refined tree-guided alignment.

```bash
bash tree_guided_ginsi_alignment.sh
```

The script uses `--treein --globalpair --maxiterate 1000 --ep 0.123`.

### Step 7 - Final phylogenetic inference

Input: the refined G-INS-i alignment.

Outputs: the final ClipKIT alignment, ModelFinder report, 100 ML searches, up to 500 bootstrap trees, the best ML tree, and TBE support.

```bash
bash final_tree.sh
```

The final alignment contains 455 sequences and 397 columns. ModelFinder selected `JTT+FO+R8`. RAxML-NG used 100 parsimony starting trees and completed all 500 bootstrap replicates under the `autoMRE{500}` limit.

### Step 8 - MRE consensus and support

Input: the 100 inferred ML trees and final bootstrap trees.

Outputs: the extended majority-rule consensus tree and its TBE-supported version.

```bash
bash mre_consensus.sh
```

The annotated family-labeled manuscript tree is included in `data/trees/` as an SVG and Newick file.

### Step 9 - Render the alignment figure

The visualization script reads the final ClipKIT alignment used by RAxML-NG.

```bash
python MSA_to_svg_and_pdf.py
```

The SVG is written to `data/visualization/`. PDF conversion is attempted when `svglib` and ReportLab are installed.

## Dependencies

- Python 3 with pandas and NumPy
- MMseqs2 for the upstream family clustering
- MAFFT
- Jalview for the documented manual 7TM boundary selection
- ClipKIT
- IQ-TREE2
- RAxML-NG
- ReportLab and svglib for optional PDF output
