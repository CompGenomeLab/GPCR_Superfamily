# Extended superfamily alignment and motif analysis

This folder contains the scripts and verified outputs used to extend the primary 455-sequence phylogeny set into the 550-sequence superfamily alignment and connect that alignment to the family-wide MSAs. These mappings support the superfamily logo shown in the manuscript and the interactive [GPCR Evolution superfamily-logo page](https://gpcrevolution.org/superfamily-logo).

The scripts are stored directly in this folder, following the organization of the other analysis stages. Processed inputs and outputs are grouped under [`data/`](data/README.md).

## Upstream inputs

The primary phylogeny contributes 455 representatives. The extended alignment adds the divergent families that were omitted from tree inference:

| Added family | Representatives |
|---|---:|
| Class C | 50 |
| Vomeronasal 2 | 20 |
| STE3 | 25 |
| **Total added** | **95** |

The resulting accession list contains 550 unique representatives. All 455 primary-phylogeny representatives occur in this extended set.

Family-wide alignments from [`family_processing/`](../family_processing/README.md) are used later to connect family positions to the extended superfamily alignment.

## Step 1 - Select the extended representative set

The representative-selection code in [`phylogeny/select_superfamily_representatives.py`](../phylogeny/select_superfamily_representatives.py) selects the complete family sampling using MMseqs2 cluster sizes and taxonomic diversity. Singleton clusters are excluded.

The full output contains:

- 100 cAMP representatives
- 100 Class A representatives
- 50 Class C representatives
- 25 representatives each for GP143, GP157, GPR1, and STE3
- 20 representatives each for Class B1, Class B2, Class F, Class T, FSLB, Mth, Nematode, Olfactory, Vomeronasal 1, and Vomeronasal 2

The saved 550-accession list is included in `data/representatives/`.

## Step 2 - Retrieve representative sequences

Script: `retrieve_superfamily_representatives.sh`

Input: the 550-accession list and the UniProtKB-derived MMseqs2 protein database.

Output: the unaligned 550-sequence representative FASTA.

```bash
bash retrieve_superfamily_representatives.sh
```

## Step 3 - Construct the initial alignment

Script: `initial_linsi_alignment.sh`

The full-length representative sequences are aligned with MAFFT L-INS-i using `--localpair --maxiterate 1000 --ep 0.123`.

```bash
bash initial_linsi_alignment.sh
```

Input: the unaligned representative FASTA.

Output: the initial 550-sequence L-INS-i alignment in `data/alignments/`.

## Step 4 - Retain the shared 7TM region

The initial alignment was inspected in Jalview. Variable terminal and poorly aligned regions outside the shared seven-transmembrane core were removed manually.

Input: the initial L-INS-i alignment.

Output: the manually trimmed 7TM alignment in `data/alignments/`.

## Step 5 - Infer the preliminary guide tree

Script: `preliminary_guide_tree.sh`

The script:

1. Trims the manually defined 7TM alignment with ClipKIT `kpic-gappy`.
2. Runs IQ-TREE2 ModelFinder using `-m MF -mfreq FO -cmin 5 -cmax 15 -T 10 --mset JTT --mrate R`.
3. Infers a preliminary RAxML-NG tree from 20 parsimony starting trees.

ModelFinder selected `JTT+FO+R8` for this alignment.

```bash
bash preliminary_guide_tree.sh
```

## Step 6 - Convert the guide tree for MAFFT

Script: `prepare_mafft_treein.py`

The script reconciles the preliminary tree and trimmed FASTA, reorders the FASTA to tree-leaf order, numbers the leaves, and writes MAFFT-compatible merge steps.

```bash
python prepare_mafft_treein.py
```

Outputs include the reordered FASTA, numbered Newick tree, MAFFT steps, and name-mismatch reports. The verified tree and FASTA contain the same 550 sequence names.

## Step 7 - Construct the refined 550-sequence alignment

Script: `tree_guided_ginsi_alignment.sh`

The manually trimmed sequences are realigned with MAFFT G-INS-i using the preliminary tree:

```bash
bash tree_guided_ginsi_alignment.sh
```

The script uses `--treein --globalpair --maxiterate 1000 --ep 0.123`. The resulting tree-guided alignment in `data/alignments/` is the extended superfamily alignment used for motif and logo analysis.

## Step 8 - Select family references for motif mapping

Script: `select_superfamily_representatives_for_logo_generation_based_on_coverage.py`

This script calculates coverage and consensus-related metrics for sequences in the 550-sequence alignment. These metrics were inspected when choosing one representative per family to link the family-wide and superfamily-wide alignments. The automatic ranking is diagnostic; final references were curated manually.

The selected reference alignments and representative tables are included under `data/alignments/` and `data/logo_outputs/`.

## Step 9 - Prepare and trim the family mappings

Scripts:

- `prepare_trim_info.py`
- `msa_trim_ref.py`

`trim_info.tsv` records:

```text
acc1    family    acc2
```

- `acc1` is the representative shared by the family-wide and superfamily-wide alignments.
- `acc2` is an optional annotated receptor used for GPCRdb numbering or region annotations.

Run:

```bash
python prepare_trim_info.py
python msa_trim_ref.py --mapping trim_info.tsv
```

The trimming script reads the family-wide MSAs produced by `family_processing/` and retains the columns required to connect their residue positions to the superfamily alignment.

## Step 10 - Generate website mappings

Script: `precompute-mappings.mjs`

The script combines the reference-trimmed family alignments, the selected superfamily references, `trim_info.tsv`, and family conservation annotations. It writes one JSON mapping per family and an index for the interactive website.

```bash
source config/paths.sh
node precompute-mappings.mjs
```

The mapping script uses the website checkout by default. Set
`GPCR_CUSTOM_MSA_DIR`, `GPCR_CONSERVATION_DIR`, `GPCR_MAPPING_OUTPUT_DIR`, or
`GPCR_SUPERFAMILY_MAPPING_MSA` when reproducing the mappings from files stored
elsewhere.

The manuscript view uses binary conservation at a 60% threshold, includes all configured families, and displays the top 36 residues in the scatter-plot view.

## Dependencies

- Python 3 with pandas and NumPy
- MMseqs2
- MAFFT
- Jalview for manual 7TM boundary selection
- ClipKIT
- IQ-TREE2
- RAxML-NG
- Node.js for website mapping generation
