# 01 · Homology searches and sequence retrieval

Methods: *Homology Searches and Sequence Retrieval*; Table S1.

Every family is searched against the eukaryotic reference proteomes of UniProtKB release 2025_02,
either with a profile HMM (Pfam or KEGG KOfam; HMMER `hmmsearch --max -E 1000`) or, for the two
families without a profile (FSL and GPR1), with an iterative MMseqs2 profile search from a seed
sequence (`-e 1000 --e-profile 0.1 -s 7.5`, 5 iterations). A family-specific number of top-ranked
hits is retrieved (the `max_hits` array; 10,000 to 500,000), on purpose well beyond the size of the
family, so the search tree also contains neighbouring families.

Two searches are shared: `classT` (TAS2R, Pfam PF05296) yields both Taste-2 and Vomeronasal 1, and
`classB` (7tm_2, PF00002) yields Class B1 and Class B2. The heterotrimeric G protein subunits
(G alpha, G beta, G gamma) and the housekeeping proteins TBP and RPB1 (`RBP1`) are searched the same
way as controls.

Each retrieved set is aligned with FAMSA, columns with 90% or more gaps are removed with trimAl
(`-gt 0.1`), and a gene tree is inferred with FastTree. HMM searches produce two alignments and trees,
one of the complete sequences and one of the HMM-matched regions (`trim_sequences_using_hmm.py`).
The Class A hits are first clustered at 80% identity (`mmseqs cluster --min-seq-id 0.8 -s 7.5
--cluster-reassign`) and the cluster representatives are aligned; stage 02 expands the extracted
representatives back to all cluster members.

| File | Role |
|---|---|
| `hmm_search.sh` | SLURM array, one task per HMM search: hmmsearch, top hits, retrieval, Class A clustering, FAMSA, trimAl, FastTree |
| `mmseqs2_search.sh` | SLURM array for the two seed searches (GPR1_YEAST, FSLB_DICDI) |
| `hmm_header_get.py` | top-ranked accessions of an hmmsearch output |
| `trim_sequences_using_hmm.py` | HMM-matched region of every hit |
| `input_sequences/<search>/` | the HMM or seed FASTA of each search |

```bash
source ../../config/paths.sh          # sets UNIPROT_FASTA, UNIPROT_MMSEQS_DB and the archive
sbatch hmm_search.sh
sbatch mmseqs2_search.sh
```

Software: HMMER 3.4, MMseqs2, FAMSA 2.2.2, trimAl 1.4, FastTree 2.1.11.

**Archive** (`01_homology_search/<search>/`): the hmmsearch output and the retrieved headers, the
retrieved sequences (`*_db.fasta.gz`), the alignments, trimmed alignments, trees and, from stage
02, the score-ordered trees, alignments and score line plots. `classA/` also holds the 80% cluster
table (`classA_hmm_E1000_db_clust0.8.tsv`).

`input_sequences/TM116/TM116.hmm` is Pfam PF30689.1 (TMEM116), copied from the cluster on 22 Oct.
