# 05 · Comparison with published classifications

Methods: *Comparison with Published Classifications*; Tables S3a-d. The full description of every
rule, written for a reviewer who wants to check it, is in `README_LITERATURE_VERIFICATION.md`.

Receptors reported in 15 published GPCR surveys (45,590 receptors) are read from the studies'
supplementary files (one reader per study in `literature_readers.py`). The files themselves are not
distributed (copyright); what each reader takes from them is, as one table per study in `literature/parsed/`,
which the comparison reads (`write_parsed_literature.py` rebuilds them from downloaded files; see
`literature/SOURCES.md`). Claims stated only in a paper's text are in `literature/curated_claims.tsv`, and
hand-resolved genome names in `literature/genome_name_decisions.tsv`. `compare_literature.py` matches every receptor to our
data in the species the study names, by identifier or by sequence (identical, contained, or
near-identical with identical residues covering at least 90% of the longer sequence; studies' genomes
outside our set are compared in a genome of the same genus when most of their receptors match it), and
follows each matched protein through our stages. Receptors not matched this way were searched with
MMseqs2 against the searched proteomes (`literature_mmseqs2.sh`; query sets written by `write_unmatched_fasta.py`,
`build_remaining_fasta.py` and `write_unsearched_fasta.py`; sequences listed only by identifier fetched
by `fetch_idonly.py` and `fetch_idonly_round2.py`; their small curated inputs are in `query_inputs/`);
`classify_all_literature.py` takes a same-species or
same-genus hit at 90% identity as a match and otherwise the best hit as an indicator only.
`build_supplementary_table.py` writes the receptor table (S3a) and `build_agreement_tables.py` the
agreement per study and per family (S3b, S3c).

`literature_deeper_origins.py` lists literature presences that would place an origin deeper than ours
(Table S3d). `write_deeper_origin_fasta.py` writes the sequences it tests, each once, as
`deeper_origin_tested_sequences.fasta` (365), the input of DeepTMHMM and of an MMseqs2 search;
`deeper_origin_best_hits.py` gives each its closest family, the family of its best MMseqs2 hit (E <= 1e-5,
query coverage >= 0.5) among our family members (`literature_closest_family.tsv`; 346 with such a hit).
`literature_7tm_gate.py` applies the 7TM pattern of stage 03 to their DeepTMHMM
predictions, and `select_literature_representatives.py` writes the 18 sequences added to the trees
(up to two per lineage outside the family's origin, from genomes outside our dataset; the choice is
recorded in the script) trimmed to their 7TM core.

```bash
python compare_literature.py
sbatch literature_mmseqs2.sh                 # MMseqs2 search of the query sets -> *_literature_sequences.m8
python classify_all_literature.py
python build_supplementary_table.py
python build_agreement_tables.py
python verify_literature_comparison.py      # checks against the supplements (PASS/FAIL); needs them downloaded
python literature_deeper_origins.py         # needs stage 06
python write_deeper_origin_fasta.py         # deeper_origin_tested_sequences.fasta
sbatch literature_mmseqs2.sh                 # fourth task: MMseqs2 search of the tested sequences
python deeper_origin_best_hits.py           # literature_closest_family.tsv
sbatch deeptmhmm_literature.sh              # DeepTMHMM on the tested sequences -> predicted_topologies.3line
python literature_7tm_gate.py
python select_literature_representatives.py
```

Verified: rerun from the archive, S3b and S3c are identical to the published tables and S3a differs only
where the earlier version used data the paper does not describe (the SLT_TRNS background family and two
extra seed searches); no receptor changes its final family.

**Archive** (`05_literature_comparison/`): the comparison tables, the receptor table and the agreement
tables, the three MMseqs2 query sets with their hits (`.fasta`, `.tsv`, `.m8`), the best-hit tables used
by the agreement tables, the fetched identifier-only sequences, the 365 sequences tested for deeper origins with their MMseqs2
hits (`deeper_origin_tested_sequences.fasta`, `.m8`), the translated TAS2R loci, the DeepTMHMM
predictions of the literature sequences, the 7TM gate (`deeper_origin_7TM_gate.tsv`) and the 18
literature sequences of the trees (`literature_reps_7TM_flank5.fasta`).

Verified: rerun from the archive, `write_deeper_origin_fasta.py` gives the same 365 sequences in the same
order; their headers differ in labels only (12 read "Liu et al. 2026" for "J Fungi 2026", and the
best_hit=<family> label, from an earlier local-alignment comparison, is no longer written). The archive
keeps the file that was used. The 7TM gate and the 18 literature sequences written from the archive are identical to those
used in the trees.
