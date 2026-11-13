# 03 · Topology filters

Methods: *Controls and Filtering* (first three paragraphs); Figure 1c-d; Table S1.

**Single GPCR-like 7TM core.** The topology of every clade member is predicted with DeepTMHMM
(`deeptmhmm_predict.sh`; the largest families split with `split_fasta_shards.py` and run with
`deeptmhmm_predict_shards.sh`). `build_deeptmhmm_tables.py` merges the per-family DeepTMHMM outputs and writes every prediction with
its helices split by direction (W, outside to inside; M, inside to outside). `write_core_membership.py`
collapses consecutive identical states and keeps a sequence when the pattern `OWIMOWIMOWIMOWI`
occurs exactly once.

**Topology score.** The 7TM cores (five flanking residues) of each family's single-core sequences are
aligned with FAMSA (`family_topology_alignment.sh`), every residue is
replaced by its topology state (`project_topology_onto_alignment.py`) and columns with 90% or more gaps
are removed. Over the columns occupied in at least 30% of the sequences, runs of at least eight
consecutive helix-consensus columns define the consensus helices. A sequence's score for one helix is
the Jaccard similarity between the consensus helix columns and the columns where the sequence has the
same helix state; its topology score is the mean over the seven helices (`per_helix_features.py`,
`score_all_family_topologies.py`, column `mean_jaccard`). Sequences flagged below 0.8 were inspected in
the topology viewer (`build_all_topology_viewers.py`), and `apply_topology_cutoff.py` removes those
scoring below 0.8. The nested-family rule of stage 02 is applied again to the passing sequences.

```bash
sbatch deeptmhmm_predict.sh                     # optional, see below (large families: deeptmhmm_predict_shards.sh)
python build_deeptmhmm_tables.py --root <DeepTMHMM output, one folder per family>   # optional
python write_core_membership.py
sbatch family_topology_alignment.sh             # FAMSA, topology projection, -gt 0.1 -> topology_alignments/
python score_all_family_topologies.py
python apply_topology_cutoff.py
python build_all_topology_viewers.py            # topology_viewer.html (inspection; website page)
```

The archive holds the merged DeepTMHMM predictions, so `build_deeptmhmm_tables.py` is only needed when
DeepTMHMM is rerun. DeepTMHMM is used under its academic license and is not included; set `DEEPTMHMM_DIR`.

Verified: from the archive, `write_core_membership.py` reproduces the single-core membership of the
paper (408,489 sequences) and `apply_topology_cutoff.py` the topology-filtered membership (391,856)
exactly, family by family.

**Archive** (`03_topology_filter/`):

| File | Contents |
|---|---|
| `deeptmhmm_predictions_oriented.3line` | DeepTMHMM prediction of every clade member (sequence, oriented states) |
| `core_membership.tsv`, `core_sequences/` | the single-core sequences per family |
| `topology_alignments/` | the topology-converted, gap-trimmed alignment of each family |
| `topology_helix_scores.tsv` | per sequence: mean Jaccard (topology score), per-helix Jaccards |
| `topology_membership.tsv` | membership after the topology filter (input of stage 04) |
| `topology_viewer.html` | the inspection viewer (also a website page) |
