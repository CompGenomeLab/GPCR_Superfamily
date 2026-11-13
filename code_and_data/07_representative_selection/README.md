# 07 · Representative sequences

Methods: *Representative Protein Sequence Selection*; Table S5.

`select_representatives.py` chooses up to a set number of representatives per family (`TARGETS`: 100
for cAMP and Class A, 50 for Class C, 25 for GP143, GP157, GPR1 and STE3, 20 for the others; 590 in
total) from the family's complete final membership. Selection walks the taxonomic ranks from the highest
rank holding more than one group down to species, round-robin across the groups of each rank, and within
a group takes the sequence with the highest topology score (`mean_jaccard`, stage 03). At most two
representatives come from one species of a family, and no minimum topology score is applied. The three
pooled sets are family filters over the same selection: all families (590), without STE2 (570), and the
core set without STE2, STE3, Class C and Vomeronasal 2 (475).

```bash
python select_representatives.py
```

Verified: rerun from the archive, the three sets are identical to those the trees were built from.

**Archive** (`07_representative_selection/`): `representatives_all_families.tsv` (every representative
with its family, taxonomy and topology score), `representatives_<set>.txt` for the three sets, and the
taxonomy table the ranks are read from (`taxonomy_eukaryotes_filtered.json`).
