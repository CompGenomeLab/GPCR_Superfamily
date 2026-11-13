# 06 · Taxonomic distribution and family origins

Methods: *Taxonomic Distribution Visualization*, *Inference of Family Origins*; Figure 2; Table S4.

`build_taxonomy_trees.py` builds the NCBI taxonomy of all 2,760 proteomes of the searched database
(`ncbi_taxonomy_create.py`, `generate-ncbi-taxonomy-newick.mjs`), adding Unikonta, Archaeplastida,
Holozoa and Holomycota, which NCBI does not have. `build_ancestral_inputs.py` scores each genome as
carrying or lacking every family and control. `dollo_origins.py` places the origin of each family by
Dollo parsimony at the last common ancestor of the carrying genomes, after setting aside, at that
ancestor, any lineage that carries the family in fewer than two genera while other sampled genera of
the lineage lack it and another lineage carries it in at least two genera (repeated until nothing more
is set aside; a lineage carrying the family in its only sampled genus is never set aside).
`build_lineage_origin_figure.py` draws Figure 2: the share of genomes carrying each family in each
lineage, the presences judged sporadic (`ORANGE`) or due to contamination or horizontal transfer
(`RED`), the Dollo origin and our interpretation (the common ancestor of the remaining presences).

```bash
python build_taxonomy_trees.py --nodes .../nodes.dmp --names .../names.dmp   # needs node; optional
python build_ancestral_inputs.py
python dollo_origins.py
python build_lineage_origin_figure.py
```

Verified: rerun from the archive, the presence table, the tree, every Dollo origin and the Figure 2 page
are identical to those of the paper.

**Archive** (`06_family_origins/`): `taxonomy/` (searched taxIDs, NCBI taxonomy JSON and trees, the
NCBI taxdump of 7 September 2025), `taxonomy_ncbi_all.nwk`, `family_presence_absence_ncbi_all.csv`,
`dollo_origins.tsv` and `lineage_origin_figure.html`. The per-proteome distribution (Fig. S1) is the
gpcrevolution.org/taxonomic-distribution page.
