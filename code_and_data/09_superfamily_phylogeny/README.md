# 09 · Phylogenetic inference of the superfamily hierarchy

Methods: *Phylogenetic Inference of the GPCR Superfamily Hierarchy*; Figure 3; Tables S6 and S7.

For each of the eight trimmed alignments, ModelFinder (IQ-TREE 2.3.6; `-m MF -mfreq FO -cmin 5 -cmax 15
--mset JTT --mrate R`) selects the model, and RAxML-NG 1.2.2 infers the tree from 100 parsimony starting
trees with autoMRE bootstrapping (up to 500 replicates) and transfer bootstrap expectation (TBE)
support (`raxml-ng --all --tree pars{100} --model <best> --bs-trees autoMRE{500} --bs-metric TBE`; the
exact command of every run is in its `.raxml.log`). booster gives the transfer index of every leaf in
each best tree (`superfamily_tree.sh`, one array task per alignment).

| Script | Writes (in `readouts/` unless noted) |
|---|---|
| `tindex_across_trees.py` | `booster_stats/tindex_mean_across_final_trees.tsv`: mean transfer index of every sequence over the trees it is in |
| `label_tree_families.py` | `trees/<tree>.families.nwk` |
| `tree_metrics.py` | per family and run: clade, stem, sister; `tree_unit_metrics.tsv`, grouping recovery |
| `core_tree_readout.py` | TBE, recovery and placement of the family groupings (`core_groupings_tbe_recovery.tsv`) |
| `split_dating.py` | the taxonomic content on either side of the splits (`split_dating.tsv`) |
| `stem_ranks.py` | stem length of every family over the 100 ML searches (`stem_ranks.tsv`) |
| `nonanimal_classA.py`, `internal_order.py` | placement of named and non-animal sequences, order within Class B and Class F |
| `literature_placement.py` | the 18 literature sequences in or outside the clade of their claimed family |
| `sCF_calculate.sh`, `run_scf_gather.sh` (+ `gather_scf.py`, `compare_scf_trees.py`, `discover_groupings.py`) | site concordance factors of every ML search tree (IQ-TREE `--scfl 100`) and every family-level split with its recovery and median sCF (`sCF/<alignment>/`) |
| `tree_runs.py` | Table S6a: sequences, columns, model, log-likelihood, searches and bootstraps of each run |
| `build_tree_support_tables.py` | the other source tables of Tables S6 and S7 (`tree_support_tables/`) |
| `website/plot_circular_trees.py`, `website/build_tree_viewer_data.py`, `website/build_superfamily_trees_data.py` | the data of the interactive trees at gpcrevolution.org/superfamily-trees, through a local `website/build/` folder (not in the archive) |

A grouping is recovered in a search tree when a branch separates its families from the rest, each family
counted on the side holding more than half of its sequences. Sequences whose mean transfer index is at
or above that of the most stable STE2 sequence (0.692) are left out of the interpretation.

```bash
sbatch superfamily_tree.sh         # 8 tasks: ModelFinder, RAxML-NG, booster -> runs/, booster_stats/
python tindex_across_trees.py && python label_tree_families.py
python tree_metrics.py && python core_tree_readout.py && python split_dating.py && python stem_ranks.py
python nonanimal_classA.py && python internal_order.py && python literature_placement.py
sbatch sCF_calculate.sh            # 800 tasks, then: sbatch run_scf_gather.sh
python tree_runs.py && python build_tree_support_tables.py
cd website && python plot_circular_trees.py && python build_tree_viewer_data.py
python build_superfamily_trees_data.py --out <GPCRevolution>/public/superfamily-trees
```

Verified: rerun from the archive, every source table of Tables S6 and S7 is identical to those of the paper,
and the website data is identical to the published page's.

**Archive** (`09_superfamily_phylogeny/`): `runs/<set>/` (RAxML-NG logs, best models, best trees, the 100
ML search trees, TBE support, IQ-TREE ModelFinder reports and booster output), `trees/` (the eight best
trees with TBE and family labels), `booster_stats/`, `readouts/`, `sCF_analysis/` (the 800 per-tree IQ-TREE `.cf.tree` files, `scf_<set>_<alignment>/scf_<i>.cf.tree`;
rebuilt from them, `scf_per_tree_summary.tsv` and the grouping tables are identical to `sCF/`, apart from
the dataset names, which carried `_Oct08`), `sCF/`, `tree_support_tables/` and
`representatives_family_map.tsv`. Each run folder also holds its 500 bootstrap replicate trees (`.raxml.bootstraps.gz`).

**Figure 3** is drawn at gpcrevolution.org/superfamily-trees with Sequence set "Core families",
Alignment "Higher gap penalty" (MAFFT --ep 0.123), View "Circular (Simplified)", Family labels "Final assignment" and
"Fade unstable sequences" checked; "Export SVG" saves the drawing (the figure itself is in `figure3/` in the
archive). The figure names the GPR1 clade "GPR1 / Class A" where the page names it GPR1. All eight trees can
be explored on the same page.
