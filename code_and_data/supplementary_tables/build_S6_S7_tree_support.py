"""Supplementary Tables S5 (tree support) and S6 (per-sequence stability), copied from step 9, which builds them
(step9_ml_hierarchy/supplementary_Oct08/build_supplementary_tables.py, on the Oct08 trees).

  S5a  the eight tree runs: sequences, columns, model, log-likelihood
  S5b  family stability: tIndex, clade F1, stem and sister in the all-family trees
  S5c  stems of each family across the 100 ML searches of each run
  S5d  groupings: TBE in each best tree, recovery in its 100 searches, median sCF
  S5e  named placements (non-animal members, Class B / F order, split dating)
  S5f  literature sequences: placement in the family clade
  S5g  every family-level split: recovery and median sCF, main alignment (G-INS-i --ep 0.123)
  S5h  same, plain G-INS-i
  S6   every representative: tIndex in each tree, mean, and the STE2-based instability flag

    python build_S5_S6_tree_support.py
"""
import shutil
from pathlib import Path

import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import stage  # noqa: E402
CODE = Path(__file__).resolve().parent
HERE = stage("supplementary_tables")
HERE.mkdir(parents=True, exist_ok=True)
SRC = stage("09_superfamily_phylogeny", "tree_support_tables")
FILES = {"Table_S6a_tree_runs.tsv": "S_tree_runs.tsv",
         "Table_S6b_family_stability.tsv": "S_family_stability.tsv",
         "Table_S6c_family_stems.tsv": "S_family_stems.tsv",
         "Table_S6d_groupings.tsv": "S_groupings.tsv",
         "Table_S6e_placements.tsv": "S_placements.tsv",
         "Table_S6f_literature_placement.tsv": "S_literature.tsv",
         "Table_S6g_family_splits_sCF_ep0.123.tsv": "S_family_splits_sCF_ginsi_ep0.123.tsv",
         "Table_S6h_family_splits_sCF_ginsi.tsv": "S_family_splits_sCF_ginsi.tsv",
         "Table_S7_sequence_tindex.tsv": "S_sequence_tindex.tsv"}
for out, src in FILES.items():
    shutil.copyfile(SRC / src, HERE / out)
    print("copied", src, "->", out)
