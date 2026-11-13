#!/usr/bin/env bash
# Family-column mapping and logo data for the two tree-guided superfamily alignments (Figure 4,
# Table S7): every superfamily column is mapped to a family-alignment column by the votes of the
# family's top 5 representatives (rank-sum topology score, anchors/bridge_representatives_ranksum.tsv),
# assigned when at least a third agree, one-to-one; the logo JSONs are rebuilt on that mapping and
# labelled with the GPCRdb number voted by each family's human receptors.
# Run from this folder: bash run_treein_logos.sh
set -euo pipefail
cd "$(dirname "$0")"

python - <<'PY'
from consensus_column_mapping import TREEIN, anchors, load_families, map_alignment
from consensus_mapping_top_reps import top_reps
keep = top_reps(5)
fams = {f: {a: v for a, v in reps.items() if (f, a) in keep} for f, reps in load_families().items()}
for name, d in TREEIN.items():
    map_alignment(name, d, fams, anchors(), 0.33, tag="_top5")
PY
python consensus_logo_mappings.py --treein --top5 0.33

python - <<'PY'
from consensus_column_mapping import TREEIN
from label_by_human_receptors import main as label
for name in TREEIN:
    label(name, "_top5")
PY
python - <<'PY'
from apply_voted_labels import main as apply_labels
from consensus_column_mapping import DATA, TREEIN
for name, d in TREEIN.items():
    apply_labels(DATA / "consensus_mapping" / "logo_conf0.33_top5" / d,
                 DATA / "gpcrdb_labels" / f"{name}_top5.tsv")
PY
