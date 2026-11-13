"""Mean booster transfer index (tIndex) of each sequence across the 8 final trees.

A sequence's tIndex is read from booster_stats/<tree>.tbe_stats.txt of every final tree it is a
leaf of, and averaged over those trees. Sequences are then ranked from the highest mean (least
stable placement) down.

Writes booster_stats/tindex_mean_across_final_trees.tsv.

Usage:
    python tindex_across_trees.py
"""
import csv
from pathlib import Path

import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import stage  # noqa: E402
S9 = stage("09_superfamily_phylogeny")
STATS = S9 / "booster_stats"
FAMILY_MAP = S9 / "representatives_family_map.tsv"

# the 8 final trees (README.md), with a short column name for each
TREES = {f"{tag}_{short}": f"representatives_{s}_{a}_kpicgappy"
         for tag, s in (("all", "all"), ("noSTE2", "no_ste2"),
                        ("core", "no_ste2ste3classCV2"), ("coreLit", "and_literature"))
         for short, a in (("ep0.123", "ginsi_ep0.123"), ("ginsi", "ginsi"))}


def tindex(path):
    out, on = {}, False
    for line in path.read_text().splitlines():
        if line.startswith("Taxon\ttIndex"):
            on = True
        elif on and line.startswith("Edge\t"):
            break
        elif on and line.strip():
            name, val = line.rsplit("\t", 1)
            out[name] = float(val)
    return out


def accession(taxon):
    """'tr|A0A123|A0A123_HUMAN/3-277' -> A0A123; 'LIT_GP157|LIT0024|...' -> LIT0024."""
    fields = taxon.split("|")
    return fields[1] if len(fields) >= 3 else taxon.split("/")[0]


family = dict(line.rstrip("\n").split("\t") for line in FAMILY_MAP.open() if line.strip())
values, label = {}, {}
for col, stem in TREES.items():
    for taxon, v in tindex(STATS / f"{stem}.tbe_stats.txt").items():
        acc = accession(taxon)
        values.setdefault(acc, {})[col] = v
        label.setdefault(acc, taxon)

rows = []
for acc, per in values.items():
    mean = sum(per.values()) / len(per)
    rows.append([acc, family.get(acc, ""), len(per), round(mean, 4), round(max(per.values()), 4),
                 label[acc]] + [per.get(c, "") for c in TREES])
rows.sort(key=lambda r: -r[3])

out = STATS / "tindex_mean_across_final_trees.tsv"
with out.open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["accession", "family", "n_trees", "mean_tIndex", "max_tIndex", "leaf"] + list(TREES))
    w.writerows(rows)
print(f"{len(rows)} sequences -> {out.name}")
