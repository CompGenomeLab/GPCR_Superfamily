"""Prefix each leaf of the 8 best ML trees with its GPCR family from the representative table.

Leaf "tr|F0ZVL1|F0ZVL1_DICPU/3-277" becomes "FSLB|tr|F0ZVL1|F0ZVL1_DICPU/3-277"; a literature
sequence keeps its candidate family ("LIT_GP157|..."). Writes <tree>.families.nwk next to each
best tree in 09_superfamily_phylogeny/trees/.

Usage:
    python label_tree_families.py
"""
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import stage  # noqa: E402

TABLE = stage("07_representative_selection", "representatives_all_families.tsv")
TREE_DIR = stage("09_superfamily_phylogeny", "trees")
LEAF = re.compile(r"(?<=[(,])([^(),:;]+)(?=[:,)])")

with open(TABLE, newline="") as fh:
    FAMILY_OF = {row["accession"]: row["gpcrfamily"] for row in csv.DictReader(fh, delimiter="	")}
trees = sorted(p for p in TREE_DIR.glob("*.raxml.support"))
if not trees:
    raise SystemExit(f"No best trees in {TREE_DIR}")

for tree_path in trees:
    family_of = FAMILY_OF
    text = tree_path.read_text()
    missing = []

    def relabel(m):
        name = m.group(1)
        # A literature sequence ("LIT_GP157|LIT0024|Arcella_intermedia|claimed=...") is labelled
        # with its candidate family, the family the study claimed, kept apart from our families.
        if name.startswith("LIT_"):
            return f"{name.split('|')[0]}|{name}"
        acc = name.split("|")[1] if "|" in name else name.split("/")[0]
        fam = family_of.get(acc)
        if fam is None:
            missing.append(name)
            return name
        return f"{fam}|{name}"

    out = LEAF.sub(relabel, text)
    out_path = tree_path.with_suffix(".families.nwk")
    out_path.write_text(out)
    n = len(LEAF.findall(text))
    print(f"{tree_path.name}: {n} leaves, {len(missing)} unmatched "
          f"-> {out_path.name}")
    for name in missing:
        print("  unmatched:", name)
