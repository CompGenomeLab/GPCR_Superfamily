"""Build the taxonomy tree and the presence/absence table the family origins are inferred on.

Stage 06. The tree is the NCBI taxonomy of every proteome in the searched database (UniProt
2025_02), built by build_taxonomy_trees.py, with Unikonta, Archaeplastida, Holozoa and
Holomycota added. The table has one row per tree leaf (genome) and one column per family and
control, 1 when the genome carries at least one final member. A searched genome with no row in
the count table is written as absent everywhere, so a genome in which nothing was found counts
as evidence of absence.

Reads  06_family_origins/taxonomy/taxonomy_ncbi_searched.rankdepth.nwk
       04_overlap_assignment/tax_counts_per_family.tsv
Writes 06_family_origins/taxonomy_ncbi_all.nwk
       06_family_origins/family_presence_absence_ncbi_all.csv

Usage:
    python build_ancestral_inputs.py
"""

import csv
import os

from ete3 import Tree

from key_nodes import FAMILIES, OUT_DIR, TAXONOMY_DIR
from gpcr_paths import stage

SOURCE_TREE = os.path.join(TAXONOMY_DIR, "taxonomy_ncbi_searched.rankdepth.nwk")
COUNT_TABLE = str(stage("04_overlap_assignment", "tax_counts_per_family.tsv"))


def read_counts():
    with open(COUNT_TABLE, newline="", encoding="utf-8") as handle:
        return {row["taxID"]: {f: int(row[f + "_count"]) for f in FAMILIES}
                for row in csv.DictReader(handle, delimiter="\t")}


def main():
    counts = read_counts()
    tree = Tree(SOURCE_TREE, format=1)
    leaf_order = [leaf.name for leaf in tree.get_leaves()]
    tree.write(format=1, outfile=os.path.join(OUT_DIR, "taxonomy_ncbi_all.nwk"))

    absent = 0
    table = os.path.join(OUT_DIR, "family_presence_absence_ncbi_all.csv")
    with open(table, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["taxID"] + FAMILIES)
        for taxid in leaf_order:
            row = counts.get(taxid)
            if row is None:
                absent += 1
                row = {}
            writer.writerow([taxid] + [1 if row.get(f, 0) > 0 else 0 for f in FAMILIES])
    print("%d genomes, %d with no counts row (written as absent) -> %s"
          % (len(leaf_order), absent, os.path.basename(table)))


if __name__ == "__main__":
    main()
