"""Dollo parsimony origin of each family on the NCBI taxonomy tree of all searched proteomes.

Stage 06. Under Dollo parsimony (Farris 1977; as in Count, Csuros 2010) a family is gained
once and may be lost any number of times but never regained. On a fixed rooted tree the
reconstruction with the fewest losses places the single gain on the branch above the last
common ancestor (LCA) of all genomes carrying the family. Inside the LCA's clade a node is
present when any genome below it carries the family, and every branch from a present node
into a clade with no carrier is one loss.

Before the reconstruction, stray hits are set aside with a genus threshold k. At the
current LCA, a child lineage is set aside when its carriers span fewer than k genera, while
other genera of that lineage were sampled without the family, and another child lineage
carries it in k genera or more. The LCA is then recomputed, until nothing more is set aside.
A single-genus lineage that carries the family in its only sampled genus (Capsaspora, for
instance) is never set aside. k = 1 would set nothing aside (plain Dollo parsimony), which
lets a single stray genome decide the origin, so the sweep starts at 2.

The paper uses k = 2. Writes 06_family_origins/dollo_origins.tsv with one row per family: the
origin node, the key node it falls in (the most recent key node containing it), a label for
plotting (the key node, or the origin node inside it when the two differ), the number of
losses, and the genera whose setting aside moved the origin there.

Usage:
    python dollo_origins.py
"""

import csv
import json
import os

from ete3 import Tree

from key_nodes import FAMILIES, KEY_NODES, OUT_DIR, TAXONOMY_DIR

THRESHOLDS = (2,)


def genus_of_taxid():
    with open(os.path.join(TAXONOMY_DIR, "taxonomy_ncbi_searched.json"), encoding="utf-8") as handle:
        records = json.load(handle)
    genus = {}
    for r in records:
        g = r.get("Genus")
        genus[str(r["taxID"])] = g if g and g != "NA" else "taxid:%s" % r["taxID"]
    return genus


def lca(tree, leaves):
    return leaves[0] if len(leaves) == 1 else tree.get_common_ancestor(leaves)


def dollo(tree, carriers, genus, k):
    """Origin node, losses and set-aside leaf names for one family."""
    kept = set(carriers)
    set_aside = []
    while True:
        origin = lca(tree, [tree.leaf_by_name[n] for n in kept])
        if origin.is_leaf():
            break
        lineages = []
        for child in origin.children:
            names = child.get_leaf_names()
            hit = [n for n in names if n in kept]
            if hit:
                lineages.append((child, hit, {genus[n] for n in hit}, {genus[n] for n in names}))
        if not any(len(g) >= k for _, _, g, _ in lineages):
            break
        stray = [hit for _, hit, g, sampled in lineages if len(g) < k and len(sampled) > len(g)]
        if not stray:
            break
        for hit in stray:
            kept.difference_update(hit)
            set_aside.extend(hit)
    # Losses: branches from a present node into a clade without a kept carrier.
    present = {}
    for node in origin.traverse("postorder"):
        present[node] = (node.name in kept) if node.is_leaf() else any(present[c] for c in node.children)
    losses = sum(1 for node in present if present[node] and not node.is_leaf()
                 for child in node.children if not present[child])
    return origin, losses, set_aside


def main():
    genus = genus_of_taxid()
    out = os.path.join(OUT_DIR, "dollo_origins.tsv")
    with open(out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["family", "taxonomy", "k_from", "k_to", "origin_node", "key_node",
                         "origin_label", "carriers", "losses", "newly_set_aside"])
        for taxonomy in ("ncbi",):
            tree = Tree(os.path.join(OUT_DIR, "taxonomy_ncbi_all.nwk"), format=1)
            tree.name = tree.name or "root"  # the tree file leaves the root unnamed
            tree.leaf_by_name = {leaf.name: leaf for leaf in tree}
            key_by_name = {ncbi: label for label, ncbi in KEY_NODES}
            with open(os.path.join(OUT_DIR, "family_presence_absence_ncbi_all.csv")) as table:
                rows = list(csv.DictReader(table))
            for family in FAMILIES:
                carriers = [r["taxID"] for r in rows if r[family] == "1"]
                steps = []
                for k in THRESHOLDS:
                    origin, losses, set_aside = dollo(tree, carriers, genus, k)
                    if steps and steps[-1][1] is origin:
                        steps[-1][0][1] = k
                        continue
                    node = origin
                    while node.name not in key_by_name:
                        node = node.up
                    label = key_by_name[node.name]
                    if len(node) != len(origin):
                        label = "%s (in %s)" % (
                            origin.name.split("__")[-1].replace("%20", " "), label)
                    before = steps[-1][2] if steps else set()
                    aside = {genus[n] for n in set_aside}
                    steps.append(([k, k], origin, aside, key_by_name[node.name], label,
                                  len(carriers) - len(set_aside), losses, sorted(aside - before)))
                for (k_from, k_to), origin, _, key, label, kept, losses, new in steps:
                    writer.writerow([family, taxonomy, k_from, k_to, origin.name, key, label, kept,
                                     losses, ",".join(new)])
    print("wrote %s" % out)


if __name__ == "__main__":
    main()
