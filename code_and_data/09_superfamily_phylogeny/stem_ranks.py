"""Stem length of every family unit over the 100 ML searches of the four Oct08 runs that contain the
divergent families (all and no STE2, ginsi_ep0.123 and ginsi): mean and median stem, and how often
the unit's stem is among the three longest of that tree, per run and pooled over the four runs.
Stem as in tree_metrics.py (the branch above the unit's best-F1 clade; nested lineages merged), with
GP143 and its sister TM116 (100/100 searches) also merged into one unit, since each one's own stem is
mostly the branch between the two. Writes results/stem_ranks.tsv.
Usage: python stem_ranks_Oct08.py
"""
import statistics as st
from collections import defaultdict

import tree_metrics
from tree_metrics import OUT, RUNS, S9, load, per_tree_unit_metrics

tree_metrics.UNIT["GP143"] = tree_metrics.UNIT["TM116"] = "GP143 + TM116"

rows = []
pooled = defaultdict(list)
for label, best, ml in RUNS:
    if not label.startswith(("all", "no STE2")):
        continue
    stems, top3, rank1 = defaultdict(list), defaultdict(int), defaultdict(int)
    lines = [l for l in (S9 / ml).read_text().splitlines() if l.strip()]
    for nwk in lines:
        t, leaves, fam, mask = load(nwk)
        m = per_tree_unit_metrics(t, leaves, fam, mask, False)
        order = sorted(m, key=lambda u: -m[u]["stem"])
        for i, u in enumerate(order):
            stems[u].append(m[u]["stem"])
            pooled[u].append(m[u]["stem"])
            top3[u] += i < 3
            rank1[u] += i == 0
    for u in sorted(stems, key=lambda u: -st.mean(stems[u])):
        rows.append([label, u, f"{st.mean(stems[u]):.3f}", f"{st.median(stems[u]):.3f}",
                     top3[u], rank1[u], len(lines)])
    print("done", label)
for u in sorted(pooled, key=lambda u: -st.mean(pooled[u])):
    rows.append(["pooled, 4 runs", u, f"{st.mean(pooled[u]):.3f}", f"{st.median(pooled[u]):.3f}",
                 "", "", len(pooled[u])])
with open(OUT / "stem_ranks.tsv", "w") as fh:
    fh.write("run\tunit\tmean_stem\tmedian_stem\tin_top3\tlongest\tof\n")
    fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)
