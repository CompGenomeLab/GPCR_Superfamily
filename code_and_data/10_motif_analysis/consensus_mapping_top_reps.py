"""Consensus column mapping (consensus_column_mapping.py) voted by each family's top N
representatives only, ranked by the rank-sum topology score that chose the anchors
(bridge_representatives_Oct08_topo0.8_ranksum.tsv, rank_in_family; rank 1 = acc1), and its
comparison with the all-representative mapping.

python consensus_mapping_top_reps.py [N] [threshold]   (defaults 5, 0.33)

Writes consensus_mapping/<alignment>_conf<threshold>_top<N>.tsv for the 8 alignments and
prints, per family on the main alignment, how many superfamily columns keep the same family
column, change it, or are mapped by only one of the two.
"""
import csv
import sys

from consensus_column_mapping import ALIGNMENTS, DATA, HERE, OUT, anchors, load_families, map_alignment

N = int(sys.argv[1]) if len(sys.argv) > 1 else 5
THRESHOLD = float(sys.argv[2]) if len(sys.argv) > 2 else 0.33
MAIN = "representatives_all_Oct08_ginsi_ep0.123"


def top_reps(n):
    rows = csv.DictReader(open(DATA / "anchors" / "bridge_representatives_ranksum.tsv"), delimiter="\t")
    return {(r["family"], r["accession"]) for r in rows if int(r["rank_in_family"]) <= n}


def read_map(path):
    out = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["consensus_family_column"]:
            out[(r["family"], int(r["sup_column"]))] = int(r["consensus_family_column"])
    return out


def main():
    keep = top_reps(N)
    fams = {f: {a: v for a, v in reps.items() if (f, a) in keep} for f, reps in load_families().items()}
    acc1 = anchors()
    for name, mapping_dir in ALIGNMENTS.items():
        map_alignment(name, mapping_dir, fams, acc1, THRESHOLD, tag=f"_top{N}")
    every = read_map(OUT / f"{MAIN}_conf{THRESHOLD}.tsv")
    top = read_map(OUT / f"{MAIN}_conf{THRESHOLD}_top{N}.tsv")
    print(f"{MAIN}: top {N} vs all representatives (threshold {THRESHOLD})")
    print("  family        reps  same  changed  only_all  only_top")
    for f in sorted({k[0] for k in every} | {k[0] for k in top}):
        a = {c: v for (g, c), v in every.items() if g == f}
        t = {c: v for (g, c), v in top.items() if g == f}
        both = a.keys() & t.keys()
        same = sum(a[c] == t[c] for c in both)
        print(f"  {f:13} {len(fams.get(f, {})):4}  {same:4}  {len(both) - same:7}  {len(a.keys() - t.keys()):8}  "
              f"{len(t.keys() - a.keys()):8}")


if __name__ == "__main__":
    main()
