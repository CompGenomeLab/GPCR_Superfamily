"""Side of the GPR1 | Class A-like split for every non-animal Class A representative.

The reassignment was made on the trees of the September 22 representative set (8 best trees: all,
no STE2 and core sets, and the core + literature set, each with MAFFT G-INS-i --ep 0.123 and with the
default penalty; 100 ML search trees for the five runs whose search trees were kept). These trees are in
the archive under 10_motif_analysis/phylogeny_reassignment/inputs_Sep22/. The trees only define the two
reference groups of map_fungal_classA_to_clades.py.

In each tree the GPR1 + Class A-like clade is the smallest split holding the majority of every family
involved, and its GPR1 side is the split within it that best matches the GPR1 leaves (F1). Each
non-animal Class A leaf is counted on the GPR1 side, the Class A side, or outside the clade.

Reads  10_motif_analysis/phylogeny_reassignment/inputs_Sep22/
Writes 10_motif_analysis/phylogeny_reassignment/tree_side_per_sequence.tsv

Usage:
    python tree_side_per_sequence.py      (needs ete3)
"""
import sys
from collections import defaultdict
from pathlib import Path

from ete3 import Tree

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from gpcr_paths import open_text, stage  # noqa: E402

IN = stage("10_motif_analysis", "phylogeny_reassignment", "inputs_Sep22")
OUT = stage("10_motif_analysis", "phylogeny_reassignment")
A_LIKE = {"classA", "Olfactory", "Vomeronasal1", "classT", "Nematode"}
# (best tree, its 100 ML search trees or None)
RUNS = [(f"representatives_{s}_{a}_kpicgappy", f"representatives_{s}_{a}_kpicgappy" if a.endswith("ginsi")
         or s.startswith("and_literature") else None)
        for s in ("all_Sep22_topo0.8", "no_ste2_Sep22_topo0.8", "no_ste2ste3classCV2_Sep22_topo0.8",
                  "and_literature_Oct4")
        for a in ("ginsi_ep0.123", "ginsi")]


def acc(name):
    return name.split("|")[1] if "|" in name else name


fam_of = {}
with open_text(IN / "representatives_Sep22_topo0.8_literature_family_map.tsv") as fh:
    for line in fh:
        if line.strip():
            a, f = line.rstrip("\n").split("\t")[:2]
            fam_of[a] = None if f.startswith("LIT_") else f
rows = {}
with open_text(IN / "representatives_all_families_Sep22_topo0.8.tsv") as fh:
    for line in list(fh)[1:]:
        f = line.rstrip("\n").split("\t")
        rows[f[1]] = f
kingdom = {a: r[5] for a, r in rows.items()}


def load(newick):
    t = Tree(newick.strip())
    leaves = t.get_leaves()
    idx = {l: i for i, l in enumerate(leaves)}
    fam = []
    for l in leaves:        # literature leaves keep their LIT_ name
        f = fam_of[acc(l.name)]
        fam.append(l.name.split("|")[0] if f is None and l.name.startswith("LIT_") else f)
    mask = {}
    for n in t.traverse("postorder"):
        mask[n] = (1 << idx[n]) if n.is_leaf() else sum(mask[c] for c in n.children)
    return t, leaves, fam, mask


def masks(fam, wanted):
    return sum(1 << i for i, f in enumerate(fam) if f in wanted)


def maj(side, m):
    return 2 * (side & m).bit_count() > m.bit_count()


cnt = defaultdict(lambda: [0, 0, 0, 0, 0, 0])  # best trees: GPR1, Class A, outside; ML searches: same


def side(nwk, k):
    t, leaves, fam, mask = load(nwk)
    full = mask[t]
    sides = []
    for n in t.traverse():
        if n is not t:
            sides += [mask[n], full ^ mask[n]]
    fx = [masks(fam, {"GPR1"})]
    fy = [masks(fam, {f}) for f in A_LIKE if masks(fam, {f})]
    S = min((s for s in sides if all(maj(s, m) for m in fx + fy)), key=lambda s: s.bit_count())
    best, bf = None, -1
    for s in sides:
        if s & ~S:
            continue
        f1 = 2 * (s & fx[0]).bit_count() / (s.bit_count() + fx[0].bit_count())
        if f1 > bf:
            best, bf = s, f1
    for i, (l, f) in enumerate(zip(leaves, fam)):
        a = acc(l.name)
        if f != "classA" or kingdom.get(a) == "Metazoa":
            continue
        b = 1 << i
        cnt[a][k + (0 if b & best else (1 if b & S else 2))] += 1


for best, ml in RUNS:
    with open_text(IN / "best_trees" / f"{best}.raxml.support") as fh:
        side(fh.read(), 0)
    if ml:
        with open_text(IN / "ml_trees" / f"{ml}.raxml.mlTrees") as fh:
            for line in fh:
                if line.strip():
                    side(line, 3)

with open(OUT / "tree_side_per_sequence.tsv", "w", encoding="utf-8", newline="\n") as out:
    out.write("acc\tkingdom\tphylum\tclass\tspecies\tbest_GPR1\tbest_ClassA\tbest_out\tml_GPR1\tml_ClassA\tml_out\n")
    for a, c in sorted(cnt.items(), key=lambda x: (rows[x[0]][5], rows[x[0]][6], -x[1][0])):
        r = rows[a]
        out.write("\t".join([a, r[5], r[6], r[9], r[19]] + list(map(str, c))) + "\n")
print("wrote", OUT / "tree_side_per_sequence.tsv")
