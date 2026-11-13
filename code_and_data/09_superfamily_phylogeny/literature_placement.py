"""Placement of the 18 literature sequences in the two Oct08 literature runs (core + literature),
over the best tree and the 100 ML search trees of each alignment.

Per literature sequence and tree:
  in clade  the sequence is inside its candidate family's clade, the best-F1 side of any branch for
            that family's own (non-literature) leaves; nested lineages merged as in tree_metrics.py
            (Class A-like = Class A, Olfactory, Vomeronasal 1, Taste-2, Nematode). For the
            Acanthamoeba "Rhodopsin" sequences (LIT_GPR1like) both Class A-like and GPR1 are tested.
  nearest   the families of the non-literature leaves in the smallest clade that holds the sequence
            and at least one non-literature leaf (trees rooted on a Class B1 leaf).
  foreign   for a sequence outside the family clade: the leaves of other families (literature
            leaves not counted) in the smallest clade holding both the family clade and the
            sequence, i.e. how far outside it sits (0 = inside; trees rooted on a Class B1 leaf).

Writes results/literature_placement.tsv.
Usage: python literature_placement_Oct08.py
"""
from collections import Counter
from pathlib import Path

from ete3 import Tree

import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import stage  # noqa: E402
S9 = stage("09_superfamily_phylogeny")
OUT = S9 / "readouts"
FAMMAP = S9 / "representatives_family_map.tsv"
TINDEX = S9 / "booster_stats" / "tindex_mean_across_final_trees.tsv"
RUNS = {a: (S9 / "trees" / f"representatives_and_literature_{a}_kpicgappy.raxml.support",
            S9 / "runs" / "representatives_and_literature"
            / f"representatives_and_literature_{a}_kpicgappy.raxml.mlTrees")
        for a in ("ginsi_ep0.123", "ginsi")}

A_LIKE = {"classA", "Olfactory", "Vomeronasal1", "classT", "Nematode"}
UNIT = {f: "Class A-like" for f in A_LIKE}
# the family clades each literature group is tested against
TARGETS = {"LIT_GP157": {"GP157": {"GP157"}},
           "LIT_classT": {"Taste-2": {"classT"}},
           "LIT_Vomeronasal1": {"Vomeronasal 1": {"Vomeronasal1"}},
           "LIT_GPR1like": {"Class A-like": A_LIKE, "GPR1": {"GPR1"}}}

fam_of = dict(l.split("\t") for l in FAMMAP.read_text().splitlines() if l.strip())
mean_ti = {r.split("\t")[0]: r.split("\t")[3] for r in TINDEX.read_text().splitlines()[1:]}


def acc(name):
    return name.split("|")[1] if "|" in name else name


def placements(newick):
    """-> {literature accession: ({target: inside?}, nearest-families string)} for one tree."""
    t = Tree(newick.strip())
    leaves = t.get_leaves()
    t.set_outgroup(next(l for l in leaves if fam_of[acc(l.name)] == "classB1"))
    idx = {l: i for i, l in enumerate(leaves)}
    fam = [fam_of[acc(l.name)] for l in leaves]
    mask = {}
    for n in t.traverse("postorder"):
        mask[n] = (1 << idx[n]) if n.is_leaf() else sum(mask[c] for c in n.children)
    full = (1 << len(leaves)) - 1
    sides = [s for n in t.traverse() if n is not t for s in (mask[n], full ^ mask[n])]

    clade = {}
    for targets in TARGETS.values():
        for name, fams in targets.items():
            if name in clade:
                continue
            target = sum(1 << i for i, f in enumerate(fam) if f in fams)
            nt = target.bit_count()
            clade[name] = max(sides, key=lambda s: 2 * (s & target).bit_count() / (s.bit_count() + nt))

    lit = sum(1 << i for i, f in enumerate(fam) if f.startswith("LIT_"))
    nodes = sorted(mask.values(), key=int.bit_count)

    def foreign(c, bit):
        """leaves of other families in the smallest clade holding clade c and the sequence"""
        hold = next((m for m in nodes if m & (c | bit) == c | bit), full)
        extra = hold & ~c & ~lit
        return Counter(UNIT.get(fam[i], fam[i]) for i in range(len(fam)) if extra >> i & 1)

    out = {}
    for l in leaves:
        f = fam[idx[l]]
        if not f.startswith("LIT_"):
            continue
        bit = 1 << idx[l]
        inside = {name: (bool(clade[name] & bit), foreign(clade[name], bit)) for name in TARGETS[f]}
        n = l.up
        while all(fam[idx[x]].startswith("LIT_") for x in n.get_leaves()):
            n = n.up
        near = Counter(UNIT.get(fam[idx[x]], fam[idx[x]]) for x in n.get_leaves()
                       if not fam[idx[x]].startswith("LIT_"))
        out[acc(l.name)] = (l.name, inside, " + ".join(sorted(near)))
    return out


def main():
    OUT.mkdir(exist_ok=True)
    rows = []
    for aln, (best, ml) in RUNS.items():
        best_p = placements(best.read_text())
        trees = [x for x in ml.read_text().splitlines() if x.strip()]
        ins, near, far = {}, {}, {}
        for nwk in trees:
            for a, (_, inside, nb) in placements(nwk).items():
                ins.setdefault(a, Counter()).update(k for k, v in inside.items() if v[0])
                for k, v in inside.items():
                    far.setdefault((a, k), []).append(sum(v[1].values()))
                near.setdefault(a, Counter())[nb] += 1
        for a, (name, inside, nb) in sorted(best_p.items(), key=lambda x: x[1][0]):
            fields = name.split("|")
            for target, (v, extra) in inside.items():
                ml_far = sorted(far[(a, target)])
                rows.append([aln, a, fields[0], fields[2], target, "yes" if v else "no",
                             ins[a][target], len(trees),
                             ", ".join(f"{k} {c}" for k, c in extra.most_common()) or "0",
                             ml_far[len(ml_far) // 2], nb,
                             "; ".join(f"{k} ({c})" for k, c in near[a].most_common(3)),
                             mean_ti.get(a, "")])
    with open(OUT / "literature_placement.tsv", "w") as fh:
        fh.write("alignment\taccession\tcandidate\tspecies\tfamily_clade\tin_clade_best_tree\t"
                 "in_clade_ml_trees\tof\tforeign_best_tree\tforeign_ml_median\tnearest_best_tree\tnearest_ml_trees_top3\tmean_tIndex\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    print(f"{len(rows)} rows -> {OUT / 'literature_placement.tsv'}")


if __name__ == "__main__":
    main()
