"""Per non-animal Class A representative (early-diverging fungi, Dikarya, other eukaryotes, plants):
in each tree, is it inside the animal Class A-like clade, inside the GPR1 clade, or in neither?

Animal Class A-like clade = best-F1 side for the animal leaves of Class A, Olfactory, Vomeronasal 1,
Taste-2 and Nematode (the sequence itself added to the target, so it is not penalised for being
tested). GPR1 clade = best-F1 side for the GPR1-labelled leaves (plus the sequence). Over the best
tree and the 100 ML searches of every Oct08 run; mean tIndex and the STE2 cutoff flag alongside.
Writes results/nonanimal_classA.tsv.
Usage: python nonanimal_classA_Oct08.py
"""
import csv
from collections import Counter
from pathlib import Path

from ete3 import Tree

import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import stage  # noqa: E402
S9 = stage("09_superfamily_phylogeny")
OUT = S9 / "readouts"
TAXA = stage("07_representative_selection", "representatives_all_families.tsv")
FAMMAP = S9 / "representatives_family_map.tsv"
TINDEX = S9 / "booster_stats" / "tindex_mean_across_final_trees.tsv"
SETS = [("all", "all"), ("no STE2", "no_ste2"), ("core", "no_ste2ste3classCV2"), ("core+lit", "and_literature")]
ALNS = [("ep0.123", "ginsi_ep0.123"), ("ginsi", "ginsi")]
A_LIKE = {"classA", "Olfactory", "Vomeronasal1", "classT", "Nematode"}
CUTOFF = 0.6919
DIKARYA = {"Ascomycota", "Basidiomycota"}

tax = {r["accession"]: r for r in csv.DictReader(TAXA.open(), delimiter="\t")}
fam_of = dict(l.split("\t") for l in FAMMAP.read_text().splitlines() if l.strip())
mean_ti = {r["accession"]: float(r["mean_tIndex"]) for r in csv.DictReader(TINDEX.open(), delimiter="\t")}
acc = lambda name: name.split("|")[1]


def lineage(r):
    if r["kingdom"] == "Fungi":
        return "Dikarya" if r["phylum"] in DIKARYA else "early-diverging fungi"
    return "plants" if r["kingdom"] == "Viridiplantae" else "other eukaryotes"


TESTED = {a: lineage(r) for a, r in tax.items() if r["gpcrfamily"] == "classA" and r["kingdom"] != "Metazoa"}


def masks(t):
    L = t.get_leaves()
    idx = {l: i for i, l in enumerate(L)}
    m = {}
    for n in t.traverse("postorder"):
        m[n] = (1 << idx[n]) if n.is_leaf() else sum(m[c] for c in n.children)
    full = (1 << len(L)) - 1
    sides = [s for n in t.traverse() if n is not t for s in (m[n], full ^ m[n])]
    return L, idx, sides


def best(sides, target):
    nt = target.bit_count()
    return max(sides, key=lambda s: 2 * (s & target).bit_count() / (s.bit_count() + nt))


def placements(nwk):
    t = Tree(nwk.strip())
    L, idx, sides = masks(t)
    bit = lambda pred: sum(1 << idx[l] for l in L if pred(acc(l.name)))
    animal_a = bit(lambda a: fam_of.get(a) in A_LIKE and tax.get(a, {}).get("kingdom") == "Metazoa")
    gpr1 = bit(lambda a: fam_of.get(a) == "GPR1")
    out = {}
    for l in L:
        a = acc(l.name)
        if a not in TESTED:
            continue
        b = 1 << idx[l]
        in_a = bool(best(sides, animal_a | b) & b)
        in_g = bool(best(sides, gpr1 | b) & b)
        out[a] = "Class A-like" if in_a and not in_g else "GPR1" if in_g and not in_a else \
            "both" if in_a else "neither"
    return out


def main():
    OUT.mkdir(exist_ok=True)
    rows = []
    for slab, s in SETS:
        for alab, al in ALNS:
            stem = f"representatives_{s}_{al}_kpicgappy"
            bp = placements((S9 / "trees" / f"{stem}.raxml.support").read_text())
            c = {}
            ml = [x for x in (S9 / "runs" / f"representatives_{s}" / f"{stem}.raxml.mlTrees").read_text().splitlines() if x.strip()]
            for nwk in ml:
                for a, p in placements(nwk).items():
                    c.setdefault(a, Counter())[p] += 1
            for a, p in bp.items():
                r = tax[a]
                rows.append([f"{slab}, {alab}", a, TESTED[a], r["phylum"], r["species"], p,
                             c[a]["Class A-like"], c[a]["GPR1"], c[a]["both"], c[a]["neither"], len(ml),
                             mean_ti[a], "yes" if mean_ti[a] >= CUTOFF else "no"])
            print("done", slab, alab)
    with open(OUT / "nonanimal_classA.tsv", "w") as fh:
        fh.write("run\taccession\tlineage\tphylum\tspecies\tbest_tree\tml_in_ClassA_like\tml_in_GPR1\t"
                 "ml_both\tml_neither\tof\tmean_tIndex\tunstable\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)


if __name__ == "__main__":
    main()
