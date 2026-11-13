"""Minimum ages of family splits from where the deep lineages fall, and paraphyly of cAMP and Class A.

1. Split dating. For a pair of sister groups X and Y, the X+Y clade is the smallest side of any
   branch holding more than half of the leaves of every family in X and in Y; X's clade is its
   best-F1 side inside it. For each taxonomic test group (e.g. early-diverging fungal Class A), the
   leaves are counted on X's side, on Y's side, or outside X+Y. A test group that falls on one side
   of the split means the split is older than that lineage's divergence from the lineages on the
   other side (no transfer assumed).
2. Paraphyly. In each tree, the smallest side holding more than half of each named part of a family
   (cAMP: fungal and non-fungal; Class A: animal and fungal), and the other families holding more
   than half of their leaves inside it. Divergent families (STE2, STE3, Class C, V2) are ignored.

Counted over the 100 ML search trees of the runs on disk and in the 8 best trees.
Writes results/split_dating.tsv and results/paraphyly.tsv.
"""
from collections import Counter

from core_tree_readout import kingdom
from tree_metrics import A_LIKE, DIVERGENT, OUT, RUNS, S9, acc, load
from gpcr_paths import stage  # noqa: E402  (path set by tree_metrics)

TAXA = stage("07_representative_selection", "representatives_all_families.tsv")
phylum = {}
for line in TAXA.read_text().splitlines()[1:]:
    f = line.split("\t")
    phylum[f[1]] = f[6] if f[6] not in ("NA", "") else f[9]   # unicellular holozoans: class column
DIKARYA = {"Ascomycota", "Basidiomycota"}
UNICELL_HOLOZOA = {"Choanoflagellata", "Filasterea", "Ichthyosporea"}


def tag(a, fam):
    """Taxonomic test group of one leaf, or None."""
    k, p = kingdom.get(a, ""), phylum.get(a, "")
    if fam == "classA" and k == "Fungi":
        return "Class A, Dikarya" if p in DIKARYA else "Class A, early-diverging fungi"
    if fam == "classA" and k not in ("Metazoa", "Fungi"):
        return "Class A, non-opisthokonts"
    if fam == "GP143" and p in UNICELL_HOLOZOA:
        return "GPR143, unicellular holozoans"
    if fam == "classB2" and p in UNICELL_HOLOZOA:
        return "Class B2, unicellular holozoans"
    if fam == "Mth" and k not in ("Metazoa",):
        return "Mth, amoebozoans"
    if fam == "FSLB":
        return "FSL, fungi" if k == "Fungi" else "FSL, amoebozoans"
    if fam == "LIT_Vomeronasal1":
        return "lamprey Vomeronasal 1 (literature)"
    if fam == "LIT_classT":
        return "shark Taste-2 (literature)"
    return None


PAIRS = [  # (name, X families, Y families, test groups)
    ("GPR1 | Class A-like", {"GPR1"}, A_LIKE,
     ["Class A, Dikarya", "Class A, early-diverging fungi", "Class A, non-opisthokonts"]),
    ("GPR143 | TMEM116", {"GP143"}, {"TM116"}, ["GPR143, unicellular holozoans"]),
    ("Mth | Class B2 + B1", {"Mth"}, {"classB2", "classB1"},
     ["Mth, amoebozoans", "Class B2, unicellular holozoans"]),
    ("Class B | Class F + FSL", {"Mth", "classB2", "classB1"}, {"classF", "FSLB"},
     ["Mth, amoebozoans", "FSL, amoebozoans", "FSL, fungi"]),
    ("Taste-2 | Vomeronasal 1", {"classT"}, {"Vomeronasal1"},
     ["shark Taste-2 (literature)", "lamprey Vomeronasal 1 (literature)"]),
]


def family_of(leaves, fam):
    """Family per leaf with literature leaves kept under their LIT_ name."""
    out = []
    for l, f in zip(leaves, fam):
        name = l.name
        out.append(name.split("|")[0] if f is None and name.startswith("LIT_") else f)
    return out


def masks(fam2, wanted):
    return sum(1 << i for i, f in enumerate(fam2) if f in wanted)


def sides_of(t, mask):
    full = mask[t]
    out = []
    for n in t.traverse():
        if n is not t:
            out += [mask[n], full ^ mask[n]]
    return out


def maj(side, m):
    return 2 * (side & m).bit_count() > m.bit_count()


def dating(newick):
    t, leaves, fam, mask = load(newick)
    fam2 = family_of(leaves, fam)
    sides = sides_of(t, mask)
    res = {}
    for name, X, Y, tests in PAIRS:
        fx = [masks(fam2, {f}) for f in X if masks(fam2, {f})]
        fy = [masks(fam2, {f}) for f in Y if masks(fam2, {f})]
        if not fx or not fy:
            continue
        cands = [s for s in sides if all(maj(s, m) for m in fx + fy)]
        if not cands:
            continue
        S = min(cands, key=lambda s: s.bit_count())
        mx = sum(fx)
        best, bf = None, -1
        for s in sides:
            if s & ~S:
                continue
            f1 = 2 * (s & mx).bit_count() / (s.bit_count() + mx.bit_count())
            if f1 > bf:
                best, bf = s, f1
        for g in tests:
            gm = sum(1 << i for i, (l, f) in enumerate(zip(leaves, fam2)) if tag(acc(l.name), f) == g)
            if not gm:
                continue
            res[(name, g)] = ((gm & best).bit_count(), (gm & S & ~best).bit_count(), (gm & ~S).bit_count())
    return res


def paraphyly(newick):
    t, leaves, fam, mask = load(newick)
    fam2 = []
    for l, f in zip(leaves, fam):
        a = acc(l.name)
        if f == "cAMP":
            f = "cAMP (fungi)" if kingdom.get(a) == "Fungi" else "cAMP (others)"
        elif f == "classA":
            f = {"Metazoa": "Class A (animals)", "Fungi": "Class A (fungi)"}.get(kingdom.get(a), "Class A (others)")
        fam2.append(f)
    ftot = Counter(f for f in fam2 if f and f not in DIVERGENT)
    fm = {f: masks(fam2, {f}) for f in ftot}
    sides = sides_of(t, mask)
    out = {}
    for name, parts in (("cAMP", ["cAMP (fungi)", "cAMP (others)"]),
                        ("Class A", ["Class A (animals)", "Class A (fungi)"])):
        cands = [s for s in sides if all(maj(s, fm[p]) for p in parts if p in fm)]
        S = min(cands, key=lambda s: s.bit_count())
        inside = sorted(f for f in ftot if f not in parts and maj(S, fm[f]))
        out[name] = " + ".join(inside) if inside else "(none: monophyletic)"
    return out


def main():
    date_rows, para_rows = [], []
    for label, best, ml in RUNS:
        trees = [("best tree", (S9 / "trees" / f"{best}.raxml.support").read_text())]
        if ml:
            trees += [("ML searches", l) for l in (S9 / ml).read_text().splitlines() if l.strip()]
        dsum, psum = {}, {}
        for kind, nwk in trees:
            for key, v in dating(nwk).items():
                dsum.setdefault((kind,) + key, []).append(v)
            for key, v in paraphyly(nwk).items():
                psum.setdefault((kind, key), Counter())[v] += 1
        for (kind, pair, g), vals in dsum.items():
            n = len(vals)
            tot = [sum(v[i] for v in vals) for i in range(3)]
            date_rows.append([label, kind, pair, g, n, *(round(x / n, 1) for x in tot)])
        for (kind, fam), c in psum.items():
            for nested, k in c.most_common(4):
                para_rows.append([label, kind, fam, nested, k])
    with open(OUT / "split_dating.tsv", "w") as fh:
        fh.write("run\ttrees\tsplit\ttest_group\tn_trees\tmean_on_first_side\tmean_on_second_side\tmean_outside\n")
        for r in date_rows:
            fh.write("\t".join(map(str, r)) + "\n")
    with open(OUT / "paraphyly.tsv", "w") as fh:
        fh.write("run\ttrees\tfamily\tfamilies_nested_in_its_smallest_clade\tn_trees\n")
        for r in para_rows:
            fh.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main()
