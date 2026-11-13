"""Internal order of Class B and Class F/FSL, and the placement of named non-animal sequences,
over the best tree and the 100 ML searches of every Oct08 run.

Only stable leaves enter the Class B / F tests (mean tIndex below the STE2 cutoff), as in the text.
Trees are rooted on a cAMP leaf for the Class B / F questions. For a family X and a leaf group Y,
"X | Y" asks whether the smallest clade holding all of X excludes Y (Y outside) and whether the
smallest clade holding X and Y holds anything else.

Tests
  B1 in B2       smallest clade of all Class B2 leaves contains Class B1 leaves (B1 nested in B2)
  B2 in B1       the reverse
  B1+B2 excl Mth smallest clade of B1 + B2 holds no animal Mth leaf
  Mth+B2 excl B1 smallest clade of animal Mth + B2 holds no B1 leaf
  amoebozoan FSL with Class F   smallest clade of Class F + amoebozoan FSL holds no fungal FSL
  fungal FSL with Class F       smallest clade of Class F + fungal FSL holds no amoebozoan FSL
Named sequences: inside their family's best-F1 clade (nested lineages merged for Class A-like).
Writes results/internal_order.tsv and results/named_sequences.tsv.
Usage: python internal_order_Oct08.py
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
FAMMAP = S9 / "representatives_family_map.tsv"
TAXA = stage("07_representative_selection", "representatives_all_families.tsv")
TINDEX = S9 / "booster_stats" / "tindex_mean_across_final_trees.tsv"
SETS = [("all", "all"), ("no STE2", "no_ste2"), ("core", "no_ste2ste3classCV2"), ("core+lit", "and_literature")]
ALNS = [("ep0.123", "ginsi_ep0.123"), ("ginsi", "ginsi")]
A_LIKE = {"classA", "Olfactory", "Vomeronasal1", "classT", "Nematode"}
CUTOFF = 0.6919      # most stable STE2 sequence (lowest STE2 mean tIndex)

fam_of = dict(l.split("\t") for l in FAMMAP.read_text().splitlines() if l.strip())
tax = {r["accession"]: r for r in csv.DictReader(TAXA.open(), delimiter="\t")}
mean_ti = {r["accession"]: float(r["mean_tIndex"]) for r in csv.DictReader(TINDEX.open(), delimiter="\t")}
acc = lambda name: name.split("|")[1]


def group(a):
    f = fam_of.get(a, "")
    if f in A_LIKE:
        return "Class A-like"
    return f


# named non-animal sequences of the text, by family and genus
NAMED = {}
for a, r in tax.items():
    g = r["genus"]
    f = r["gpcrfamily"]
    if (f, g) in {("GP143", "Capsaspora"), ("GP157", "Capsaspora"), ("GP157", "Acanthamoeba"),
                  ("GP157", "Triparma"), ("GP157", "Anaeramoeba"), ("GPR1", "Planoprotostelium"),
                  ("classA", "Triparma"), ("Mth", "Dictyostelium"), ("Mth", "Polysphondylium"),
                  ("Mth", "Heterostelium"), ("Mth", "Cavenderia"), ("Mth", "Tieghemostelium")}:
        NAMED[a] = (f, g + " " + r["species"].split()[-1])


def smallest(t, leaves):
    return t.get_common_ancestor(leaves) if len(leaves) > 1 else leaves[0]


def best_clade(t, members):
    L = set(t.get_leaves())
    best = None
    for n in t.traverse():
        if n.is_root():
            continue
        inside = set(n.get_leaves())
        for s in (inside, L - inside):
            f1 = 2 * len(s & members) / (len(s) + len(members))
            if best is None or f1 > best[0]:
                best = (f1, s)
    return best[1]


def tests(nwk):
    t = Tree(nwk.strip())
    L = t.get_leaves()
    by = lambda pred: [l for l in L if pred(acc(l.name)) and mean_ti.get(acc(l.name), 0) < CUTOFF]
    t.set_outgroup(by(lambda a: fam_of.get(a) == "cAMP" and tax.get(a, {}).get("kingdom") == "Metazoa")[0])
    b2_an = by(lambda a: fam_of.get(a) == "classB2" and tax.get(a, {}).get("kingdom") == "Metazoa")
    b1, b2 = by(lambda a: fam_of.get(a) == "classB1"), by(lambda a: fam_of.get(a) == "classB2")
    mth_an = by(lambda a: fam_of.get(a) == "Mth" and tax.get(a, {}).get("kingdom") == "Metazoa")
    F = by(lambda a: fam_of.get(a) == "classF")
    fsl_am = by(lambda a: fam_of.get(a) == "FSLB" and tax.get(a, {}).get("kingdom") != "Fungi")
    fsl_fu = by(lambda a: fam_of.get(a) == "FSLB" and tax.get(a, {}).get("kingdom") == "Fungi")
    has = lambda node, leaves: any(l in set(node.get_leaves()) for l in leaves)
    out = {
        "B1 nested in B2": has(smallest(t, b2), b1),
        "B2 nested in B1": has(smallest(t, b1), b2),
        "B1+B2 exclude animal Mth": not has(smallest(t, b1 + b2), mth_an),
        "animal Mth+B2 exclude B1": not has(smallest(t, mth_an + b2), b1),
        "animal B2 monophyletic (no B1, no Mth)": not has(smallest(t, b2_an), b1 + mth_an),
        "B1 monophyletic": not has(smallest(t, b1), b2 + mth_an),
        "animal Mth monophyletic": not has(smallest(t, mth_an), b1 + b2),
        "Class F monophyletic": not has(smallest(t, F), fsl_am + fsl_fu),
        "fungal FSL monophyletic": not has(smallest(t, fsl_fu), F + fsl_am),
        "Class F + amoebozoan FSL exclude fungal FSL": not has(smallest(t, F + fsl_am), fsl_fu),
        "Class F + fungal FSL exclude amoebozoan FSL": not has(smallest(t, F + fsl_fu), fsl_am),
        "amoebozoan FSL monophyletic": not has(smallest(t, fsl_am), F + fsl_fu + b1 + b2),
    }
    named = {}
    for a, (f, _) in NAMED.items():
        leaf = next((l for l in L if acc(l.name) == a), None)
        if leaf is None:
            continue
        g = group(a)
        members = {l for l in L if group(acc(l.name)) == g and l is not leaf}
        named[a] = leaf in best_clade(t, members | {leaf}) if g != "Mth" else \
            leaf in best_clade(t, set(mth_an))   # dictyostelid Mth: inside the animal Mth clade?
    return out, named


def main():
    OUT.mkdir(exist_ok=True)
    order_rows, named_rows = [], []
    for slab, s in SETS:
        for alab, a in ALNS:
            run = f"{slab}, {alab}"
            stem = f"representatives_{s}_{a}_kpicgappy"
            best = (S9 / "trees" / f"{stem}.raxml.support").read_text()
            ml = (S9 / "runs" / f"representatives_{s}" / f"{stem}.raxml.mlTrees").read_text().splitlines()
            bo, bn = tests(best)
            oc, nc = Counter(), Counter()
            for line in ml:
                if line.strip():
                    o, n = tests(line)
                    oc.update(k for k, v in o.items() if v)
                    nc.update(k for k, v in n.items() if v)
            for k, v in bo.items():
                order_rows.append([run, k, "yes" if v else "no", oc[k], len(ml)])
            for k, v in bn.items():
                f, sp = NAMED[k]
                named_rows.append([run, k, f, sp, "yes" if v else "no", nc[k], len(ml),
                                   mean_ti.get(k, ""), "yes" if mean_ti.get(k, 0) >= CUTOFF else "no"])
            print("done", run)
    with open(OUT / "internal_order.tsv", "w") as fh:
        fh.write("run\ttest\tbest_tree\tml_trees\tof\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in order_rows)
    with open(OUT / "named_sequences.tsv", "w") as fh:
        fh.write("run\taccession\tfamily\tspecies\tin_family_clade_best_tree\tml_trees\tof\tmean_tIndex\tunstable\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in named_rows)


if __name__ == "__main__":
    main()
