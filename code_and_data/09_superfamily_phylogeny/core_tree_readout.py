"""Readout of the multi-family groupings in the 8 final trees, for the hierarchy section.

Groupings are family-level splits named in "core" mode (divergent families ignored, step10 rule).
Fungal cAMP is split off as its own pseudo-family ("cAMP_fungi") so the fungal-cAMP block that
sits next to Class A-like + GPR1 in the best trees can be tested.

For each best tree: TBE of the branch carrying each grouping (highest TBE if several branches give
the same family-level split; "-" if absent). For every run (100 ML search trees each): number
of trees recovering it. Writes results/core_groupings_tbe_recovery.tsv.
"""
from collections import Counter

from tree_metrics import A_LIKE, B, DIVERGENT, F, OUT, RUNS, S9, acc, family_split, load
from gpcr_paths import stage  # noqa: E402  (path set by tree_metrics)

TAXA = stage("07_representative_selection", "representatives_all_families.tsv")
kingdom = {}
for line in TAXA.read_text().splitlines()[1:]:
    f = line.split("\t")
    kingdom[f[1]] = f[5]

GROUPS = {
    "Taste-2 + V1": {"classT", "Vomeronasal1"},
    "GP143 + TM116": {"GP143", "TM116"},
    "Class B1 + B2": {"classB1", "classB2"},
    "Class B (B1+B2+Mth)": B,
    "Class F + FSL": F,
    "Class A-like": A_LIKE,
    "Class A-like + GPR1": A_LIKE | {"GPR1"},
    "Class B + Class F/FSL": B | F,
    "Class B + F/FSL + GP143 + TM116": B | F | {"GP143", "TM116"},
    "cAMP (non-fungal) + GP157": {"cAMP", "GP157"},
    "fungal cAMP": {"cAMP_fungi"},
    "Class A-like + GPR1 + fungal cAMP": A_LIKE | {"GPR1", "cAMP_fungi"},
}


def relabel(leaves, fam):
    return [("cAMP_fungi" if f == "cAMP" and kingdom.get(acc(l.name)) == "Fungi" else f)
            for l, f in zip(leaves, fam)]


def groupings(newick):
    """-> {grouping: best TBE among branches carrying it}"""
    t, leaves, fam, mask = load(newick)
    fam = relabel(leaves, fam)
    ftot = Counter(f for f in fam if f)
    present = frozenset(f for f in ftot if f not in DIVERGENT)
    found = {}
    for n in t.traverse():
        if n is t:
            continue
        a = family_split(mask[n], fam, ftot, DIVERGENT)
        for side in (a, present - a):
            for g, fams in GROUPS.items():
                if side == frozenset(fams) and fams <= present:
                    found[g] = max(found.get(g, 0), n.support)
    return found


def main():
    rows = []
    for label, best, ml in RUNS:
        found = groupings((S9 / "trees" / f"{best}.raxml.support").read_text())
        rec = Counter()
        if ml:
            for line in (S9 / ml).read_text().splitlines():
                if line.strip():
                    rec.update(groupings(line).keys())
        for g in GROUPS:
            rows.append([label, g, f"{found[g]:.2f}" if g in found else "-",
                         rec[g] if ml else ""])
    with open(OUT / "core_groupings_tbe_recovery.tsv", "w") as fh:
        fh.write("run\tgrouping\tTBE_best_tree\tml_trees_recovered_of_100\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main()
