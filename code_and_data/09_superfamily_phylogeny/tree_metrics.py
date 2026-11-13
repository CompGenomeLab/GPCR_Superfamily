"""Per-family tree metrics on the 8 final Oct08 trees, for deciding which families can
enter the superfamily tree.

Families are merged into units (nested lineages with their parent class, see UNIT); literature
sequences are ignored. For every tree:

  * the unit's clade is its best-F1 side of any branch (tolerates a few stray leaves);
  * stem = length of that branch; reach = stem + median distance from the unit's leaves to the
    clade's own end of the branch (how far a typical tip is from where the unit attaches);
  * isolation = per leaf, patristic distance to the nearest leaf of another unit (median);
  * sister = the smaller of the two subtrees that meet the unit's stem at its far end, named by
    the units holding >= half of their leaves in it (rooting-free).

Grouping recovery follows step10 (family-level, each family on the side holding > half of its
leaves), counted over the 100 ML search trees of each run and in the best trees, in two modes:
"strict" names every family, "core" ignores the divergent families when naming, so it asks
whether the relationships among the core families themselves changed.

Writes results/tree_unit_metrics.tsv, results/sister_counts.tsv, results/grouping_recovery.tsv.
"""
import statistics as st
from collections import Counter, defaultdict
from pathlib import Path

from ete3 import Tree

HERE = Path(__file__).resolve().parent
import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import stage  # noqa: E402
S9 = stage("09_superfamily_phylogeny")
OUT = S9 / "readouts"
FAMMAP = S9 / "representatives_family_map.tsv"

UNIT = {"classA": "Class A-like", "Olfactory": "Class A-like", "Vomeronasal1": "Class A-like",
        "classT": "Class A-like", "Nematode": "Class A-like", "GPR1": "GPR1",
        "classB1": "Class B", "classB2": "Class B", "Mth": "Class B",
        "classF": "Class F/FSL", "FSLB": "Class F/FSL", "GP143": "GP143", "TM116": "TM116",
        "GP157": "GP157", "cAMP": "cAMP", "classC": "Class C/V2", "Vomeronasal2": "Class C/V2",
        "STE2": "STE2", "STE3": "STE3"}
DIVERGENT = {"STE2", "STE3", "classC", "Vomeronasal2"}

# (label, best tree stem in trees_Oct08, 100 ML search trees in runs_Oct08 or None)
_SETS = [("all", "all"), ("no STE2", "no_ste2"), ("core", "no_ste2ste3classCV2"),
         ("core+lit", "and_literature")]
RUNS = [(f"{label}, {short}", f"representatives_{s}_{a}_kpicgappy",
         f"runs/representatives_{s}/representatives_{s}_{a}_kpicgappy.raxml.mlTrees")
        for label, s in _SETS for short, a in (("ep0.123", "ginsi_ep0.123"), ("ginsi", "ginsi"))]

A_LIKE = {"classA", "Olfactory", "Vomeronasal1", "classT", "Nematode"}
B = {"classB1", "classB2", "Mth"}
F = {"classF", "FSLB"}
GROUPINGS = {
    "Taste-2 + V1": {"classT", "Vomeronasal1"},
    "GP143 + TM116": {"GP143", "TM116"},
    "Class B (B1+B2+Mth)": B,
    "Class F + FSL": F,
    "Class A-like": A_LIKE,
    "Class A-like + GPR1": A_LIKE | {"GPR1"},
    "Class B + Class F/FSL": B | F,
    "Class C/V2 + Nematode": {"classC", "Vomeronasal2", "Nematode"},
    "Class A-like without Nematode": A_LIKE - {"Nematode"},
}

fam_of = {}
for line in FAMMAP.read_text().splitlines():
    if line.strip():
        acc, fam = line.split("\t")[:2]
        fam_of[acc] = None if fam.startswith("LIT_") else fam


def acc(name):
    return name.split("|")[1] if "|" in name else name


def load(newick):
    t = Tree(newick.strip())
    leaves = t.get_leaves()
    idx = {l: i for i, l in enumerate(leaves)}
    fam = [fam_of[acc(l.name)] for l in leaves]
    mask = {}
    for n in t.traverse("postorder"):
        mask[n] = (1 << idx[n]) if n.is_leaf() else sum(mask[c] for c in n.children)
    return t, leaves, fam, mask


def members(fam, wanted):
    return sum(1 << i for i, f in enumerate(fam) if f is not None and wanted(f))


def popcount(x):
    return x.bit_count()


def best_side(t, mask, full, target):
    """Best-F1 side of any branch for the leaf set `target` -> (f1, node, side_is_below)."""
    nt, best = popcount(target), (0, None, True)
    for n in t.traverse():
        if n is t:
            continue
        for below, side in ((True, mask[n]), (False, full ^ mask[n])):
            inter = popcount(side & target)
            f1 = 2 * inter / (popcount(side) + nt)
            if f1 > best[0]:
                best = (f1, n, below)
    return best


def adjacent(t, mask, full, node, below):
    """Leaf masks of the two subtrees meeting the unit's stem at its far end."""
    if below:
        v = node.up
        parts = [mask[c] for c in v.children if c is not node]
        if v is not t:
            parts.append(full ^ mask[v])
    else:
        parts = [mask[c] for c in node.children]
    return parts


def unit_counts(fam):
    return Counter(UNIT[f] for f in fam if f)


def name_side(side, fam, totals):
    inside = Counter(UNIT[fam[i]] for i in range(len(fam)) if side >> i & 1 and fam[i])
    named = sorted(u for u, k in inside.items() if 2 * k >= totals[u])
    if named:
        return " + ".join(named)
    u, k = inside.most_common(1)[0] if inside else ("-", 0)
    return f"part of {u}"


def family_split(side, fam, ftot, ignore):
    """Families (not units) holding > half of their leaves on `side`."""
    inside = Counter(fam[i] for i in range(len(fam)) if side >> i & 1 and fam[i])
    return frozenset(f for f, k in inside.items() if 2 * k > ftot[f] and f not in ignore)


def recovered(t, mask, full, fam, ignore):
    ftot = Counter(f for f in fam if f)
    present = frozenset(f for f in ftot if f not in ignore)
    splits = set()
    for n in t.traverse():
        if n is t:
            continue
        a = family_split(mask[n], fam, ftot, ignore)
        splits.add(a)
        splits.add(present - a)
    return {g for g, fams in GROUPINGS.items() if fams <= present and frozenset(fams) in splits}


def leaf_distances(t, leaves):
    """Patristic distance matrix between leaves (list of lists)."""
    nodes = list(t.traverse())
    adj = defaultdict(list)
    for n in nodes:
        if n.up is not None:
            adj[n].append((n.up, n.dist))
            adj[n.up].append((n, n.dist))
    out = []
    for l in leaves:
        d, stack = {l: 0.0}, [l]
        while stack:
            x = stack.pop()
            for y, w in adj[x]:
                if y not in d:
                    d[y] = d[x] + w
                    stack.append(y)
        out.append([d[m] for m in leaves])
    return out


def per_tree_unit_metrics(t, leaves, fam, mask, with_distances):
    full = (1 << len(leaves)) - 1
    totals = unit_counts(fam)
    unit_of_leaf = [UNIT[f] if f else None for f in fam]
    dist = leaf_distances(t, leaves) if with_distances else None
    rows = {}
    for u in totals:
        target = members(fam, lambda f: UNIT[f] == u)
        f1, node, below = best_side(t, mask, full, target)
        side = mask[node] if below else full ^ mask[node]
        parts = adjacent(t, mask, full, node, below)
        sister = min(parts, key=popcount)
        row = {"f1": f1, "stem": node.dist, "sister": name_side(sister, fam, totals)}
        if dist is not None:
            # distance from each member leaf on the unit's side to the unit's end of the stem
            end = node if below else node.up
            inside = [i for i in range(len(leaves)) if side >> i & 1 and unit_of_leaf[i] == u]
            row["reach"] = node.dist + st.median(leaves[i].get_distance(end) for i in inside)
            iso = [min(dist[i][j] for j in range(len(leaves))
                       if unit_of_leaf[j] not in (None, u))
                   for i in range(len(leaves)) if unit_of_leaf[i] == u]
            row["isolation"] = st.median(iso)
        rows[u] = row
    return rows


def main():
    OUT.mkdir(exist_ok=True)
    metric_rows, sister_rows, recovery_rows = [], [], []
    for label, best_stem, ml in RUNS:
        best_path = S9 / "trees" / f"{best_stem}.raxml.support"
        t, leaves, fam, mask = load(best_path.read_text())
        full = (1 << len(leaves)) - 1
        for u, r in per_tree_unit_metrics(t, leaves, fam, mask, True).items():
            metric_rows.append([label, u, f"{r['f1']:.2f}", f"{r['stem']:.3f}",
                                f"{r['reach']:.3f}", f"{r['isolation']:.3f}", r["sister"]])
        best_rec = {mode: recovered(t, mask, full, fam, ign)
                    for mode, ign in (("strict", set()), ("core", DIVERGENT))}

        sisters = defaultdict(Counter)
        rec = {"strict": Counter(), "core": Counter()}
        n_ml = 0
        if ml:
            for line in (S9 / ml).read_text().splitlines():
                if not line.strip():
                    continue
                n_ml += 1
                t2, l2, f2, m2 = load(line)
                full2 = (1 << len(l2)) - 1
                for u, r in per_tree_unit_metrics(t2, l2, f2, m2, False).items():
                    sisters[u][r["sister"]] += 1
                for mode, ign in (("strict", set()), ("core", DIVERGENT)):
                    for g in recovered(t2, m2, full2, f2, ign):
                        rec[mode][g] += 1
            for u, c in sisters.items():
                for s, k in c.most_common():
                    sister_rows.append([label, u, s, k, n_ml])
        for g in GROUPINGS:
            for mode in ("strict", "core"):
                recovery_rows.append([label, g, mode, "yes" if g in best_rec[mode] else "no",
                                      rec[mode][g] if ml else "", n_ml if ml else ""])
        print("done", label, n_ml)

    def write(name, header, rows):
        with open(OUT / name, "w") as fh:
            fh.write("\t".join(header) + "\n")
            for r in rows:
                fh.write("\t".join(map(str, r)) + "\n")

    write("tree_unit_metrics.tsv",
          ["run", "unit", "F1", "stem", "reach", "isolation_median", "sister_in_best_tree"],
          metric_rows)
    write("sister_counts.tsv", ["run", "unit", "sister", "n_trees", "of"], sister_rows)
    write("grouping_recovery.tsv",
          ["run", "grouping", "mode", "in_best_tree", "ml_trees_recovered", "of"], recovery_rows)


if __name__ == "__main__":
    main()
