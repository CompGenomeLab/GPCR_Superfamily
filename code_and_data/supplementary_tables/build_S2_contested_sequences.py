"""Supplementary Table S2: sequences claimed by two families after the topology filter, and how each was resolved.

Claiming families: step3 accessions_by_family_7TM_topology0.8_Sep24.tsv (accessions present in more than one
family column). Resolution: step5 mrca_selections.tsv (Sep30), the record written by the assignment tool:
  sweep      assigned by its clade in the assignment tree (the clade of the family's top-scoring members)
  manual     decided individually during review (clustering not decisive; median-normalized scores used)
  ambiguous  left unresolved and removed from every claiming family
Final family: step5 final_families/accessions_by_family_final_Sep24.tsv.

    python build_S2_contested_sequences.py
"""
import csv
from pathlib import Path

import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import stage  # noqa: E402
CODE = Path(__file__).resolve().parent
HERE = stage("supplementary_tables")
HERE.mkdir(parents=True, exist_ok=True)
NICE = {"classA": "Class A", "classB1": "Class B1", "classB2": "Class B2", "classC": "Class C", "classF": "Class F",
        "classT": "Taste-2", "FSLB": "FSL", "Vomeronasal1": "Vomeronasal 1", "Vomeronasal2": "Vomeronasal 2"}
HOW = {"sweep": "clade in assignment tree", "manual": "individual review", "ambiguous": "unresolved, removed"}


def columns(path):
    rows = list(csv.reader(open(path, encoding="utf-8"), delimiter="\t"))
    cols = {c: set() for c in rows[0]}
    for r in rows[1:]:
        for c, v in zip(rows[0], r):
            if v.strip():
                cols[c].add(v.strip())
    return cols


def main():
    topo = columns(stage("03_topology_filter", "topology_membership.tsv"))
    topo.pop("SLT_TRNS", None)
    final = columns(stage("04_overlap_assignment", "final_membership.tsv"))
    claims = {}
    for fam, accs in topo.items():
        for a in accs:
            claims.setdefault(a, []).append(fam)
    contested = {a: f for a, f in claims.items() if len(f) > 1}
    sel = {r["accession"]: r for r in csv.DictReader(
        open(CODE.parent / "04_overlap_assignment" / "mrca_selections.tsv", encoding="utf-8"), delimiter="\t")
        if r["kind"] != "rep"}
    in_final = {a: f for f, accs in final.items() for a in accs}

    rows, missing = [], []
    for a in sorted(contested, key=lambda a: (sorted(contested[a]), a)):
        r = sel.get(a)
        if r is None:
            missing.append(a)
        rows.append({
            "Accession": a,
            "Claiming families": ", ".join(NICE.get(f, f) for f in sorted(contested[a])),
            "Assigned family": NICE.get(in_final.get(a, ""), in_final.get(a, "")) or "none",
            "How resolved": HOW.get(r["kind"], r["kind"]) if r else "not in assignment record",
        })
    out = HERE / "Table_S2_contested_sequences.tsv"
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(contested)} contested; {len(missing)} without an assignment row; "
          f"{len(set(sel) - set(contested))} assignment rows not contested at this stage")
    from collections import Counter
    print(Counter(r["How resolved"] for r in rows))
    print(Counter(r["Claiming families"] for r in rows).most_common(12))
    print("wrote", out)


if __name__ == "__main__":
    main()
