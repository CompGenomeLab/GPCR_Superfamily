"""Apply the phylogeny-based family reassignments to the final membership and counts.

The reassignments are listed in phylogeny_reassignments.tsv (accession, from_family, to_family, ...):
the 136 Dikarya "Class A" sequences that fall in the GPR1 clade of the superfamily trees. The stage 04
membership is left unchanged; the corrected copies are what the family-wide alignments of this stage
(and the choice of structures in stage 11) are built from.

Reads  04_overlap_assignment/final_membership.tsv, acc_taxid_map.tsv, tax_counts_per_family.tsv
Writes 10_motif_analysis/phylogeny_reassignment/final_membership_reassigned.tsv,
       acc_taxid_map_reassigned.tsv and tax_counts_per_family_reassigned.tsv

Usage:
    python apply_phylogeny_reassignments.py
"""
import csv
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from gpcr_paths import stage  # noqa: E402

HERE = Path(__file__).resolve().parent
IN = stage("04_overlap_assignment")
OUT = stage("10_motif_analysis", "phylogeny_reassignment")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    with open(HERE / "phylogeny_reassignments.tsv", encoding="utf-8") as fh:
        moves = {r["accession"]: (r["from_family"], r["to_family"]) for r in csv.DictReader(fh, delimiter="\t")}

    # wide membership table
    with open(IN / "final_membership.tsv", encoding="utf-8") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    header = rows[0]
    members = {f: [r[i] for r in rows[1:] if i < len(r) and r[i]] for i, f in enumerate(header)}
    for acc, (src, dst) in moves.items():
        members[src].remove(acc)          # fails loudly if the accession is not where the list says
        members[dst].append(acc)
    height = max(len(v) for v in members.values())
    ordered = {f: sorted(members[f]) for f in header}
    with open(OUT / "final_membership_reassigned.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for i in range(height):
            w.writerow([ordered[f][i] if i < len(ordered[f]) else "" for f in header])

    # accession -> taxID map, and the two affected count columns
    moved = Counter()
    with open(IN / "acc_taxid_map.tsv", encoding="utf-8") as src, \
            open(OUT / "acc_taxid_map_reassigned.tsv", "w", newline="", encoding="utf-8") as dst:
        rd = csv.DictReader(src, delimiter="\t")
        w = csv.DictWriter(dst, fieldnames=rd.fieldnames, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rd:
            mv = moves.get(r["accession"])
            if mv and r["family"] == mv[0]:
                r["family"] = mv[1]
                moved[(r["taxID"], mv[0], mv[1])] += 1
            w.writerow(r)
    with open(IN / "tax_counts_per_family.tsv", encoding="utf-8") as src, \
            open(OUT / "tax_counts_per_family_reassigned.tsv", "w", newline="", encoding="utf-8") as dst:
        rd = csv.DictReader(src, delimiter="\t")
        w = csv.DictWriter(dst, fieldnames=rd.fieldnames, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rd:
            for (taxid, a, b), n in moved.items():
                if r["taxID"] == taxid:
                    r[a + "_count"] = str(int(r[a + "_count"]) - n)
                    r[b + "_count"] = str(int(r[b + "_count"]) + n)
            w.writerow(r)
    print("%d sequences moved" % sum(moved.values()))


if __name__ == "__main__":
    main()
