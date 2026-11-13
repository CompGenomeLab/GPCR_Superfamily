"""Write the final per-family accession lists, one file per family.

The end of the pipeline. Each list is a family's membership after the 7TM topology
filter and after every surviving cross-family overlap has been resolved, which is the
set to align and cluster.

Three things happen to the 7TM membership here.

A contested accession that was assigned to one family is removed from the families
that lost it. A contested accession left ambiguous is removed from every family
claiming it, because the evidence did not favour one over another and leaving it in
all of them would assert a membership no analysis supports.

A row naming a family that never claimed the sequence is refused rather than applied.
Applying it would strip the sequence from the families that did claim it and add it to
none, losing it silently. Such a row means the assignment was made against a stale
table, so it is reported and the sequence is left contested.

Uncontested sequences, which are the overwhelming majority, pass through untouched.

Usage (from this folder, after `source ../../config/paths.sh`):
    python write_final_members.py         --after-tsv $GPCR_ZENODO_DIR/03_topology_filter/topology_membership.tsv         --selections mrca_selections.tsv         --output-dir $GPCR_ZENODO_DIR/04_overlap_assignment/final_families         --merged-tsv $GPCR_ZENODO_DIR/04_overlap_assignment/final_membership.tsv
"""

import argparse
import csv
import os


def read_wide_table(path):
    """Column order and {family: set(accessions)} from a family-per-column table."""
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    header = [name.strip() for name in rows[0]]
    table = {name: set() for name in header}
    for row in rows[1:]:
        for name, value in zip(header, row):
            value = value.strip()
            if value:
                table[name].add(value)
    return header, table


def main():
    parser = argparse.ArgumentParser(
        description="Per-family accession lists after 7TM and overlap resolution")
    parser.add_argument("--after-tsv", required=True,
                        help="wide family-by-accession table after 7TM filtering")
    parser.add_argument("--selections", required=True,
                        help="mrca_selections.tsv from the assignment tool")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--suffix", default="_final_members.txt")
    parser.add_argument("--merged-tsv", default=None,
                        help="also write the result as one wide table, for any step "
                             "that reads the family-per-column form")
    args = parser.parse_args()

    families, members = read_wide_table(args.after_tsv)
    before = {family: len(ids) for family, ids in members.items()}

    assigned = ambiguous = refused = 0
    refused_rows = []
    with open(args.selections, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            kind = (row.get("kind") or "").strip()
            accession = row["accession"].strip()
            if kind == "rep":
                continue
            # `ambiguous` is a sequence no evidence resolved, `excluded` one set aside
            # on purpose. They are recorded apart because only the second has to
            # survive a rebuild of the picker, but neither was awarded to a family, so
            # both leave every family that claimed them.
            if kind in ("ambiguous", "excluded"):
                for ids in members.values():
                    ids.discard(accession)
                ambiguous += 1
                continue
            if kind not in ("clade", "sweep", "manual"):
                continue
            winner = row["family"].strip()
            if accession not in members.get(winner, ()):
                refused += 1
                refused_rows.append((accession, winner, row.get("label", "")))
                continue
            for family, ids in members.items():
                if family != winner:
                    ids.discard(accession)
            assigned += 1

    os.makedirs(args.output_dir, exist_ok=True)
    for family in families:
        path = os.path.join(args.output_dir, family + args.suffix)
        with open(path, "w", newline="\n", encoding="utf-8") as handle:
            handle.write("\n".join(sorted(members[family])) + "\n")

    if args.merged_tsv:
        height = max(len(ids) for ids in members.values())
        ordered = {f: sorted(members[f]) for f in families}
        with open(args.merged_tsv, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(families)
            for index in range(height):
                writer.writerow([ordered[f][index] if index < len(ordered[f]) else ""
                                 for f in families])

    if refused:
        print("WARNING: %d rows name a family that never claimed the sequence and were "
              "ignored, so those sequences are still contested:" % refused)
        for accession, winner, label in refused_rows[:10]:
            print("    %-12s -> %-12s but claimed by %s" % (accession, winner, label))
        print()

    print("%d assigned, %d left ambiguous\n" % (assigned, ambiguous))
    print("%-14s %10s %10s %8s" % ("family", "after 7TM", "final", "removed"))
    total_before = total_after = 0
    for family in families:
        after = len(members[family])
        total_before += before[family]
        total_after += after
        print("  %-13s %9d %10d %8d" % (family, before[family], after,
                                        before[family] - after))
    print("  %-13s %9d %10d %8d" % ("TOTAL", total_before, total_after,
                                    total_before - total_after))

    overlapping = sum(1 for family in families for other in families
                      if family < other and members[family] & members[other])
    print("\nfamily pairs still sharing an accession: %d" % overlapping)
    print("written -> %s" % args.output_dir)


if __name__ == "__main__":
    main()
