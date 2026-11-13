#!/usr/bin/env python3
"""Write the anchor table the superfamily logo is built from.

Each family is placed onto the superfamily alignment through two sequences:

    acc1  the anchor, a representative present in both the family alignment and the
          superfamily alignment. Family columns are mapped to superfamily columns
          through it. Taken from `rank_bridge_representatives.py`, rank 1 per family.
    acc2  optional, a human reference whose GPCRdb numbering labels the columns. It is
          fixed per family and read from `gpcrdb_references.tsv`, because the numbering
          must not change when the anchor does.

The output has the columns `acc1`, `family`, `acc2`, the format the website's
`scripts/precompute-mappings.mjs` reads as `public/superfamily_logo_source_files/trim_info.tsv`,
and the one `trim_family_alignments.py` trims by.

Usage:
    python write_logo_anchors.py \\
        --ranking bridge_representatives_Sep24_topo0.8.tsv \\
        --out superfamily_logo_anchors_Sep24_topo0.8.tsv
"""

import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ranking", required=True,
                        help="bridge_representatives_<tag>.tsv from rank_bridge_representatives.py")
    parser.add_argument("--references", default=os.path.join(HERE, "gpcrdb_references.tsv"),
                        help="family -> acc2 table of the GPCRdb-numbered human references")
    parser.add_argument("--exclude", nargs="*", default=[],
                        help="families to leave out of the table")
    parser.add_argument("--add", nargs="*", default=[], metavar="FAMILY=ACC",
                        help="anchors chosen outside the ranking, such as a background family "
                             "placed by place_by_topology.py")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    with open(args.references, newline="") as handle:
        references = {r["family"]: r["acc2"].strip()
                      for r in csv.DictReader(handle, delimiter="\t")}

    anchors = {}
    with open(args.ranking, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["rank_in_family"] == "1":
                anchors[row["family"]] = row["accession"]

    for item in args.add:
        family, _, accession = item.partition("=")
        if not family or not accession:
            sys.exit(f"--add expects FAMILY=ACC, got {item!r}")
        if family in anchors:
            sys.exit(f"--add names {family}, which the ranking already anchors")
        anchors[family] = accession

    unknown = set(references) - set(anchors)
    if unknown:
        sys.exit(f"references name families the ranking has no anchor for: {sorted(unknown)}")

    excluded = {f.lower() for f in args.exclude}
    families = [f for f in sorted(anchors, key=str.lower) if f.lower() not in excluded]
    with open(args.out, "w", newline="") as out:
        out.write("acc1\tfamily\tacc2\n")
        for family in families:
            out.write(f"{anchors[family]}\t{family}\t{references.get(family, '')}\n")
            print(f"{family:<14}{anchors[family]:<13}{references.get(family, '')}")
    print(f"Wrote {args.out} ({len(families)} families)")


if __name__ == "__main__":
    main()
