"""Label every leaf of the assignment tree with the families that claim it.

The assignment tree holds two kinds of leaf, the contested
accessions and the per-family representatives, and neither is distinguishable from a
bare accession. This appends the claiming families to each leaf label so the tree can
be read, and so the picker can group by them.

A contested accession is labelled with every family that claimed it, in the column
order of the wide table, joined by a hyphen, for example `classA-Olfactory`. A
representative is labelled with the single family it was picked to represent, with a
trailing asterisk, for example `classA*`. The asterisk is what makes representatives
searchable on their own in a tree viewer.

A representative that is also contested is labelled as a representative. It was chosen
to anchor its family and that is the role it plays in the tree, so labelling it with
its full overlap would hide the anchor among the sequences it is there to resolve.

Only the label token changes. Topology, branch lengths and support values are written
back exactly as they were read.

Usage (from this folder, after `source ../../config/paths.sh`):
    python annotate_tree_with_families.py         --tree $GPCR_ZENODO_DIR/04_overlap_assignment/assignment_tree/overlaps_and_reps_accession_db_FAMSA_gt0.1.tree         --after-tsv $GPCR_ZENODO_DIR/03_topology_filter/topology_membership.tsv         --rep-set $GPCR_ZENODO_DIR/04_overlap_assignment/assignment_tree/overlaps_and_reps_accession.txt         --output $GPCR_ZENODO_DIR/04_overlap_assignment/assignment_tree/overlaps_and_reps.with_families.nwk
"""

import argparse
import collections
import csv
import re


def accession_of(label):
    """tr|A0A8H5TWS9|A0A8H5TWS9_FUSHE -> A0A8H5TWS9; anything else passes through."""
    parts = label.split("|") if "|" in label else label.split("_")
    return parts[1].strip() if len(parts) >= 2 else label.strip()


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
    parser = argparse.ArgumentParser(description="Append family labels to tree leaves")
    parser.add_argument("--tree", required=True)
    parser.add_argument("--diagnostics", default=None,
                        help="optional. Without it the claiming families "
                             "are read off --after-tsv directly.")
    parser.add_argument("--after-tsv", required=True)
    parser.add_argument("--rep-set", required=False,
                        help="the accession list from build_overlap_rep_set.py. Accessions in it that "
                             "are not contested are the representatives. Without it "
                             "every non-contested leaf is treated as one.")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    order, members = read_wide_table(args.after_tsv)
    rank = {family: index for index, family in enumerate(order)}

    # Contested accessions carry every family that claimed them. Without a
    # diagnostics file the claims are read off the wide table itself, which is the
    # same information and cannot fall out of step with the table the tree was
    # built from.
    claimed = {}
    if args.diagnostics:
        with open(args.diagnostics, newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                families = [f.strip() for f in row["families"].split(",") if f.strip()]
                if families:
                    claimed[row["accession"].strip()] = sorted(
                        set(families), key=lambda f: rank.get(f, len(order)))
    else:
        claims = {}
        for family in order:
            for accession in members[family]:
                claims.setdefault(accession, []).append(family)
        claimed = {accession: families for accession, families in claims.items()
                   if len(families) > 1}

    # Anything else in the step-4 set is a representative, labelled by the single
    # family whose member list it belongs to.
    single = {}
    if args.rep_set:
        with open(args.rep_set, encoding="utf-8") as handle:
            for line in handle:
                accession = line.strip()
                if not accession or accession in claimed:
                    continue
                owners = [f for f in order if accession in members[f]]
                if len(owners) == 1:
                    single[accession] = owners[0]

    text = open(args.tree, encoding="utf-8").read()
    stats = collections.Counter()
    unmatched = []

    def relabel(match):
        prefix, raw = match.group(1), match.group(2).strip()
        if not raw:
            return match.group(0)
        accession = accession_of(raw)
        if accession in single:
            tag = single[accession] + "*"
            stats["representative"] += 1
        elif accession in claimed:
            tag = "-".join(claimed[accession])
            stats["n%d" % min(len(claimed[accession]), 4)] += 1
        else:
            stats["unmatched"] += 1
            unmatched.append(raw)
            return match.group(0)
        stats["labelled"] += 1
        return prefix + raw + "|" + tag

    # A leaf label sits after '(' or ',' and runs to ':' ',' ')' or ';'.
    out = re.sub(r"([(,])([^(),:;]+?)(?=:|[,);])", relabel, text)
    with open(args.output, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(out)

    print("leaves labelled  : %d" % stats["labelled"])
    print("  representatives: %d" % stats["representative"])
    print("  2 families     : %d" % stats["n2"])
    print("  3 families     : %d" % stats["n3"])
    print("  4 or more      : %d" % stats["n4"])
    print("unmatched        : %d %s" % (stats["unmatched"], unmatched[:5]))
    print("written -> %s" % args.output)


if __name__ == "__main__":
    main()
