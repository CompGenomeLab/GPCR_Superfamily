"""Collect the contested accessions and each family's top-scoring representatives.

The overlap analysis says which accessions more than one
family claims, but a list of contested sequences on its own cannot be resolved,
because there is nothing in it to resolve them against. This writes the set that
the assignment tree is built from, which is every contested accession plus the
highest-scoring members of each family, so that the contested sequences have
something unambiguous to cluster with.

Representatives are taken from each family's own search, ranked by the full-sequence
bit score for a HMMER search or the bit score column for an MMseqs2 seed search, and
restricted to accessions that are members of that family after the topology filter. The
routing that sends a family to another family's search is the same as in
analyze_family_overlaps.py and build_picker.py, and is repeated
here rather than imported so this script stays runnable on its own.

Contested accessions can also be representatives. A sequence picked as one of a
family's best-scoring members is not removed from the set when it is also claimed by
another family, it simply appears once.

Usage (from this folder, after `source ../../config/paths.sh`):
    python build_overlap_rep_set.py         --diagnostics $GPCR_ZENODO_DIR/04_overlap_assignment/overlap_diagnostics.tsv         --after-tsv $GPCR_ZENODO_DIR/03_topology_filter/topology_membership.tsv         --search-dir $GPCR_ZENODO_DIR/01_homology_search         --output $GPCR_ZENODO_DIR/04_overlap_assignment/assignment_tree/overlaps_and_reps_accession.txt
"""

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import exists, open_text  # noqa: E402

# Families ranked against a broader shared search rather than one of their own.
SHARED_SEARCH = {
    "classB1": "classB",
    "classB2": "classB",
    "Vomeronasal1": "classT",
}

# Families whose search is an MMseqs2 seed search, given as the directory holding it.
M8_SEARCH = {
    "GPR1": "GPR1_YEAST",
    "FSLB": "FSLB_DICDI",
}


def accession_of(token):
    """tr|A0A8H5TWS9|A0A8H5TWS9_FUSHE -> A0A8H5TWS9; anything else passes through."""
    fields = token.split("|")
    return fields[1].strip() if len(fields) >= 3 else token.strip()


def read_wide_table(path):
    """Read a family-per-column accession table into {family: set(accessions)}."""
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


def hmm_ranked(path):
    """Accessions in descending full-sequence bit score order."""
    hits = []
    in_scores = header_seen = reading = False
    with open_text(path) as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("Domain annotation for each sequence:"):
                break
            if not in_scores:
                if line.startswith("Scores for complete sequences"):
                    in_scores = True
                continue
            if not header_seen:
                if line.lstrip().startswith("E-value"):
                    header_seen = True
                continue
            if not reading:
                if set(line.strip()) <= {"-", " "} and "-" in line:
                    reading = True
                continue
            parts = line.split()
            if len(parts) < 9:
                continue
            try:
                score = float(parts[1])
            except ValueError:
                continue
            hits.append((score, accession_of(parts[8])))
    hits.sort(key=lambda hit: -hit[0])
    return [accession for _, accession in hits]


def m8_ranked(path):
    """Accessions in descending bit score order.

    The layout is query, target, theader, pident, fident, nident, alnlen, qlen, tlen,
    mismatch, raw, bits, qcov, tcov, evalue, so the bit score is column 11.
    """
    hits = []
    with open_text(path) as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            try:
                score = float(parts[11])
            except ValueError:
                continue
            hits.append((score, accession_of(parts[1])))
    hits.sort(key=lambda hit: -hit[0])
    return [accession for _, accession in hits]


def ranked_for(family, search_dir):
    """Accessions of the search this family's members were found by, best first."""
    search_dir = Path(search_dir)
    if family in M8_SEARCH:
        folder = search_dir / M8_SEARCH[family]
        paths = sorted(list(folder.glob("*.m8")) + list(folder.glob("*.m8.gz")))
        if not paths:
            sys.exit(f"no .m8 file for {family} in {folder}")
        return m8_ranked(paths[0])
    search = SHARED_SEARCH.get(family, family)
    path = search_dir / search / f"{search}_hmmsearch_E1000.txt"
    if not exists(path):
        sys.exit(f"missing search output for {family}: {path}")
    return hmm_ranked(path)


def main():
    parser = argparse.ArgumentParser(
        description="Contested accessions plus each family's top representatives")
    parser.add_argument("--diagnostics", default=None,
                        help="overlap_diagnostics_7TM_*.tsv from "
                             "analyze_family_overlaps.py. Optional: without it the "
                             "contested accessions are derived from --after-tsv "
                             "directly, which needs no scoring run and cannot fall "
                             "out of step with the table the tree is built from.")
    parser.add_argument("--after-tsv", required=True,
                        help="wide family-by-accession table after 7TM filtering")
    parser.add_argument("--search-dir", required=True,
                        help="search_results tree, one directory per search")
    parser.add_argument("--output", required=True)
    parser.add_argument("--top", type=int, default=10,
                        help="representatives per family (default 10)")
    parser.add_argument("--reps-for", choices=("overlapping", "all"),
                        default="overlapping",
                        help="which families get representatives. 'overlapping' is "
                             "only those still sharing an accession with another "
                             "family, which is the default because a family with "
                             "nothing contested has nothing for the tree to resolve.")
    parser.add_argument("--drop-pair", action="append", default=[],
                        metavar="FAM_A,FAM_B",
                        help="exclude accessions whose claiming families are exactly "
                             "this pair, for an overlap already resolved by other "
                             "means. Repeatable. Rows naming a third family are kept, "
                             "because they still carry that family's signal.")
    args = parser.parse_args()

    drop = [set(pair.split(",")) for pair in args.drop_pair]

    families, members = read_wide_table(args.after_tsv)

    contested, dropped = [], 0
    if args.diagnostics:
        with open(args.diagnostics, newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                claiming = {f.strip() for f in row["families"].split(",") if f.strip()}
                if claiming in drop:
                    dropped += 1
                    continue
                contested.append(row["accession"].strip())
    else:
        claims = {}
        for family in families:
            for accession in members[family]:
                claims.setdefault(accession, set()).add(family)
        for accession in sorted(a for a, seen in claims.items() if len(seen) > 1):
            if claims[accession] in drop:
                dropped += 1
                continue
            contested.append(accession)

    # A family with nothing contested has nothing for the tree to resolve, so
    # anchoring it would add ten sequences the tree never has to arbitrate.
    involved = {family for family in families
                if any(a in members[family] for a in contested)}
    if args.reps_for == "overlapping":
        rep_families = [f for f in families if f in involved]
        skipped = [f for f in families if f not in involved]
        if skipped:
            print("no surviving overlap, so no representatives: %s"
                  % ", ".join(skipped))
    else:
        rep_families = list(families)

    reps, per_family = [], {}
    for family in rep_families:
        seen, picked = set(), []
        for accession in ranked_for(family, args.search_dir):
            if accession in seen or accession not in members[family]:
                continue
            seen.add(accession)
            picked.append(accession)
            if len(picked) == args.top:
                break
        per_family[family] = picked
        reps.extend(picked)

    final, seen = [], set()
    for accession in contested + reps:
        if accession not in seen:
            seen.add(accession)
            final.append(accession)

    with open(args.output, "w", newline="\n", encoding="utf-8") as handle:
        handle.write("\n".join(final) + "\n")

    contested_set = set(contested)
    for family in rep_families:
        new = len([a for a in per_family[family] if a not in contested_set])
        print(f"  {family:<14} {len(per_family[family]):2d} representatives, {new:2d} new")
    if dropped:
        print(f"\ncontested rows dropped by --drop-pair: {dropped:,}")
    print(f"contested accessions kept : {len(contested):,}")
    print(f"representatives           : {len(reps):,}")
    print(f"total written             : {len(final):,} -> {args.output}")


if __name__ == "__main__":
    main()
