#!/usr/bin/env python3
"""Per-genome family counts at each filtering step, for the phyletic distribution page of gpcrevolution.org.

One count table per step, in the layout of the page's final table (taxID, then one <family>_count column per
family and control):

    clade        members of the extracted clades (02_clade_extraction/clade_membership.tsv)
    core         with a single GPCR-like 7TM core (03_topology_filter/core_membership.tsv)
    topology     with a topology score of at least 0.8 (03_topology_filter/topology_membership.tsv)
    overlap      after the overlap assignment (04_overlap_assignment/final_membership.tsv)
    final        after the Dikarya Class A sequences moved to GPR1
                 (10_motif_analysis/phylogeny_reassignment/final_membership_reassigned.tsv)

A sequence's genome is the OX= taxID of its UniProt header (02_clade_extraction/clade_sequences/). The controls
are not filtered, so their columns are the same at every step. Rows and columns follow the page's current final
table (--template), so the page reads every step the same way; genomes that have sequences only at an earlier step
are added as rows of that step. The final and overlap tables are checked against the published count tables.

Usage:
    python build_step_count_tables.py --template <public/phyletic-distribution/tax_counts_per_family_final.tsv> \
        --out <public/phyletic-distribution>
"""
import argparse
import csv
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from gpcr_paths import CONTROLS, FAMILIES, open_text, stage  # noqa: E402

STEPS = [
    ("clade", ("02_clade_extraction", "clade_membership.tsv")),
    ("core", ("03_topology_filter", "core_membership.tsv")),
    ("topology", ("03_topology_filter", "topology_membership.tsv")),
    ("overlap", ("04_overlap_assignment", "final_membership.tsv")),
    ("final", ("10_motif_analysis", "phylogeny_reassignment", "final_membership_reassigned.tsv")),
]
# the published tables the last two steps must reproduce
CHECKS = {
    "overlap": ("04_overlap_assignment", "tax_counts_per_family.tsv"),
    "final": ("10_motif_analysis", "phylogeny_reassignment", "tax_counts_per_family_reassigned.tsv"),
}
OX = re.compile(r"^>(?:sp|tr)\|([^|]+)\|.*?\sOX=(\d+)")


def wide(path):
    with open_text(path) as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    table = {name: set() for name in rows[0]}
    for row in rows[1:]:
        for name, value in zip(rows[0], row):
            if value.strip():
                table[name].add(value.strip())
    return table


def read_counts(path):
    with open_text(path) as handle:
        return {r["taxID"]: r for r in csv.DictReader(handle, delimiter="\t")}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--template", required=True, help="the page's current final count table")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    with open(args.template, encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    template = read_counts(args.template)
    taxids = list(template)

    taxid_of = {}
    for fam in FAMILIES:
        with open_text(stage("02_clade_extraction", "clade_sequences", f"{fam}.fasta")) as handle:
            for line in handle:
                m = OX.match(line)
                if m:
                    taxid_of[m.group(1)] = m.group(2)

    for step, parts in STEPS:
        members = wide(stage(*parts))
        counts = defaultdict(lambda: defaultdict(int))
        for fam in FAMILIES:
            for acc in members[fam]:
                counts[taxid_of[acc]][fam] += 1
        # The final table leaves out genomes without a single final member or control; at the earlier steps
        # some of them have sequences, and they are added at the end (their controls are 0).
        extra = sorted(set(counts) - set(template), key=int)
        if extra:
            print(f"{step}: {len(extra)} genomes not in the final table added: {', '.join(extra)}")
        if step in CHECKS:
            published = read_counts(stage(*CHECKS[step]))
            bad = [(t, f) for t in taxids for f in FAMILIES
                   if counts[t][f] != int(published.get(t, {}).get(f + "_count", 0) or 0)]
            if bad:
                sys.exit(f"{step}: {len(bad)} cells differ from {CHECKS[step][-1]}, e.g. {bad[:3]}")
            print(f"{step}: equals {CHECKS[step][-1]}")
        path = os.path.join(args.out, f"tax_counts_per_family_{step}.tsv")
        with open(path, "w", encoding="utf-8", newline="") as handle:
            handle.write("\t".join(header) + "\n")
            for t in taxids + extra:
                row = [t]
                for col in header[1:]:
                    key = col[:-len("_count")]
                    row.append(str(counts[t][key]) if key in FAMILIES else template.get(t, {}).get(col, "0"))
                handle.write("\t".join(row) + "\n")
        total = sum(sum(c.values()) for c in counts.values())
        print(f"{step}: {total} family memberships -> {path}")
    assert set(CONTROLS) <= {c[:-len("_count")] for c in header}, "control columns missing from the template"


if __name__ == "__main__":
    main()
