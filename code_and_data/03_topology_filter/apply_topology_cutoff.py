"""Remove the sequences whose topology score is below the cutoff (0.8).

The topology score of a sequence (score_all_family_topologies.py, column `mean_jaccard`) is
the mean of its seven per-helix Jaccard similarities to the consensus helices of its family's
topology alignment. Sequences scoring below 0.8 were inspected family by family in the
topology viewer (build_all_topology_viewers.py) before the cutoff was applied.

Each family keeps the sequences of its own alignment that pass. The nested-family rule of
stage 02 is then applied to the passing sequences: a sequence that passes in Vomeronasal 2 is
removed from Class C, and one that passes in GPR1 or Vomeronasal 1 from Class A.

Reads  03_topology_filter/topology_helix_scores.tsv
Writes 03_topology_filter/topology_membership.tsv

Usage:
    python apply_topology_cutoff.py [--cutoff 0.8]
"""
import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import FAMILIES, stage  # noqa: E402

NESTED_MODEL_EXCLUSIONS = {"classA": ["GPR1", "Vomeronasal1"], "classC": ["Vomeronasal2"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cutoff", type=float, default=0.8)
    args = parser.parse_args()

    aligned = defaultdict(set)
    passed = defaultdict(set)
    with open(stage("03_topology_filter", "topology_helix_scores.tsv"), encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            aligned[row["family"]].add(row["accession"])
            if float(row["mean_jaccard"]) >= args.cutoff:
                passed[row["family"]].add(row["accession"])

    kept = {f: set(passed[f]) for f in FAMILIES}
    for broad, specific in NESTED_MODEL_EXCLUSIONS.items():
        for family in specific:
            kept[broad] -= passed[family]
    ordered = {f: sorted(kept[f]) for f in FAMILIES}
    height = max(len(v) for v in ordered.values())
    with open(stage("03_topology_filter", "topology_membership.tsv"), "w", newline="",
              encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(FAMILIES)
        for i in range(height):
            writer.writerow([ordered[f][i] if i < len(ordered[f]) else "" for f in FAMILIES])
    for family in FAMILIES:
        print("%-14s single core %8d  topology >= %.2f %8d" % (family, len(aligned[family]),
                                                              args.cutoff, len(kept[family])))
    print("unique sequences: %d" % len(set().union(*kept.values())))


if __name__ == "__main__":
    main()
