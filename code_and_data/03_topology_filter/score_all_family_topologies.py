#!/usr/bin/env python3
"""Score every sequence of every family's topology alignment and write one table.

The score is the helix distance score the topology viewer shows, computed here
offline so that later steps can read it from a file instead of a browser. It is
the same quantity `build_topology_viewer.py` computes in `computeScores()` and
the same one `per_helix_features.py` plots, and this script calls that module's
`helix_features` rather than reimplementing it, so the three cannot drift apart.

For each family alignment the consensus topology is taken over the columns that
clear `--min-occupancy`, consensus helix runs of at least `--min-block` columns
become blocks, and each block is given a window reaching to the midpoint of the
loop on either side. A sequence scores a Jaccard against each block inside that
block's own window,

    both / (both + helix off register + consensus helix not covered)

A column counts towards `both` only when the sequence carries the same membrane state as
the consensus. M runs inside to outside and W runs outside to inside, so a helix crossing
the wrong way is charged as off register rather than credited as a match.

and its overall score is the Euclidean distance of those per-helix Jaccards from
the origin, `sqrt(sum of squares)`, on a scale of 0 to `sqrt(K)` for a family
with K helix blocks. Squaring lets the helices that match well dominate the sum, so
one helix missing entirely costs less on this measure than the same deficit spread
thinly across the bundle -- the opposite of what a flagging threshold wants, which is
why `mean_jaccard`, not this distance, is the number to read.

`mean_jaccard` is the mean of those per-helix Jaccards, on 0 to 1, and is the number the
topology viewer shows. It is the one to read and to quote, because it is a similarity in
the same units as the per-helix values rather than a distance on a scale that moves with
the helix count.

`score` and `score_norm` are the Euclidean distance and that distance divided by `sqrt(K)`,
kept for reference; the topology filter (apply_topology_cutoff.py) and the representative
selection (stage 07) both read `mean_jaccard`.

Usage:
    python score_all_family_topologies.py
"""

import argparse
import csv
import os
import sys
import glob

import numpy as np

from per_helix_features import parse_fasta, accession_of, helix_features
from gpcr_paths import stage  # noqa: E402  (path set by per_helix_features)

HERE = os.path.dirname(os.path.abspath(__file__))
SUFFIX = "_genes_7TMfilter_db_FAMSA_topology_converted_gt0.1.fasta"

FIELDS = ["family", "accession", "mean_jaccard", "score", "score_norm", "overall_jaccard",
          "coverage", "residues", "n_helices", "helix_jaccards"]


def family_of(path):
    name = os.path.basename(path)
    if name.endswith(".gz"):
        name = name[:-3]
    return name[:-len(SUFFIX)] if name.endswith(SUFFIX) else os.path.splitext(name)[0]


def score_one(path, min_occupancy, min_block):
    records = list(parse_fasta(path))
    if not records:
        sys.exit(f"No sequences in {path}")
    features, overall, blocks, windows, n_core = helix_features(
        records, min_occupancy, min_block)

    # Coverage needs the core columns again; recomputing the occupancy mask is
    # cheaper than threading it back out of helix_features.
    width = len(records[0][1])
    matrix = np.frombuffer("".join(seq for _, seq in records).encode(),
                           dtype="S1").reshape(len(records), width)
    occupancy = (matrix != b"-").sum(axis=0) / len(records)
    core = np.where(occupancy >= min_occupancy)[0]
    residues = (matrix[:, core] != b"-").sum(axis=1)

    mean_jaccard = features.mean(axis=1)
    distance = np.sqrt((features ** 2).sum(axis=1))
    n_helices = features.shape[1]
    rows = []
    for i, (header, _) in enumerate(records):
        rows.append({
            "family": family_of(path),
            "accession": accession_of(header),
            "mean_jaccard": f"{mean_jaccard[i]:.5f}",
            "score": f"{distance[i]:.5f}",
            "score_norm": f"{distance[i] / np.sqrt(n_helices):.5f}",
            "overall_jaccard": f"{overall[i]:.5f}",
            "coverage": f"{residues[i] / n_core:.5f}",
            "residues": int(residues[i]),
            "n_helices": n_helices,
            "helix_jaccards": ";".join(f"{v:.4f}" for v in features[i]),
        })
    return rows, n_helices, n_core


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--folder", default=str(stage("03_topology_filter", "topology_alignments")),
                        help="folder of topology-converted per-family alignments")
    parser.add_argument("--out", default=str(stage("03_topology_filter", "topology_helix_scores.tsv")))
    parser.add_argument("--min-occupancy", type=float, default=0.30,
                        help="the viewer's default, and what the exclusions were made at")
    parser.add_argument("--min-block", type=int, default=8,
                        help="shortest consensus helix run that counts as a block")
    args = parser.parse_args()

    alignments = sorted(glob.glob(os.path.join(args.folder, "*.fasta")) +
                        glob.glob(os.path.join(args.folder, "*.fasta.gz")),
                        key=lambda p: os.path.basename(p).lower())
    if not alignments:
        sys.exit(f"No .fasta files in {args.folder}")

    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        total = 0
        for path in alignments:
            rows, n_helices, n_core = score_one(path, args.min_occupancy, args.min_block)
            writer.writerows(rows)
            total += len(rows)
            values = np.array([float(r["score"]) for r in rows])
            print(f"{family_of(path):<14}{len(rows):>8} sequences "
                  f"{n_helices:>2} helices {n_core:>5} core columns  "
                  f"score median {np.median(values):.3f} of {np.sqrt(n_helices):.3f}")
    print(f"\nwrote {total} rows to {args.out}")


if __name__ == "__main__":
    main()
