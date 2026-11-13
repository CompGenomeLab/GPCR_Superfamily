#!/usr/bin/env python3
"""Build the data of the topology-filter page of gpcrevolution.org (transparency/topology-filter).

For every topology alignment of the archive (03_topology_filter/topology_alignments/, the single-core
sequences of each family, with the nested families already assigned: Class C without Vomeronasal 2,
Class A without GPR1 and Vomeronasal 1) this writes the payload of the topology viewer, built by
../build_topology_viewer.py exactly as for the archive's topology_viewer.html, as <out>/<family>.json.gz,
and writes <out>/index.json with each family's counts before and after the 0.8 topology-score filter
and the distribution of its scores (topology_helix_scores.tsv, column mean_jaccard).

The page scores every sequence over the columns occupied in at least 30% of the sequences, as the
filter did, so the sequences it shows as removed are the ones the filter removed.

Usage:
    python build_topology_filter_data.py --out <GPCRevolution>/public/topology-filter
"""
import argparse
import base64
import csv
import glob
import json
import os
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[1]))
import build_topology_viewer as viewer  # noqa: E402
from gpcr_paths import stage  # noqa: E402

SUFFIX = "_genes_7TMfilter_db_FAMSA_topology_converted_gt0.1.fasta"
CUTOFF = 0.8
BIN = 0.02
MAX_GAP_PAYLOAD = 0.90
DISPLAY = {
    "classA": "Class A", "classB1": "Class B1", "classB2": "Class B2", "classC": "Class C",
    "classF": "Class F", "classT": "Taste-2", "FSLB": "FSL", "Vomeronasal1": "Vomeronasal 1",
    "Vomeronasal2": "Vomeronasal 2",
}


def family_of(path):
    name = os.path.basename(path)
    name = name[:-3] if name.endswith(".gz") else name
    return name[:-len(SUFFIX)] if name.endswith(SUFFIX) else os.path.splitext(name)[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    folder = stage("03_topology_filter", "topology_alignments")
    paths = sorted(glob.glob(os.path.join(str(folder), "*.fasta.gz")) +
                   glob.glob(os.path.join(str(folder), "*.fasta")), key=lambda p: family_of(p).lower())

    scores = defaultdict(list)
    with open(stage("03_topology_filter", "topology_helix_scores.tsv"), encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            scores[row["family"]].append(float(row["mean_jaccard"]))

    wanted = {viewer.short_id(h) for p in paths for h, _ in viewer.parse_fasta(p)}
    print(f"reading descriptions for {len(wanted)} accessions ...", flush=True)
    descriptions = viewer.load_descriptions([str(stage("01_homology_search"))], wanted)

    index = []
    for path in paths:
        family = family_of(path)
        print(f"[{family}]", flush=True)
        blob = viewer.build_payload(path, MAX_GAP_PAYLOAD, None, descriptions)
        with open(os.path.join(args.out, family + ".json.gz"), "wb") as fh:
            fh.write(base64.b64decode(blob))
        values = scores[family]
        bins = [0] * int(round(1 / BIN))
        for v in values:
            bins[min(len(bins) - 1, int(v / BIN))] += 1
        kept = sum(1 for v in values if v >= CUTOFF)
        index.append({
            "family": family,
            "display": DISPLAY.get(family, family),
            "single_core": len(values),
            "kept": kept,
            "removed": len(values) - kept,
            "bin_width": BIN,
            "histogram": bins,
        })
    with open(os.path.join(args.out, "index.json"), "w", encoding="utf-8") as fh:
        json.dump({"cutoff": CUTOFF, "families": index}, fh, separators=(",", ":"))
    print(f"wrote {len(index)} families to {args.out}")


if __name__ == "__main__":
    main()
