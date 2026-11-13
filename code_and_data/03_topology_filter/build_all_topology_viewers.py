#!/usr/bin/env python3
"""Build one topology viewer holding every alignment in a folder.

Passes every *.fasta in the folder to build_topology_viewer.py, which puts them
all in a single page with a selector to switch between families, so a family
added later (class A, say) is picked up by simply re-running this. Rows are
ordered by the helix distance score, best first, as in the class C viewer.
"""

import argparse
import glob
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import stage  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SUFFIX = "_genes_7TMfilter_db_FAMSA_topology_converted_gt0.1.fasta"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--folder", default=str(stage("03_topology_filter", "topology_alignments")),
                        help="folder of topology-converted alignments")
    parser.add_argument("--out", default=None,
                        help="output HTML (default: <folder>/topology_viewer.html)")
    parser.add_argument("--search-outputs", default=str(stage("01_homology_search")),
                        help="hmmsearch/m8 outputs the UniProt descriptions are read from, since "
                             "the gt0.1 alignments carry only the IDs")
    parser.add_argument("--if-stale", action="store_true",
                        help="skip the build when the page is newer than every alignment and "
                             "than the builder itself")
    args = parser.parse_args()

    alignments = sorted(glob.glob(os.path.join(args.folder, "*.fasta")) +
                        glob.glob(os.path.join(args.folder, "*.fasta.gz")),
                        key=lambda path: os.path.basename(path).lower())
    if not alignments:
        sys.exit(f"No .fasta files in {args.folder}")
    names = []
    for path in alignments:
        name = os.path.basename(path)
        name = name[:-3] if name.endswith(".gz") else name
        names.append(name[:-len(SUFFIX)] if name.endswith(SUFFIX) else os.path.splitext(name)[0])

    out = args.out or str(stage("03_topology_filter", "topology_viewer.html"))
    if args.if_stale and os.path.exists(out):
        inputs = alignments + [os.path.join(HERE, "build_topology_viewer.py"), os.path.abspath(__file__)]
        if os.path.getmtime(out) > max(os.path.getmtime(path) for path in inputs):
            print(f"{os.path.basename(out)} is up to date with {len(alignments)} alignments")
            return
    subprocess.run([sys.executable, os.path.join(HERE, "build_topology_viewer.py"),
                    "--topology-alignment", *alignments, "--names", *names,
                    "--descriptions", args.search_outputs, "--out", out], check=True)


if __name__ == "__main__":
    main()
