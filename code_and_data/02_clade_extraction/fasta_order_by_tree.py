#!/usr/bin/env python3
"""Rewrite a FASTA file in the leaf order of a Newick tree."""

from __future__ import annotations

import argparse
from pathlib import Path

from Bio import SeqIO
from ete3 import Tree


def reorder_fasta(tree_file: Path, fasta_file: Path, output: Path) -> None:
    tree = Tree(str(tree_file), format=1)
    records = {record.id: record for record in SeqIO.parse(str(fasta_file), "fasta")}
    leaf_order = [leaf.name for leaf in tree.iter_leaves()]
    missing = [identifier for identifier in leaf_order if identifier not in records]
    if missing:
        preview = ", ".join(missing[:5])
        raise ValueError(f"{len(missing)} tree leaves are absent from the FASTA: {preview}")
    SeqIO.write((records[identifier] for identifier in leaf_order), str(output), "fasta")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output = args.output or args.fasta.with_name(f"{args.fasta.stem}_treeordered.fasta")
    reorder_fasta(args.tree, args.fasta, output)
    print(f"Reordered FASTA written to {output}")


if __name__ == "__main__":
    main()

