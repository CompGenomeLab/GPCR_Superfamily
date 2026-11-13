#!/usr/bin/env python3
"""Root and orient an exploratory tree using HMMER or MMseqs2 hit scores."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from ete3 import Tree


M8_COLUMNS = [
    "query",
    "target",
    "theader",
    "pident",
    "fident",
    "nident",
    "alnlen",
    "qlen",
    "tlen",
    "mismatch",
    "raw",
    "bits",
    "qcov",
    "tcov",
    "evalue",
]


def import_hmm_data(path: Path, limit: int) -> pd.DataFrame:
    rows: list[tuple[str, str, str]] = []
    in_hit_table = False
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if "E-value  score  bias    E-value  score  bias    exp  N  Sequence" in line:
                in_hit_table = True
                continue
            if not in_hit_table:
                continue
            if "Domain annotation for each sequence:" in line:
                break
            fields = line.split()
            if len(fields) < 9 or fields[0].startswith("-") or "threshold" in line:
                continue
            try:
                float(fields[0])
                float(fields[1])
            except ValueError:
                continue
            rows.append((fields[0], fields[1], fields[8]))
            if len(rows) >= limit:
                break
    frame = pd.DataFrame(rows, columns=["E-value", "Score", "ID"])
    return numeric_scores(frame)


def import_m8_data(path: Path, limit: int) -> pd.DataFrame:
    raw = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=M8_COLUMNS,
        usecols=["theader", "bits", "evalue"],
        nrows=limit,
        dtype=str,
    )
    frame = pd.DataFrame(
        {
            "E-value": raw["evalue"],
            "Score": raw["bits"],
            "ID": raw["theader"].str.split().str[0],
        }
    )
    return numeric_scores(frame)


def numeric_scores(frame: pd.DataFrame) -> pd.DataFrame:
    frame["Score"] = pd.to_numeric(frame["Score"], errors="coerce")
    frame["E-value"] = pd.to_numeric(frame["E-value"], errors="coerce")
    frame = frame.dropna(subset=["Score", "E-value", "ID"]).copy()
    frame["E-value"] = frame["E-value"].replace(0, np.finfo(float).tiny)
    frame["-logeval"] = -np.log10(frame["E-value"])
    return frame.drop_duplicates(subset="ID", keep="first")


def orient_by_mean_score(tree: Tree, scores: dict[str, float]) -> int:
    swaps = 0
    cached_mean: dict[int, float] = {}

    def mean_score(node: Tree) -> float:
        key = id(node)
        if key not in cached_mean:
            values = [scores[name] for name in node.get_leaf_names() if name in scores]
            cached_mean[key] = float(np.mean(values)) if values else 0.0
        return cached_mean[key]

    for node in tree.traverse("postorder"):
        if len(node.children) == 2:
            left, right = node.children
            if mean_score(right) > mean_score(left):
                node.swap_children()
                swaps += 1
    return swaps


def default_output(tree_path: Path, source_type: str) -> Path:
    suffix = f"_{source_type}ordered.tree"
    if tree_path.suffix == ".tree":
        return tree_path.with_name(f"{tree_path.stem}{suffix}")
    return Path(f"{tree_path}{suffix}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--hmm", type=Path)
    source.add_argument("--m8", type=Path)
    parser.add_argument("--tree", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--seq-limit", type=int, default=500000)
    args = parser.parse_args()

    source_type = "hmm" if args.hmm else "m8"
    frame = (
        import_hmm_data(args.hmm, args.seq_limit)
        if args.hmm
        else import_m8_data(args.m8, args.seq_limit)
    )
    if frame.empty:
        sys.exit("No usable search-score records were found")

    tree = Tree(str(args.tree), format=1)
    leaves = set(tree.get_leaf_names())
    frame = frame[frame["ID"].isin(leaves)].copy()
    if frame.empty:
        sys.exit("No search-result identifiers match the tree leaves")

    outgroup = frame.sort_values("Score", ascending=True).iloc[0]["ID"]
    tree.set_outgroup(tree & outgroup)
    scores = dict(zip(frame["ID"], frame["Score"]))
    swaps = orient_by_mean_score(tree, scores)
    output = args.output or default_output(args.tree, source_type)
    tree.write(format=1, outfile=str(output))
    print(f"Matched {len(frame)} leaves; outgroup={outgroup}; swaps={swaps}")
    print(f"Ordered tree written to {output}")


if __name__ == "__main__":
    main()

