#!/usr/bin/env python3
"""Plot normalized search score and search-result rank in ordered-tree order."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
from ete3 import Tree

from tree_order_by_hmm import import_hmm_data, import_m8_data


def normalize(values: pd.Series) -> pd.Series:
    minimum, maximum = values.min(), values.max()
    if maximum == minimum:
        return pd.Series(np.zeros(len(values)), index=values.index)
    return (values - minimum) / (maximum - minimum)


def smooth(values: pd.Series, window: int, weighting: str) -> pd.Series:
    if window < 1:
        return values
    weights = (
        np.ones(window + 1)
        if weighting == "none"
        else np.arange(1, window + 2, dtype=float)
    )
    if weighting == "reversed":
        weights = weights[::-1]

    result = values.astype(float).copy()
    array = values.to_numpy(dtype=float)
    for index in range(window, len(values) - window):
        left = array[index - window : index + 1]
        right = array[index : index + window + 1]
        left_mean = np.average(left, weights=weights)
        right_mean = np.average(right, weights=weights[::-1])
        result.iloc[index] = (
            left_mean
            if abs(left_mean - array[index]) <= abs(right_mean - array[index])
            else right_mean
        )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--hmm", type=Path)
    source.add_argument("--m8", type=Path)
    parser.add_argument("--tree", required=True, type=Path)
    parser.add_argument("--weighting-type", "--weighting_type", choices=["normal", "reversed", "none"], default="none")
    parser.add_argument("--window-size", "--window_size", type=int, default=50)
    parser.add_argument("--html-out", "--html_out", required=True, type=Path)
    parser.add_argument("--seq-limit", type=int, default=500000)
    args = parser.parse_args()

    frame = (
        import_hmm_data(args.hmm, args.seq_limit)
        if args.hmm
        else import_m8_data(args.m8, args.seq_limit)
    )
    if frame.empty:
        parser.error("No usable search-result records were found")

    # The retained DataFrame index is the zero-based row in the HMMER/MMseqs2
    # hit list. This reproduces the index calculation in the original script.
    frame["Search result index"] = frame.index
    frame["Normalized score"] = normalize(frame["Score"])
    frame["Normalized search index"] = 1.0 - normalize(frame["Search result index"])

    tree_order = Tree(str(args.tree), format=1).get_leaf_names()
    plot_data = frame.drop_duplicates("ID").set_index("ID").reindex(tree_order).reset_index()
    plot_data["Tree position"] = np.arange(len(plot_data))
    plot_data[["Normalized score", "Normalized search index"]] = plot_data[
        ["Normalized score", "Normalized search index"]
    ].fillna(0.0)
    plot_data["Smoothed normalized score"] = smooth(
        plot_data["Normalized score"], args.window_size, args.weighting_type
    )
    plot_data["Smoothed normalized search index"] = smooth(
        plot_data["Normalized search index"], args.window_size, args.weighting_type
    )

    # Match the presentation used by the retained line plots.  The underlying
    # values and tree-leaf ordering above are unchanged; these names and the
    # sequence-ID x-axis are the historical display convention.
    metric_labels = {
        "Smoothed normalized score": "Score_smoothed",
        "Smoothed normalized search index": "Normalized_Index_smoothed",
        "Normalized score": "Score",
        "Normalized search index": "Normalized_Index",
    }
    long_data = plot_data.melt(
        id_vars=["ID"],
        value_vars=[
            "Smoothed normalized score",
            "Smoothed normalized search index",
            "Normalized score",
            "Normalized search index",
        ],
        var_name="Metric",
        value_name="Smoothed Value",
    )
    long_data["Metric"] = long_data["Metric"].map(metric_labels)
    long_data = long_data.rename(columns={"ID": "Sequence IDs"})

    figure = px.line(
        long_data,
        x="Sequence IDs",
        y="Smoothed Value",
        color="Metric",
        markers=True,
        title="Smoothed Line Plot of Multiple Metrics",
    )
    figure.update_xaxes(showticklabels=False)
    figure.update_layout(legend_title_text="Metrics")
    figure.write_html(args.html_out)
    print(f"Interactive plot written to {args.html_out}")


if __name__ == "__main__":
    main()
