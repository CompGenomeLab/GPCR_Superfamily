#!/usr/bin/env python3
"""Plot each transmembrane helix's Jaccard score along the overall ranking.

Sequences are ordered once, by how far they sit from the origin of the space
whose axes are the per-helix Jaccard scores, best to worst. Each
helix is then drawn as its own line in that fixed order, so the seven can be read
against each other: a fall that appears in every line at the same place is a
general loss of agreement, while a fall in one line alone is that helix slipping
out of register while the rest of the bundle holds.

Every value is plotted as it is. The per-sequence scores are quantised and jump
about, so the lines are dense; --window applies a rolling median instead when a
smoothed shape is wanted.

A helix's window runs to the midpoint of the loop on either side, so a helix that
has slipped is charged inside its own window rather than its neighbour's.
"""

import argparse
import csv
import os
import re
import sys

import numpy as np
import plotly.graph_objects as go
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import open_text, stage  # noqa: E402


# The reference categorical order, which clears every adjacent-pair gate.
HELIX_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
                "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
AVERAGE_COLOR = "#111111"
RATIO_COLOR = "#e34948"
DISTANCE_COLOR = "#4a3aa7"
RANGE_RE = re.compile(r"[/:]\d+-\d+$")


def parse_fasta(path):
    header, chunks = None, []
    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:].strip(), []
            elif header is not None:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def accession_of(header):
    token = RANGE_RE.sub("", header.split()[0] if header.split() else header)
    fields = token.split("|")
    return fields[1] if len(fields) >= 3 else token


def load_ratios(path, headers):
    """The pipeline's median normalized search score, one per sequence."""
    table = {}
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="	"):
            try:
                table[row["accession"]] = float(row["median_ratio"])
            except (KeyError, ValueError):
                continue
    if not table:
        sys.exit(f"No median_ratio column found in {path}")
    values = np.array([table.get(accession_of(h), np.nan) for h in headers])
    found = int(np.isfinite(values).sum())
    print(f"median ratio found for {found} of {len(headers)} sequences")
    return values


def helix_features(records, min_occupancy, min_block):
    """Per-sequence Jaccard against each consensus helix, plus the overall score."""
    n = len(records)
    widths = {len(seq) for _, seq in records}
    if len(widths) != 1:
        sys.exit(f"Sequences are not all the same width: {sorted(widths)[:5]}...")

    matrix = np.frombuffer("".join(seq for _, seq in records).encode(),
                           dtype="S1").reshape(n, widths.pop())
    occupancy = (matrix != b"-").sum(axis=0) / n
    core = np.where(occupancy >= min_occupancy)[0]
    if core.size == 0:
        sys.exit(f"No columns reach {min_occupancy:.0%} occupancy")

    sub = matrix[:, core]
    consensus = np.empty(sub.shape[1], dtype="S1")
    for j in range(sub.shape[1]):
        column = sub[:, j]
        column = column[column != b"-"]
        values, counts = np.unique(column, return_counts=True)
        consensus[j] = values[counts.argmax()]

    helix_consensus = (consensus == b"M") | (consensus == b"W")
    blocks, start = [], None
    for index, value in enumerate(helix_consensus):
        if value and start is None:
            start = index
        elif not value and start is not None:
            blocks.append((start, index))
            start = None
    if start is not None:
        blocks.append((start, len(helix_consensus)))
    blocks = [b for b in blocks if b[1] - b[0] >= min_block]
    if not blocks:
        sys.exit("The consensus has no helix block wide enough to score")

    edges = [0] + [(blocks[j][1] + blocks[j + 1][0]) // 2
                   for j in range(len(blocks) - 1)] + [len(helix_consensus)]
    windows = [(edges[j], edges[j + 1]) for j in range(len(blocks))]

    is_helix = (sub == b"M") | (sub == b"W")
    present = sub != b"-"

    # M and W are not the same state. M crosses the membrane inside to outside and W
    # crosses it outside to inside, so the consensus bundle alternates M, W, M, W ...
    # and the direction of a helix is part of what it means to match that consensus. A
    # column therefore counts as agreement only when the sequence carries the same
    # membrane state as the consensus; a helix running the wrong way is charged as
    # off-register, exactly like a helix sitting in a loop, and is counted once.
    agrees = is_helix & (sub == consensus[None, :])

    def jaccard(a, b):
        both = (agrees[:, a:b] & helix_consensus[None, a:b]).sum(axis=1)
        off = (is_helix[:, a:b]
               & ~(agrees[:, a:b] & helix_consensus[None, a:b])).sum(axis=1)
        missed = ((~is_helix[:, a:b]) & present[:, a:b]
                  & helix_consensus[None, a:b]).sum(axis=1)
        return both / np.maximum(1, both + off + missed)

    features = np.column_stack([jaccard(a, b) for a, b in windows])
    overall = jaccard(0, sub.shape[1])
    return features, overall, blocks, windows, core.size


def rolling_median(values, window):
    if window <= 1:
        return values
    half = window // 2
    padded = np.pad(values, (half, half), mode="edge")
    strided = np.lib.stride_tricks.sliding_window_view(padded, 2 * half + 1)
    return np.median(strided, axis=1)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--topology-alignment", required=True)
    parser.add_argument("--out", required=True, help="output HTML file")
    parser.add_argument("--min-occupancy", type=float, default=0.30)
    parser.add_argument("--min-block", type=int, default=8)
    parser.add_argument("--window", type=int, default=1,
                        help="rolling median width in sequences; 1 (the default) plots "
                             "every value as it is, with no smoothing")
    parser.add_argument("--ratios", default=None,
                        help="acc_taxid_map_final_with_ratio.tsv, to draw the median "
                             "normalized search score alongside the helices")
    parser.add_argument("--title", default=None)
    args = parser.parse_args()

    records = list(parse_fasta(args.topology_alignment))
    if not records:
        sys.exit(f"No sequences found in {args.topology_alignment}")

    features, overall, blocks, windows, n_core = helix_features(
        records, args.min_occupancy, args.min_block)
    n_helix = features.shape[1]
    print(f"{len(records)} sequences, {n_core} core columns, {n_helix} consensus helices")
    for j, ((a, b), (wa, wb)) in enumerate(zip(blocks, windows), 1):
        print(f"  TM{j}: helix columns {a}-{b - 1} ({b - a} wide), window {wa}-{wb - 1}")

    headers = [h for h, _ in records]
    ratios = load_ratios(args.ratios, headers) if args.ratios else None

    # Each sequence is a point in as many dimensions as there are helices, its
    # coordinates being its Jaccard against each one. The worst a sequence can
    # do is sit at the origin, matching no helix at all; the best is to sit at
    # (1, 1, ... 1). Distance from the origin therefore measures the whole
    # bundle at once, on its own scale from 0 to sqrt(K).
    distance = np.sqrt((features ** 2).sum(axis=1))
    print(f"\ndistance from the origin in {n_helix} dimensions: "
          f"median {np.median(distance):.3f} of a possible {np.sqrt(n_helix):.3f}, "
          f"min {distance.min():.3f}")

    order = np.argsort(-distance, kind="stable")
    x = np.arange(1, len(order) + 1)

    figure = go.Figure()
    for j in range(n_helix):
        figure.add_trace(go.Scattergl(
            x=x, y=rolling_median(features[order, j], args.window),
            mode="lines", name=f"TM{j + 1}",
            line=dict(color=HELIX_COLORS[j % len(HELIX_COLORS)], width=2),
            hovertemplate=f"TM{j + 1}" " %{y:.3f} at rank %{x}<extra></extra>",
        ))
    figure.add_trace(go.Scattergl(
        x=x, y=rolling_median(distance[order], args.window),
        mode="lines", name="distance from the origin",
        line=dict(color=DISTANCE_COLOR, width=2),
        hovertemplate="distance %{y:.3f} of "
                      f"{np.sqrt(n_helix):.3f}" " at rank %{x}<extra></extra>",
    ))
    figure.add_trace(go.Scattergl(
        x=x, y=rolling_median(overall[order], args.window),
        mode="lines", name="single Jaccard over all helix columns",
        line=dict(color=AVERAGE_COLOR, width=2, dash="dash"),
        hovertemplate="overall %{y:.3f} at rank %{x}<extra></extra>",
    ))
    if ratios is not None:
        # Same axis: both are ratios on the same scale, though the search score
        # can sit a little above 1. Plotted in the same order as the helices, so
        # a sequence sits at the same x in every line.
        figure.add_trace(go.Scattergl(
            x=x, y=ratios[order], mode="lines",
            name="median normalized search score",
            line=dict(color=RATIO_COLOR, width=2, dash="dot"),
            hovertemplate="median ratio %{y:.3f} at rank %{x}<extra></extra>",
        ))

    figure.update_layout(
        title=dict(text=f"{args.title or os.path.basename(args.topology_alignment)}"
                        f"<br><sub>each helix scored on its own, sequences ordered once "
                        f"by the overall score; "
                        + ("every value plotted, no smoothing"
                           if args.window <= 1 else
                           f"rolling median over {args.window} sequences")
                        + "</sub>"),
        xaxis=dict(title="sequences, ordered by distance from the origin, best to worst",
                   gridcolor="#eceae5"),
        yaxis=dict(title="Jaccard against the consensus helix, and the search score",
                   gridcolor="#eceae5", range=[0, np.sqrt(n_helix) * 1.03]),
        template="plotly_white", height=650,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0),
        margin=dict(l=80, r=40, t=120, b=70), hoverlabel=dict(bgcolor="#ffffff"),
        hovermode="x unified",
    )
    figure.write_html(args.out, include_plotlyjs="cdn")
    print(f"Wrote {args.out} ({os.path.getsize(args.out) / 1e6:.1f} MB)")

    # A terminal readout of the same curves, so the shape is legible without the plot.
    print("\nvalues at deciles of the ranking"
          + (" (raw)" if args.window <= 1
             else f" (rolling median, window {args.window})") + ":")
    print(f"{'rank':>8s} {'dist':>8s} {'overall':>8s} " +
          " ".join(f"{'TM' + str(j + 1):>7s}" for j in range(n_helix)))
    smooth = [rolling_median(features[order, j], args.window) for j in range(n_helix)]
    smooth_overall = rolling_median(overall[order], args.window)
    for q in list(np.linspace(0, 0.9, 10)) + [0.95, 0.98, 0.99, 0.999]:
        i = min(len(order) - 1, int(q * len(order)))
        print(f"{i + 1:8d} {rolling_median(distance[order], args.window)[i]:8.3f} "
              f"{smooth_overall[i]:8.3f} " +
              " ".join(f"{smooth[j][i]:7.3f}" for j in range(n_helix)))


if __name__ == "__main__":
    main()
