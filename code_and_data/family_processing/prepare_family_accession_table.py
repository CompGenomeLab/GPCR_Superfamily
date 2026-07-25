#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Prepare a family-by-accession TSV from extracted family FASTA files.

This script is intended to run after `identify_gene_clade.py` extracts family
alignments into `extracted_clades/alignments`. Its main output is the accession
TSV consumed by `create_taxID_family_count_table.py`.
"""

from __future__ import annotations

from itertools import combinations
import os
from pathlib import Path
from typing import Dict, Iterable, Optional, Sequence, Set

import pandas as pd


FASTA_EXTENSIONS = {".fa", ".fasta", ".fna", ".faa", ".ffn"}


def read_short_headers(fasta_path: Path) -> Set[str]:
    """Return accession-like headers, using the second pipe-delimited field when present."""
    headers: Set[str] = set()
    with fasta_path.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith(">"):
                full_header = line[1:].strip()
                parts = full_header.split("|")
                headers.add((parts[1] if len(parts) > 1 else full_header).strip())
    return headers


def load_clusters(cluster_tsv_path: Optional[Path]) -> Dict[str, Set[str]]:
    """
    Load representative-to-member mappings from a cluster TSV.

    Expected format: column 1 is representative accession, column 2 contains
    member accessions separated by commas, semicolons, or whitespace.
    """
    if not cluster_tsv_path:
        return {}

    df = pd.read_csv(
        cluster_tsv_path,
        sep="\t",
        header=None,
        usecols=[0, 1],
        names=["rep", "members"],
        dtype=str,
    )
    df["members"] = df["members"].fillna("")

    exploded = df.assign(member=df["members"].str.split(r"[,\s;]+", regex=True)).explode("member")
    exploded["member"] = exploded["member"].astype(str).str.strip()
    exploded = exploded[exploded["member"] != ""]

    reps = df[["rep"]].copy()
    reps["member"] = reps["rep"]

    all_rows = pd.concat([exploded[["rep", "member"]], reps], ignore_index=True)
    return all_rows.groupby("rep")["member"].apply(set).to_dict()


def expand_representatives(headers: Set[str], family: str, clusters: Dict[str, Set[str]], expand_families: Set[str]) -> Set[str]:
    """Expand selected representative accessions to their full cluster members."""
    if family not in expand_families or not clusters:
        return headers

    expanded: Set[str] = set()
    for header in headers:
        expanded.update(clusters.get(header, {header}))
    return expanded


def collect_family_accessions(
    fasta_dir: Path,
    cluster_tsv_path: Optional[Path] = None,
    expand_families: Iterable[str] = ("classA",),
) -> Dict[str, Set[str]]:
    """Collect accessions from extracted family FASTA files."""
    clusters = load_clusters(cluster_tsv_path)
    expand_family_set = set(expand_families)
    family_to_accessions: Dict[str, Set[str]] = {}

    for fasta_path in sorted(fasta_dir.iterdir()):
        if not fasta_path.is_file() or fasta_path.suffix.lower() not in FASTA_EXTENSIONS:
            continue

        family = fasta_path.stem.split("_", 1)[0]
        headers = read_short_headers(fasta_path)
        headers = expand_representatives(headers, family, clusters, expand_family_set)
        family_to_accessions.setdefault(family, set()).update(headers)

    return family_to_accessions


def apply_exclusions(family_to_accessions: Dict[str, Set[str]], exclude_map: Dict[str, Sequence[str]]) -> Dict[str, Set[str]]:
    """Remove accessions from selected families when they also belong to excluded families."""
    cleaned = {family: set(accessions) for family, accessions in family_to_accessions.items()}

    for family, excluded_families in exclude_map.items():
        if family not in cleaned:
            continue

        to_remove: Set[str] = set()
        for excluded_family in excluded_families:
            to_remove.update(cleaned.get(excluded_family, set()))
        cleaned[family] -= to_remove

    return cleaned


def write_accession_tsv(family_to_accessions: Dict[str, Set[str]], output_tsv: Path) -> None:
    """Write one TSV column per family, with accessions as rows."""
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    family_lists = {
        family: sorted(accessions)
        for family, accessions in sorted(
            family_to_accessions.items(), key=lambda item: item[0].casefold()
        )
        if accessions
    }
    if not family_lists:
        raise ValueError("No accessions found to write.")

    max_len = max(len(accessions) for accessions in family_lists.values())
    for accessions in family_lists.values():
        accessions.extend([""] * (max_len - len(accessions)))

    pd.DataFrame(family_lists).to_csv(output_tsv, sep="\t", index=False)


def write_family_accession_txts(family_to_accessions: Dict[str, Set[str]], output_dir: Optional[Path]) -> None:
    """Optionally write one accession text file per family."""
    if output_dir is None:
        return

    output_dir.mkdir(parents=True, exist_ok=True)
    for family, accessions in sorted(family_to_accessions.items()):
        out_path = output_dir / f"{family}_accessions.txt"
        out_path.write_text("\n".join(sorted(accessions)) + "\n", encoding="utf-8")


def compute_pairwise_overlaps(family_to_accessions: Dict[str, Set[str]]) -> list[dict]:
    """Return non-empty pairwise overlaps for optional reporting/plotting."""
    overlaps = []
    families = sorted(family_to_accessions)
    for family_a, family_b in combinations(families, 2):
        shared = family_to_accessions[family_a] & family_to_accessions[family_b]
        if shared:
            overlaps.append(
                {
                    "family_a": family_a,
                    "family_b": family_b,
                    "count": len(shared),
                    "accessions": shared,
                    "size_a": len(family_to_accessions[family_a]),
                    "size_b": len(family_to_accessions[family_b]),
                }
            )
    return overlaps


def write_overlap_plot(
    overlaps: list[dict],
    output_html: Optional[Path],
    title: str = "Pairwise FASTA Header Overlaps",
) -> None:
    """Optionally write a Plotly HTML overlap plot."""
    if output_html is None:
        return

    import plotly.graph_objects as go

    output_html.parent.mkdir(parents=True, exist_ok=True)

    if not overlaps:
        fig = go.Figure()
        fig.update_layout(
            title=f"{title}: no non-zero pairwise overlaps found",
            xaxis_title="Family Pair",
            yaxis_title="# Overlaps",
            template="plotly_white",
        )
        fig.write_html(output_html, include_plotlyjs="cdn")
        return

    labels = [f"{row['family_a']} / {row['family_b']}" for row in overlaps]
    counts = [row["count"] for row in overlaps]
    hover = []
    for row in overlaps:
        preview = "<br>".join(sorted(row["accessions"])[:10])
        more = "" if row["count"] <= 10 else f"<br>... and {row['count'] - 10} more"
        hover.append(
            f"{row['family_a']} / {row['family_b']}<br>"
            f"Overlap: {row['count']}<br>"
            f"Coverage: {row['family_a']} {row['count']}/{row['size_a']} | "
            f"{row['family_b']} {row['count']}/{row['size_b']}<br>"
            f"{preview}{more}"
        )

    fig = go.Figure(
        [go.Bar(x=labels, y=counts, text=counts, textposition="auto", hovertext=hover, hovertemplate="%{hovertext}<extra></extra>")]
    )
    fig.update_layout(
        title=title,
        xaxis_title="Family Pair",
        yaxis_title="# Overlapping Accessions",
        xaxis_tickangle=45,
        template="plotly_white",
        margin=dict(l=60, r=30, t=70, b=80),
    )
    fig.write_html(output_html, include_plotlyjs="cdn")


def parse_exclude_map(values: Sequence[str]) -> Dict[str, list[str]]:
    """
    Parse exclusions passed as FAMILY:EXCLUDED1,EXCLUDED2.

    Example items:
      "classA:GPR1,Vomeronasal1"
      "classC:Vomeronasal2"
    """
    exclude_map: Dict[str, list[str]] = {}
    for value in values:
        family, _, excluded = value.partition(":")
        family = family.strip()
        excluded_families = [item.strip() for item in excluded.split(",") if item.strip()]
        if family and excluded_families:
            exclude_map[family] = excluded_families
    return exclude_map


def prepare_family_accessions(
    fasta_dir: Path,
    cluster_tsv: Optional[Path],
    expand_families: Sequence[str],
    exclude_items: Sequence[str],
) -> Dict[str, Set[str]]:
    """Load family accessions and reproduce the study's documented curation."""
    family_to_accessions = collect_family_accessions(
        fasta_dir=fasta_dir,
        cluster_tsv_path=cluster_tsv,
        expand_families=expand_families,
    )
    family_to_accessions = apply_exclusions(
        family_to_accessions, parse_exclude_map(exclude_items)
    )

    # LRLA was retained as a search/control dataset but was not propagated into
    # the 23-family accession table or the final 18 family databases.
    family_to_accessions.pop("LRLA", None)

    additions = Path(__file__).with_name("curated_accession_additions.tsv")
    with additions.open(encoding="utf-8") as handle:
        next(handle)
        for line in handle:
            family, accession, _note = line.rstrip("\n").split("\t", 2)
            family_to_accessions.setdefault(family, set()).add(accession)

    return family_to_accessions


def generate_accession_outputs(
    family_to_accessions: Dict[str, Set[str]],
    output_tsv: Path,
    accession_txt_dir: Optional[Path],
) -> None:
    """Generate TSV and optional per-family accession text outputs."""
    write_accession_tsv(family_to_accessions, output_tsv)
    write_family_accession_txts(family_to_accessions, accession_txt_dir)


def generate_overlap_outputs(
    family_to_accessions: Dict[str, Set[str]],
    plot_html: Optional[Path],
) -> int:
    """Calculate overlaps and optionally write overlap plot. Returns overlap count."""
    if plot_html is None:
        return 0

    overlaps = compute_pairwise_overlaps(family_to_accessions)
    write_overlap_plot(overlaps, plot_html)
    return len(overlaps)


def main(
    fasta_dir: Path,
    cluster_tsv: Optional[Path],
    output_tsv: Path,
    accession_txt_dir: Optional[Path],
    plot_html: Optional[Path],
    expand_families: Sequence[str],
    exclude_items: Sequence[str],
    generate_tsv_outputs: bool = True,
    generate_overlap_plot: bool = True,
) -> None:
    family_to_accessions = prepare_family_accessions(
        fasta_dir=fasta_dir,
        cluster_tsv=cluster_tsv,
        expand_families=expand_families,
        exclude_items=exclude_items,
    )

    if generate_tsv_outputs:
        generate_accession_outputs(
            family_to_accessions=family_to_accessions,
            output_tsv=output_tsv,
            accession_txt_dir=accession_txt_dir,
        )

    overlap_count = generate_overlap_outputs(
        family_to_accessions=family_to_accessions,
        plot_html=plot_html if generate_overlap_plot else None,
    )

    if generate_tsv_outputs:
        print(f"Wrote accession TSV: {output_tsv}")
    print(f"Families: {len(family_to_accessions)}")
    if generate_overlap_plot:
        print(f"Pairwise overlaps after exclusions: {overlap_count}")


if __name__ == "__main__":
    # ================= USER VARIABLES ================= #
    BASE_DIR = Path(os.environ.get("GPCR_ANALYSIS_DIR", os.environ.get("GPCR_ZENODO_DIR", "/path/to/GPCR_Superfamily")))
    SEARCH_RESULTS_DIR = Path(
        os.environ.get("GPCR_SEARCH_RESULTS_DIR", BASE_DIR / "search_results")
    )

    FAMILY_DIR = Path(os.environ.get("GPCR_FAMILY_DIR", BASE_DIR / "families"))
    FASTA_DIR = Path(os.environ.get(
        "GPCR_EXTRACTED_CLADE_DIR", BASE_DIR / "extracted_clades"
    )) / "alignments"
    TABLE_DIR = Path(os.environ.get("GPCR_FAMILY_TABLE_DIR", FAMILY_DIR / "tables"))
    CLUSTER_TSV = Path(os.environ.get(
        "GPCR_CLASSA_CLUSTER_TSV",
        SEARCH_RESULTS_DIR / "classA" / "classA_hmm_E1000_db_clust0.8.tsv",
    ))
    OUTPUT_TSV = TABLE_DIR / "accessions_by_family_prefilter_May22.tsv"

    # Set to None to disable optional outputs
    ACCESSION_TXT_DIR = Path(os.environ.get("GPCR_ACCESSION_TXT_DIR", FAMILY_DIR / "accessions"))
    PLOT_HTML = TABLE_DIR / "overlap_plot_prefilter_May22.html"

    # Only Class A was reduced to 80% representatives during the search stage.
    # Olfactory was clustered later, independently, during family processing.
    EXPAND_FAMILIES = ["classA"]
    EXCLUDE_ITEMS = [
        "classA:GPR1,Vomeronasal1",
        "classC:Vomeronasal2",
    ]

    # Modular step toggles
    GENERATE_TSV_OUTPUTS = True
    GENERATE_OVERLAP_PLOT = True
    # ================================================== #

    main(
        fasta_dir=FASTA_DIR,
        cluster_tsv=CLUSTER_TSV,
        output_tsv=OUTPUT_TSV,
        accession_txt_dir=ACCESSION_TXT_DIR,
        plot_html=PLOT_HTML,
        expand_families=EXPAND_FAMILIES,
        exclude_items=EXCLUDE_ITEMS,
        generate_tsv_outputs=GENERATE_TSV_OUTPUTS,
        generate_overlap_plot=GENERATE_OVERLAP_PLOT,
    )
