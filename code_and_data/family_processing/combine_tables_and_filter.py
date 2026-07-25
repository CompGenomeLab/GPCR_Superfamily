# -*- coding: utf-8 -*-
"""
Created on Sun Oct 26 15:20:33 2025

@author: selcuk.1
"""

import os

PROJECT_DIR = os.environ.get("GPCR_ANALYSIS_DIR", os.environ.get("GPCR_ZENODO_DIR", "/path/to/GPCR_Superfamily"))

# --- Set your two inputs here (unchanged) ---
TABLE_DIR = os.environ.get("GPCR_FAMILY_TABLE_DIR", os.path.join(PROJECT_DIR, "families", "tables"))
BASE_PATH = os.path.join(TABLE_DIR, "tax_counts_per_family_filtered_Oct24_e0.1_BestEval.tsv")
SECOND_PATH = os.path.join(TABLE_DIR, "tax_counts_per_family_prefilter_Oct24_e0.1_all.tsv")
# --------------------------------------------

# --- Your single output path + optional filtering ---
# This is a derived tax-count table, so it belongs beside the two input tables
# rather than in families/accessions (which contains one-accession-per-line
# inputs for sequence retrieval).
OUT_PATH = os.environ.get(
    "GPCR_CONTROL_FILTERED_TABLE",
    os.path.join(TABLE_DIR, "tax_counts_per_family_filtered_Oct24_e0.1_BestEval_withControls_Filtered.tsv"),
)
FILTER_ALL_FIVE = True   # True -> write filtered rows (all five counts >= MIN_COUNT) to OUT_PATH; False -> write merged to OUT_PATH
MIN_COUNT = 1          # Optional; if None and FILTER_ALL_FIVE is True, defaults to 1 (i.e., >0). If set (e.g., 2), uses that threshold.
# ----------------------------------------------------

import re
from pathlib import Path
import pandas as pd

FAMILY_PATTERN = re.compile(r"^(Galpha|Gbeta|Ggamma|RBP1|TBP).*_count$")
REQUIRED_EXACT = ["Galpha_count", "Gbeta_count", "Ggamma_count", "RBP1_count", "TBP_count"]

def find_taxid_column(df: pd.DataFrame) -> str:
    if "taxID" in df.columns:
        return "taxID"
    for c in df.columns:
        if c.lower() == "taxid":
            return c
    for c in df.columns:
        if c.lower() in {"tax_id","taxonomyid","ncbi_taxid","taxonid","taxon_id"}:
            return c
    raise KeyError("Could not find a taxID column (case-insensitive).")

def select_family_count_cols(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if FAMILY_PATTERN.match(c)]

def to_int_series(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce").fillna(0).astype(int)

def merge_and_optionally_filter(
    base_path: str,
    second_path: str,
    out_path: str,
    filter_all_five: bool = False,
    min_count: int | None = None
):
    base_p = Path(base_path).expanduser()
    second_p = Path(second_path).expanduser()
    out_p = Path(out_path).expanduser()

    if not base_p.exists():
        raise FileNotFoundError(f"Base TSV not found: {base_p}")
    if not second_p.exists():
        raise FileNotFoundError(f"Second TSV not found: {second_p}")

    # Load data
    base_df = pd.read_csv(base_p, sep="\t", dtype=str, low_memory=False)
    second_df = pd.read_csv(second_p, sep="\t", dtype=str, low_memory=False)

    base_tax = find_taxid_column(base_df)
    second_tax = find_taxid_column(second_df)

    # Pick requested *_count columns from the second file
    count_cols = select_family_count_cols(second_df)
    if not count_cols:
        raise ValueError("No columns matched ^(Galpha|Gbeta|Ggamma|RBP1|TBP).*_count$ in the second TSV.")

    # Ensure numeric for those counts, aggregate duplicates by taxID
    for c in count_cols:
        second_df[c] = to_int_series(second_df[c])

    reduced = second_df[[second_tax] + count_cols].copy()
    agg = reduced.groupby(second_tax, as_index=False)[count_cols].sum()

    # Align taxID column name to base
    if second_tax != base_tax:
        agg = agg.rename(columns={second_tax: base_tax})

    # Merge (left join keeps ALL original base columns)
    merged = base_df.merge(agg, on=base_tax, how="left")

    if not filter_all_five:
        # Write the merged file directly to the single OUT_PATH
        merged.to_csv(out_p, sep="\t", index=False)
        subset = None
    else:
        # Determine threshold: default to 1 if not provided
        threshold = 1 if min_count is None else int(min_count)

        # Ensure the five exact columns exist for filtering
        for col in REQUIRED_EXACT:
            if col not in merged.columns:
                merged[col] = 0
            merged[col] = to_int_series(merged[col])

        subset = merged[
            (merged["Galpha_count"] >= threshold) &
            (merged["Gbeta_count"]  >= threshold) &
            (merged["Ggamma_count"] >= threshold) &
            (merged["RBP1_count"]   >= threshold) &
            (merged["TBP_count"]    >= threshold)
        ].copy()

        # Write the filtered subset to the same OUT_PATH
        subset.to_csv(out_p, sep="\t", index=False)

    report = {
        "base_file": str(base_p),
        "second_file": str(second_p),
        "out_path": str(out_p),
        "taxid_base_col": base_tax,
        "taxid_second_col": second_tax,
        "selected_columns": count_cols,
        "filter_all_five": filter_all_five,
        "min_count_used": (1 if (filter_all_five and min_count is None) else min_count),
        "rows_in_base": len(base_df),
        "rows_in_merged": len(merged),
        "rows_in_output": (len(subset) if subset is not None else len(merged))
    }
    return merged, (subset if filter_all_five else None), report

if __name__ == "__main__":
    merged_df, subset_df, info = merge_and_optionally_filter(
        BASE_PATH,
        SECOND_PATH,
        OUT_PATH,
        filter_all_five=FILTER_ALL_FIVE,
        min_count=MIN_COUNT,  # optional; None -> defaults to 1 when filtering
    )
    print(info)
