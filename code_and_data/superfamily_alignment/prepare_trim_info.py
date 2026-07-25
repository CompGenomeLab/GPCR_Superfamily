#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 7: prepare trim_info.tsv from logo-representative outputs.

Output format is compatible with msa_trim_ref.py:
  acc1    family    acc2
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, Optional

import pandas as pd


def extract_accession_from_header(header: str) -> str:
    """Extract accession token from a FASTA-like header."""
    h = str(header).strip()
    if "|" in h:
        parts = h.split("|")
        if len(parts) > 1 and parts[1].strip():
            return parts[1].strip()
    if "_" in h:
        parts = h.split("_")
        if len(parts) > 1 and parts[1].strip():
            return parts[1].strip()
    return h


def load_logo_representative_table(path: Path) -> pd.DataFrame:
    """
    Load table from CSV/TSV using pandas auto-detection.
    Must contain a family column and either accession or sequence/header column.
    """
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    df = pd.read_csv(path, sep=None, engine="python", dtype=str, keep_default_na=False)
    if df.empty:
        raise ValueError(f"Input table is empty: {path}")
    return df


def canonicalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = [str(c).strip().lower() for c in out.columns]
    return out


def build_trim_info(
    reps_df: pd.DataFrame,
    family_col: str = "family",
    accession_col: Optional[str] = "accession",
    sequence_col: Optional[str] = "sequence",
    score_col: Optional[str] = None,
    acc2_by_family: Optional[Dict[str, str]] = None,
) -> pd.DataFrame:
    """
    Build trim-info table with columns: acc1, family, acc2.
    If multiple rows per family exist, choose one row per family:
      - by highest score_col when provided and present
      - otherwise by alphabetical accession
    """
    df = canonicalize_columns(reps_df)
    fam_col = family_col.lower()
    acc_col = accession_col.lower() if accession_col else None
    seq_col = sequence_col.lower() if sequence_col else None
    scr_col = score_col.lower() if score_col else None

    if fam_col not in df.columns:
        raise ValueError(f"Missing family column '{family_col}'. Columns: {list(df.columns)}")

    def _col_as_series(frame: pd.DataFrame, col_name: str) -> pd.Series:
        """Return a single Series even when duplicate column names exist."""
        col = frame[col_name]
        if isinstance(col, pd.DataFrame):
            return col.iloc[:, 0]
        return col

    if acc_col and acc_col in df.columns:
        df["__acc1"] = _col_as_series(df, acc_col).astype(str).str.strip()
    elif seq_col and seq_col in df.columns:
        df["__acc1"] = _col_as_series(df, seq_col).astype(str).map(extract_accession_from_header).str.strip()
    else:
        raise ValueError(
            f"Need either accession column '{accession_col}' or sequence/header column '{sequence_col}'. "
            f"Columns: {list(df.columns)}"
        )

    df["__family"] = _col_as_series(df, fam_col).astype(str).str.strip()
    df = df[(df["__family"] != "") & (df["__acc1"] != "")]
    if df.empty:
        raise ValueError("No valid family/accession rows after filtering.")

    if scr_col and scr_col in df.columns:
        df["__score"] = pd.to_numeric(df[scr_col], errors="coerce")
        df = df.sort_values(["__family", "__score", "__acc1"], ascending=[True, False, True])
    else:
        df = df.sort_values(["__family", "__acc1"], ascending=[True, True])

    one_per_family = df.drop_duplicates(subset=["__family"], keep="first").copy()
    out = one_per_family.loc[:, ["__acc1", "__family"]].copy()
    out.columns = ["acc1", "family"]

    out["acc2"] = ""
    if acc2_by_family:
        lower_map = {str(k).lower(): str(v) for k, v in acc2_by_family.items()}
        out["acc2"] = out["family"].map(lambda fam: acc2_by_family.get(fam, lower_map.get(str(fam).lower(), "")))

    return out


def write_trim_info(df: pd.DataFrame, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)


if __name__ == "__main__":
    # ================= USER VARIABLES ================= #
    # Input table exported by Step 5 representative selection
    STAGE_DIR = Path(__file__).resolve().parent
    LOGO_REPS_TABLE = Path(os.environ.get(
        "GPCR_LOGO_REPS_TABLE",
        STAGE_DIR / "data" / "logo_outputs" / "logo_representatives_May22_custom.tsv",
    ))

    # Output naming (match your existing style like *_May22.tsv)
    TRIM_INFO_OUT = Path(os.environ.get("GPCR_TRIM_INFO_OUT", STAGE_DIR / "trim_info.tsv"))

    # Column mapping for your Step 6 output table
    FAMILY_COL = "family"
    ACCESSION_COL = "accession"   # set to None to derive from sequence/header column
    SEQUENCE_COL = "sequence"     # used when ACCESSION_COL is None
    SCORE_COL = None              # e.g., "consensus_identity" if multiple rows per family

    # Default acc2 references: auto-filled only when known from prior trim-info files
    ACC2_BY_FAMILY = {
        "GP143": "P51810",
        "Olfactory": "Q9H255",
        "classA": "P25021",
        "classB1": "P43220",
        "classB2": "Q9HAR2",
        "classC": "P41180",
        "classF": "O75084",
        "classT": "P59534",
        "GP157": "Q5UAW9",
    }
    # ================================================== #

    reps = load_logo_representative_table(LOGO_REPS_TABLE)
    trim_df = build_trim_info(
        reps_df=reps,
        family_col=FAMILY_COL,
        accession_col=ACCESSION_COL,
        sequence_col=SEQUENCE_COL,
        score_col=SCORE_COL,
        acc2_by_family=ACC2_BY_FAMILY,
    )
    write_trim_info(trim_df, TRIM_INFO_OUT)

    print(f"Wrote trim info: {TRIM_INFO_OUT}")
    print(f"Families in trim info: {trim_df['family'].nunique()}")
