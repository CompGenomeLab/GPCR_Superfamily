#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Prepare family accession text files after e-value filtering.

Supported input formats:
1) Mapping table (preferred; direct output from create_taxID_family_count_table.py):
   columns include: family, accession
2) Family-column table:
   one family per column, accession values under each column

Output:
  - <OUT_DIR>/<family>_accessions.txt
"""

from __future__ import annotations

import os
from pathlib import Path
import pandas as pd


def write_family_accession_files(input_tsv: Path, out_dir: Path) -> dict:
    if not input_tsv.exists():
        raise FileNotFoundError(f"Input TSV not found: {input_tsv}")

    df = pd.read_csv(input_tsv, sep="\t", dtype=str)
    if df.empty or len(df.columns) == 0:
        raise ValueError(f"Input TSV has no usable columns: {input_tsv}")

    out_dir.mkdir(parents=True, exist_ok=True)

    lower_cols = {str(c).strip().lower(): c for c in df.columns}
    written = {}

    # Preferred format: mapping table with family/accession columns.
    if "family" in lower_cols and "accession" in lower_cols:
        fam_col = lower_cols["family"]
        acc_col = lower_cols["accession"]

        tmp = df[[fam_col, acc_col]].copy()
        tmp[fam_col] = tmp[fam_col].fillna("").astype(str).map(str.strip)
        tmp[acc_col] = tmp[acc_col].fillna("").astype(str).map(str.strip)
        tmp = tmp[(tmp[fam_col] != "") & (tmp[acc_col] != "")]

        grouped = tmp.groupby(fam_col)[acc_col].apply(lambda s: sorted(set(s))).to_dict()
        items = grouped.items()
    else:
        # Fallback format: one family per column.
        items = []
        for family in df.columns:
            vals = (
                df[family]
                .fillna("")
                .astype(str)
                .map(str.strip)
            )
            items.append((family, sorted({v for v in vals if v})))

    for family, accessions in items:
        if not accessions:
            continue

        out_path = out_dir / f"{family}_accessions.txt"
        out_path.write_text("\n".join(accessions) + "\n", encoding="utf-8")
        written[family] = {"n_accessions": len(accessions), "path": str(out_path)}

    return written


if __name__ == "__main__":
    # ================= USER VARIABLES ================= #
    BASE_DIR = Path(os.environ.get("GPCR_ANALYSIS_DIR", os.environ.get("GPCR_ZENODO_DIR", "/path/to/GPCR_Superfamily")))

    # Exact broad-family mapping used to construct the archived family DBs.
    # The E<=0.1 mapping remains a downstream taxonomy/representative filter.
    FILTERED_FAMILY_ACCESSION_TSV = Path(os.environ.get(
        "GPCR_FAMILY_DB_ACCESSION_TSV",
        BASE_DIR / "families" / "tables" / "acc_taxid_map_prefilter_Oct24_e1000_all.tsv",
    ))

    # Output directory expected by familydb_cluster_align.sh
    OUT_DIR = Path(os.environ.get("GPCR_ACCESSION_TXT_DIR", BASE_DIR / "accessions"))
    # ================================================== #

    summary = write_family_accession_files(
        input_tsv=FILTERED_FAMILY_ACCESSION_TSV,
        out_dir=OUT_DIR,
    )

    print(f"Wrote family accession files for {len(summary)} families.")
    for fam, info in sorted(summary.items()):
        print(f"{fam}: {info['n_accessions']} -> {info['path']}")
