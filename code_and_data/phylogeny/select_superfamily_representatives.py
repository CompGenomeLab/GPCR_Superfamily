# -*- coding: utf-8 -*-
from __future__ import annotations
import os
import json
import glob
import time
from typing import Dict, List, Optional, Iterable, Set, Tuple
import pandas as pd
import numpy as np
import logging
from contextlib import contextmanager
from collections import OrderedDict


# Final representatives used in the study. These were selected manually after
# comparing automated metrics and within-/between-family alignment quality.
MANUAL_REPRESENTATIVES = OrderedDict([
    ("GP143", "A0A835NM65"),
    ("GPR1", "A0A3M2SGC9"),
    ("Mth", "A0A8J4PRC8"),
    ("Nematode", "A0A6V7VX31"),
    ("Olfactory", "A0A8C0QWZ4"),
    ("STE3", "A0A395RUA3"),
    ("Vomeronasal1", "A0A2I4B258"),
    ("Vomeronasal2", "A0A7N6AN26"),
    ("cAMP", "A0A7J6UXT4"),
    ("classA", "A0A556TZG5"),
    ("classB1", "A0A6J2X0B2"),
    ("classB2", "A0A0M9A0X0"),
    ("classC", "A0A8C7JL03"),
    ("classF", "A0AA97IXW4"),
    ("classT", "Q646E2"),
    ("GP157", "A0A6P7ZWE1"),
    ("FSLB", "A0A9P6IC23"),
])

# =========================
# Logging / Verbosity knobs
# =========================
LOG_LEVEL = os.environ.get("GPCR_LOG_LEVEL", "INFO").upper()  # DEBUG | INFO | WARNING
TRACE_FAMILIES: Set[str] = {"classA", "cAMP"}  # Families included in detailed diagnostic output.
TRACE_SAMPLES = 10

logging.basicConfig(
    level=getattr(logging, LOG_LEVEL, logging.INFO),
    format="%(asctime)s | %(levelname)-7s | %(message)s",
)

@contextmanager
def stage(name: str):
    t0 = time.time()
    logging.info(f"â–¶ {name} ...")
    try:
        yield
    finally:
        dt = time.time() - t0
        logging.info(f"âœ“ {name} done in {dt:.2f}s")

def _head(df: pd.DataFrame, n=5) -> pd.DataFrame:
    try:
        return df.head(n)
    except Exception:
        return df

def _require_file(path: str, kind: str) -> None:
    if not isinstance(path, str) or not path:
        raise ValueError(f"{kind}: expected non-empty string path.")
    if not os.path.exists(path):
        raise FileNotFoundError(f"{kind}: file not found -> {path}")

def _normalize_str(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip()

def infer_gpcrfamily_from_filename(path: str) -> str:
    base = os.path.basename(path)
    return base.split("_", 1)[0] if "_" in base else os.path.splitext(base)[0]

# ---------- 0) Canonicalization ----------
UNK = "__unk__"

def _canonicalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Lowercase all column names and coalesce duplicates that collide after lowercasing.
    Coalescing rule: leftâ†’right, take the first non-null / non-empty value across duplicates.
    """
    if df is None:
        return df
    g = df.copy()
    g.columns = [str(c).strip().lower() for c in g.columns]
    # handle collisions
    lower = list(g.columns)
    unique_order = list(OrderedDict.fromkeys(lower).keys())
    newcols = OrderedDict()
    for name in unique_order:
        idxs = [i for i, c in enumerate(lower) if c == name]
        block = g.iloc[:, idxs]
        if block.shape[1] == 1:
            col = block.iloc[:, 0]
            newcols[name] = col
        else:
            tmp = block.copy()
            for c in tmp.columns:
                tmp[c] = tmp[c].replace("", np.nan)
            combined = tmp.bfill(axis=1).iloc[:, 0]
            newcols[name] = combined
    out = pd.DataFrame(newcols, index=g.index)
    return out

def _nunique_nonunk(series: pd.Series) -> int:
    return (
        series.fillna(UNK)
        .astype(str).str.strip()
        .replace({"": UNK})
        .pipe(lambda s: s[s != UNK])
        .nunique()
    )

def _prep_tax_ranks(df: pd.DataFrame, ranks_priority: Iterable[str]) -> pd.DataFrame:
    """
    Lowercase + dedupe columns, normalize unknowns on listed ranks (and 'kingdom' if present).
    """
    g = _canonicalize_columns(df)
    ranks = [r.lower() for r in ranks_priority]
    for r in ranks:
        if r in g.columns:
            g[r] = g[r].astype(str).str.strip()
            g[r] = g[r].replace({"": UNK, "na": UNK, "nan": UNK, "none": UNK, "None": UNK}).fillna(UNK)
    if "kingdom" in g.columns:
        g["kingdom"] = g["kingdom"].astype(str).str.strip().replace({"": UNK, "na": UNK, "nan": UNK}).fillna(UNK)
    return g

# ---------- 1) Load mappings ----------
def load_accession_taxid_tsv(tsv_path: str) -> pd.DataFrame:
    with stage("Load accessionâ†’taxID TSV"):
        _require_file(tsv_path, "accession_tsv")
        df = pd.read_csv(tsv_path, sep="\t", dtype=str, keep_default_na=False)
        df = _canonicalize_columns(df)
        # Accept 'taxid' (lowercased from TaxID/taxID), 'family', 'accession'
        required = ["family", "accession", "taxid"]
        missing = [c for c in required if c not in df.columns]
        if missing:
            raise ValueError(f"accession_tsv: missing required columns {missing}; got {list(df.columns)}")

        df["family"]    = _normalize_str(df["family"])
        df["accession"] = _normalize_str(df["accession"])
        df["taxid"]     = _normalize_str(df["taxid"])
        before = len(df)
        df = df[(df["accession"] != "") & (df["taxid"] != "")]
        logging.info(f"accession_tsv rows: {before} â†’ {len(df)} after empty accession/taxid drop")
        fam_counts = df["family"].value_counts().to_dict()
        logging.info(f"Families in TSV: {len(fam_counts)} (top 10 counts shown)")
        for i, (fam, cnt) in enumerate(sorted(fam_counts.items(), key=lambda x: (-x[1], x[0]))[:10]):
            logging.info(f"  - {fam}: {cnt}")
        return df[["family", "accession", "taxid"]]

def load_taxonomy_json(json_path: str) -> pd.DataFrame:
    with stage("Load taxonomy JSON"):
        _require_file(json_path, "taxonomy_json")
        with open(json_path, "r", encoding="utf-8") as fh:
            data = json.load(fh)
        if not isinstance(data, list) or (len(data) and not isinstance(data[0], dict)):
            raise ValueError("taxonomy_json: must be a JSON list of objects.")
        tax = pd.DataFrame(data)
        tax = _canonicalize_columns(tax)
        if tax.empty:
            raise ValueError("taxonomy_json: contains zero rows.")
        if "taxid" not in tax.columns:
            raise ValueError("taxonomy_json: missing 'taxid' column (case-insensitive).")
        tax["taxid"] = _normalize_str(tax["taxid"])
        for c in tax.columns:
            if c != "taxid" and pd.api.types.is_string_dtype(tax[c]):
                tax[c] = tax[c].astype(str).str.strip()
        if "family" in tax.columns:
            tax = tax.rename(columns={"family": "family_tax"})
        logging.info(f"taxonomy rows: {len(tax)}; columns: {list(tax.columns)}")
        return tax

# ---------- 2) Parse cluster files ----------
def parse_cluster_file_multiline(path: str, gpcrfamily: Optional[str] = None) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, comment="#", dtype=str, engine="python")
    if df.shape[1] < 1:
        return pd.DataFrame(columns=["gpcrfamily","accession","cluster_size"])
    reps = _normalize_str(df.iloc[:, 0]).replace("", pd.NA).dropna()
    if reps.empty:
        return pd.DataFrame(columns=["gpcrfamily","accession","cluster_size"])
    out = reps.value_counts().rename_axis("accession").reset_index(name="cluster_size")
    out["cluster_size"] = out["cluster_size"].astype(int)
    fam = gpcrfamily or infer_gpcrfamily_from_filename(path)
    out["gpcrfamily"] = str(fam)
    out = _canonicalize_columns(out)
    return out[["gpcrfamily","accession","cluster_size"]]

def parse_all_clusters(cluster_glob: str) -> pd.DataFrame:
    with stage(f"Parse cluster files from glob: {cluster_glob}"):
        if not isinstance(cluster_glob, str) or not cluster_glob:
            raise ValueError("cluster_glob: expected non-empty glob string.")
        files = glob.glob(cluster_glob)
        if not files:
            raise FileNotFoundError(f"cluster_glob matched 0 files: {cluster_glob}")
        logging.info(f"Matched {len(files)} cluster files")
        rows = [parse_cluster_file_multiline(fp) for fp in files]
        out = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["gpcrfamily","accession","cluster_size"])
        logging.info(f"Parsed reps: {len(out)} across {out['gpcrfamily'].nunique()} families")
        fam_counts = out["gpcrfamily"].value_counts().to_dict()
        logging.info("Families in clusters (top 10):")
        for i, (fam, cnt) in enumerate(sorted(fam_counts.items(), key=lambda x: (-x[1], x[0]))[:10]):
            logging.info(f"  - {fam}: {cnt}")
        logging.debug("Cluster sample:\n%s", _head(out, 5).to_string(index=False))
        return out

# ---------- 3) Hierarchical, recursive diversity selection ----------
def _first_diverse_rank(df: pd.DataFrame, ranks: List[str]) -> Optional[str]:
    for r in ranks:
        if r in df.columns and _nunique_nonunk(df[r]) > 1:
            return r
    return None

def _sort_global(df: pd.DataFrame) -> pd.DataFrame:
    return df.sort_values(["cluster_size", "accession"], ascending=[False, True]).reset_index(drop=True)

def _tuple_prefix(row: pd.Series, ranks: List[str], upto: int) -> Tuple[str, ...]:
    vals = []
    for i, r in enumerate(ranks[:upto]):
        vals.append(str(row.get(r, UNK)))
    return tuple(vals)

def hierarchical_diverse_topk(
    family_df: pd.DataFrame,
    k: int,
    ranks_priority: Iterable[str],
    min_cluster_size: int = 0,
) -> pd.DataFrame:
    """
    Select up to k rows maximizing diversity along ranks_priority:
      - Choose PRIMARY = first rank in ranks_priority with >1 non-UNK values.
      - Round-robin across PRIMARY bins (least-filled first).
      - Within a PRIMARY bin, greedily open new deeper prefixes (primary+phylum, then +class, ...).
      - Falls back naturally if a rank has only one non-UNK value or bins are exhausted.
    """
    fam_name = str(family_df["gpcrfamily"].iloc[0]) if not family_df.empty and "gpcrfamily" in family_df.columns else "?"
    logging.debug(f"[{fam_name}] selection start: k={k}, min_cluster_size={min_cluster_size}, ranks={list(ranks_priority)}")
    if k <= 0 or family_df.empty:
        logging.warning(f"[{fam_name}] k<=0 or empty family_df; returning empty.")
        return family_df.head(0)

    ranks = [r.lower() for r in ranks_priority]
    df = _prep_tax_ranks(family_df, ranks)

    need = {"accession", "cluster_size"}
    if not need.issubset(df.columns):
        raise ValueError(f"[{fam_name}] need columns {sorted(list(need - set(df.columns)))}.")

    df["accession"] = df["accession"].astype(str).str.strip()
    df["cluster_size"] = pd.to_numeric(df["cluster_size"], errors="coerce").fillna(0).astype(int)

    before = len(df)
    df = df[df["cluster_size"] >= min_cluster_size].copy()
    logging.debug(f"[{fam_name}] filter min_cluster_size: {before} â†’ {len(df)}")
    if df.empty:
        return df.head(0)

    primary = _first_diverse_rank(df, ranks)
    if primary is None:
        logging.info(f"[{fam_name}] no rank with diversity; falling back to global top-k.")
        return _sort_global(df).head(k).reset_index(drop=True)

    ranks_ordered = [primary] + [r for r in ranks if r != primary]
    df = _sort_global(df)

    picked_ids: Set[str] = set()
    picks: List[pd.Series] = []

    primary_bins = sorted(set(df[primary].fillna(UNK).astype(str).tolist()))
    counts_by_primary: Dict[str, int] = {p: 0 for p in primary_bins}

    max_depth = len(ranks_ordered)
    used_prefixes_by_primary: Dict[str, List[Set[Tuple[str, ...]]]] = {
        pv: [set() for _ in range(max_depth + 1)] for pv in primary_bins
    }

    def remaining_in_primary(pv: str) -> pd.DataFrame:
        rem = df[~df["accession"].isin(picked_ids)]
        return rem[rem[primary] == pv]

    def choose_candidate(sub: pd.DataFrame, pv: str) -> Optional[pd.Series]:
        if sub.empty:
            return None
        # Try to open a new prefix at increasing depth (primary+next ... full chain)
        for depth in range(2, max_depth + 1):
            seen = used_prefixes_by_primary[pv][depth]
            for _, row in sub.iterrows():
                prefix = _tuple_prefix(row, ranks_ordered, depth)
                if prefix not in seen:
                    return row
        # fallback: best remaining in this primary
        return sub.iloc[0]

    made_progress = True
    while len(picks) < k and made_progress:
        made_progress = False
        order = sorted(primary_bins, key=lambda pv: (counts_by_primary.get(pv, 0), pv))
        for pv in order:
            if len(picks) >= k:
                break
            sub = remaining_in_primary(pv)
            if sub.empty:
                continue
            row = choose_candidate(sub, pv)
            if row is None:
                continue
            acc = str(row["accession"])
            if acc in picked_ids:
                continue
            picked_ids.add(acc)
            picks.append(row)
            counts_by_primary[pv] = counts_by_primary.get(pv, 0) + 1
            # mark all prefixes used for this pick
            for depth in range(1, max_depth + 1):
                prefix = _tuple_prefix(row, ranks_ordered, depth)
                used_prefixes_by_primary[pv][depth].add(prefix)
            made_progress = True

    if len(picks) < k:
        rem = df[~df["accession"].isin(picked_ids)]
        if not rem.empty:
            extra = rem.head(k - len(picks))
            for _, row in extra.iterrows():
                picks.append(row)

    if not picks:
        return df.head(0)

    out = pd.DataFrame(picks).drop_duplicates(subset=["accession"], keep="first").head(k).reset_index(drop=True)
    logging.info(f"[{fam_name}] selected {len(out)}/{k} (balanced by {primary}; recursive diversification)")
    return out

# ---------- 4) Utilities: reporting & diagnostics ----------
def _compare_targets_vs_clusters(rep_sizes: pd.DataFrame, targets_per_family: Dict[str, int]) -> None:
    cluster_fams = set(rep_sizes["gpcrfamily"].unique())
    target_fams = set([f.lower() for f in targets_per_family.keys()])
    missing_in_clusters = [f for f in target_fams if f not in cluster_fams]
    if missing_in_clusters:
        logging.warning("TARGETS contains families not present in cluster files: %s", missing_in_clusters)

    lower_clusters = {f.lower(): f for f in cluster_fams}
    suspicious_case = []
    for t in target_fams:
        if t not in cluster_fams and t.lower() in lower_clusters:
            suspicious_case.append((t, lower_clusters[t.lower()]))
    if suspicious_case:
        logging.warning("Case/typo mismatch candidates (TARGETS -> cluster family): %s", suspicious_case)

def _diagnose_family_pipeline(
    fam: str,
    rep_sizes: pd.DataFrame,
    acc_tax_df: pd.DataFrame,
    tax_df: pd.DataFrame,
    min_cluster_size: int
) -> None:
    logging.info(f"--- Diagnostics for family '{fam}' ---")
    fam_rep = rep_sizes[rep_sizes["gpcrfamily"] == fam].copy()
    logging.info(f"[{fam}] from clusters: {len(fam_rep)} rows")
    if fam_rep.empty:
        logging.info(f"[{fam}] No representatives parsed from cluster files. Check filename prefixes/glob.")
        return

    acc_tax = _canonicalize_columns(acc_tax_df[["family", "accession", "taxid"]].drop_duplicates().copy())
    acc_tax["family"]    = _normalize_str(acc_tax["family"])
    acc_tax["accession"] = _normalize_str(acc_tax["accession"])
    acc_tax["taxid"]     = _normalize_str(acc_tax["taxid"])

    fam_rep["gpcrfamily"] = _normalize_str(fam_rep["gpcrfamily"])
    fam_rep["accession"]  = _normalize_str(fam_rep["accession"])

    joined = fam_rep.merge(
        acc_tax,
        left_on=["accession", "gpcrfamily"],
        right_on=["accession", "family"],
        how="inner"
    )
    logging.info(f"[{fam}] after accession TSV (family-matched) join: {len(joined)} rows (from {len(fam_rep)})")
    if len(joined) == 0:
        in_tsv = acc_tax_df[acc_tax_df["accession"].isin(fam_rep["accession"])].copy()
        logging.info(f"[{fam}] accessions present in TSV (any family): {len(in_tsv)}")
        missing_acc = fam_rep[~fam_rep["accession"].isin(acc_tax_df["accession"])]
        if not missing_acc.empty:
            logging.info(f"[{fam}] sample accessions missing in TSV (up to {TRACE_SAMPLES}):\n%s", _head(missing_acc, TRACE_SAMPLES).to_string(index=False))
        return

    tax = _canonicalize_columns(tax_df.copy())
    tax["taxid"] = _normalize_str(tax["taxid"])

    joined_tax = joined.merge(tax, on="taxid", how="inner", suffixes=("", "_tax"))
    logging.info(f"[{fam}] after taxonomy join: {len(joined_tax)} rows (from {len(joined)})")
    if len(joined_tax) == 0:
        missing_tax = joined[~joined["taxid"].isin(tax["taxid"])]
        logging.info(f"[{fam}] sample taxids missing from taxonomy JSON (up to {TRACE_SAMPLES}):\n%s",
                     _head(missing_tax[['accession','taxid']], TRACE_SAMPLES).to_string(index=False))
        return

    before = len(joined_tax)
    after_df = joined_tax[joined_tax["cluster_size"] >= min_cluster_size]
    logging.info(f"[{fam}] after min_cluster_size={min_cluster_size}: {before} â†’ {len(after_df)}")
    logging.info(f"[{fam}] sample candidates (up to {TRACE_SAMPLES}):\n%s",
                 _head(after_df[['accession','cluster_size','taxid']], TRACE_SAMPLES).to_string(index=False))

def _keep_case_insensitive(df: pd.DataFrame, names: Iterable[str]) -> pd.DataFrame:
    present = {c.lower(): c for c in df.columns}
    keep = [present[n.lower()] for n in names if n.lower() in present]
    return df.loc[:, keep] if keep else df.iloc[:, 0:0]

# ---------- 5) Kingdom summaries ----------
def print_kingdom_summary(df: pd.DataFrame, label: str = "pre-selection") -> None:
    g = _canonicalize_columns(df)
    if "gpcrfamily" not in g.columns:
        logging.warning(f"[{label}] Need 'gpcrfamily' column to summarize kingdoms.")
        return
    if "kingdom" not in g.columns:
        logging.warning(f"[{label}] No 'kingdom' column found.")
        return

    tmp = g.copy()
    tmp["gpcrfamily"] = tmp["gpcrfamily"].astype(str).str.strip()
    tmp["kingdom"] = tmp["kingdom"].astype(str).str.strip().replace({"": UNK, "na": UNK, "nan": UNK}).fillna(UNK)

    per_fam_counts = tmp.groupby(["gpcrfamily", "kingdom"]).size().rename("count").reset_index()

    print(f"=== {label}: number of kingdoms represented per family ===")
    for fam, fam_grp in per_fam_counts.groupby("gpcrfamily", sort=True):
        sub = fam_grp.sort_values("kingdom")
        distinct_n = sub[sub["kingdom"] != UNK]["kingdom"].nunique()
        parts = [f"{k}:{int(c)}" for k, c in zip(sub["kingdom"], sub["count"])]
        print(f"{fam}\t{distinct_n}\t[{', '.join(parts)}]")

def print_selected_kingdom_counts(df: pd.DataFrame) -> None:
    print_kingdom_summary(df, label="post-selection")

# ---------- 6) Orchestration ----------
def sample_representatives_for_all_families(
    cluster_glob: str,
    accession_tsv: str,
    taxonomy_json: str,
    targets_per_family: Dict[str, int],
    ranks_priority: Iterable[str],
    min_cluster_size: int = 0,
    accession_filter_file: Optional[str] = None,
) -> pd.DataFrame:

    BASE_COLS = ["gpcrfamily", "accession", "cluster_size", "taxid"]

    # 1) Load mapping and taxonomy
    acc_tax_df = load_accession_taxid_tsv(accession_tsv)  # family, accession, taxid (lower)
    tax_df = load_taxonomy_json(taxonomy_json)            # taxid, taxonomy columns (lower)

    # 2) Optional whitelist
    whitelist: Optional[set] = None
    if accession_filter_file:
        with stage("Load accession whitelist"):
            _require_file(accession_filter_file, "accession_filter_file")
            with open(accession_filter_file, "r", encoding="utf-8") as fh:
                whitelist = {line.strip() for line in fh if line.strip()}
            logging.info(f"Whitelist accessions: {len(whitelist)}")

    # 3) Parse cluster files
    rep_sizes = parse_all_clusters(cluster_glob)
    if whitelist is not None:
        before = len(rep_sizes)
        rep_sizes = rep_sizes[rep_sizes["accession"].isin(whitelist)].copy()
        logging.info(f"Apply whitelist: {before} â†’ {len(rep_sizes)}")
        if rep_sizes.empty:
            logging.warning("After whitelist, no representatives remain.")
            desired = BASE_COLS + [r.lower() for r in ranks_priority]
            return pd.DataFrame(columns=desired)

    _compare_targets_vs_clusters(rep_sizes, targets_per_family)

    # 4) STRICT mapping with FAMILY MATCH (lowercase everything)
    acc_tax = acc_tax_df[["family", "accession", "taxid"]].drop_duplicates().copy()
    acc_tax = _canonicalize_columns(acc_tax)
    acc_tax["family"]    = _normalize_str(acc_tax["family"])
    acc_tax["accession"] = _normalize_str(acc_tax["accession"])
    acc_tax["taxid"]     = _normalize_str(acc_tax["taxid"])

    rep_sizes = rep_sizes.copy()
    rep_sizes["gpcrfamily"] = _normalize_str(rep_sizes["gpcrfamily"])
    rep_sizes["accession"]  = _normalize_str(rep_sizes["accession"])

    with stage("Join clusters â†” TSV (strict family match)"):
        reps_with_taxid = rep_sizes.merge(
            acc_tax,
            left_on=["accession", "gpcrfamily"],
            right_on=["accession", "family"],
            how="inner"
        ).drop(columns=["family"])
        logging.info(f"Strict family match join: {len(rep_sizes)} â†’ {len(reps_with_taxid)} rows")

    if reps_with_taxid.empty:
        logging.warning("No rows after strict family match. Run diagnostics for suspect families:")
        for fam in sorted(set([f.lower() for f in targets_per_family.keys()]) | set(rep_sizes["gpcrfamily"].unique())):
            if fam in TRACE_FAMILIES:
                _diagnose_family_pipeline(fam, rep_sizes, acc_tax_df, tax_df, min_cluster_size)
        desired = BASE_COLS + [r.lower() for r in ranks_priority]
        return pd.DataFrame(columns=desired)

    # 5) Taxonomy join
    with stage("Join with taxonomy JSON"):
        tax = _canonicalize_columns(tax_df.copy())
        tax["taxid"] = _normalize_str(tax["taxid"])
        reps_full = reps_with_taxid.merge(tax, on="taxid", how="inner", suffixes=("", "_tax"))
        logging.info(f"Taxonomy join: {len(reps_with_taxid)} â†’ {len(reps_full)} rows")

    if reps_full.empty:
        logging.warning("No rows after taxonomy join. Running diagnostics for traced families.")
        for fam in TRACE_FAMILIES:
            _diagnose_family_pipeline(fam, rep_sizes, acc_tax_df, tax_df, min_cluster_size)
        desired = BASE_COLS + [r.lower() for r in ranks_priority]
        return pd.DataFrame(columns=desired)

    # 5.5) Pre-selection kingdom possibilities
    with stage("Pre-selection kingdom possibilities"):
        print_kingdom_summary(reps_full, label="pre-selection")

    # 6) Per-family selection with recursive diversification
    selections: List[pd.DataFrame] = []
    with stage("Hierarchical selection per family"):
        fam_list = sorted(reps_full["gpcrfamily"].unique())
        logging.info(f"Families available post-joins: {fam_list}")
        for fam, fam_df in reps_full.groupby("gpcrfamily", dropna=False):
            k = int(targets_per_family.get(fam, 0) or targets_per_family.get(str(fam).lower(), 0))
            if k <= 0:
                logging.info(f"[{fam}] TARGETS requests k={k}; skipping.")
                continue

            if fam in TRACE_FAMILIES:
                _diagnose_family_pipeline(fam, rep_sizes, acc_tax_df, tax_df, min_cluster_size)

            sel = hierarchical_diverse_topk(
                family_df=fam_df.assign(gpcrfamily=fam),
                k=k,
                ranks_priority=[r.lower() for r in ranks_priority],
                min_cluster_size=min_cluster_size,
            )
            if sel.empty:
                logging.warning(f"[{fam}] selection returned 0 (k={k}). Consider checking ranks/filters.")
                continue

            desired_cols = BASE_COLS + [r.lower() for r in ranks_priority]
            if "kingdom" not in desired_cols:
                desired_cols.append("kingdom")
            sel = _keep_case_insensitive(sel, desired_cols)

            selections.append(sel.assign(gpcrfamily=fam) if "gpcrfamily" not in sel.columns else sel)

    if not selections:
        logging.warning("No selections produced for any family.")
        desired = BASE_COLS + [r.lower() for r in ranks_priority]
        return pd.DataFrame(columns=desired)

    result = pd.concat(selections, ignore_index=True)

    # final trim
    desired = BASE_COLS + [r.lower() for r in ranks_priority]
    if "kingdom" not in desired:
        desired.append("kingdom")
    result = _keep_case_insensitive(result, desired)

    logging.info("Total selected: %d unique accessions", result["accession"].nunique() if "accession" in result.columns else 0)
    return result

def write_accession_list(df: pd.DataFrame, out_path: str) -> None:
    if df.empty or "accession" not in df.columns:
        raise ValueError("write_accession_list: DataFrame empty or missing 'accession'.")
    df["accession"].dropna().astype(str).str.strip().replace("", pd.NA).dropna().drop_duplicates() \
        .to_csv(out_path, index=False, header=False)


def write_manual_representatives(out_path: str) -> None:
    pd.DataFrame(
        MANUAL_REPRESENTATIVES.items(), columns=["family", "accession"]
    ).to_csv(out_path, sep="\t", index=False)

#%% ---------- main ----------
if __name__ == "__main__":
    # EDIT THESE:
    PROJECT_DIR = os.environ.get("GPCR_ANALYSIS_DIR", os.environ.get("GPCR_ZENODO_DIR", "/path/to/GPCR_Superfamily"))
    CLUSTER_GLOB = os.path.join(
        os.environ.get("GPCR_CLUSTER_DIR", os.path.join(PROJECT_DIR, "families", "clusters")),
        "*_genes_filtered_db_clust0.7_clust0_cov1.tsv",
    )
    ACCESSION_TSV = os.environ.get(
        "GPCR_ACCESSION_TSV",
        os.path.join(PROJECT_DIR, "families", "tables", "acc_taxid_map_BestEval_May22_e0.1.tsv"),
    )
    TAXONOMY_JSON  = os.environ.get(
        "GPCR_TAXONOMY_JSON",
        os.path.join(PROJECT_DIR, "families", "reference_data", "taxonomy_eukaryotes_filtered.json"),
    )

    # Lowercase keys to match canonical family names.
    TARGETS = {
        "camp": 100, "classa": 100, "classb1": 20, "classb2": 20, "classf": 20,
        "classt": 20, "gp143": 25, "mth": 20, "nematode": 20, "olfactory": 20,
        "vomeronasal1": 20, "gpr1": 25, "fslb": 20, "gp157": 25,"vomeronasal2":20,
        "classc":50,"ste3":25
    }
    


    RANKS_PRIORITY = [
        "kingdom","phylum","subphylum","superclass","class","superorder","order",
        "suborder","infraorder","parvorder","superfamily","family_tax","subfamily","genus","species"
    ]

    MIN_CLUSTER_SIZE = 2
    ACCESSION_FILTER_FILE = None

    os.makedirs(PROJECT_DIR, exist_ok=True)
    manual_out = os.path.join(PROJECT_DIR, "manual_superfamily_representatives.tsv")
    write_manual_representatives(manual_out)
    logging.info("[OK] wrote manually selected study representatives -> %s", manual_out)

    try:
        chosen = sample_representatives_for_all_families(
            cluster_glob=CLUSTER_GLOB,
            accession_tsv=ACCESSION_TSV,
            taxonomy_json=TAXONOMY_JSON,
            targets_per_family=TARGETS,
            ranks_priority=RANKS_PRIORITY,
            min_cluster_size=MIN_CLUSTER_SIZE,
            accession_filter_file=ACCESSION_FILTER_FILE,
        )

        out_ids = os.path.join(PROJECT_DIR, "candidate_superfamily_representatives_May22.txt")
        if not chosen.empty:
            # Post-selection kingdom counts with names
            print_selected_kingdom_counts(chosen)

            write_accession_list(chosen, out_ids)
            logging.info(f"[OK] wrote {chosen['accession'].nunique()} IDs -> {out_ids}")
            logging.debug("Chosen sample:\n%s", _head(chosen, 20).to_string(index=False))
        else:
            logging.warning("[WARN] No selections produced â€” see logs above for exact stage where rows disappeared.")

    except Exception as e:
        logging.exception("Fatal error: %s", e)
