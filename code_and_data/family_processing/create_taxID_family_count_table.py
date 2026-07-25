#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build per-family taxonomic count tables LIMITED to accessions provided in a TSV,
with optional special-cases for one or more families that should use MMseqs m8 files
(3rd column contains a header with OX=... taxID) instead of HMMs.

Strictness: only accessions present in the TSV are ever mapped or counted.
"""

import os
import re
from pathlib import Path
import pandas as pd
from typing import Dict, Optional, Set, List, Tuple
from collections import defaultdict

# --------------------- Patterns --------------------- #
OX_PAT = re.compile(r"\bOX=(\d+)\b")

# --------------------- Helpers --------------------- #
def normalize_accession_token(tok: Optional[str]) -> Optional[str]:
    if tok is None:
        return None
    tok = str(tok).strip()
    if not tok:
        return None
    if "|" in tok:
        parts = tok.split("|")
        if len(parts) >= 2 and parts[1]:
            return parts[1].strip()
    return tok

def extract_taxid(text: str) -> Optional[str]:
    m = OX_PAT.search(text)
    return m.group(1) if m else None

# --------------------- IO: TSV --------------------- #
def read_family_accession_tsv(path: str) -> Dict[str, Set[str]]:
    df = pd.read_csv(path, sep="\t", dtype=str)
    fam2acc: Dict[str, Set[str]] = {}
    for fam in df.columns:
        col = df[fam].dropna().map(str).map(str.strip)
        col = col[col.astype(bool)]
        accs = {normalize_accession_token(x) for x in col}
        accs.discard(None)
        fam2acc[fam] = accs
    return fam2acc

# --------------------- IO: HMMER --------------------- #
def parse_all_hmm_for_acc_taxid(hmm_dir: str, restrict_to: Set[str]) -> Dict[str, Dict[str, Tuple[str, float]]]:
    """
    Parse all HMM result files in `hmm_dir`.

    Returns:
      { accession : { prefix : (taxid, best_evalue) } }

    - `prefix` = filename before first underscore (lowercased)
    - Accession extracted from "sp|ACC|..." / "tr|ACC|..." if present, else token
    - Keeps the BEST (smallest) e-value per (accession, prefix)
    - Only records accessions present in `restrict_to`
    """
    start_pat  = re.compile(r"^Scores for complete sequences")
    stop_pat   = re.compile(r"^Domain annotation for each sequence:")
    header_pat = re.compile(r"^\s*E-value\s+score\s+bias", re.I)
    sep_pat    = re.compile(r"^\s*(?:-+\s+)+-+\s*$")
    data_pat   = re.compile(r"^\s*[0-9.eE+-]+\s")
    # numeric like 0, 3.1, .5, 1e-10, 3.2E+5
    EV_NUM     = re.compile(r"^(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$")

    acc2by_prefix: Dict[str, Dict[str, Tuple[str, float]]] = defaultdict(dict)

    root = Path(hmm_dir)
    if not root.exists():
        raise FileNotFoundError(f"Search-results directory not found: {root}")

    # Zenodo preserves one directory per family. Restrict discovery to the
    # primary study output and avoid parsing FASTAs, trees, databases, or
    # exploratory *_new_model files stored beside it.
    for path_obj in sorted(root.glob("*/*_hmmsearch_E1000.txt")):
        path = str(path_obj)
        prefix = os.path.basename(path).split("_")[0].strip().lower()
        if path_obj.parent.name.lower() not in {prefix, f"{prefix}_hmm"}:
            continue
        print(f"[hmm] Scanning: {path} (prefix={prefix})")

        with open(path, "r", errors="ignore") as fh:
            in_scores = header_seen = reading_rows = False
            for raw in fh:
                line = raw.rstrip("\n")

                if stop_pat.match(line):
                    break
                if not in_scores:
                    if start_pat.match(line):
                        in_scores = True
                    continue
                if in_scores and not header_seen and header_pat.match(line):
                    header_seen = True
                    continue
                if in_scores and header_seen and not reading_rows:
                    if sep_pat.match(line):
                        reading_rows = True
                    continue
                if not reading_rows or not line.strip() or not data_pat.match(line):
                    continue

                # Robust split: works with single or multiple spaces
                parts = line.strip().split()
                if len(parts) < 9:
                    continue

                evalue_str = parts[0]
                if not EV_NUM.match(evalue_str):
                    continue
                evalue = float(evalue_str)

                full_id = parts[8]                               # e.g., "sp|Q96RI0|PAR4_HUMAN"
                desc    = " ".join(parts[9:]) if len(parts) > 9 else ""

                # Accession
                if full_id.startswith(("sp|","tr|")) and full_id.count("|") >= 2:
                    acc = full_id.split("|")[1].strip()
                else:
                    acc = full_id.strip()

                if not acc or acc not in restrict_to:
                    continue

                m_ox = OX_PAT.search(desc)
                if not m_ox:
                    continue

                taxid = m_ox.group(1)
                prev = acc2by_prefix[acc].get(prefix)
                if (prev is None) or (evalue < prev[1]):
                    acc2by_prefix[acc][prefix] = (taxid, evalue)

    print(f"[hmm] Done. Accessions with hits: {len(acc2by_prefix)}")
    return acc2by_prefix


# --------------------- IO: MMseqs m8 (family special) --------------------- #
def parse_mmseqs_m8(m8_path: str, restrict_to: Set[str]) -> Dict[str, Tuple[str, float, str]]:
    """
    Parse an MMseqs m8 file where the 3rd column contains a header with OX=.
    Returns accession -> (taxID, evalue, prefix) only for accessions in `restrict_to`.

    Accession extraction:
      - If sseqid looks like "sp|ACC|..." or "tr|ACC|...", take split('|')[1].
      - Else if header looks like that, take split('|')[1].
      - Else use sseqid as-is.
    """
    acc2info: Dict[str, Tuple[str, float, str]] = {}
    if not m8_path or not os.path.isfile(m8_path):
        return acc2info

    prefix = os.path.basename(m8_path).split("_")[0].strip().lower()
    EV_NUM = re.compile(r"^(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$")

    with open(m8_path, "r", errors="ignore") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 3:
                continue

            sseqid = parts[1].strip()
            header = parts[2].strip()
            evalue_str = parts[-1].strip()

            # Accession extraction
            def _extract_acc(token: str) -> str:
                token = token.strip()
                if token.startswith(("sp|", "tr|")) and token.count("|") >= 2:
                    return token.split("|")[1].strip()
                return token

            acc = _extract_acc(sseqid) if sseqid else None
            if (not acc or acc == sseqid) and header:
                maybe = _extract_acc(header)
                if maybe != header:
                    acc = maybe

            if not acc or acc in acc2info or acc not in restrict_to:
                continue

            tax = extract_taxid(header) or extract_taxid(line)
            if not tax:
                continue

            if not EV_NUM.match(evalue_str):
                continue
            evalue = float(evalue_str)

            acc2info[acc] = (tax, evalue, prefix)

    return acc2info


# --------------------- IO: MMseqs m8 discovery --------------------- #
def discover_m8_by_prefix(m8_dir: str) -> Dict[str, str]:
    """
    Scan `m8_dir` for *.m8 files and return { prefix_lower : m8_path }.
    If multiple files share the same prefix, the most-recently-modified file wins.
    Prefix is the filename part BEFORE the first underscore, lowercased.
      e.g.,  FSLB_DICDI_....m8  ->  prefix 'fslb'
    """
    mapping: Dict[str, str] = {}
    if not m8_dir or not os.path.isdir(m8_dir):
        return mapping

    for path_obj in Path(m8_dir).glob("*/*.m8"):
        p = str(path_obj)
        fname = path_obj.name
        prefix = fname.split("_")[0].strip().lower()
        if not path_obj.parent.name.lower().startswith(prefix):
            continue
        prev = mapping.get(prefix)
        if (prev is None) or (os.path.getmtime(p) > os.path.getmtime(prev)):
            mapping[prefix] = p
    return mapping


# --------------------- Core --------------------- #
def build_tables_from_accession_tsv(
    hmm_dir: str,
    accession_tsv: str,
    output_dir: str,
    special_m8_family_names: Optional[List[str]] = None,   # list of families to source from m8
    special_m8_search_dir: Optional[str] = None,           # directory to auto-detect *.m8
    explicit_m8_paths: Optional[Dict[str, str]] = None,    # {family_name: "/path/to/file.m8"} overrides auto-detect; data stored under family name as prefix
    accession_filter_file: Optional[str] = None,
    max_evalue: Optional[float] = None,
    output_suffix: Optional[str] = None,
    dedupe_by_best_evalue: bool = False,                   # keep one row per accession (lowest e)
) -> None:
    """
    Build two TSVs:
      1) acc_taxid_map[_{suffix}].tsv: family, accession, taxID, evalue, prefix
      2) tax_counts_per_family[_{suffix}].tsv: per-taxID counts by family

    NEW BEHAVIOR: If any (family, accession) has hits but none under the family's
    allowed prefixes (i.e., only other families' files), raise an error and stop.
    """

    # Some final families were defined as separate clades in a shared broad
    # search. Their family labels are retained in the output, while E-values
    # and taxIDs are read from the actual source search.
    shared_hmm_search_prefixes = {
        "classb1": ["classb"],
        "classb2": ["classb"],
        "vomeronasal1": ["classt"],
    }

    def resolve_allowed_prefixes(family: str) -> List[str]:
        fam = str(family).strip().lower()
        return shared_hmm_search_prefixes.get(fam, [fam])

    # --- suffix helper ---
    def _sanitize_suffix(s: Optional[str]) -> str:
        if not s:
            return ""
        cleaned = "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in s.strip())
        return cleaned

    def _with_suffix(basename: str, suffix: str) -> str:
        return f"{basename}_{suffix}.tsv" if suffix else f"{basename}.tsv"

    suffix = _sanitize_suffix(output_suffix)
    os.makedirs(output_dir, exist_ok=True)

    # TSV families and whitelist
    fam2acc = read_family_accession_tsv(accession_tsv)
    whitelist: Set[str] = set().union(*fam2acc.values()) if fam2acc else set()

    # Optional accession filter
    if accession_filter_file:
        with open(accession_filter_file) as fh:
            filter_ids = {line.strip() for line in fh if line.strip()}
        whitelist &= filter_ids
        fam2acc = {
            fam: sorted(set(accs) & filter_ids)
            for fam, accs in fam2acc.items()
            if set(accs) & filter_ids
        }

    # Parse HMMs â†’ {acc: {prefix: (tax, e)}}
    acc2info: Dict[str, Dict[str, Tuple[str, float]]] = parse_all_hmm_for_acc_taxid(
        hmm_dir, restrict_to=whitelist
    )

    # --- Auto-detect & merge m8s for the listed families ---
    added_from_m8_total = 0
    m8_root = special_m8_search_dir or hmm_dir
    if special_m8_family_names:
        m8_index = discover_m8_by_prefix(m8_root)  # {prefix_lower: path}
        for fam_name in special_m8_family_names:
            fam_key = str(fam_name).strip()
            fam_key_lc = fam_key.lower()
            if fam_key not in fam2acc:
                continue

            # Explicit path overrides auto-discovery
            if explicit_m8_paths and fam_key in explicit_m8_paths:
                m8_path = explicit_m8_paths[fam_key]
            else:
                m8_path = m8_index.get(fam_key_lc)
            if not m8_path or not os.path.isfile(m8_path):
                print(f"[m8] WARNING: No .m8 detected for family '{fam_key}' in {m8_root}")
                continue

            # For explicit-path families, store data under the family name so that
            # resolve_allowed_prefixes() finds it without any extra overrides.
            if explicit_m8_paths and fam_key in explicit_m8_paths:
                m8_prefix = fam_key_lc
            else:
                m8_prefix = os.path.basename(m8_path).split("_")[0].strip().lower()
            print(f"[m8] Scanning: {m8_path} (family={fam_key}, prefix={m8_prefix})")
            m8_map = parse_mmseqs_m8(m8_path, restrict_to=set(fam2acc[fam_key]))

            def _merge_one(acc: str, tax: str, e: float, pfx: str):
                nonlocal added_from_m8_total
                by_prefix = acc2info.setdefault(acc, {})
                prev = by_prefix.get(pfx)
                if (prev is None) or (e < prev[1]):
                    by_prefix[pfx] = (tax, e)
                    if prev is None:
                        added_from_m8_total += 1

            if isinstance(m8_map, dict):
                for acc, val in m8_map.items():
                    if isinstance(val, dict):
                        for pfx, (tax, e) in val.items():
                            _merge_one(acc, tax, float(e), str(pfx).strip().lower())
                    elif isinstance(val, tuple):
                        if len(val) == 2:
                            tax, e = val
                        elif len(val) == 3:
                            tax, e, _ = val  # ignore embedded prefix; use m8_prefix set above
                        else:
                            continue
                        _merge_one(acc, tax, float(e), m8_prefix)

    # --------------------- Build mapping (STOP on cross-family-only hits) ---------------------
    records: List[Dict[str, object]] = []
    missing_details: List[Dict[str, object]] = []

    def _classify_reason(prefmap: Optional[Dict[str, Tuple[str, float]]],
                         allowed_prefixes: Set[str],
                         threshold: Optional[float]) -> Dict[str, object]:
        if not prefmap:
            return {"reason": "no_hit_for_accession"}
        allowed_hits = [(p, tax, e) for p, (tax, e) in prefmap.items()
                        if str(p).strip().lower() in allowed_prefixes]
        if not allowed_hits:
            seen_prefixes = ",".join(sorted(str(p).strip().lower() for p in prefmap.keys()))
            return {
                "reason": "no_allowed_prefix_hit",
                "seen_prefixes": seen_prefixes,
                "allowed_prefixes": ",".join(sorted(allowed_prefixes)),
            }
        if threshold is not None and all(e > threshold for _, _, e in allowed_hits):
            return {
                "reason": "evalue_above_threshold",
                "min_evalue_allowed": min(e for _, _, e in allowed_hits),
                "threshold": threshold,
            }
        return {"reason": "unknown"}

    for fam, accs in fam2acc.items():
        allowed = set(resolve_allowed_prefixes(fam))
        for acc in accs:
            prefmap = acc2info.get(acc)

            # choose best allowed-prefix hit (respecting max_evalue if set)
            candidate = None  # (prefix, tax, e)
            if prefmap:
                for pfx, (tax, e) in prefmap.items():
                    p = str(pfx).strip().lower()
                    if p in allowed:
                        if (max_evalue is not None) and (e > max_evalue):
                            continue
                        if (candidate is None) or (e < candidate[2]):
                            candidate = (pfx, tax, e)

            if candidate is None:
                info = _classify_reason(prefmap, allowed, max_evalue)
                # NEW: hard-stop if hits exist but only under other prefixes
                if info.get("reason") == "no_allowed_prefix_hit":
                    raise RuntimeError(
                        f"Found hits only under other families for {fam}/{acc}. "
                        f"Seen prefixes: [{info.get('seen_prefixes','')}], "
                        f"allowed: [{info.get('allowed_prefixes','')}]. "
                        f"Fix filenames/prefix routing or TSV family assignment."
                    )
                # Otherwise, record as not-mapped (true missing) and continue
                info.update({"family": fam, "accession": acc})
                missing_details.append(info)
                continue

            sel_prefix, tax, evalue = candidate
            # Preserve the family-assignment prefix used by the archived table
            # even when the underlying HMMER source is a shared search.
            output_prefix = str(fam).strip().lower()
            records.append({
                "family": fam,
                "accession": acc,
                "taxID": tax,
                "evalue": evalue,
                "prefix": output_prefix,
            })

    mapping_df = pd.DataFrame.from_records(
        records, columns=["family", "accession", "taxID", "evalue", "prefix"]
    ).drop_duplicates()

    # --- Optional global dedupe ---
    if dedupe_by_best_evalue and not mapping_df.empty:
        mapping_df = (
            mapping_df.sort_values(["accession", "evalue"], ascending=[True, True])
                      .drop_duplicates(subset=["accession"], keep="first")
                      .reset_index(drop=True)
        )

    # Counts by taxID/family
    if mapping_df.empty:
        counts_df = pd.DataFrame(columns=["taxID"])
    else:
        families = sorted(mapping_df["family"].unique())
        gb = mapping_df.groupby(["taxID", "family"]).size().unstack(fill_value=0)
        gb = gb.reindex(columns=families, fill_value=0)
        gb.columns = [f"{c}_count" for c in gb.columns]
        counts_df = gb.reset_index()

    # Write outputs
    map_out = os.path.join(output_dir, _with_suffix("acc_taxid_map", suffix))
    counts_out = os.path.join(output_dir, _with_suffix("tax_counts_per_family", suffix))
    mapping_df.to_csv(map_out, sep="\t", index=False)
    counts_df.to_csv(counts_out, sep="\t", index=False)

    # Missing audit (excludes wrong-prefix casesâ€”those now raise)
    missing_out = os.path.join(output_dir, _with_suffix("missing_acc_taxid", suffix))
    true_missing_count = 0
    if missing_details:
        md = pd.DataFrame(missing_details)
        mask = md["reason"].isin(["no_hit_for_accession", "evalue_above_threshold", "unknown"])
        md_true = md.loc[mask, ["family", "accession"]].drop_duplicates()
        true_missing_count = len(md_true)
        if true_missing_count:
            md_true.to_csv(missing_out, sep="\t", index=False)

    # Report
    tsv_total = sum(len(s) for s in fam2acc.values())
    matched = len(mapping_df)
    unique_tax = mapping_df["taxID"].nunique() if not mapping_df.empty else 0

    print("\nâœ“ Finished (STRICT: only TSV-listed families & accessions).")
    print(f"Families in TSV: {len(fam2acc)}")
    print(f"Accessions in TSV (non-empty): {tsv_total}")
    print(f"Accessions mapped (unique rows): {matched}")
    if special_m8_family_names:
        print(f"M8 families requested: {len(special_m8_family_names)} | Added prefix entries: {added_from_m8_total}")
    print(f"Unique taxIDs: {unique_tax}")
    if not mapping_df.empty:
        num_with_e = mapping_df['evalue'].notna().sum()
        print(f"Rows with e-values: {num_with_e}")
    if true_missing_count:
        print(f"Missing mappings (written to): {missing_out}  [count={true_missing_count}]")
    else:
        print("Missing mappings: none.")
    print("Outputs:")
    print(f"  - {map_out}")
    print(f"  - {counts_out}")


def filter_tax_count_table_for_controls(
    input_path: str,
    output_path: str,
    min_count: int = 1,
) -> None:
    """Retain taxa with at least ``min_count`` hits for every study control."""
    counts = pd.read_csv(input_path, sep="\t")
    controls = [
        "Galpha_count", "Gbeta_count", "Ggamma_count", "RBP1_count", "TBP_count"
    ]
    missing = [column for column in controls if column not in counts.columns]
    if missing:
        raise ValueError(f"Missing required control columns: {', '.join(missing)}")

    for column in controls:
        counts[column] = pd.to_numeric(counts[column], errors="coerce").fillna(0).astype(int)

    counts = counts[counts[controls].ge(min_count).all(axis=1)].copy()
    counts.to_csv(output_path, sep="\t", index=False)
    print(f"Control-complete count table: {output_path} ({len(counts)} taxa)")


def write_mapping_overlap_plot(
    mapping_path: str,
    output_path: str,
    title: str,
    tax_count_path: Optional[str] = None,
) -> None:
    """Plot pairwise family overlaps among accessions in a mapping view."""
    from prepare_family_accession_table import compute_pairwise_overlaps, write_overlap_plot

    mapping = pd.read_csv(mapping_path, sep="\t", dtype=str)
    if tax_count_path:
        retained_taxa = set(pd.read_csv(tax_count_path, sep="\t", dtype=str)["taxID"])
        mapping = mapping[mapping["taxID"].isin(retained_taxa)]

    family_to_accessions = {
        family: set(rows["accession"])
        for family, rows in mapping.groupby("family")
    }
    overlaps = compute_pairwise_overlaps(family_to_accessions)
    write_overlap_plot(overlaps, Path(output_path), title=title)
    print(f"Accession-overlap plot: {output_path} ({len(overlaps)} non-zero pairs)")


# ================= USER VARIABLES ================= #
PROJECT_DIR = os.environ.get("GPCR_ANALYSIS_DIR", os.environ.get("GPCR_ZENODO_DIR", "/path/to/GPCR_Superfamily"))
HMM_DIR = os.environ.get(
    "GPCR_SEARCH_RESULTS_DIR", os.path.join(PROJECT_DIR, "search_results")
)
OUTPUT_DIR = os.environ.get("GPCR_FAMILY_TABLE_DIR", os.path.join(PROJECT_DIR, "families", "tables"))
ACCESSION_TSV = os.environ.get(
    "GPCR_PREFILTER_ACCESSION_TSV",
    os.path.join(OUTPUT_DIR, "accessions_by_family_prefilter_May22.tsv"),
)

# Optional: file with one accession per line to filter analysis
ACCESSION_FILTER_FILE = None
# ================================================== #

if __name__ == "__main__":
    special_m8_family_names = ["GPR1", "FSLB"]

    # View 1: all mapped extracted-clade assignments, without an E-value cutoff
    # or global best-family deduplication.
    build_tables_from_accession_tsv(
        hmm_dir=HMM_DIR,
        accession_tsv=ACCESSION_TSV,
        output_dir=OUTPUT_DIR,
        special_m8_family_names=special_m8_family_names,
        special_m8_search_dir=HMM_DIR,
        accession_filter_file=ACCESSION_FILTER_FILE,
        max_evalue=None,
        output_suffix="unfiltered",
        dedupe_by_best_evalue=False,
    )

    # View 2: retain E <= 0.1 and assign each accession only to its
    # lowest-E-value family when cross-family assignments remain.
    build_tables_from_accession_tsv(
        hmm_dir=HMM_DIR,
        accession_tsv=ACCESSION_TSV,
        output_dir=OUTPUT_DIR,
        special_m8_family_names=special_m8_family_names,
        special_m8_search_dir=HMM_DIR,
        accession_filter_file=ACCESSION_FILTER_FILE,
        max_evalue=0.1,
        output_suffix="e0.1_BestEval",
        dedupe_by_best_evalue=True,
    )

    # View 3: apply the all-five-controls-present taxon filter to View 2.
    filter_tax_count_table_for_controls(
        os.path.join(OUTPUT_DIR, "tax_counts_per_family_e0.1_BestEval.tsv"),
        os.path.join(
            OUTPUT_DIR,
            "tax_counts_per_family_e0.1_BestEval_allControls.tsv",
        ),
        min_count=1,
    )

    write_mapping_overlap_plot(
        os.path.join(OUTPUT_DIR, "acc_taxid_map_unfiltered.tsv"),
        os.path.join(OUTPUT_DIR, "overlap_plot_unfiltered.html"),
        "Pairwise accession overlaps: unfiltered",
    )
    write_mapping_overlap_plot(
        os.path.join(OUTPUT_DIR, "acc_taxid_map_e0.1_BestEval.tsv"),
        os.path.join(OUTPUT_DIR, "overlap_plot_e0.1_BestEval.html"),
        "Pairwise accession overlaps: E≤0.1 BestEval",
    )
    write_mapping_overlap_plot(
        os.path.join(OUTPUT_DIR, "acc_taxid_map_e0.1_BestEval.tsv"),
        os.path.join(OUTPUT_DIR, "overlap_plot_e0.1_BestEval_allControls.html"),
        "Pairwise accession overlaps: E≤0.1 BestEval, all controls present",
        tax_count_path=os.path.join(
            OUTPUT_DIR, "tax_counts_per_family_e0.1_BestEval_allControls.tsv"
        ),
    )
