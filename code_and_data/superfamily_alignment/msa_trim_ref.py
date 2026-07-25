#!/usr/bin/env python3
"""
Memory-safe, two-pass streaming trimmer for huge MSAs.

Input TSV (tab-delimited) must have columns: acc1, family, (optional) acc2.
For each row, the script reads the family's FAMSA alignment through the
configured family-data directory, then writes the reference-trimmed alignment
to the configured trimmed-MSA directory.
"""

import os
import sys
import csv
import gzip
import argparse
from typing import Optional, Tuple, Iterable

PROJECT_DIR = os.environ.get("GPCR_ANALYSIS_DIR", os.environ.get("GPCR_ZENODO_DIR", "/path/to/GPCR_Superfamily"))
BASE_DIR = os.environ.get("GPCR_FAMILY_DIR", os.path.join(PROJECT_DIR, "families"))
OUT_DIR  = os.environ.get("GPCR_TRIMMED_MSA_DIR", os.path.join(BASE_DIR, "trimmed_msas"))
GAP = set("-.")  # gap chars


def parse_args():
    ap = argparse.ArgumentParser(description="Stream-trim MSA columns using acc1/acc2 from TSV (two-pass, low RAM).")
    ap.add_argument("--mapping", required=True, help="TSV with header: acc1, family, [acc2]")
    ap.add_argument("--strict", action="store_true", help="Fail on first error.")
    ap.add_argument("--wrap", type=int, default=80, help="FASTA line wrap (0 = no wrap). Default: 80")
    return ap.parse_args()


def iter_fasta_records(path: str) -> Iterable[Tuple[str, str]]:
    """Yield (header_without_gt, sequence_str) streaming from file."""
    opener = gzip.open if path.endswith(".gz") else open
    kwargs = {"encoding": "utf-8"} if path.endswith(".gz") else {"encoding": "utf-8", "buffering": 1024*1024}
    with opener(path, "rt", **kwargs) as fh:
        header = None
        chunks = []
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].rstrip("\n").strip()
                chunks = []
            else:
                s = line.strip()
                if s:
                    chunks.append(s)
        if header is not None:
            yield header, "".join(chunks)


def accession_from_header(h: str) -> str:
    parts = h.split("|")
    return parts[1].strip() if len(parts) > 1 else ""


def collect_ref_sequences(path: str, acc1: str, acc2: Optional[str]) -> Tuple[Optional[str], Optional[str]]:
    """Pass 1: stream file and capture only the reference sequence(s)."""
    ref1 = ref2 = None
    for hdr, seq in iter_fasta_records(path):
        acc = accession_from_header(hdr)
        if acc1 and acc == acc1:
            ref1 = seq
        elif acc2 and acc == acc2:
            ref2 = seq
        # Stop scanning after both reference sequences are found.
        if ref1 is not None and ((acc2 and ref2 is not None) or not acc2):
            break
    return ref1, ref2


def build_keep_mask(ref1: str, ref2: Optional[str]) -> bytearray:
    """Return bytearray mask with 1=keep, 0=drop based on gaps in provided refs."""
    if ref2 is not None and len(ref1) != len(ref2):
        raise ValueError(f"Reference sequence lengths differ: {len(ref1)} vs {len(ref2)}")
    L = len(ref1)
    mask = bytearray(L)
    if ref2 is None:
        # keep where ref1 is non-gap
        for i, c in enumerate(ref1):
            mask[i] = 1 if c not in GAP else 0
    else:
        for i, (a, b) in enumerate(zip(ref1, ref2)):
            mask[i] = 1 if (a not in GAP or b not in GAP) else 0
    return mask


def write_wrapped(outfh, seq: str, width: int) -> None:
    if width and width > 0:
        for i in range(0, len(seq), width):
            outfh.write(seq[i:i+width] + "\n")
    else:
        outfh.write(seq + "\n")


def trim_sequence_streaming(seq: str, mask: bytearray) -> str:
    if len(seq) != len(mask):
        raise ValueError("Sequence length does not match mask length.")
    # Build trimmed seq by filtering with mask (fast join on list)
    return "".join(ch for ch, keep in zip(seq, mask) if keep)


def trim_alignment_two_pass(in_path: str, out_path: str, acc1: str, acc2: Optional[str], wrap: int) -> Tuple[int, int]:
    """
    Returns (n_kept_columns, n_dropped_columns).
    """
    # PASS 1: collect refs
    ref1, ref2 = collect_ref_sequences(in_path, acc1, acc2)
    if ref1 is None:
        raise KeyError(f"acc1 '{acc1}' not found in {in_path}")
    if acc2 and ref2 is None:
        raise KeyError(f"acc2 '{acc2}' not found in {in_path}")

    mask = build_keep_mask(ref1, ref2)
    kept = int(sum(mask))
    dropped = len(mask) - kept

    # PASS 2: stream-trim all sequences and write
    with open(out_path, "w", encoding="utf-8", buffering=1024*1024) as out:
        for hdr, seq in iter_fasta_records(in_path):
            out.write(f">{hdr}\n")
            trimmed = trim_sequence_streaming(seq, mask)
            write_wrapped(out, trimmed, wrap)

    return kept, dropped


def process_row(family: str, acc1: str, acc2: Optional[str], wrap: int) -> bool:
    in_path = os.path.join(BASE_DIR, family, f"{family}_genes_filtered_db_FAMSA.fasta")
    if not os.path.exists(in_path) and os.path.exists(in_path + ".gz"):
        in_path += ".gz"
    if not os.path.exists(in_path):
        print(f"[WARN] Missing alignment for {family}: {in_path}", file=sys.stderr)
        return False

    os.makedirs(OUT_DIR, exist_ok=True)
    out_path = os.path.join(OUT_DIR, f"{family}_genes_filtered_db_FAMSA.ref_trimmed.fasta")

    kept, dropped = trim_alignment_two_pass(in_path, out_path, acc1, acc2 or None, wrap)
    print(f"[OK] {family}: wrote {out_path}  (kept {kept}, dropped {dropped} columns)")
    return True


def main():
    args = parse_args()

    with open(args.mapping, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            print("[ERROR] TSV appears empty or lacks a header.", file=sys.stderr)
            sys.exit(2)

        required = {"acc1", "family"}
        header = {h.lower(): h for h in reader.fieldnames if h}
        have = set(header.keys())
        if not required.issubset(have):
            print(f"[ERROR] TSV must contain columns: {sorted(required)}; found: {sorted(have)}", file=sys.stderr)
            sys.exit(2)

        ok = 0
        for row in reader:
            r = {k.lower(): (v or "").strip() for k, v in row.items() if k}
            family = r.get("family", "")
            acc1   = r.get("acc1", "")
            acc2   = r.get("acc2", "")

            if not family or not acc1:
                msg = f"[WARN] Skipping row with missing family/acc1: family='{family}' acc1='{acc1}'"
                if args.strict:
                    print(msg, file=sys.stderr)
                    sys.exit(1)
                print(msg, file=sys.stderr)
                continue

            try:
                if process_row(family, acc1, acc2, args.wrap):
                    ok += 1
            except Exception as e:
                print(f"[ERROR] family='{family}' acc1='{acc1}' acc2='{acc2}': {e}", file=sys.stderr)
                if args.strict:
                    sys.exit(1)

    if ok == 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
