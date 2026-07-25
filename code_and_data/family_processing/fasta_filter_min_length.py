#!/usr/bin/env python3
"""
Filter an alignment by per-sequence ungapped length.

Outputs (auto-named from the input file):
  - <input>_filtered.fasta    (sequences kept)
  - <input>_removed.fasta     (sequences removed)
  - <input>_filtered_ids.txt  (accession IDs of kept sequences; header split by '|' -> second token)

A residue is NON-GAP iff it is a nucleotide or amino-acid IUPAC code (case-insensitive).
Everything else (e.g., '-', '.', '?', '*', digits, spaces) is treated as a gap.

Example:
    python filter_alignment_by_ungapped_length.py \
        --input aln.fasta \
        --min-ungapped-length 300
"""

import argparse
import os
import sys
from typing import List, Tuple

# IUPAC nucleotide codes (DNA/RNA), including ambiguity codes.
_IUPAC_NT = set("ACGTURYSWKMBDHVN")
# Amino-acid codes including common extensions (B, Z, J) and U (Sec), O (Pyl), X (unknown).
_IUPAC_AA = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")
_VALID_RESIDUES = _IUPAC_NT | _IUPAC_AA

def _is_residue(c: str) -> bool:
    """True if c is a nucleotide or amino-acid code (case-insensitive)."""
    return c.upper() in _VALID_RESIDUES

def read_fasta(path: str) -> List[Tuple[str, str]]:
    records = []
    with open(path, "r") as f:
        header = None
        seq_chunks = []
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_chunks)))
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            records.append((header, "".join(seq_chunks)))
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records

def write_fasta(path: str, records: List[Tuple[str, str]]) -> None:
    with open(path, "w") as out:
        for header, seq in records:
            out.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")

def derive_output_paths(input_path: str):
    """
    Return (filtered_fasta, removed_fasta, filtered_ids_txt)
    - If input ends with '.fasta' (case-insensitive), strip it for the stem.
    - Else, use the full basename as stem.
    """
    d = os.path.dirname(input_path)
    base = os.path.basename(input_path)
    stem = base[:-6] if base.lower().endswith(".fasta") else base  # remove ".fasta"
    filtered = os.path.join(d, f"{stem}_filtered.fasta")
    removed = os.path.join(d, f"{stem}_removed.fasta")
    ids_txt = os.path.join(d, f"{stem}_filtered_ids.txt")
    return filtered, removed, ids_txt

def parse_args():
    p = argparse.ArgumentParser(
        description="Keep sequences whose ungapped length >= threshold."
    )
    p.add_argument("-i", "--input", required=True, help="Input FASTA alignment path")
    p.add_argument("--min-ungapped-length", type=int, required=True,
                   help="Minimum ungapped length required to KEEP a sequence (integer = 0)")
    return p.parse_args()

def main():
    args = parse_args()
    if args.min_ungapped_length < 0:
        raise ValueError("--min-ungapped-length must be = 0")

    filtered_fa, removed_fa, ids_txt = derive_output_paths(args.input)

    records = read_fasta(args.input)
    headers = [h for h, _ in records]
    seqs = [s for _, s in records]

    # Compute ungapped lengths and split into kept/removed.
    kept_idx = []
    removed_idx = []
    for i, s in enumerate(seqs):
        ungapped_len = sum(1 for c in s if _is_residue(c))
        if ungapped_len >= args.min_ungapped_length:
            kept_idx.append(i)
        else:
            removed_idx.append(i)

    kept_records = [(headers[i], seqs[i]) for i in kept_idx]
    removed_records = [(headers[i], seqs[i]) for i in removed_idx]

    # Write outputs
    write_fasta(filtered_fa, kept_records)
    write_fasta(removed_fa, removed_records)

    with open(ids_txt, "w") as ids_out:
        for i in kept_idx:
            header = headers[i]
            parts = header.split("|")
            accession = parts[1].strip() if len(parts) >= 2 else header.strip()
            ids_out.write(accession + "\n")

    # Summary to stderr
    print(f"Input sequences: {len(seqs)}", file=sys.stderr)
    print(f"Min ungapped length: {args.min_ungapped_length}", file=sys.stderr)
    print(f"Kept sequences: {len(kept_idx)}", file=sys.stderr)
    print(f"Removed sequences: {len(removed_idx)}", file=sys.stderr)
    print(f"Wrote kept sequences to: {filtered_fa}", file=sys.stderr)
    print(f"Wrote removed sequences to: {removed_fa}", file=sys.stderr)
    print(f"Wrote kept accessions to: {ids_txt}", file=sys.stderr)

if __name__ == "__main__":
    main()

