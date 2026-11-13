#!/usr/bin/env python3
"""Extract ranked target accessions from standard hmmsearch text output."""

from __future__ import annotations

import argparse
from pathlib import Path


def normalize_identifier(token: str) -> str:
    """Return the accession from a UniProt-style token, or the token itself."""
    parts = token.split("|")
    return parts[1] if len(parts) >= 3 and parts[1] else token


def extract_accessions(hmm_file: Path, limit: int, output: Path) -> int:
    count = 0
    with hmm_file.open(encoding="utf-8", errors="replace") as source, output.open(
        "w", encoding="utf-8"
    ) as destination:
        for line in source:
            if not line.startswith(">>"):
                continue
            token = line[2:].strip().split(maxsplit=1)[0]
            destination.write(normalize_identifier(token) + "\n")
            count += 1
            if count >= limit:
                break
    return count


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("hmm_output", type=Path)
    parser.add_argument("limit", type=int)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    if args.limit < 1:
        parser.error("limit must be a positive integer")
    output = args.output or args.hmm_output.with_name(
        f"{args.hmm_output.stem}_headers.txt"
    )
    count = extract_accessions(args.hmm_output, args.limit, output)
    print(f"Wrote {count} accessions to {output}")


if __name__ == "__main__":
    main()

