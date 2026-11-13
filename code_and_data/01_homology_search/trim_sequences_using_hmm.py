#!/usr/bin/env python3
"""Extract the combined HMM-aligned span from every retrieved target sequence."""

from __future__ import annotations

import argparse
from pathlib import Path

from Bio import SeqIO


NO_DOMAIN_MESSAGE = (
    "[No individual domains that satisfy reporting thresholds "
    "(although complete target did)]"
)


def parse_domain_spans(hmm_output: Path) -> dict[str, tuple[int, int]]:
    """Parse 1-based inclusive ali-from/ali-to coordinates from HMMER output."""
    spans: dict[str, tuple[int, int]] = {}
    current_id: str | None = None
    domain_table_row = 0

    with hmm_output.open(encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line.startswith(">>"):
                current_id = line[2:].strip().split(maxsplit=1)[0]
                domain_table_row = 0
                continue
            if current_id is None:
                continue
            if NO_DOMAIN_MESSAGE in line:
                current_id = None
                continue

            domain_table_row += 1
            if domain_table_row < 3:
                continue
            fields = line.split()
            if len(fields) <= 13:
                continue
            try:
                # HMMER domain-table columns 12/13 are envfrom/envto. The
                # envelope includes the aligned domain plus plausible flanking
                # residues and is the boundary used in the original analysis.
                env_from, env_to = int(fields[12]), int(fields[13])
            except ValueError:
                continue

            old_from, old_to = spans.get(current_id, (env_from, env_to))
            spans[current_id] = (min(old_from, env_from), max(old_to, env_to))

    # Add accession-only aliases for UniProt-style identifiers. This supports
    # MMseqs2 databases whose FASTA conversion shortens `tr|ACC|NAME` to `ACC`.
    for identifier, span in list(spans.items()):
        parts = identifier.split("|")
        if len(parts) >= 3 and parts[1]:
            spans.setdefault(parts[1], span)
    return spans


def candidate_ids(record_id: str) -> tuple[str, ...]:
    """Support both complete UniProt tokens and accession-only FASTA headers."""
    parts = record_id.split("|")
    if len(parts) >= 3 and parts[1]:
        return record_id, parts[1]
    return (record_id,)


def write_regions(
    fasta_file: Path,
    spans: dict[str, tuple[int, int]],
    output: Path,
    minimum_length: int,
) -> tuple[int, int]:
    matched = written = 0
    with output.open("w", encoding="utf-8") as destination:
        for record in SeqIO.parse(str(fasta_file), "fasta"):
            span = next((spans[key] for key in candidate_ids(record.id) if key in spans), None)
            if span is None:
                continue
            matched += 1
            ali_from, ali_to = span
            region = record.seq[ali_from - 1 : ali_to]
            if len(region) < minimum_length:
                continue
            record.seq = region
            SeqIO.write(record, destination, "fasta")
            written += 1
    return matched, written


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("hmm_output", type=Path)
    parser.add_argument("fasta", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--minimum-length", type=int, default=20)
    args = parser.parse_args()

    output = args.output or args.fasta.with_name(f"{args.fasta.stem}_hmmregions.fasta")
    spans = parse_domain_spans(args.hmm_output)
    matched, written = write_regions(args.fasta, spans, output, args.minimum_length)
    print(f"Parsed spans for {len(spans)} targets")
    print(f"Matched {matched} FASTA records; wrote {written} regions to {output}")


if __name__ == "__main__":
    main()
