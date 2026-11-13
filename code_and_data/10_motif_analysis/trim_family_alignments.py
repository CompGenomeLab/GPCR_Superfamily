#!/usr/bin/env python3
"""Trim each family's residue alignment to the columns its logo anchors occupy.

For every row of the anchor table (`acc1`, `family`, optional `acc2`, as written by
`write_logo_anchors.py`) the family alignment is read twice. The first pass captures
the two reference rows and builds a column mask keeping every column where acc1 or
acc2 has a residue; the second streams every sequence through that mask. Sequences
are never dropped and headers are never changed, and because a kept column is one a
reference occupies, the references lose no residues. Other sequences lose only the
insertions relative to both references, which no superfamily column could be mapped
to anyway. Two passes keep memory flat on class A and Olfactory.

The references are checked before anything is written. The logo precompute locates a
reference residue by counting residues from the start of the `/start-end` range in
its header, so the reference must hold exactly that range's residues. An alignment
whose columns were already gap-trimmed, such as a `_gt0.1` file, has lost residues
from its references wherever a gappy column held one, and every position after the
first loss would be mapped to the wrong residue without any error. Such a family is
refused, and `--check` runs only this test.

Output is `<out-dir>/<family>_genes_filtered_db_FAMSA.ref_trimmed.fasta`, the name
the website's `precompute-mappings.mjs` looks for.

Usage:
    python trim_family_alignments.py \\
        --anchors superfamily_logo_anchors_Sep24_topo0.8.tsv \\
        --out-dir logo_trimmed_msas

The family alignments default to the TM-only FAMSA alignments written by family_alignment.sh,
`$FAMILIES_DIR/<family>/<family>_filtered_genes_db_TMonly_FAMSA.fasta` (default `families/`); pass
`--alignment-pattern` with `{family}` in it to read them from anywhere else.
"""

import argparse
import csv
import os
import re
import sys

GAP = set("-.")
RANGE_RE = re.compile(r"/(\d+)-(\d+)")
OUT_SUFFIX = "_genes_filtered_db_FAMSA.ref_trimmed.fasta"
DEFAULT_PATTERN = os.path.join(os.environ.get("FAMILIES_DIR", "families"),
                               "{family}", "{family}_filtered_genes_db_TMonly_FAMSA.fasta")


def iter_fasta(path):
    """Yield (header without '>', sequence) from a FASTA file, streaming."""
    with open(path, "r", encoding="utf-8", buffering=1024 * 1024) as handle:
        header, chunks = None, []
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:].strip(), []
            else:
                line = line.strip()
                if line:
                    chunks.append(line)
        if header is not None:
            yield header, "".join(chunks)


def accession_of(header):
    """`tr|ACC|NAME/start-end ...` -> ACC."""
    token = RANGE_RE.sub("", header.split()[0]) if header.split() else header
    fields = token.split("|")
    return fields[1].strip() if len(fields) >= 3 else token


def check_reference(label, accession, header, seq):
    """Return a problem description, or None when the row holds its whole range."""
    residues = sum(1 for c in seq if c not in GAP)
    match = RANGE_RE.search(header.split()[0])
    if not match:
        return None  # full length; residue numbering starts at 1 and cannot be checked
    span = int(match.group(2)) - int(match.group(1)) + 1
    if residues != span:
        return (f"{label} {accession} holds {residues} residues but its header range "
                f"{match.group(0)} spans {span}; the alignment has lost reference residues")
    return None


def collect_references(path, acc1, acc2):
    ref1 = ref2 = None
    for header, seq in iter_fasta(path):
        accession = accession_of(header)
        if accession == acc1 and ref1 is None:
            ref1 = (header, seq)
        elif acc2 and accession == acc2 and ref2 is None:
            ref2 = (header, seq)
        if ref1 is not None and (ref2 is not None or not acc2):
            break
    return ref1, ref2


def build_mask(seq1, seq2):
    if seq2 is not None and len(seq1) != len(seq2):
        raise ValueError(f"reference rows differ in width: {len(seq1)} vs {len(seq2)}")
    if seq2 is None:
        return [c not in GAP for c in seq1]
    return [a not in GAP or b not in GAP for a, b in zip(seq1, seq2)]


def write_wrapped(handle, seq, width):
    if width > 0:
        for i in range(0, len(seq), width):
            handle.write(seq[i:i + width] + "\n")
    else:
        handle.write(seq + "\n")


def process(family, acc1, acc2, args):
    path = args.alignment_pattern.format(family=family)
    if not os.path.exists(path):
        return f"no alignment at {path}"

    ref1, ref2 = collect_references(path, acc1, acc2)
    if ref1 is None:
        return f"acc1 {acc1} not in {os.path.basename(path)}"
    if acc2 and ref2 is None:
        return f"acc2 {acc2} not in {os.path.basename(path)}"

    problems = [p for p in (check_reference("acc1", acc1, *ref1),
                            check_reference("acc2", acc2, *ref2) if ref2 else None) if p]
    if problems:
        return "; ".join(problems)
    if args.check:
        print(f"[OK] {family}: references intact")
        return None

    mask = build_mask(ref1[1], ref2[1] if ref2 else None)
    keep = [i for i, k in enumerate(mask) if k]
    os.makedirs(args.out_dir, exist_ok=True)
    out_path = os.path.join(args.out_dir, family + OUT_SUFFIX)
    tmp_path = out_path + ".partial"
    n = 0
    problem = None
    # Written under a temporary name first, so a failure never leaves a partial file
    # under the name the website reads.
    with open(tmp_path, "w", encoding="utf-8", newline="\n", buffering=1024 * 1024) as out:
        for header, seq in iter_fasta(path):
            if len(seq) != len(mask):
                problem = (f"{accession_of(header)} is {len(seq)} columns wide, "
                           f"the references {len(mask)}")
                break
            out.write(f">{header}\n")
            write_wrapped(out, "".join(seq[i] for i in keep), args.wrap)
            n += 1
    if problem:
        os.remove(tmp_path)
        return problem
    os.replace(tmp_path, out_path)
    print(f"[OK] {family}: {n} sequences, kept {len(keep)} of {len(mask)} columns -> {out_path}")
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--anchors", required=True, help="acc1/family/acc2 table")
    parser.add_argument("--alignment-pattern", default=DEFAULT_PATTERN,
                        help="path of each residue family alignment, not gap-trimmed, with "
                             "{family} where the family name goes (default: %(default)s)")
    parser.add_argument("--out-dir", default="logo_trimmed_msas")
    parser.add_argument("--wrap", type=int, default=80, help="FASTA line width, 0 for none")
    parser.add_argument("--check", action="store_true",
                        help="only verify that every reference holds its full range")
    args = parser.parse_args()

    with open(args.anchors, newline="") as handle:
        rows = [{k.strip().lower(): (v or "").strip() for k, v in r.items()}
                for r in csv.DictReader(handle, delimiter="\t")]

    failures = 0
    for row in rows:
        family, acc1, acc2 = row.get("family", ""), row.get("acc1", ""), row.get("acc2", "")
        if not family or not acc1:
            continue
        problem = process(family, acc1, acc2 or None, args)
        if problem:
            failures += 1
            print(f"[FAIL] {family}: {problem}", file=sys.stderr)
    print(f"\n{len(rows) - failures} of {len(rows)} families passed")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
