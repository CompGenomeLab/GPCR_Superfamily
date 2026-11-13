#!/usr/bin/env python3
"""Trim each family's residue alignment to the columns of its GPCRdb numbering reference
plus the family columns the consensus mapping assigned (consensus_column_mapping.py).

The column table (family, reference, columns) is written locally from the consensus mappings:
`reference` is acc2, or acc1 where the family has no GPCRdb reference, and `columns` lists the
1-based family alignment columns assigned in any of the 8 superfamily alignments. A column is
kept where the reference has a residue or where it is listed, so the reference loses no
residues; this is checked as in trim_family_alignments.py, whose helpers are reused.
Sequences are never dropped and headers never changed.

Writes <out-dir>/<family>_genes_filtered_db_FAMSA.consensus_trimmed.fasta and
<out-dir>/<family>_kept_columns.txt (the 1-based family alignment column of each kept column,
in order), with which the logo builder maps family columns to trimmed ones.

Usage:
    python3.9 trim_family_alignments_consensus.py \\
        --columns consensus_trim_columns_conf0.33_Oct10.tsv --out-dir consensus_trimmed_msas
"""

import argparse
import csv
import os
import sys

from trim_family_alignments import (DEFAULT_PATTERN, GAP, accession_of, check_reference, iter_fasta,
                                    write_wrapped)

OUT_SUFFIX = "_genes_filtered_db_FAMSA.consensus_trimmed.fasta"


def process(family, reference, listed, args):
    path = args.alignment_pattern.format(family=family)
    if not os.path.exists(path):
        return f"no alignment at {path}"
    ref = next(((h, s) for h, s in iter_fasta(path) if accession_of(h) == reference), None)
    if ref is None:
        return f"reference {reference} not in {os.path.basename(path)}"
    problem = check_reference("reference", reference, *ref)
    if problem:
        return problem
    width = len(ref[1])
    if listed and max(listed) > width:
        return f"listed column {max(listed)} beyond the alignment's {width} columns"
    keep = [i for i, c in enumerate(ref[1]) if c not in GAP or i + 1 in listed]

    os.makedirs(args.out_dir, exist_ok=True)
    out_path = os.path.join(args.out_dir, family + OUT_SUFFIX)
    tmp_path = out_path + ".partial"
    n = 0
    with open(tmp_path, "w", encoding="utf-8", newline="\n", buffering=1024 * 1024) as out:
        for header, seq in iter_fasta(path):
            if len(seq) != width:
                out.close()
                os.remove(tmp_path)
                return f"{accession_of(header)} is {len(seq)} columns wide, the reference {width}"
            out.write(f">{header}\n")
            write_wrapped(out, "".join(seq[i] for i in keep), args.wrap)
            n += 1
    os.replace(tmp_path, out_path)
    with open(os.path.join(args.out_dir, f"{family}_kept_columns.txt"), "w", newline="\n") as out:
        out.write("\n".join(str(i + 1) for i in keep) + "\n")
    n_ref = sum(c not in GAP for c in ref[1])
    print(f"[OK] {family}: {n} sequences, kept {len(keep)} of {width} columns "
          f"({n_ref} reference, {len(listed)} listed) -> {out_path}")
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--columns", required=True, help="family / reference / columns table")
    parser.add_argument("--alignment-pattern", default=DEFAULT_PATTERN,
                        help="path of each residue family alignment, not gap-trimmed, with "
                             "{family} where the family name goes (default: %(default)s)")
    parser.add_argument("--out-dir", default="consensus_trimmed_msas")
    parser.add_argument("--wrap", type=int, default=80, help="FASTA line width, 0 for none")
    args = parser.parse_args()

    with open(args.columns, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    failures = 0
    for row in rows:
        listed = {int(c) for c in row["columns"].split(",") if c}
        problem = process(row["family"], row["reference"], listed, args)
        if problem:
            failures += 1
            print(f"[FAIL] {row['family']}: {problem}", file=sys.stderr)
    print(f"\n{len(rows) - failures} of {len(rows)} families passed")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
