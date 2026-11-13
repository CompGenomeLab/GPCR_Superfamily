#!/usr/bin/env python3
"""Map motif columns to real residue numbers using the full (untrimmed) family alignments.

The reference-trimmed alignments behind the superfamily logo drop columns, so counting residues in
their rows gives wrong residue numbers. The full family alignments (before reference-based
trimming) keep every residue of each row's /start-end range, so counting from the start of the
range gives the exact UniProt number at every column.

For every row of the request table (one per structure and motif position), the member's residue
at the requested full family alignment column (family_column, 1-based, from the consensus column
mapping) is read off with its exact number.

Standard library only, so it runs on the cluster as is.

Inputs:
  --request        residue_number_request.tsv (written by residue_number_request.py)
  --families-dir   directory holding <family>/<family>_filtered_genes_db_TMonly_FAMSA.fasta
Outputs:
  --out            motif_residue_numbers.tsv: the request with the member residue filled in
  --sequences-out  member_sequences.tsv: each member's numbered sequence from the full alignment,
                   used by build_superposition.py to translate UniProt numbers to structure numbers

Usage:
  python3 map_motif_residues.py --request residue_number_request.tsv \
      --families-dir $GPCR_ZENODO_DIR/10_motif_analysis/family_alignments/families \
      --out motif_residue_numbers.tsv --sequences-out member_sequences.tsv
"""
import argparse
import csv
import os
import re
import sys

PATTERN = "{fam}/{fam}_filtered_genes_db_TMonly_FAMSA.fasta"


def read_fasta(path):
    name, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:], []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def accession(header):
    parts = header.split("|")
    return parts[1] if len(parts) > 2 else header.split()[0].split("/")[0]


def numbered(header, row):
    """Row -> (start, end, intact, {column: (letter, number)})."""
    m = re.search(r"/\s*(\d+)\s*-\s*(\d+)", header)
    start, end = (int(m.group(1)), int(m.group(2))) if m else (1, None)
    cols, k = {}, start - 1
    for i, a in enumerate(row):
        if a not in "-.":
            k += 1
            cols[i] = (a.upper(), k)
    intact = end is None or k == end
    return start, end if end is not None else k, intact, cols


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--request", required=True)
    ap.add_argument("--families-dir", required=True)
    ap.add_argument("--pattern", default=PATTERN,
                    help=f"alignment path under --families-dir (default {PATTERN})")
    ap.add_argument("--out", default="motif_residue_numbers.tsv")
    ap.add_argument("--sequences-out", default="member_sequences.tsv")
    args = ap.parse_args()

    request = list(csv.DictReader(open(args.request), delimiter="\t"))
    wanted = {}
    for r in request:
        wanted.setdefault(r["family"], set()).add(r["member_accession"])

    rows = {}  # (family, accession) -> (start, end, intact, columns, sequence)
    for fam, accs in sorted(wanted.items()):
        path = os.path.join(args.families_dir, args.pattern.format(fam=fam))
        if not os.path.exists(path):
            print(f"missing alignment for {fam}: {path}", file=sys.stderr)
            continue
        found = 0
        for h, s in read_fasta(path):
            acc = accession(h)
            if acc in accs and (fam, acc) not in rows:
                start, end, intact, cols = numbered(h, s)
                rows[(fam, acc)] = (start, end, intact, cols,
                                    "".join(a.upper() for a in s if a not in "-."))
                found += 1
                if found == len(accs):
                    break
        print(f"{fam}: found {found}/{len(accs)} sequences in {path}", file=sys.stderr)

    with open(args.out, "w", newline="\n") as f:
        fields = list(request[0].keys())
        fields = [c for c in fields if c != "member_residue"] + [
            "member_residue", "member_row_intact", "status"]
        f.write("\t".join(fields) + "\n")
        for r in request:
            fam, mem = r["family"], r["member_accession"]
            out = dict(r)
            m = rows.get((fam, mem))
            if m is None:
                out.update(member_residue="", member_row_intact="",
                           status="member not in full alignment")
            elif r["family_column"] == "unmapped":
                out.update(member_residue="", member_row_intact=str(m[2]),
                           status="superfamily column not mapped in this family")
            else:
                res = m[3].get(int(r["family_column"]) - 1)
                out.update(member_residue=f"{res[0]}{res[1]}" if res else "gap",
                           member_row_intact=str(m[2]),
                           status="ok" if m[2] else "member row not intact")
            f.write("\t".join(str(out.get(c, "")) for c in fields) + "\n")

    with open(args.sequences_out, "w", newline="\n") as f:
        f.write("family\taccession\tstart\tend\trow_intact\tsequence\n")
        for (fam, acc), (start, end, intact, _, seq) in sorted(rows.items()):
            f.write(f"{fam}\t{acc}\t{start}\t{end}\t{intact}\t{seq}\n")
    print(f"wrote {args.out} and {args.sequences_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
