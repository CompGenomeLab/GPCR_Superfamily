"""For every family alignment, write the rows whose accessions are in a superfamily
alignment, with every column kept (no trimming), for mapping family columns to
superfamily columns through all representatives instead of one anchor.

Family alignments are <families_dir>/<family>/<family>_filtered_genes_db_TMonly_FAMSA.fasta;
output is <out_dir>/<family>_family_reps_TMonly_FAMSA.fasta, rows copied unchanged so the
column count equals the family alignment's (every row is checked against it). Accessions
are the second '|' field of the header.

python extract_family_representatives.py <superfamily.fasta> <families_dir> <out_dir>
"""
import sys
from pathlib import Path

PATTERN = "*/*_filtered_genes_db_TMonly_FAMSA.fasta"


def read_fasta(path):
    header, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header, seq = line[1:], []
            elif line:
                seq.append(line)
    if header is not None:
        yield header, "".join(seq)


def accession(header):
    parts = header.split("|")
    return parts[1].strip() if len(parts) > 1 else header.split()[0]


def main():
    superfamily, families_dir, out_dir = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
    reps = {accession(h) for h, _ in read_fasta(superfamily)}
    print(f"{len(reps)} representatives in {superfamily}")
    found_all = set()
    for aln in sorted(families_dir.glob(PATTERN)):
        family = aln.parent.name
        found, ncol = set(), None
        with open(out_dir / f"{family}_family_reps_TMonly_FAMSA.fasta", "w") as out:
            for h, s in read_fasta(aln):
                if ncol is None:
                    ncol = len(s)
                if len(s) != ncol:
                    sys.exit(f"{aln}: row {h.split()[0]} has {len(s)} columns, alignment has {ncol}")
                acc = accession(h)
                if acc in reps and acc not in found:
                    out.write(f">{h}\n{s}\n")
                    found.add(acc)
        found_all |= found
        print(f"{family}\t{len(found)} reps\t{ncol} columns")
    print(f"found {len(found_all)} of {len(reps)}; missing: {' '.join(sorted(reps - found_all))}")


if __name__ == "__main__":
    main()
