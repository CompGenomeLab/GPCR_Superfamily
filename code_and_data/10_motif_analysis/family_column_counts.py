"""Residue counts in every column of every family alignment, so the logo can be built on any
family column (the consensus mapping's), not only on the anchor's.

For each <families_dir>/<family>/<family>_filtered_genes_db_TMonly_FAMSA.fasta writes
<out_dir>/<family>_column_counts.tsv: one row per column and group with at least one
residue, giving the group's sequence count (gaps included) and the count of each residue.
Group 'all' is the whole family; --groups adds the groups of a TSV (accession, family,
group), e.g. the Class C and cAMP kingdoms. --refs writes the rows of the given accessions
(the GPCRdb references) to <out_dir>/<family>_refs_TMonly_FAMSA.fasta, all columns kept.

python family_column_counts.py <families_dir> <out_dir> [--groups kingdom_groups.tsv] [--refs P25021 ...]
"""
import argparse
from collections import Counter, defaultdict
from pathlib import Path

PATTERN = "*/*_filtered_genes_db_TMonly_FAMSA.fasta"
AA = "ACDEFGHIKLMNPQRSTVWY"


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
    ap = argparse.ArgumentParser()
    ap.add_argument("families_dir", type=Path)
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("--groups")
    ap.add_argument("--refs", nargs="*", default=[])
    a = ap.parse_args()

    groups = defaultdict(dict)
    if a.groups:
        for line in open(a.groups):
            acc, family, group = line.rstrip("\n").split("\t")
            groups[family][acc] = group
    refs = set(a.refs)

    for aln in sorted(a.families_dir.glob(PATTERN)):
        family = aln.parent.name
        counts = defaultdict(lambda: defaultdict(Counter))  # group -> column -> residue counts
        nseq, ncol, ref_rows = Counter(), None, []
        for h, s in read_fasta(aln):
            if ncol is None:
                ncol = len(s)
            if len(s) != ncol:
                raise SystemExit(f"{aln}: row {h.split()[0]} has {len(s)} columns, alignment has {ncol}")
            acc = accession(h)
            if acc in refs:
                ref_rows.append(f">{h}\n{s}\n")
            for g in ("all", groups[family].get(acc)):
                if g is None:
                    continue
                nseq[g] += 1
                for i, aa in enumerate(s.upper()):
                    if aa in AA:
                        counts[g][i][aa] += 1
        with open(a.out_dir / f"{family}_column_counts.tsv", "w") as out:
            out.write("column\tgroup\tn_sequences\t" + "\t".join(AA) + "\n")
            for g in counts:
                for i in sorted(counts[g]):
                    out.write(f"{i + 1}\t{g}\t{nseq[g]}\t" + "\t".join(str(counts[g][i][x]) for x in AA) + "\n")
        if ref_rows:
            (a.out_dir / f"{family}_refs_TMonly_FAMSA.fasta").write_text("".join(ref_rows))
        print(f"{family}\t{ncol} columns\t" + "\t".join(f"{g}={n}" for g, n in nseq.items())
              + f"\trefs={len(ref_rows)}")


if __name__ == "__main__":
    main()
