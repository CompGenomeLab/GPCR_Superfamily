"""GPCRdb labels for superfamily columns by majority vote of every human receptor with GPCRdb
numbering in each family alignment, instead of the single reference (acc2).

A family's voters are the human rows of its trimmed alignment (family_alignments/consensus_trimmed_msas/)
whose gene has a GPCRdb numbering file (website_inputs/conservation_files/<GENE>_conservation.txt,
column 6) and whose row kept every residue of its /start-end range, so residue numbers can be
counted from the trimmed row. A superfamily column takes the family column assigned by the
consensus mapping (consensus_mapping/<alignment>_conf<threshold><tag>.tsv); every voter
with a residue there votes with its generic number (NxNN, loops NNxNN; plain residue numbers
outside the 7TM are ignored). The family label is the majority, its agreement the majority's
share of the voters. Each family is labelled in its own GPCRdb scheme.

python label_by_human_receptors.py <alignment> [mapping tag, e.g. _top5]
Writes gpcrdb_labels_by_human_receptors/<alignment><tag>.tsv: sup_column, family, label,
agreement, n_voters, top 3 labels.
"""
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

from consensus_column_mapping import DATA, HERE, OUT, PUBLIC, read_fasta

TRIMMED = DATA / "family_alignments" / "consensus_trimmed_msas"
NUMBERS = PUBLIC / "conservation_files"
GENERIC = re.compile(r"^\d+x\d+$")
RANGE = re.compile(r"/(\d+)-(\d+)")


def gpcrdb_numbers(gene):
    out = {}
    for line in (NUMBERS / f"{gene}_conservation.txt").read_text().splitlines()[1:]:
        cols = line.split("\t")
        if len(cols) >= 6 and GENERIC.match(cols[5].strip()):
            out[int(cols[0])] = cols[5].strip()
    return out


def family_votes(family):
    """{full family column (0-based): Counter(label)} over the family's usable human receptors."""
    genes = {p.name.split("_conservation")[0] for p in NUMBERS.glob("*_conservation.txt")}
    kept = [int(x) - 1 for x in (TRIMMED / f"{family}_kept_columns.txt").read_text().split()]
    votes, used = defaultdict(Counter), 0
    for h, s in read_fasta(TRIMMED / f"{family}_genes_filtered_db_FAMSA.consensus_trimmed.fasta"):
        parts = h.split("|")
        if len(parts) < 3 or not parts[2].split("_")[-1].startswith("HUMAN"):
            continue
        gene = parts[2].split("_")[0]
        m = RANGE.search(h)
        if gene not in genes or not m:
            continue
        start, end = int(m.group(1)), int(m.group(2))
        if sum(c not in "-." for c in s) != end - start + 1:
            continue  # lost residues to trimming: numbering would drift
        numbers, n = gpcrdb_numbers(gene), start - 1
        used += 1
        for t, c in enumerate(s):
            if c not in "-.":
                n += 1
                if n in numbers:
                    votes[kept[t]][numbers[n]] += 1
    return votes, used


def main(alignment, tag=""):
    mapping = OUT / f"{alignment}_conf0.33{tag}.tsv"
    assigned = defaultdict(dict)
    for r in csv.DictReader(open(mapping), delimiter="\t"):
        if r["consensus_family_column"]:
            assigned[r["family"]][int(r["sup_column"])] = int(r["consensus_family_column"]) - 1
    out_dir = DATA / "gpcrdb_labels"
    out_dir.mkdir(exist_ok=True)
    rows = []
    for family in sorted(assigned):
        if not (TRIMMED / f"{family}_kept_columns.txt").exists():
            continue
        votes, used = family_votes(family)
        if not used:
            continue
        print(f"{family}: {used} human receptors vote")
        for c, f in sorted(assigned[family].items()):
            v = votes.get(f)
            if v:
                lab, k = v.most_common(1)[0]
                n = sum(v.values())
                rows.append((c, family, lab, k / n, n, v.most_common(3)))
    with open(out_dir / f"{alignment}{tag}.tsv", "w") as fh:
        fh.write("sup_column\tfamily\tlabel\tagreement\tn_voters\ttop_labels\n")
        for c, family, lab, a, n, top in rows:
            fh.write(f"{c}\t{family}\t{lab}\t{a:.2f}\t{n}\t{top}\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "")
