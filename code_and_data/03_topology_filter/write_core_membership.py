"""Keep the clade members whose predicted topology holds exactly one GPCR-like 7TM core.

Each DeepTMHMM prediction is written as a string of topology states (O and I for the
extracellular and cytoplasmic sides, W and M for helices crossing from outside to inside and
from inside to outside) and consecutive identical states are collapsed. A sequence is kept when
the pattern OWIMOWIMOWIMOWI (seven helices, N terminus outside, C terminus inside) occurs
exactly once, counting overlapping matches. A signal peptide predicted as an extra helix
outside the core is therefore kept, while truncated sequences and those whose core cannot be
located unambiguously are not.

Reads  02_clade_extraction/clade_membership.tsv
       03_topology_filter/deeptmhmm_predictions_oriented.3line (build_deeptmhmm_tables.py)
Writes 03_topology_filter/core_membership.tsv             one column per family
       03_topology_filter/core_sequences/<family>.fasta    full-length sequences, the input of
                                                          the per-family alignments

Usage:
    python write_core_membership.py
"""
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import FAMILIES, open_text, stage  # noqa: E402

GPCR_PATTERN = "OWIMOWIMOWIMOWI"


def accession(header):
    parts = header.split("|")
    return parts[1] if len(parts) >= 3 else header


def core_count(oriented_labels):
    collapsed = "".join(m.group()[0] for m in re.finditer(r"(.)\1*", oriented_labels))
    count, index = 0, collapsed.find(GPCR_PATTERN)
    while index != -1:
        count += 1
        index = collapsed.find(GPCR_PATTERN, index + 1)
    return count


def read_wide(path):
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        members = {f: set() for f in header}
        for row in reader:
            for family, value in zip(header, row):
                if value.strip():
                    members[family].add(value.strip())
    return members


def main():
    clades = read_wide(stage("02_clade_extraction", "clade_membership.tsv"))
    single = {}
    with open_text(stage("03_topology_filter", "deeptmhmm_predictions_oriented.3line")) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().strip()
            labels = handle.readline().strip()
            name = header[1:].split()[0]
            if core_count(labels) == 1:
                single[accession(name)] = (name, sequence)

    kept = {f: clades[f] & single.keys() for f in FAMILIES}
    out_dir = stage("03_topology_filter", "core_sequences")
    out_dir.mkdir(parents=True, exist_ok=True)
    for family in FAMILIES:
        with open(out_dir / f"{family}.fasta", "w", encoding="utf-8") as handle:
            for acc in sorted(kept[family]):
                name, sequence = single[acc]
                handle.write(f">{name}\n{sequence}\n")
        print("%-14s clade %8d  single core %8d" % (family, len(clades[family]), len(kept[family])))

    ordered = {f: sorted(kept[f]) for f in FAMILIES}
    height = max(len(v) for v in ordered.values())
    with open(stage("03_topology_filter", "core_membership.tsv"), "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(FAMILIES)
        for i in range(height):
            writer.writerow([ordered[f][i] if i < len(ordered[f]) else "" for f in FAMILIES])
    print("unique sequences: %d" % len(set().union(*kept.values())))


if __name__ == "__main__":
    main()
