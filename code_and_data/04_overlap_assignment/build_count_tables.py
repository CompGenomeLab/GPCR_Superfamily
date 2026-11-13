"""Map every final family member to its genome and count members per genome and family.

The taxID of an accession is the OX= field of its UniProt header in the search outputs
(hmmsearch hit lines and the target header column of the MMseqs2 m8 files). The controls
(G alpha, G beta, G gamma, TBP and RPB1) are counted from their clade membership, since the
topology filter does not apply to them.

Reads  04_overlap_assignment/final_membership.tsv           (write_final_members.py)
       02_clade_extraction/clade_membership_controls.tsv
       01_homology_search/<search>/*_hmmsearch_E1000.txt[.gz], *.m8
Writes 04_overlap_assignment/acc_taxid_map.tsv              family, accession, taxID
       04_overlap_assignment/tax_counts_per_family.tsv      one row per genome (taxID),
                                                           <family>_count columns, controls last

Usage:
    python build_count_tables.py
"""
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import CONTROLS, FAMILIES, open_text, stage  # noqa: E402

OX = re.compile(r"\bOX=(\d+)\b")


def accession(token):
    parts = token.split("|")
    return parts[1] if len(parts) >= 3 else token


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


def taxid_map(wanted):
    """accession -> taxID for every wanted accession found in the search outputs."""
    found = {}
    searches = stage("01_homology_search")
    paths = sorted(searches.glob("*/*_hmmsearch_E1000.txt*")) + sorted(searches.glob("*/*.m8*"))
    for path in paths:
        m8 = ".m8" in path.name
        with open_text(path) as handle:
            for line in handle:
                if m8:
                    fields = line.split("\t", 3)
                    if len(fields) < 3:
                        continue
                    acc, text = accession(fields[1]), fields[2]
                elif line.startswith(">> "):
                    token, _, text = line[3:].strip().partition(" ")
                    acc = accession(token)
                else:
                    continue
                if acc in wanted and acc not in found:
                    match = OX.search(text)
                    if match:
                        found[acc] = match.group(1)
        if len(found) == len(wanted):
            break
    return found


def main():
    final = read_wide(stage("04_overlap_assignment", "final_membership.tsv"))
    controls = read_wide(stage("02_clade_extraction", "clade_membership_controls.tsv"))
    members = {**{f: final[f] for f in FAMILIES}, **{c: controls[c] for c in CONTROLS}}
    wanted = set().union(*members.values())
    taxid = taxid_map(wanted)
    missing = wanted - taxid.keys()
    if missing:
        print("WARNING: no taxID for %d accessions, e.g. %s" % (len(missing), sorted(missing)[:5]))

    counts = defaultdict(lambda: defaultdict(int))
    with open(stage("04_overlap_assignment", "acc_taxid_map.tsv"), "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["family", "accession", "taxID"])
        for family in FAMILIES + CONTROLS:
            for acc in sorted(members[family]):
                if acc in taxid:
                    counts[taxid[acc]][family] += 1
                    if family in FAMILIES:
                        writer.writerow([family, acc, taxid[acc]])

    columns = FAMILIES + CONTROLS
    with open(stage("04_overlap_assignment", "tax_counts_per_family.tsv"), "w", newline="",
              encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["taxID"] + [f + "_count" for f in columns])
        for tid in sorted(counts, key=int):
            writer.writerow([tid] + [counts[tid][f] for f in columns])
    print("%d genomes, %d accessions mapped" % (len(counts), len(wanted) - len(missing)))


if __name__ == "__main__":
    main()
