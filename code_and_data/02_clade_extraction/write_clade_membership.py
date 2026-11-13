"""Write the clade membership of every family as one wide table, the input of stage 03.

Three steps turn the extracted clades into family membership.

1. Class A was aligned and its tree built from the representatives of an 80% identity
   clustering of its search hits, so every extracted representative is expanded back to all
   members of its cluster (01_homology_search/classA/classA_hmm_E1000_db_clust0.8.tsv).
2. Where the search of a broader family also retrieved the members of a family defined by a
   more specific search, the shared sequences are given to the specific family: GPR1 and
   Vomeronasal 1 are removed from Class A, and Vomeronasal 2 from Class C.
3. Five Olfactory accessions recorded in curated_accession_additions.tsv are added.

Reads  02_clade_extraction/extracted_clades/alignments/*_extracted_clade.fasta (identify_gene_clade.py)
Writes 02_clade_extraction/clade_membership.tsv           one column per GPCR family
       02_clade_extraction/clade_membership_controls.tsv  the G protein and housekeeping controls

Usage:
    python write_clade_membership.py
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import CONTROLS, FAMILIES, open_text, stage  # noqa: E402

HERE = Path(__file__).resolve().parent
CLADES = stage("02_clade_extraction", "extracted_clades", "alignments")
CLASSA_CLUSTERS = stage("01_homology_search", "classA", "classA_hmm_E1000_db_clust0.8.tsv")
NESTED_MODEL_EXCLUSIONS = {"classA": ["GPR1", "Vomeronasal1"], "classC": ["Vomeronasal2"]}
# Output name of identify_gene_clade.py -> family name used from here on.
RENAME = {"GPR1_YEAST": "GPR1", "FSLB_DICDI": "FSLB"}


def accession(header):
    parts = header.split("|")
    return parts[1] if len(parts) >= 3 else header


def read_clades():
    members = {}
    for path in sorted(CLADES.glob("*_extracted_clade.fasta*")):
        name = path.name.split("_hmm_E1000_db")[0].split("_db_FAMSA")[0]
        name = RENAME.get(name, name)
        with open_text(path) as handle:
            members[name] = {accession(line[1:].split()[0]) for line in handle if line.startswith(">")}
    return members


def expand_classA(representatives):
    clusters = defaultdict(set)
    with open_text(CLASSA_CLUSTERS) as handle:
        for line in handle:
            rep, member = line.split()
            clusters[accession(rep)].add(accession(member))
    expanded = set()
    for rep in representatives:
        expanded |= clusters.get(rep, set()) | {rep}
    return expanded


def write_wide(path, columns, members):
    ordered = {f: sorted(members[f]) for f in columns}
    height = max(len(v) for v in ordered.values())
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for i in range(height):
            writer.writerow([ordered[f][i] if i < len(ordered[f]) else "" for f in columns])


def main():
    members = read_clades()
    missing = [f for f in FAMILIES + CONTROLS if f not in members]
    if missing:
        sys.exit("no extracted clade for: " + ", ".join(missing))
    members["classA"] = expand_classA(members["classA"])
    for broad, specific in NESTED_MODEL_EXCLUSIONS.items():
        for family in specific:
            members[broad] -= members[family]
    with open(HERE / "curated_accession_additions.tsv", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            members[row["family"]].add(row["accession"])

    write_wide(stage("02_clade_extraction", "clade_membership.tsv"), FAMILIES, members)
    write_wide(stage("02_clade_extraction", "clade_membership_controls.tsv"), CONTROLS, members)
    for family in FAMILIES:
        print("%-14s %8d" % (family, len(members[family])))
    print("unique GPCR sequences: %d" % len(set().union(*(members[f] for f in FAMILIES))))


if __name__ == "__main__":
    main()
