"""Column mapping for the motif superposition, and the residue of every family member at each
motif column.

The motif positions are columns of the ML-tree-guided superfamily alignment
(representatives_all_Oct08_ginsi_ep0.123_treein). Each superfamily column is tied to a family
alignment column by the consensus column mapping (consensus_column_mapping.py; the table in
COLUMN_MAP), as in the website logos (superfamily_logo_mappings_treein). Family rows are read from
the consensus-trimmed family alignments (consensus_trimmed_msas_Oct11), whose
<family>_kept_columns.txt gives the full family alignment column of each kept column.

main() writes member_motif_residues_<family>.tsv (letter and number from the /start-end header
range at every motif position) and checks the residue counts against the website JSONs.

Usage: python motif_columns.py [family ...]
"""
import csv
import json
import re
import sys
from collections import Counter
from pathlib import Path

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import open_text, stage  # noqa: E402

# Every script of this stage reads and writes its data in 11_structure_superposition/ of the
# archive (HERE); the motif columns come from stage 10.
HERE = stage("11_structure_superposition")
MOTIFS = stage("10_motif_analysis")
SRC = MOTIFS / "website_inputs" / "superfamily_logo_source_files"
MAPS = MOTIFS / "consensus_mapping" / "logo_conf0.33_top5" / "superfamily_logo_mappings_treein"
SUP_FASTA = SRC / "representatives_all_Oct08_ginsi_ep0.123_treein_plus_SLT_TRNS_topology.fasta"
# superfamily column -> family column: the top-5 consensus, as in the logos
COLUMN_MAP = (MOTIFS / "consensus_mapping"
              / "representatives_all_Oct08_ginsi_ep0.123_treein_conf0.33_top5.tsv")
TRIMMED = MOTIFS / "family_alignments" / "consensus_trimmed_msas"
FAMILIES = ["classA", "Olfactory", "classT", "Vomeronasal1", "Nematode", "GPR1", "cAMP", "classF",
            "FSLB", "GP143", "TM116", "GP157", "Mth", "classB2", "classB1", "STE3", "classC",
            "Vomeronasal2", "STE2"]


def read_fasta(path):
    name, seq = None, []
    with open_text(path) as f:
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
    p = header.split("|")
    return p[1] if len(p) > 2 else header.split()[0]


def start_of(header):
    m = re.search(r"/\s*(\d+)\s*-\s*(\d+)", header)
    return int(m.group(1)) if m else 1


def anchors():
    out = {}
    for line in open(SRC / "trim_info.tsv"):
        c = line.rstrip("\n").split("\t")
        if c[0] and c[0] != "acc1":
            out[c[1]] = c[0]
    return out


def residue_numbers(aligned, start):
    """Column index -> real residue number, for columns where the sequence has a residue."""
    nums, n = {}, start - 1
    for i, aa in enumerate(aligned):
        if aa != "-":
            n += 1
            nums[i] = n
    return nums


def trimmed_alignment(family):
    """The family's consensus-trimmed alignment (the rows behind the website logo)."""
    return TRIMMED / f"{family}_genes_filtered_db_FAMSA.consensus_trimmed.fasta"


def full_family_columns(family, motifs):
    """Motif name -> full family alignment column (1-based) through COLUMN_MAP, or None."""
    sup_to_fam = {int(r["sup_column"]): int(r["consensus_family_column"])
                  for r in csv.DictReader(open(COLUMN_MAP), delimiter="\t")
                  if r["family"] == family and r["consensus_family_column"]}
    return {m: sup_to_fam.get(c) for m, c in motifs.items()}


def family_motif_columns(family, acc1, sup, motifs):
    """Motif name -> column (0-based) of the consensus-trimmed family alignment, or None.
    acc1 and sup are unused, kept from the single-anchor mapping this replaced."""
    kept = {int(c): i for i, c in
            enumerate(open(TRIMMED / f"{family}_kept_columns.txt").read().split())}
    return {m: kept.get(c) for m, c in full_family_columns(family, motifs).items()}


def main():
    from final_motif_list import POSITIONS  # imported here: final_motif_list imports this module
    motifs = {label: col for label, col, *_ in POSITIONS}
    for family in sys.argv[1:] or FAMILIES:
        cols = family_motif_columns(family, None, None, motifs)
        site = json.load(open(MAPS / f"{family}.json"))["positions"]
        counts = {m: Counter() for m in motifs}
        out = open(HERE / f"member_motif_residues_{family}.tsv", "w")
        out.write("accession\theader\t" + "\t".join(motifs) + "\n")
        for h, s in read_fasta(trimmed_alignment(family)):
            nums = residue_numbers(s, start_of(h))
            cells = []
            for m, c in cols.items():
                if c is None or s[c] == "-":
                    cells.append("-")
                    continue
                aa = s[c].upper()
                if aa in "ACDEFGHIKLMNPQRSTVWY":  # the website counts standard residues only
                    counts[m][aa] += 1
                cells.append(f"{aa}{nums[c]}")
            out.write(f"{accession(h)}\t{h.split(' ', 1)[-1]}\t" + "\t".join(cells) + "\n")
        out.close()
        bad = [m for m, c in motifs.items()
               if dict(counts[m]) != {k: v for k, v in (site[c - 1].get("residueCounts") or {}).items()
                                      if k in "ACDEFGHIKLMNPQRSTVWY"}]
        print(f"{family:14s} website counts {'match' if not bad else 'DIFFER at ' + ','.join(bad)}")


if __name__ == "__main__":
    main()
