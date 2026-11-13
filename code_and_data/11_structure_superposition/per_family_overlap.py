"""Per motif position, call each named family's structure as overlapping or not.

Reads local_backbone_check_<tag>.tsv (carriers, distances) and motif_residues_per_structure_<tag>.tsv
(residues of the named families that do not carry the motif). Calls:
  overlaps        residue within OVERLAP A of the other carriers' C-alpha centroid
  backbone        further than OVERLAP, but within EXCESS A of its local backbone offset (the region
                  is shifted, not the residue)
  residue off     more than EXCESS A beyond its local backbone offset
  no motif        family named for the position, but this structure's residue does not carry it
Writes per_family_overlap_<tag>.tsv and prints one line per position.
Usage: python per_family_overlap.py [tag]   (default gp157)
"""
import csv
import sys

from local_backbone_check import TEXT_FAMILIES, tags
from motif_columns import HERE

OVERLAP = 2.0
EXCESS = 1.0


def main():
    tag = sys.argv[1] if len(sys.argv) > 1 else "gp157"
    rows = [r for r in csv.DictReader(open(HERE / f"local_backbone_check_{tag}.tsv"), delimiter="\t")
            if r["role"] == "carrier"]
    per = list(csv.DictReader(open(HERE / f"motif_residues_per_structure_{tag}.tsv"), delimiter="\t"))
    out = []
    for pos, fams in TEXT_FAMILIES.items():
        named = set(fams[1] if isinstance(fams, tuple) else fams)
        for r in (r for r in rows if r["position"] == pos):
            dev, ex = float(r["residue_dev"]), float(r["excess"])
            call = "overlaps" if dev <= OVERLAP else "backbone" if ex <= EXCESS else "residue off"
            out.append((pos, r["structure"], r["family"], r["lineage"], r["residue"],
                        r["residue_dev"], r["local_bb"], r["excess"], call))
        if not isinstance(fams, tuple):  # data-defined carriers have no named non-carriers
            for r in per:
                if (r["position"] == pos and r["family"] != "STE2" and tags(r) & named
                        and r["carries_motif"] != "True"):
                    out.append((pos, r["structure"], r["family"], "", r["residue"], "", "", "", "no motif"))
    with open(HERE / f"per_family_overlap_{tag}.tsv", "w") as f:
        f.write("position\tstructure\tfamily\tlineage\tresidue\tresidue_dev\tlocal_bb\texcess\tcall\n")
        for o in out:
            f.write("\t".join(o) + "\n")
    for pos in TEXT_FAMILIES:
        o = [x for x in out if x[0] == pos]
        n_ok = sum(x[8] == "overlaps" for x in o)
        n_car = sum(x[8] != "no motif" for x in o)
        rest = "; ".join(f"{x[1]} {x[4]} {x[5] or '-'} ({x[8]})" for x in o if x[8] != "overlaps")
        print(f"{pos}\t{n_ok}/{n_car} overlap\t{rest}")


if __name__ == "__main__":
    main()
