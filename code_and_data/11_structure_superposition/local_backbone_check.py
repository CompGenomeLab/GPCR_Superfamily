"""Motif residue overlap next to the overlap of its local backbone, from a saved superposition.

For each shared motif position the carriers are the structures of the families the motif text names
(TEXT_FAMILIES; None = the families the data call carriers, within the named lineages). For every
structure (carriers, and STE2 at the same site):
  residue_dev   C-alpha distance of its column residue to the carriers' C-alpha centroid
                (the structure itself left out of the centroid)
  local_bb      sequence-free local backbone offset: for each C-alpha within +-WINDOW residues of
                the column residue, the distance to the nearest C-alpha of each other carrier,
                medianed over carriers and averaged over the window
  excess        residue_dev - local_bb: large when the backbone lines up but the residue does not
                (a register or facing difference), small when the residue just follows its backbone

Writes local_backbone_check_<tag>.tsv (per structure) and prints a per-position summary.
Usage: python local_backbone_check.py [tag]   (default gp157)
"""
import csv
import sys

import numpy as np
from pymol import cmd

from motif_columns import HERE

WINDOW = 3
RHO = ["classA", "Olfactory", "classT", "Vomeronasal1", "Nematode"]
CAMP = ["cAMP", "classB1", "classB2", "Mth", "classF", "FSLB", "GP143", "GP157", "GPR1", "STE3", "TM116"]
GLU = ["classC", "Vomeronasal2"]
ALL = RHO + CAMP + GLU


def but(fams, out):
    return [f for f in fams if f not in out]


# shared positions of the motif text and the families it names for each; ("lineages", ...) = data
# carriers within those families; tags cAMP_plant / cAMP_fungal / C_plant pick single structures
TEXT_FAMILIES = {
    "3x25": but(ALL, ["classT", "Nematode"]),
    # the build takes the cysteine bonded to C3x25 here (Vomeronasal 1 C205 lies off the column)
    "45x50": but(ALL, ["classT", "Nematode"]),
    "6x48": RHO + GLU + ["classB1", "classB2", "Mth", "cAMP", "GP143", "TM116"],
    "7x53": ("lineages", ALL),
    "5x65": ("lineages", ALL),
    "2x46": ("lineages", ALL),
    "1x57": ("lineages", CAMP + GLU),
    "2x42": ("lineages", CAMP + GLU),
    "12x51": ("lineages", CAMP + GLU),
    "4x53": ("lineages", CAMP + ["C_plant"]),
    "4x50": but(RHO + CAMP, ["classT", "Nematode", "TM116", "STE3"]),
    "2x50": ["classA", "Olfactory", "GP143", "GP157", "GPR1", "TM116", "cAMP_plant"],
    "3x39": ["classA", "Olfactory", "cAMP_fungal", "cAMP_plant"],
    "3x49": ["classA", "Olfactory", "GP143", "Mth"],
    "3x51": but(RHO + CAMP, ["classT", "Vomeronasal1", "STE3"]),
    "3x54": ("lineages", RHO + CAMP),
    "6x44": ["classA", "classT", "Vomeronasal1", "classF", "GP143", "GP157", "GPR1", "TM116"],
    "5x58": ["classA", "Olfactory", "GP143", "GP157", "GPR1", "classB1", "Mth", "cAMP"],
    "7x49": ["classA", "Olfactory", "cAMP", "GP157", "GPR1", "TM116"],
    "P5x50": ("lineages", RHO),
    "P5x461": ("lineages", CAMP),
}


def tags(r):
    s, f = r["structure"], r["family"]
    if s.startswith("cAMP_plant"):
        return {"cAMP", "cAMP_plant"}
    if s.startswith("cAMP_"):
        return {"cAMP", "cAMP_fungal"}
    if s.startswith("C_plant"):
        return {"classC", "C_plant"}
    return {f}


def carriers(rows, fams):
    """Structures of the named families whose residue matches the motif (non-matching residues,
    e.g. a Taste-2 cysteine at an aromatic position, never enter the centroid or the medians)."""
    fams = fams[1] if isinstance(fams, tuple) else fams
    return [r for r in rows if r["carries_motif"] == "True" and tags(r) & set(fams)]


def ca_atoms(obj):
    atoms = [a for a in cmd.get_model(f"{obj} and name CA and not alt B").atom
             if a.resi.lstrip("-").isdigit()]
    return {int(a.resi): np.array(a.coord) for a in atoms}


def main():
    tag = sys.argv[1] if len(sys.argv) > 1 else "gp157"
    cmd.load(str(HERE / f"motif_superposition_tmcore_{tag}.pse"))
    rows = list(csv.DictReader(open(HERE / f"motif_residues_per_structure_{tag}.tsv"), delimiter="\t"))
    lineage = {}
    for r in csv.DictReader(open(HERE / f"motif_residue_deviation_{tag}.tsv"), delimiter="\t"):
        lineage[r["structure"]] = r["lineage"]
    cas = {}
    out = []
    for pos, fams in TEXT_FAMILIES.items():
        here = [r for r in rows if r["position"] == pos]
        car = carriers([r for r in here if r["family"] != "STE2"], fams)
        others = [r for r in here if r["family"] == "STE2"]
        for r in car + others:
            cas.setdefault(r["structure"], ca_atoms(r["structure"]))
        site = {r["structure"]: cas[r["structure"]].get(int(r["residue"][1:])) for r in car}
        site = {k: v for k, v in site.items() if v is not None}
        for r in car + others:
            s, n = r["structure"], int(r["residue"][1:])
            me = cas[s].get(n)
            ref = [v for k, v in site.items() if k != s]
            if me is None or len(ref) < 2:
                continue
            dev = float(np.linalg.norm(me - np.mean(ref, 0)))
            arrays = {k: np.array(list(cas[k].values())) for k in site if k != s}
            per_ca = []
            for i in range(n - WINDOW, n + WINDOW + 1):
                if i in cas[s]:
                    d = [np.min(np.linalg.norm(a - cas[s][i], axis=1)) for a in arrays.values()]
                    per_ca.append(float(np.median(d)))
            bb = float(np.mean(per_ca))
            role = "carrier" if r in car else "STE2"
            out.append((pos, s, r["family"], lineage.get(s, "background"), role, r["residue"],
                        dev, bb, dev - bb))
    with open(HERE / f"local_backbone_check_{tag}.tsv", "w") as f:
        f.write("position\tstructure\tfamily\tlineage\trole\tresidue\tresidue_dev\tlocal_bb\texcess\n")
        for o in out:
            f.write("\t".join(map(str, o[:6])) + f"\t{o[6]:.2f}\t{o[7]:.2f}\t{o[8]:.2f}\n")
    print("position\tn\tresidue_dev median (max)\tlocal_bb median (max)\texcess median (max: structure)\tSTE2 dev/bb/excess")
    for pos in TEXT_FAMILIES:
        c = [o for o in out if o[0] == pos and o[4] == "carrier"]
        s2 = [o for o in out if o[0] == pos and o[4] == "STE2"]
        if not c:
            continue
        dev, bb, ex = (np.array([o[k] for o in c]) for k in (6, 7, 8))
        top = max(c, key=lambda o: o[8])
        print(f"{pos}\t{len(c)}\t{np.median(dev):.1f} ({dev.max():.1f})\t{np.median(bb):.1f} ({bb.max():.1f})\t"
              f"{np.median(ex):.1f} ({top[8]:.1f}: {top[1]} {top[5]})\t"
              + "; ".join(f"{o[5]} {o[6]:.1f}/{o[7]:.1f}/{o[8]:.1f}" for o in s2))


if __name__ == "__main__":
    main()
