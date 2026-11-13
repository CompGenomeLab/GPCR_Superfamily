"""Supplementary Table S9: structures superimposed and the spatial overlap of every shared motif residue.

From step12_structure_superposition (GP157-anchored build refined against the mean TM core, Oct 12-13):
  S9a  one row per structure: family, source (PDB / AlphaFold), PDB ID or UniProt accession, TM residues
       used, core RMSD after refinement (superposition_core_fit_gp157.tsv)
  S9b  one row per motif position x structure: residue, distance of its CA from the centroid of the other
       carrying residues, local backbone variation (+-3 residues), and the call; STE2 rows give the STE2
       residue at the site (local_backbone_check_gp157.tsv). Positions are named by Class A GPCRdb number,
       as in step 12; the text's own numbering is added where it differs.

    python build_S9_structures.py
"""
import csv
from pathlib import Path

import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import stage  # noqa: E402
CODE = Path(__file__).resolve().parent
HERE = stage("supplementary_tables")
HERE.mkdir(parents=True, exist_ok=True)
S12 = stage("11_structure_superposition")
NICE = {"classA": "Class A", "classB1": "Class B1", "classB2": "Class B2", "classC": "Class C", "classF": "Class F",
        "classT": "Taste-2", "FSLB": "FSL", "Vomeronasal1": "Vomeronasal 1", "Vomeronasal2": "Vomeronasal 2"}
TEXT_NAME = {"1x57": "1x61C", "2x42": "2x38C", "4x53": "4x53B (plant Class C: 4x43C)", "P5x461": "5x42B",
             "P5x50": "5x50A", "12x51": "12x51C"}
CALL = {"overlaps": "overlaps", "backbone": "within local backbone variation",
        "residue off": "displaced beyond local backbone variation"}


def write(name, rows):
    with open(HERE / name, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print("wrote", name, len(rows), "rows")


def main():
    a = []
    for r in csv.DictReader(open(S12 / "superposition_core_fit_gp157.tsv", encoding="utf-8"), delimiter="\t"):
        ident = r["structure"].rsplit("_", 1)[-1]
        a.append({"Structure": r["structure"], "Family": NICE.get(r["family"], r["family"]),
                  "Source": "experimental (PDB)" if r["source"] == "PDB" else "AlphaFold model",
                  "PDB ID or UniProt accession": ident, "TM residues used": r["TM_residues"],
                  "Core residue pairs": r["core_pairs"], "Core RMSD (A)": r["core_RMSD"]})
    write("Table_S10a_structures.tsv", a)

    calls = {(r["position"], r["structure"]): r["call"] for r in csv.DictReader(
        open(S12 / "per_family_overlap_gp157.tsv", encoding="utf-8"), delimiter="\t")}
    b = []
    for r in csv.DictReader(open(S12 / "local_backbone_check_gp157.tsv", encoding="utf-8"), delimiter="\t"):
        p = r["position"]
        b.append({"Position (Class A GPCRdb)": p.lstrip("P") if p.startswith("P5") else p,
                  "Name in text (if different)": TEXT_NAME.get(p, ""),
                  "Structure": r["structure"], "Family": NICE.get(r["family"], r["family"]), "Lineage": r["lineage"],
                  "Role": "motif carrier" if r["role"] == "carrier" else "STE2 (control)",
                  "Residue": r["residue"], "Distance from centroid (A)": r["residue_dev"],
                  "Local backbone variation (A)": r["local_bb"],
                  "Call": CALL.get(calls.get((p, r["structure"]), ""), "")})
    write("Table_S10b_motif_residue_distances.tsv", b)


if __name__ == "__main__":
    main()
