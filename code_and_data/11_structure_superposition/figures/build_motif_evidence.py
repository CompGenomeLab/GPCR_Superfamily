"""Shared-motif evidence per family for Fig. 6 (motif columns) and Fig. S3 : writes motif_evidence_data.js.

One set of criteria, applied the same way to every family and motif (Oct 16):
  sequence  = conserved (>= 60%, Table S8 '*') with a residue of the motif's type; at the two I/V/L columns (3x54,
              5x65) instead I + V + L together in >= 60% of the sequences (gaps counted, as in Table S8), because the
              S8 similarity groups split them (V/I and L/M); for cAMP and Class C, the whole family or any of its
              kingdoms (animals, fungi, plants) counts
  structure = a carrying residue of the motif's type in one of the family's structures (any kingdom) overlaps with
              the other carriers or lies within the local backbone variation (Table S10b); at C45x50 (ECL2, a loop)
              any carrier counts. STE2 residues were measured as the control, not as carriers, so STE2 has none.
  STE2 is the one exception: its sequence matches are not counted either, because in both STE2 structures its
              residues at these positions lie 4-17 A from the shared sites (alignment artifacts, as in the Results).
  evidence  = both | sequence (only) | structure (only) | none
Fig. S3 shows every shared motif grouped by the lineages that share it; Fig. 6 shows the same cells in sequence order,
without R/K12x51 (no common structural site).

    python build_motif_evidence.py
"""
import csv
import json
from pathlib import Path

import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from gpcr_paths import stage  # noqa: E402

SUPP = stage("supplementary_tables")

# (group, label (_X = GPCRdb class letter, drawn as a capital subscript), Class A label(s) in S8/S10b,
#  residues of the motif's type)
MOTIFS = [
    ("All three lineages", "L2x46_A", ["2x46"], "LIVM"),
    ("All three lineages", "C3x25_A", ["3x25"], "C"),
    ("All three lineages", "C45x50_A", ["45x50"], "C"),
    ("All three lineages", "Y7x53_A", ["7x53"], "YFWH"),
    ("All three lineages", "I/V/L5x65_A", ["5x65"], "IVL"),
    ("All three lineages", "W6x48_A", ["6x48"], "WFYH"),
    ("cAMP- and Rhodopsin-related", "D2x50_A", ["2x50"], "DE"),
    ("cAMP- and Rhodopsin-related", "S/D/E3x39_A", ["3x39"], "SDE"),
    ("cAMP- and Rhodopsin-related", "W4x50_A", ["4x50"], "W"),
    ("cAMP- and Rhodopsin-related", "D3x49_A", ["3x49"], "DE"),
    ("cAMP- and Rhodopsin-related", "Y3x51_A", ["3x51"], "YFWH"),
    ("cAMP- and Rhodopsin-related", "I/V/L3x54_A", ["3x54"], "IVL"),
    ("cAMP- and Rhodopsin-related", "Y5x58_A", ["5x58"], "YFWH"),
    ("cAMP- and Rhodopsin-related", "F6x44_A", ["6x44"], "FYWH"),
    ("cAMP- and Rhodopsin-related", "N7x49_A", ["7x49"], "ND"),
    ("cAMP- and Rhodopsin-related", "P5x42_B / P5x50_A", ["5x50", "5x461"], "P"),
    ("Glutamate- and cAMP-related", "F1x61_C", ["1x57"], "FYWH"),
    ("Glutamate- and cAMP-related", "L2x38_C", ["2x42"], "LIVM"),
    ("Glutamate- and cAMP-related", "P4x53_B / P4x43_C", ["4x53"], "P"),
    ("Glutamate- and cAMP-related", "R/K12x51_C", ["12x51"], "RK"),
]
# Fig. 6 columns in sequence order (Class A numbers shown without the class letter)
FIG6_ORDER = ["F1x61_C", "L2x38_C", "L2x46_A", "D2x50_A", "C3x25_A", "S/D/E3x39_A", "D3x49_A", "Y3x51_A",
              "I/V/L3x54_A", "W4x50_A", "P4x53_B / P4x43_C", "C45x50_A", "P5x42_B / P5x50_A", "Y5x58_A",
              "I/V/L5x65_A", "F6x44_A", "W6x48_A", "N7x49_A", "Y7x53_A"]
# rows: (lineage, family, placed in the trees of the 15 retained families)
ROWS = [
    ("Rho", "Class A", True), ("Rho", "Olfactory", True), ("Rho", "Vomeronasal 1", True),
    ("Rho", "Taste-2", True), ("Rho", "Nematode", True),
    ("cAMP", "cAMP", True), ("cAMP", "Mth", True), ("cAMP", "Class B2", True), ("cAMP", "Class B1", True),
    ("cAMP", "FSL", True), ("cAMP", "Class F", True), ("cAMP", "GP143", True), ("cAMP", "TM116", True),
    ("cAMP", "GP157", True), ("cAMP", "GPR1", True), ("cAMP", "STE3", False),
    ("Glu", "Class C", False), ("Glu", "Vomeronasal 2", False),
    ("none", "STE2", False),
]
KINGDOMS = ["animals", "fungi", "plants"]
# residue counts per alignment column, for the I/V/L totals: the logo mappings Table S8 is built from
MAPPINGS = stage("10_motif_analysis", "consensus_mapping", "logo_conf0.33_top5", "superfamily_logo_mappings_treein")
MAPPING_NAME = {"Class A": "classA", "Class B1": "classB1", "Class B2": "classB2", "Class C": "classC",
                "Class F": "classF", "Taste-2": "classT", "FSL": "FSLB", "Vomeronasal 1": "Vomeronasal1",
                "Vomeronasal 2": "Vomeronasal2"}


def column_share(col, column, residues, _cache={}):
    """Share (%) of a family's (or kingdom's) sequences carrying any of `residues` at 1-based alignment `column`."""
    base, _, kingdom = col.partition(" (")
    name = MAPPING_NAME.get(base, base) + (f"_{kingdom.rstrip(')').capitalize()}" if kingdom else "")
    if name not in _cache:
        _cache[name] = json.loads((MAPPINGS / f"{name}.json").read_text())["positions"]
    p = _cache[name][column - 1]
    if not p or not p.get("totalSequences"):
        return 0.0
    counts = p.get("residueCounts") or {}
    return 100 * sum(n for a, n in counts.items() if a in residues) / p["totalSequences"]


def evidence(s7, s9):
    by_label = {r["GPCRdb (Class A)"]: r for r in s7}
    out = {}
    for _, fam, _ in ROWS:
        cols = [fam] + [f"{fam} ({k})" for k in KINGDOMS if f"{fam} ({k})" in s7[0]]
        cells = []
        for _, label, pos, residues in MOTIFS:
            if label.startswith("I/V/L"):
                shares = [(c, column_share(c, int(by_label[p]["Alignment column"]), residues)) for p in pos for c in cols]
                seq = [(c, f"I+V+L {v:.0f}%") for c, v in shares if v >= 60]
            else:
                seq = [(c, by_label[p][c]) for p in pos for c in cols
                       if by_label[p][c].endswith("*") and by_label[p][c].split()[0] in residues]
            carriers = [r for r in s9 if r["Family"] == fam and r["Position (Class A GPCRdb)"] in pos
                        and r["Role"] == "motif carrier" and r["Residue"][0] in residues]
            fits = [r for r in carriers if label.startswith("C45x50")
                    or r["Call"] in ("overlaps", "within local backbone variation")]
            ev = "both" if seq and fits else "sequence" if seq else "structure" if fits else None
            if fam == "STE2":
                ev = None  # sequence matches displaced 4-17 A in structure: alignment artifacts
            parts = [f"{c}: {v}" for c, v in seq] or [f"{fam}: {by_label[pos[0]][fam]}"]
            parts += [f"{r['Structure']} {r['Residue']} {r['Distance from centroid (A)']} Å, {r['Call']}"
                      for r in carriers]
            cells.append({"ev": ev, "note": "; ".join(parts)})
        out[fam] = cells
    return out


def main():
    s7 = list(csv.DictReader(open(SUPP / "Table_S8_motif_conservation.tsv", encoding="utf-8"), delimiter="\t"))
    s9 = list(csv.DictReader(open(SUPP / "Table_S10b_motif_residue_distances.tsv", encoding="utf-8"), delimiter="\t"))
    ev = evidence(s7, s9)
    labels = [m[1] for m in MOTIFS]
    fig7 = {"motifs": [{"group": g, "label": l} for g, l, _, _ in MOTIFS],
            "rows": [{"lineage": lin, "label": fam, "placed": placed, "cells": ev[fam]} for lin, fam, placed in ROWS]}
    idx = [labels.index(l) for l in FIG6_ORDER]
    fig6 = {"motifs": [l.replace("_A", "") for l in FIG6_ORDER],
            "rows": {fam: [ev[fam][i] for i in idx] for _, fam, _ in ROWS}}
    (HERE / "motif_evidence_data.js").write_text("/* generated by build_motif_evidence.py from Tables S8 and S10b */\nconst FIG7 = "
                                       + json.dumps(fig7, ensure_ascii=False, indent=1) + ";\nconst FIG6 = "
                                       + json.dumps(fig6, ensure_ascii=False, indent=1) + ";\n", encoding="utf-8")
    sym = {"both": "B", "sequence": "s", "structure": "t", None: "."}  # B both, s sequence only, t structure only
    print(" " * 16 + " ".join(l.split("_")[0].split(" ")[0][:5] for l in labels))
    for _, fam, _ in ROWS:
        print(f"{fam:16s}", "     ".join(sym[c["ev"]] for c in ev[fam]))
    print("wrote motif_evidence_data.js")


if __name__ == "__main__":
    main()
