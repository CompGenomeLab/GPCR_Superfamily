"""Superfamily motif positions selected with the manuscript's method (website scatter plot,
binary mode), and each family's state at them.

Per family and superfamily column (logo mapping JSONs, columns of the 590-sequence --ep 0.123
alignment): conservation = (count of the most frequent residue + counts of the residues in its
similarity group) / all sequences, gaps included; groups D/E, F/W/Y/H, K/R, N/Q, V/I, L/M
(FamilyScatterPlot.tsx). A family is conserved at a column when this reaches 60%. Columns are
ranked by the number of conserved families (STE2 and SLT_TRNS never count), then by the
Shannon entropy of the families' most frequent symbols (gap included), lowest first.
Columns with at least MIN_FAMILIES conserved families are kept; MIN_FAMILIES is the smallest
value that keeps at most TARGET columns.

A family's state at a column is its top residue's group (or the residue itself when it has no
group), only where the family is conserved. Writes results/motif_positions.tsv.
"""
import json
import math
from collections import Counter
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import stage  # noqa: E402
MAPPINGS = stage("10_motif_analysis", "consensus_mapping", "logo_conf0.33_top5", "superfamily_logo_mappings_treein")
THRESHOLD = 60.0
TARGET = 31  # 31 since the Oct09 GPR1 alignment brought 3x42 to 12 families (user, 2026-10-09)
GROUPS = {"D": "DE", "E": "DE", "W": "FWYH", "Y": "FWYH", "H": "FWYH", "F": "FWYH",
          "R": "KR", "K": "KR", "Q": "NQ", "N": "NQ", "V": "IV", "I": "IV", "M": "LM", "L": "LM"}
COUNTED = ["cAMP", "classA", "classB1", "classB2", "classC", "classF", "classT", "FSLB", "GP143",
           "GP157", "GPR1", "Mth", "Nematode", "Olfactory", "STE3", "TM116", "Vomeronasal1",
           "Vomeronasal2"]
SHOWN = COUNTED + ["STE2"]


def family_columns(family):
    data = json.loads((MAPPINGS / f"{family}.json").read_text())
    cols = []
    for p in data["positions"]:
        if not p or not p.get("totalSequences"):
            cols.append(None)
            continue
        counts, total = p.get("residueCounts") or {}, p["totalSequences"]
        if not counts:
            cols.append({"label": p.get("gpcrdb"), "top": None, "cons": 0.0, "mode": "-"})
            continue
        top, n = max(counts.items(), key=lambda kv: kv[1])
        group = GROUPS.get(top, top)
        similar = sum(c for a, c in counts.items() if a in group)
        gaps = total - sum(counts.values())
        cols.append({"label": p.get("gpcrdb"), "top": top, "group": group,
                     "cons": 100 * similar / total, "mode": "-" if gaps >= n else top})
    return cols


def load_all():
    return {f: family_columns(f) for f in SHOWN}


def select(fams):
    ncol = len(fams["classA"])
    rows = []
    for c in range(ncol):
        present = [fams[f][c] for f in COUNTED if fams[f][c]]
        if not present:
            continue
        n_cons = sum(p["cons"] >= THRESHOLD for p in present)
        modes = Counter(p["mode"] for p in present)
        k = sum(modes.values())
        entropy = -sum(v / k * math.log2(v / k) for v in modes.values())
        rows.append((c, n_cons, entropy))
    for min_fam in range(len(COUNTED), 0, -1):
        if sum(r[1] >= min_fam for r in rows) > TARGET:
            min_fam += 1
            break
    chosen = sorted((r for r in rows if r[1] >= min_fam), key=lambda r: (-r[1], r[2], r[0]))
    return chosen, min_fam


def state(fams, family, c):
    p = fams[family][c]
    return p["group"] if p and p["top"] and p["cons"] >= THRESHOLD else None


LABEL_ORDER = ("classA", "classB2", "classB1", "classF", "classC", "GP143", "GP157")


def label(fams, c):
    """GPCRdb generic number (NxNN) at the column, from the first family in LABEL_ORDER whose
    reference has one there, and that family; the column number if none has one."""
    for f in LABEL_ORDER:
        p = fams[f][c]
        if p and p["label"] and "x" in str(p["label"]):
            return str(p["label"]), f
    return f"col{c + 1}", ""


def main():
    fams = load_all()
    chosen, min_fam = select(fams)
    print(f"{len(chosen)} columns with >= {min_fam} of {len(COUNTED)} families conserved at {THRESHOLD}%")
    with open(stage("10_motif_analysis", "motif_positions.tsv"), "w") as fh:
        fh.write("column\tlabel\tlabel_from\tn_conserved\tentropy\t" + "\t".join(SHOWN) + "\n")
        for c, n, e in chosen:
            cells = []
            for f in SHOWN:
                p = fams[f][c]
                cells.append(f"{p['top']} {p['cons']:.0f}" if p and p["top"] else "-")
            lab, src = label(fams, c)
            fh.write(f"{c + 1}\t{lab}\t{src}\t{n}\t{e:.2f}\t" + "\t".join(cells) + "\n")


if __name__ == "__main__":
    main()
