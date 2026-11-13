"""Supplementary Table S7: conservation at the superfamily positions discussed in the paper.

Rule and data as Fig. 3 (step11 make_fig3.py): ML-tree-guided G-INS-i --ep 0.123 superfamily alignment, logo
mapping JSONs in GPCRevolution/public/superfamily_logo_mappings_treein; conservation = share of sequences
carrying the most frequent residue or its similarity group (gaps counted); conserved at >= 60%.
Rows: every column conserved in >= 12 of the 18 counted families (the Fig. 3 set) plus every other GPCRdb
position named in the motif and structure Results sections or the figure legends (read from
second_draft/manuscript_text; schemes A, B, C, F and T).
Each family cell: most frequent residue and conservation (%), '*' when >= 60%.

    python build_S7_motif_conservation.py
"""
import csv
import re
import sys
from pathlib import Path

import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import stage  # noqa: E402
CODE = Path(__file__).resolve().parent
HERE = stage("supplementary_tables")
HERE.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(CODE.parent / "10_motif_analysis"))
import superfamily_motif_positions as smp  # noqa: E402

smp.MAPPINGS = stage("10_motif_analysis", "consensus_mapping", "logo_conf0.33_top5", "superfamily_logo_mappings_treein")
MIN_FAMILIES = 12
SCHEME = {"A": "classA", "B": "classB1", "C": "classC", "F": "classF", "T": "classT"}
ORDER = ["classA", "Olfactory", "Vomeronasal1", "classT", "Nematode", "cAMP", "cAMP_Animals", "cAMP_Fungi",
         "cAMP_Plants", "classB1", "classB2", "Mth", "classF", "FSLB", "GP143", "TM116", "GP157", "GPR1", "STE3",
         "classC", "classC_Animals", "classC_Fungi", "classC_Plants", "Vomeronasal2", "STE2"]
NICE = {"classA": "Class A", "classB1": "Class B1", "classB2": "Class B2", "classC": "Class C", "classF": "Class F",
        "classT": "Taste-2", "FSLB": "FSL", "Vomeronasal1": "Vomeronasal 1", "Vomeronasal2": "Vomeronasal 2"}


def nice(f):
    base, _, kingdom = f.partition("_")
    return f"{NICE.get(base, base)} ({kingdom.lower()})" if kingdom else NICE.get(f, f)


# GPCRdb positions (number, scheme letter) named in the motif and structure Results sections and the
# figure legends of the manuscript (Results 4-7, legends), read from the text on 21 October 2026;
# 3x39A (named as "3x39" in the Figure 5 legend) is kept from the earlier version of the table.
MENTIONED = [("1x61", "C"), ("2x38", "C"), ("12x51", "C"), ("4x43", "C"), ("4x53", "B"), ("4x50", "A"),
             ("3x49", "A"), ("3x51", "A"), ("3x54", "A"), ("5x65", "A"), ("6x44", "A"), ("5x50", "A"),
             ("5x42", "B"), ("5x58", "A"), ("7x49", "A"), ("1x50", "B"), ("3x46", "B"), ("45x51", "B"),
             ("7x50", "B"), ("6x35", "C"), ("3x58", "C"), ("2x50", "A"), ("6x53", "B"), ("3x39", "A")]


def mentioned():
    return list(MENTIONED)


def main():
    fams = {f: smp.family_columns(f) for f in ORDER}
    chosen, min_fam = smp.select({k: fams[k] for k in smp.SHOWN})
    chosen = [r for r in chosen if r[1] >= MIN_FAMILIES]
    top = {c: n for c, n, _ in chosen}
    label_col = {}
    for scheme, fam in SCHEME.items():
        for c, p in enumerate(fams[fam]):
            if p and p.get("label"):
                label_col.setdefault((str(p["label"]), scheme), c)
    extra, unplaced = {}, []
    for num, scheme in mentioned():
        c = label_col.get((num, scheme))
        if c is None:
            unplaced.append(num + scheme)
        elif c not in top:
            extra.setdefault(c, []).append(num + scheme)

    def n_cons(c):
        return sum(1 for f in smp.COUNTED if fams[f][c] and fams[f][c]["cons"] >= smp.THRESHOLD)

    rows = []
    for c in sorted(set(top) | set(extra)):
        lab = lambda f: str((fams[f][c] or {}).get("label") or "")
        row = {"Alignment column": c + 1,
               "GPCRdb (Class A)": lab("classA"), "GPCRdb (Class B)": lab("classB1"),
               "GPCRdb (Class C)": lab("classC"), "GPCRdb (Class F)": lab("classF"), "GPCRdb (Taste-2)": lab("classT"),
               "Set": "top (>= 12 families, Fig. 3)" if c in top else "named in text: " + ", ".join(extra[c]),
               "Families conserved (of 18)": n_cons(c)}
        for f in ORDER:
            p = fams[f][c]
            row[nice(f)] = (f"{p['top']} {p['cons']:.0f}%" + ("*" if p["cons"] >= smp.THRESHOLD else "")
                            if p and p.get("top") else "-")
        rows.append(row)
    out = HERE / "Table_S8_motif_conservation.tsv"
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(top)} top positions (>= {MIN_FAMILIES} families), {len(extra)} extra columns named in text:",
          {c + 1: v for c, v in extra.items()})
    print("named in text but not found in the scheme family's labels:", unplaced)
    print("wrote", out)


if __name__ == "__main__":
    main()
