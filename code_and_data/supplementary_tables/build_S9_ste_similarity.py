"""Supplementary Table S8: conserved residues STE3 and STE2 share with every other family.

From step11 ste_consensus_similarity.py, run on both MSAs built with the maximum likelihood guide tree
(G-INS-i --ep 0.123, the main alignment, and plain G-INS-i). For each family pair and threshold t, a column counts
when both families' most frequent residue reaches a share >= t (column >= 50% occupied in both); the number
reported is how many of those columns carry the identical consensus residue in both families.

    python build_S8_ste_similarity.py
"""
import csv
from pathlib import Path

import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import stage  # noqa: E402
CODE = Path(__file__).resolve().parent
HERE = stage("supplementary_tables")
HERE.mkdir(parents=True, exist_ok=True)
S11 = stage("10_motif_analysis")
RUNS = {"main (--ep 0.123)": "ste_consensus_similarity.tsv", "plain G-INS-i": "ste_consensus_similarity_lower_gap.tsv"}
THRESHOLDS = {"0.4": ">= 40%", "0.6": ">= 60%"}
LINEAGE = {"classA": "Rhodopsin", "Olfactory": "Rhodopsin", "Vomeronasal1": "Rhodopsin", "classT": "Rhodopsin",
           "Nematode": "Rhodopsin", "classC": "Glutamate", "Vomeronasal2": "Glutamate", "STE2": "-", "STE3": "-"}
NICE = {"classA": "Class A", "classB1": "Class B1", "classB2": "Class B2", "classC": "Class C", "classF": "Class F",
        "classT": "Taste-2", "FSLB": "FSL", "Vomeronasal1": "Vomeronasal 1", "Vomeronasal2": "Vomeronasal 2"}
ORDER = ["cAMP", "classB1", "classB2", "Mth", "classF", "FSLB", "GP143", "TM116", "GP157", "GPR1", "classC",
         "Vomeronasal2", "classA", "Olfactory", "Vomeronasal1", "classT", "Nematode", "STE3", "STE2"]

data = {}
for run, f in RUNS.items():
    for r in csv.DictReader(open(S11 / f, encoding="utf-8"), delimiter="\t"):
        if r["comparison"] == "identity" and r["threshold"] in THRESHOLDS:
            data[run, r["query"], r["family"], r["threshold"]] = r["matches"]
rows = []
for q in ("STE3", "STE2"):
    for fam in ORDER:
        if fam == q:
            continue
        row = {"Query": q, "Family": NICE.get(fam, fam), "Lineage": LINEAGE.get(fam, "cAMP")}
        for run in RUNS:
            for t, tl in THRESHOLDS.items():
                row[f"Shared conserved residues, {run}, {tl}"] = data.get((run, q, fam, t), "")
        rows.append(row)
out = HERE / "Table_S9_STE2_STE3_shared_residues.tsv"
with open(out, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)
print("wrote", out, len(rows), "rows")
