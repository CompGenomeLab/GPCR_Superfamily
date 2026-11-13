"""Supplementary Table S4: representative sequences used for the superfamily alignments and trees.

590 family representatives (step8 representatives_all_families_Oct08_topo0.8.tsv: family, taxonomy, topology
score) plus the 18 literature sequences added to the literature trees (headers of
step8 fasta/representatives_and_literature_Oct08; study from step6 deeper_origin_7TM_gate.tsv).
Tree-set membership from step8 representatives_<set>_Oct08_topo0.8.txt.

    python build_S4_representatives.py
"""
import csv
from pathlib import Path

import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import stage  # noqa: E402
CODE = Path(__file__).resolve().parent
HERE = stage("supplementary_tables")
HERE.mkdir(parents=True, exist_ok=True)
S8 = stage("07_representative_selection")
NICE = {"classA": "Class A", "classB1": "Class B1", "classB2": "Class B2", "classC": "Class C", "classF": "Class F",
        "classT": "Taste-2", "FSLB": "FSL", "Vomeronasal1": "Vomeronasal 1", "Vomeronasal2": "Vomeronasal 2",
        "GPR1like": "GPR1-like (claimed Rhodopsin)"}
SETS = {"all": "All families", "no_ste2": "Without STE2", "no_ste2ste3classCV2": "Without STE2, STE3, Class C, V2"}


def acc(line):
    h = line.strip().split()[0]
    return h.split("|")[1] if "|" in h else h


def main():
    members = {k: {acc(l) for l in open(S8 / f"representatives_{k}.txt", encoding="utf-8") if l.strip()}
               for k in SETS}
    rows = []
    for r in csv.DictReader(open(S8 / "representatives_all_families.tsv", encoding="utf-8"), delimiter="\t"):
        a = r["accession"]
        na = lambda k: "" if r[k] == "NA" else r[k]
        rows.append({
            "Accession or ID": a, "Source": "UniProt (our family set)", "Family": NICE.get(r["gpcrfamily"], r["gpcrfamily"]),
            "Species": r["species"], "NCBI taxID": r["taxid"], "Kingdom": na("kingdom"), "Phylum": na("phylum"),
            "Class": na("class"), "Order": na("order"),
            "Topology score": r["score"],
            **{f"In tree set: {v}": "yes" if a in members[k] else "" for k, v in SETS.items()},
            "In tree set: With literature sequences": "yes", "Literature study": "",
        })
    study = {r["id"]: r["study"].replace("_", " ") for r in csv.DictReader(
        open(stage("05_literature_comparison", "deeper_origin_7TM_gate.tsv"), encoding="utf-8"),
        delimiter="\t")}
    lit = stage("08_superfamily_alignment", "fasta", "representatives_and_literature",
                "representatives_and_literature.fasta")
    for line in open(lit, encoding="utf-8"):
        if line.startswith(">LIT_"):
            f = line[1:].strip().split("|")
            tags = dict(x.split("=", 1) for x in f[3:] if "=" in x)
            fam = f[0][4:]
            rows.append({
                "Accession or ID": f[1], "Source": "literature sequence", "Family": NICE.get(fam, fam),
                "Species": f[2].replace("_", " "), "NCBI taxID": "", "Kingdom": tags.get("layer", ""),
                "Phylum": "", "Class": "", "Order": "", "Topology score": "",
                **{f"In tree set: {v}": "" for v in SETS.values()},
                "In tree set: With literature sequences": "yes",
                "Literature study": study.get(f[1], tags.get("study", "").replace("_", " ")),
            })
    out = HERE / "Table_S5_representatives.tsv"
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    n_lit = sum(r["Source"] == "literature sequence" for r in rows)
    print(f"{len(rows) - n_lit} representatives + {n_lit} literature sequences; set sizes:",
          {v: sum(r[f'In tree set: {v}'] == 'yes' for r in rows) for v in SETS.values()},
          "; literature without study:", [r["Accession or ID"] for r in rows if r["Source"] != "UniProt (our family set)" and not r["Literature study"]])
    print("wrote", out)


if __name__ == "__main__":
    main()
