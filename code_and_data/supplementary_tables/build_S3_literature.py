"""Supplementary Table S3: literature comparison (from step 6, which builds the tables).

  S3a  every literature receptor and its fate in our pipeline   (Supplementary_Table_literature_receptors.tsv)
       plus one explicit family per receptor:
         matched receptors  -> our final family of the matched protein, or 'none (<stage that removed it>)'
         loose matches      -> final family of the best MMseqs2 hit (step5 final membership Sep24; SLT_TRNS
                               from the topology table), an indicator only, never used for assignment
  S3b  agreement by study                                         (..._agreement_by_study.tsv)
  S3c  agreement by family                                        (..._agreement_by_family.tsv), plus the studies
       using each label (synonymous labels are joined with '; '), with their receptor counts from S3a
  S3d  literature receptors that could move a family origin deeper (deeper_origin_candidates.tsv)

    python build_S3_literature.py
"""
import csv
import re
from collections import Counter
import shutil
from pathlib import Path

import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import stage  # noqa: E402
CODE = Path(__file__).resolve().parent
HERE = stage("supplementary_tables")
HERE.mkdir(parents=True, exist_ok=True)
SRC = stage("05_literature_comparison")
COPY = {"Table_S3b_literature_agreement_by_study.tsv": "Supplementary_Table_literature_agreement_by_study.tsv",
        "Table_S3d_deeper_origin_candidates.tsv": "deeper_origin_candidates.tsv"}
NO_SEQUENCE = ("no significant hit in our proteome set", "no sequence; identifier not in our data")


def membership():
    fam = {}
    for path in (stage("03_topology_filter", "topology_membership.tsv"),
                 stage("04_overlap_assignment", "final_membership.tsv")):
        rows = list(csv.reader(open(path, encoding="utf-8"), delimiter="\t"))
        for r in rows[1:]:
            for c, v in zip(rows[0], r):
                if v.strip() and (c == "SLT_TRNS" or path.name.startswith("accessions_by_family_final")):
                    fam[v.strip()] = c
    return fam


def main():
    fam = membership()
    rows = list(csv.DictReader(open(SRC / "Supplementary_Table_literature_receptors.tsv", encoding="utf-8"),
                               delimiter="\t"))
    out = []
    for r in rows:
        o, final = r["outcome"], r["step8_final_family"]
        hit = ""
        if o.startswith("not matched"):
            m = re.search(r"\(([A-Z0-9]+),", r["nearest_family_indicator"])
            hit = m.group(1) if m else ""
            f = fam.get(hit, "")
            assigned = ("none (best hit is SLT_TRNS background)" if f == "SLT_TRNS"
                        else f or "none (best hit in no family of ours)")
            via = "best MMseqs2 hit (indicator only)"
        elif o in NO_SEQUENCE:
            assigned, via = f"none ({o})", "-"
        elif final and final not in ("-", "SLT_TRNS"):
            assigned, via = final, "matched protein"
        elif final == "SLT_TRNS":
            assigned, via = "none (in the SLT_TRNS background clade)", "matched protein"
        else:
            assigned, via = f"none ({o})", "matched protein"
        row = {}
        for k, v in r.items():
            row["families corresponding to study label" if k == "our_family" else k] = v
            if k == "step8_final_family":
                row["assigned family"], row["assigned via"], row["best hit (loose match)"] = assigned, via, hit
        out.append(row)
    with open(HERE / "Table_S3a_literature_receptors.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out)
    print(Counter(r["assigned via"] for r in out))
    print(Counter(r["assigned family"] for r in out).most_common(30))
    by_label = {}
    for r in rows:
        by_label.setdefault(r["study_family"], Counter())[r["study"]] += 1
    fam_rows = list(csv.DictReader(open(SRC / "Supplementary_Table_literature_agreement_by_family.tsv",
                                        encoding="utf-8"), delimiter="\t"))
    with open(HERE / "Table_S3c_literature_agreement_by_family.tsv", "w", newline="", encoding="utf-8") as fh:
        keys = list(fam_rows[0])
        keys.insert(1, "studies (receptors)")
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t")
        w.writeheader()
        for r in fam_rows:
            tot = Counter()
            for label in r["published family"].split("; "):
                tot.update(by_label.get(label, Counter()))
            missing = [x for x in r["published family"].split("; ") if x not in by_label]
            if missing:
                print("label not in S3a:", missing)
            r["studies (receptors)"] = "; ".join(f"{k} ({v})" for k, v in tot.most_common())
            w.writerow(r)
    for dst, src in COPY.items():
        shutil.copyfile(SRC / src, HERE / dst)
    print("wrote S3a-d")


if __name__ == "__main__":
    main()
