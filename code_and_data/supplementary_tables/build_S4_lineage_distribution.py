"""Supplementary Table S4: the taxonomic distribution behind Figure 2.

One row per family (and control) and lineage of Figure 2: genomes sampled and genomes carrying the
family, genera sampled and genera carrying it, how the presence is drawn in the figure, the number of
carrying genomes of each kind, and the family's Dollo origin and our interpretation. The data are
read from the figure page itself (06_family_origins/lineage_origin_figure.html, built by
build_lineage_origin_figure.py), so the table and the figure cannot disagree.

A lineage is drawn as a presence when any of its carrying genomes is a presence our interpretation
rests on (black), as sporadic when its carriers are only sporadic ones (orange), and as contamination
or horizontal gene transfer when they are only those (red).

    python build_S4_lineage_distribution.py
"""
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import stage  # noqa: E402

FIGURE = stage("06_family_origins", "lineage_origin_figure.html")
CONTROLS = {"RPB1", "TBP", "Gα", "Gβ", "Gγ"}


def figure_data():
    page = FIGURE.read_text(encoding="utf-8")
    start = page.index("{", page.index("const DATA = "))
    data, _ = json.JSONDecoder().raw_decode(page[start:])
    return data


def presence(cell):
    if cell["k"] == 0:
        return "absent"
    if cell["black"]:
        return "present"
    return "sporadic" if cell["orange"] else "contamination or HGT"


def main():
    data = figure_data()
    lineages = data["lineages"]
    rows = []
    for r in data["rows"]:
        for lin, cell in zip(lineages, r["cells"]):
            rows.append({
                "Family": r["name"],
                "Type": "control" if r["name"] in CONTROLS else "GPCR family",
                "Lineage": lin["label"],
                "Genomes sampled": cell["n"],
                "Genomes with the family": cell["k"],
                "Genera sampled": cell["gn"],
                "Genera with the family": cell["gk"],
                "Presence in Figure 2": presence(cell),
                "Genomes counted as present": cell["black"],
                "Genomes counted as sporadic": cell["orange"],
                "Genomes counted as contamination or HGT": cell["red"],
                "Dollo origin": r["dollo"],
                "Our interpretation of the origin": r["interp"],
            })
    out = stage("supplementary_tables", "Table_S4_lineage_distribution.tsv")
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print("wrote %s: %d families x %d lineages" % (out.name, len(data["rows"]), len(lineages)))


if __name__ == "__main__":
    main()
