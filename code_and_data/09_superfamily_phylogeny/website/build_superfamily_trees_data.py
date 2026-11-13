#!/usr/bin/env python3
"""Build the data of the superfamily-trees page of gpcrevolution.org (/superfamily-trees).

Reads the eight final trees as build_tree_viewer_data.py writes them (build/tree_viewer_data/<tree>.js,
after plot_circular_trees.py) and writes one compact
JSON per tree, <out>/<tree>.json, and <out>/index.json listing them. Nothing about a tree is recomputed:
node order, parents, root-to-tip depths, TBE, family, transfer indices and family colours are copied.
Added for the page:
  - each leaf's lineage group, the five groups of the Figure 3 legend, in place of the highest rank the
    taxonomy record happened to carry (which split nested groups such as Discoba and Heterolobosea);
  - for the Class A leaves later reassigned to GPR1 (the Dikarya sequences of the GPR1 clade;
    10_motif_analysis/phylogeny_reassignment/final_membership_reassigned.tsv), `fam2` = GPR1, so the page
    can show the trees' labels or the final assignment;
  - for the 18 literature sequences, the family their study reported and our closest family, the
    family of their best MMseqs2 hit among our family members (05_literature_comparison/
    literature_closest_family.tsv).

Usage:
    python plot_circular_trees.py && python build_tree_viewer_data.py
    python build_superfamily_trees_data.py --out <GPCRevolution>/public/superfamily-trees
"""
import argparse
import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from gpcr_paths import stage  # noqa: E402

SETS = [("no_ste2ste3classCV2", "Core families", "without STE2, STE3, Class C and Vomeronasal 2"),
        ("and_literature", "Core families + literature", "the core set and 18 literature sequences"),
        ("no_ste2", "All families but STE2", "all families except STE2"),
        ("all", "All 19 families", "every family")]
ALIGNMENTS = [("ginsi_ep0.123", "Higher gap penalty", "MAFFT G-INS-i, --ep 0.123 (main)"),
              ("ginsi", "Default gap penalty", "MAFFT G-INS-i, default --ep")]
GROUPS = ["Metazoa", "Fungi", "Amoebozoa and unicellular holozoans", "Viridiplantae", "SAR and other eukaryotes"]
GROUP_OF = {"Metazoa": GROUPS[0], "Fungi": GROUPS[1], "Viridiplantae": GROUPS[3]}
for lin in ("Amoebozoa", "Evosea", "Discosea", "Filasterea", "Choanoflagellata", "Ichthyosporea"):
    GROUP_OF[lin] = GROUPS[2]
for lin in ("SAR", "Haptophyta", "Haptista", "Rhodophyta", "Heterolobosea", "Discoba", "Anaeramoebidae",
            "Parabasalia", "Fornicata", "Metamonada", "CRuMs", "Apusomonadida", "Apusozoa"):
    GROUP_OF[lin] = GROUPS[4]
TAX_KEEP = ("kingdom", "phylum", "class", "order", "species")


def load(path):
    s = path.read_text(encoding="utf-8")
    return json.loads(s[s.index("{"):s.rstrip().rstrip(";").rindex("}") + 1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    out = Path(ap.parse_args().out)
    out.mkdir(parents=True, exist_ok=True)
    src = HERE / "build" / "tree_viewer_data"
    rows = list(csv.reader(open(stage("10_motif_analysis", "phylogeny_reassignment", "final_membership_reassigned.tsv"),
                                encoding="utf-8"), delimiter="	"))
    final_fam = {acc: fam for row in rows[1:] for fam, acc in zip(rows[0], row) if acc}
    closest = {r["id"]: r for r in csv.DictReader(
        open(stage("05_literature_comparison", "literature_closest_family.tsv"), encoding="utf-8"), delimiter="\t")}

    index = []
    for set_key, set_name, set_note in SETS:
        for aln_key, aln_name, aln_note in ALIGNMENTS:
            name = f"representatives_{set_key}_{aln_key}_kpicgappy.raxml"
            T = load(src / f"{name}.js")
            colors = {}
            for d in T["D"]:
                if d.get("leaf"):
                    colors.setdefault(d["fam"], d["rc"][0])
            nodes = []
            for d in T["D"]:
                n = {"p": d["p"], "rt": d["rt"]}
                if not d.get("leaf"):
                    n["s"] = d.get("sup")
                    nodes.append(n)
                    continue
                tax = dict(d["tax"])
                n.update(leaf=1, acc=d["acc"], fam=d["fam"], grp=GROUP_OF[d["lin"]], lin=d["lin"],
                         ti=d.get("ti"), mti=d.get("mti"))
                if not d.get("lit") and final_fam.get(d["acc"], d["fam"]) != d["fam"]:
                    n["fam2"] = final_fam[d["acc"]]
                if d.get("lit"):
                    c = closest[d["acc"]]
                    n.update(lit=1, sp=tax.get("species", ""), claimed=tax.get("claimed", ""),
                             study=tax.get("study") or c["study"], closest=c["closest_family"])
                else:
                    n["tax"] = [[k, tax[k]] for k in TAX_KEEP if k in tax]
                nodes.append(n)
            leaves = sum(1 for n in nodes if n.get("leaf"))
            (out / f"{name}.json").write_text(json.dumps(
                {"name": name, "ste2Cut": T["ste2Cut"], "colors": colors, "nodes": nodes},
                separators=(",", ":")), encoding="utf-8")
            index.append({"file": name, "set": set_key, "setName": set_name, "setNote": set_note,
                          "alignment": aln_key, "alignmentName": aln_name, "alignmentNote": aln_note,
                          "leaves": leaves})
            print(name, leaves, "leaves;", sum(1 for n in nodes if n.get("fam2")), "relabelled",
                  sorted({(n["fam"], n["fam2"]) for n in nodes if n.get("fam2")}))
    (out / "index.json").write_text(json.dumps({"groups": GROUPS, "trees": index}, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
