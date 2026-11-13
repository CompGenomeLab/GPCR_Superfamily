"""One figure of family origins: Dollo parsimony and observed presence over 19 lineages.

Stage 06. The genomes are split into 20 non-overlapping lineages, the fewest that still separate
every origin the Dollo sweep calls (plants, SAR, the remaining bikonts, Amoebozoa, five fungal
groups, the unicellular holozoans and nine animal groups). For each family and lineage the page
shows two things side by side:

- Observed presence: a dot whose area is the share of the lineage's genomes (or genera) carrying
  the family, coloured by how we judge the presence (RED and ORANGE below): black for presences we
  rely on, orange for questionable ones, red for contamination or horizontal transfer.
- Our interpretation: the origin implied by the black presences alone (black_origin below).

The root is drawn as a polytomy labelled LECA. Two columns on the right give the Dollo origin at
k = 2 on the NCBI tree (from dollo_origins.tsv, which uses every presence) and our interpretation.
Hovering a row marks our origin on the lineage tree; hovering a cell gives the counts and, for
merged lineages, the breakdown by sub-lineage.

Writes 06_family_origins/lineage_origin_figure.html (Figure 2).

Usage:
    python build_lineage_origin_figure.py
"""

import csv
import json
import os
import re

from ete3 import Tree

from dollo_origins import genus_of_taxid
from key_nodes import OUT_DIR

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Lineages, left to right: (key, label, NCBI nodes or taxIDs included, excluded, sub-lineages).
FONTICULA = "691883"
LINEAGES = [
    ("chloro", "Chlorophyta", ["phylum__Chlorophyta"], [], []),
    ("strepto", "Streptophyta", ["phylum__Streptophyta"], [], []),
    ("stramen", "Stramenopiles", ["clade 1__Stramenopiles"], [], []),
    ("alveol", "Alveolata", ["clade 1__Alveolata"], [], []),
    ("other", "Other eukaryotes", ["root"],
     ["kingdom__Viridiplantae", "clade 1__Stramenopiles", "clade 1__Alveolata", "clade__Unikonta"],
     [("Rhodophyta", "phylum__Rhodophyta"), ("Rhizaria", "clade 1__Rhizaria"),
      ("Haptophyta", "phylum__Haptophyta"), ("Cryptophyta", "905079"), ("Discoba", "clade__Discoba"),
      ("Metamonada", "clade__Metamonada"), ("Apusozoa", "461836")]),
    ("amoeb", "Amoebozoa", ["clade 1__Amoebozoa"], [], []),
    ("earlyfungi", "Early-diverging fungi", ["kingdom__Fungi", FONTICULA],
     ["subkingdom__Dikarya", "phylum__Mucoromycota", "phylum__Zoopagomycota"],
     [("Chytridiomycota", "no rank 1__Chytridiomycota%20incertae%20sedis"),
      ("Blastocladiomycota", "order__Blastocladiales"), ("Microsporidia", "phylum__Microsporidia"),
      ("Cryptomycota", "no rank 1__Cryptomycota%20incertae%20sedis"), ("Fonticula", FONTICULA)]),
    ("mucoro", "Mucoromycota", ["phylum__Mucoromycota"], [], []),
    ("zoopago", "Zoopagomycota", ["phylum__Zoopagomycota"], [], []),
    ("basidio", "Basidiomycota", ["phylum__Basidiomycota"], [], []),
    ("asco", "Ascomycota", ["phylum__Ascomycota"], [], []),
    ("holoprot", "Unicellular holozoans", ["clade 1__Opisthokonta"],
     ["kingdom__Fungi", "kingdom__Metazoa", FONTICULA],
     [("Choanoflagellata", "family__Salpingoecidae"), ("Capsaspora", "595528"),
      ("Ichthyosporea", "667725")]),
    ("porifera", "Porifera", ["subclass__Heteroscleromorpha"], [], []),
    ("cnid", "Placozoa & Cnidaria", ["genus__Trichoplax", "phylum__Cnidaria"], [],
     [("Placozoa", "genus__Trichoplax"), ("Cnidaria", "phylum__Cnidaria")]),
    ("nemat", "Nematoda", ["phylum__Nematoda"], [], []),
    ("protost", "Other protostomes", ["clade 5__Protostomia"], ["phylum__Nematoda"],
     [("Spiralia", "clade 6__Spiralia"), ("Other Ecdysozoa", "clade 6__Ecdysozoa")]),
    ("invdeut", "Invertebrate deuterostomes", ["clade 5__Deuterostomia"], ["subphylum__Craniata"], []),
    ("cyclo", "Cyclostomata", ["subphylum__Craniata"], ["clade 7__Gnathostomata"], []),
    ("chondr", "Chondrichthyes", ["class__Chondrichthyes"], [], []),
    ("osteich", "Osteichthyes", ["clade 9__Euteleostomi"], [], []),
]

# The lineage tree drawn above the columns: (label, children), a child being a lineage key.
HIERARCHY = ("", [
    ("Viridiplantae", ["chloro", "strepto"]),
    ("SAR", ["stramen", "alveol"]),
    "other",
    ("Unikonta", [
        "amoeb",
        ("Opisthokonta", [
            # Fungi are a polytomy in NCBI (early-diverging fungi, Mucoromycota, Zoopagomycota and
            # Dikarya side by side), drawn as such; only Dikarya is resolved. (OTT groups
            # Glomeromycota, within our Mucoromycota, with Dikarya.)
            ("Holomycota", ["earlyfungi", "mucoro", "zoopago", ("Dikarya", ["basidio", "asco"])]),
            ("Holozoa", ["holoprot", ("Metazoa", [
                "porifera", ("Eumetazoa", ["cnid", ("Bilateria", [
                    # Protostomia is left out: no origin falls on it, so Nematoda and the other
                    # protostomes join Bilateria directly.
                    "nemat", "protost",
                    ("Deuterostomia", ["invdeut", ("Craniata", ["cyclo", ("Gnathostomata", [
                        "chondr", "osteich"])])]),
                ])]),
            ])]),
        ]),
    ]),
])

# The families are not split by lineage, since the relationships between them are what the paper
# tests. They follow the order in which the Results text discusses them: the two LECA families,
# Class A, Mth with the Class B families, FSL with Class F, the other metazoan families, the
# families that likely originated within Class A and Class C, and the fungal families last.
# Section titles: "above" draws the title horizontally over the rows, "side" as one line of
# vertical text beside them.
GROUPS = [
    ("Controls", ["RBP1", "TBP", "Galpha", "Gbeta", "Ggamma"], "side"),
    ("GPCR Families", ["cAMP", "classC",
                       "classA",
                       "Mth", "classB2", "classB1",
                       "FSLB", "classF",
                       "GP143", "TM116", "GP157",
                       "Olfactory", "Vomeronasal1", "classT", "Nematode", "Vomeronasal2",
                       "GPR1", "STE3", "STE2"], "side"),
]

DISPLAY = {"Galpha": "Gα", "Gbeta": "Gβ", "Ggamma": "Gγ", "RBP1": "RPB1", "classC": "Class C",
           "Vomeronasal2": "Vomeronasal 2", "FSLB_classF": "FSL + Class F", "FSLB": "FSL",
           "classF": "Class F", "Mth_classB2": "Mth + Class B2", "classB2": "Class B2",
           "classB1": "Class B1", "classA": "Class A", "GPR1_classA": "GPR1 + Class A",
           "Vomeronasal1": "Vomeronasal 1", "classT": "Taste-2"}

# How each presence is drawn: red for contamination or horizontal transfer, orange for presences we
# regard as questionable, black for everything else. Our interpretation of each family's origin rests
# on the black presences only; the Dollo column uses every presence. Entries are NCBI node names,
# taxIDs, or "genus:<name>". Orange is applied strictly: presences that would set a deeper origin but
# are sporadic, fit the family poorly, or sit apart from the family's own sequences in the trees.
SAR, DISCOBA, METAMONADA, HAPTO = "clade__Sar", "clade__Discoba", "clade__Metamonada", "phylum__Haptophyta"
CHYTRIDS = "no rank 1__Chytridiomycota%20incertae%20sedis"
RED = {
    "classA": ["kingdom__Viridiplantae"],  # Taxus chinensis, Olea europaea
    "STE2": ["kingdom__Viridiplantae"],    # Carpinus fangiana
}
ORANGE = {
    # Outside Unikonta: scattered single genera; most sit beside unrelated families in the trees.
    # Ascomycota is 4 of 323 genera; Conidiobolus fits weakly. Basidiomycota stay black.
    "classA": [SAR, DISCOBA, METAMONADA, HAPTO, "phylum__Ascomycota",
               "phylum__Mucoromycota", "phylum__Zoopagomycota"],
    # Romanomermis: also claimed by Class A, base of the Nematode tree. Adineta (rotifer): claimed by
    # Class A and Nematode, placed with Olfactory and Class A in the superfamily trees.
    "Nematode": ["genus:Romanomermis", "genus:Adineta"],
    "classB2": ["genus:Achlya", "genus:Auricularia"],
    "GP143": ["genus:Capsaspora"],         # single genome, weak fit
    "GP157": ["genus:Triparma", "genus:Anaeramoeba", "genus:Acanthamoeba", "genus:Capsaspora"],
    "GPR1": ["genus:Chaetoceros", "genus:Planoprotostelium", CHYTRIDS],
}

# Origins are named by their NCBI node, so both columns use the same taxonomy; only the root is
# relabelled, as LECA. (Class B1's NCBI "Eumetazoa" includes Placozoa, i.e. ParaHoxozoa, which NCBI
# does not have as a node; the text should say so.)
def node_label(name):
    return "LECA" if name == "root" else name.split("__")[-1].replace("%20", " ")


def black_origin(tree, black):
    """Common ancestor of the black carriers on the NCBI tree."""
    if not black:
        return "", set()
    origin = tree.get_common_ancestor([tree.leaf_by_name[t] for t in black]) if len(black) > 1         else tree.leaf_by_name[next(iter(black))]
    return node_label(origin.name), set(origin.get_leaf_names())


def leaves_of(tree, nodes, spec, genus=None):
    out = set()
    for name in spec:
        if name.startswith("genus:"):
            out |= {t for t, g in genus.items() if g == name[6:]}
        else:
            out |= set(nodes[name].get_leaf_names()) if name in nodes else {name}
    return out


def origin_text(rows, family):
    """The Dollo origin at k = 2 on the NCBI tree, as a plain node name; the root is LECA."""
    first = next(r for r in rows if r["family"] == family and r["taxonomy"] == "ncbi")
    return re.sub(r" \(in .*\)$", "", first["origin_label"])


def main():
    tree = Tree(os.path.join(OUT_DIR, "taxonomy_ncbi_all.nwk"), format=1)
    tree.name = tree.name or "root"  # the input tree leaves the root unnamed
    tree.leaf_by_name = {leaf.name: leaf for leaf in tree}
    nodes = {n.name: n for n in tree.traverse() if not n.is_leaf()}
    genus = genus_of_taxid()
    with open(os.path.join(OUT_DIR, "family_presence_absence_ncbi_all.csv")) as handle:
        table = {r["taxID"]: r for r in csv.DictReader(handle)}
    with open(os.path.join(OUT_DIR, "dollo_origins.tsv")) as handle:
        dollo_rows = list(csv.DictReader(handle, delimiter="\t"))

    lineages, seen = [], set()
    for key, label, include, exclude, subs in LINEAGES:
        members = leaves_of(tree, nodes, include) - leaves_of(tree, nodes, exclude)
        assert not members & seen, key
        seen |= members
        lineages.append({"key": key, "label": label, "members": members,
                         "subs": [(name, leaves_of(tree, nodes, [spec]) & members)
                                  for name, spec in subs]})
    assert seen == set(tree.get_leaf_names()), "lineages must cover every genome"

    families = [f for _, fams, _ in GROUPS for f in fams]
    rows = []
    for family in families:
        carriers = {t for t, r in table.items() if r[family] == "1"}
        red = carriers & leaves_of(tree, nodes, RED.get(family, []), genus)
        orange = (carriers & leaves_of(tree, nodes, ORANGE.get(family, []), genus)) - red
        black = carriers - red - orange
        interp, call = black_origin(tree, black)
        cells = []
        for lin in lineages:
            m = lin["members"]
            hit = m & carriers
            cell = {
                "n": len(m), "k": len(hit),
                "gn": len({genus[t] for t in m}), "gk": len({genus[t] for t in hit}),
                "black": len(hit & black), "orange": len(hit & orange), "red": len(hit & red),
                "subs": [[name, len(s), sum(table[t][family] == "1" for t in s)]
                         for name, s in lin["subs"] if s],
            }
            cells.append(cell)
        rows.append({"family": family, "name": DISPLAY.get(family, family),
                     "dollo": origin_text(dollo_rows, family), "interp": interp, "cells": cells,
                     "call": sorted(i for i, lin in enumerate(lineages) if lin["members"] & call)})

    data = {
        "lineages": [{"key": l["key"], "label": l["label"], "n": len(l["members"]),
                      "g": len({genus[t] for t in l["members"]})} for l in lineages],
        "hierarchy": HIERARCHY,
        "groups": [[name, [DISPLAY.get(f, f) for f in fams], how] for name, fams, how in GROUPS],
        "rows": rows,
    }
    with open(os.path.join(SCRIPT_DIR, "lineage_origin_template.html"), encoding="utf-8") as handle:
        page = handle.read().replace("__DATA__", json.dumps(data, ensure_ascii=False))
    out = os.path.join(OUT_DIR, "lineage_origin_figure.html")
    with open(out, "w", encoding="utf-8") as handle:
        handle.write(page)
    print("wrote %s" % out)


if __name__ == "__main__":
    main()
