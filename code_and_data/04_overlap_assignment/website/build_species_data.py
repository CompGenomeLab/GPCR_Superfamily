#!/usr/bin/env python3
"""Build the data of the per-species page of gpcrevolution.org (species).

For every searched genome (06_family_origins/taxonomy/searched_taxids_2025_02.tsv) this lists every
potential GPCR, i.e. every member of an extracted clade (stage 02), with its search E-value (stage 01)
and what became of it in each family whose clade it was in:

    kept        a final member of the family
    core        removed: no single GPCR-like 7TM core (03_topology_filter/core_membership.tsv)
    topology    removed: topology score below 0.8 (03_topology_filter/topology_membership.tsv)
    overlap     claimed by two families and assigned to the other one (04_overlap_assignment); every
                row of a claimed sequence lists the other families that claimed it
    dikarya     Dikarya Class A moved to GPR1 (10_motif_analysis/phylogeny_reassignment)

and the members of the control clades (G alpha, G beta, G gamma, TBP, RPB1), which are not filtered.
Counts per family, per lineage and in total are the final members after the GPR1 reassignment, the
counts of the website's phyletic distribution page. Each genome is also placed in its Figure 2 lineage
(LINEAGES of 06_family_origins/build_lineage_origin_figure.py).

Writes <out>/index.json (one entry per genome), <out>/species/<taxID>.json and <out>/search.json (the
names, with the NCBI common names, that the site-wide species search completes).

Usage:
    python build_species_data.py --out <GPCRevolution>/public/species
"""
import argparse
import csv
import json
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parents[1] / "06_family_origins"))
from gpcr_paths import CONTROLS, DISPLAY, FAMILIES, open_text, stage  # noqa: E402
from build_taxonomy_trees import MERGED_TAXIDS  # noqa: E402  UniProt taxIDs NCBI has since merged

LINEAGE_OF_FAMILY = {
    "Glutamate-related": ["classC", "Vomeronasal2"],
    "cAMP-related": ["cAMP", "GPR1", "GP157", "GP143", "STE3", "TM116",
                     "classB1", "classB2", "Mth", "classF", "FSLB"],
    "Rhodopsin-related": ["classA", "Olfactory", "Vomeronasal1", "classT", "Nematode"],
    "STE2 (Class D)": ["STE2"],
}

# STE2 is shown with its name in the A-F classification, on this page only.
PAGE_DISPLAY = {**DISPLAY, "STE2": "STE2 (Class D)"}

# The search each family's clade was extracted from: (folder, kind). Vomeronasal 1 and Taste-2 share
# the TAS2R search, Class B1 and B2 the 7tm_2 search.
SEARCH = {f: (f, "hmm") for f in FAMILIES + CONTROLS}
SEARCH.update({"Vomeronasal1": ("classT", "hmm"), "classB1": ("classB", "hmm"),
               "classB2": ("classB", "hmm"), "FSLB": ("FSLB_DICDI", "m8"),
               "GPR1": ("GPR1_YEAST", "m8")})

HEADER_RE = re.compile(r"^(?:sp|tr)\|([^|]+)\|\S+\s*(.*?)\s+OS=(.*?)\s+OX=(\d+)(?:\s+GN=(\S+))?")


def wide(path):
    with open_text(path) as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    table = {name: set() for name in rows[0]}
    for row in rows[1:]:
        for name, value in zip(rows[0], row):
            if value.strip():
                table[name].add(value.strip())
    return table


def read_search(folder, kind, wanted):
    """{accession: (E-value as printed, description, taxID, gene)} for the wanted accessions."""
    out = {}
    base = stage("01_homology_search", folder)
    if kind == "hmm":
        with open_text(base / f"{folder}_hmmsearch_E1000.txt") as handle:
            in_table = False
            for line in handle:
                if line.lstrip().startswith("-------"):
                    in_table = True
                    continue
                if not in_table or "inclusion threshold" in line:
                    continue
                if not line.strip():
                    break
                parts = line.split(None, 9)
                m = HEADER_RE.match(parts[8] + " " + (parts[9] if len(parts) > 9 else ""))
                if m and m.group(1) in wanted and m.group(1) not in out:
                    out[m.group(1)] = (parts[0], m.group(2), m.group(4), m.group(5) or "")
    else:
        with open_text(base / f"{folder}_UniProt_s7.5_filter0_eprofile0.1_db.m8") as handle:
            for line in handle:
                cols = line.rstrip("\n").split("\t")
                m = HEADER_RE.match(cols[2])
                if not m or m.group(1) not in wanted:
                    continue
                acc, ev = m.group(1), cols[-1]
                if acc not in out or float(ev) < float(out[acc][0]):
                    out[acc] = (ev, m.group(2), m.group(4), m.group(5) or "")
    return out


def taxid_names():
    """{taxID: scientific name} and {taxID: common name}; the GenBank common name is preferred."""
    names, common = {}, {}
    with open_text(stage("06_family_origins", "taxonomy", "ncbi_taxdump", "names.dmp")) as handle:
        for line in handle:
            cols = line.split("\t|\t")
            if cols[3].startswith("scientific name"):
                names[cols[0]] = cols[1]
            elif cols[3].startswith("genbank common name"):
                common[cols[0]] = cols[1]
            elif cols[3].startswith("common name"):
                common.setdefault(cols[0], cols[1])
    return names, common


def taxonomy_paths(taxids):
    """{taxID: [[taxID, name, rank], ...]} from Eukaryota down to the genome, in the taxonomy of the
    paper: the NCBI taxdump with the synthetic groupings of ncbi_taxonomy_create.CUSTOM_GROUPS
    (Unikonta, Archaeplastida, Holozoa, Holomycota) inserted. Synthetic nodes have negative taxIDs."""
    import ncbi_taxonomy_create as taxo
    parent, rank = {}, {}
    with open_text(stage("06_family_origins", "taxonomy", "ncbi_taxdump", "nodes.dmp")) as handle:
        for line in handle:
            p = taxo.split_dmp_line(line)
            parent[int(p[0])], rank[int(p[0])] = int(p[1]), p[2]
    names = {}
    with open_text(stage("06_family_origins", "taxonomy", "ncbi_taxdump", "names.dmp")) as handle:
        for line in handle:
            p = taxo.split_dmp_line(line)
            if p[3] == "scientific name":
                names[int(p[0])] = p[1]
    taxo.apply_custom_groups(parent, rank, names)
    return {t: [[n, names.get(n, str(n)), rank.get(n, "")]
                for n in taxo.lineage_to_ancestor(int(MERGED_TAXIDS.get(t, t)), parent)] for t in taxids}


def figure2_lineages():
    """{taxID: Figure 2 lineage label}, with the sub-lineage where Figure 2 gives one."""
    from ete3 import Tree
    import build_lineage_origin_figure as fig2
    from key_nodes import OUT_DIR
    tree = Tree(os.path.join(OUT_DIR, "taxonomy_ncbi_all.nwk"), format=1)
    tree.name = tree.name or "root"
    nodes = {n.name: n for n in tree.traverse() if not n.is_leaf()}
    out, order = {}, []
    for key, label, include, exclude, subs in fig2.LINEAGES:
        members = fig2.leaves_of(tree, nodes, include) - fig2.leaves_of(tree, nodes, exclude)
        order.append(label)
        for t in members:
            out[t] = [label, ""]
        for name, spec in subs:
            for t in fig2.leaves_of(tree, nodes, [spec]) & members:
                out[t][1] = name
    return out, order


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(os.path.join(args.out, "species"), exist_ok=True)

    clade = wide(stage("02_clade_extraction", "clade_membership.tsv"))
    controls = wide(stage("02_clade_extraction", "clade_membership_controls.tsv"))
    core = wide(stage("03_topology_filter", "core_membership.tsv"))
    topo = wide(stage("03_topology_filter", "topology_membership.tsv"))
    final = wide(stage("04_overlap_assignment", "final_membership.tsv"))
    reassigned = wide(stage("10_motif_analysis", "phylogeny_reassignment",
                            "final_membership_reassigned.tsv"))
    final_of = {a: f for f, accs in reassigned.items() for a in accs}

    scores = {}
    with open_text(stage("03_topology_filter", "topology_helix_scores.tsv")) as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            scores[(r["family"], r["accession"])] = round(float(r["mean_jaccard"]), 3)

    wanted = defaultdict(set)
    for fam, accs in list(clade.items()) + list(controls.items()):
        wanted[SEARCH[fam]] |= accs
    info = {}
    for (folder, kind), accs in wanted.items():
        print("reading search", folder, flush=True)
        found = read_search(folder, kind, accs)
        missing = accs - set(found)
        if missing:
            print(f"  WARNING: {len(missing)} accessions not in the {folder} search output")
        for acc, rec in found.items():
            info[(folder, acc)] = rec

    per_species = defaultdict(lambda: {"gpcr": [], "controls": []})
    status_counts = defaultdict(int)
    for fam in FAMILIES:
        folder = SEARCH[fam][0]
        for acc in sorted(clade[fam]):
            ev, desc, taxid, gene = info[(folder, acc)]
            if acc not in core[fam]:
                status = "core"
            elif acc not in topo[fam]:
                status = "topology"
            elif acc not in final[fam]:
                status = "overlap"
            elif final_of.get(acc) != fam:
                status = "dikarya"
            else:
                status = "kept"
            if status in ("overlap", "dikarya") and acc not in final_of:
                raise SystemExit(f"{acc} ({fam}) left {status} without a final family")
            status_counts[status] += 1
            # the other families whose topology-filtered clade also held it (the overlap step)
            partners = [f for f in FAMILIES if f != fam and acc in topo[f]] if acc in topo[fam] else []
            per_species[taxid]["gpcr"].append(
                [acc, fam, ev, scores.get((fam, acc)), status, final_of.get(acc, ""), desc, gene, partners])
    for fam in CONTROLS:
        folder = SEARCH[fam][0]
        for acc in sorted(controls[fam]):
            ev, desc, taxid, gene = info[(folder, acc)]
            per_species[taxid]["controls"].append([acc, fam, ev, desc, gene])
    print("rows by status:", dict(status_counts))

    with open_text(stage("06_family_origins", "taxonomy", "searched_taxids_2025_02.tsv")) as handle:
        searched = [line.split("\t")[0].strip() for line in handle if line.strip()]
    searched = [t for t in searched if t.isdigit()]
    extra = set(per_species) - set(searched)
    if extra:
        raise SystemExit(f"{len(extra)} taxIDs with sequences are not in the searched list: {sorted(extra)[:5]}")
    names, common = taxid_names()
    paths = taxonomy_paths(searched)
    lineage, lineage_order = figure2_lineages()

    index = []
    for taxid in searched:
        data = per_species.get(taxid, {"gpcr": [], "controls": []})
        counts = defaultdict(int)
        for acc in {row[0] for row in data["gpcr"]}:
            if acc in final_of:
                counts[final_of[acc]] += 1
        ctrl = defaultdict(int)
        for row in data["controls"]:
            ctrl[row[1]] += 1
        lin, sub = lineage.get(taxid, ["", ""])
        entry = {"taxID": taxid, "name": names.get(MERGED_TAXIDS.get(taxid, taxid), taxid), "lineage": lin, "sublineage": sub,
                 "candidates": len({r[0] for r in data["gpcr"]}),
                 "total": sum(counts.values()),
                 "counts": {f: counts[f] for f in FAMILIES if counts[f]},
                 "controls": {f: ctrl[f] for f in CONTROLS if ctrl[f]}}
        index.append(entry)
        with open(os.path.join(args.out, "species", f"{taxid}.json"), "w", encoding="utf-8") as handle:
            json.dump({**entry, "taxonomy": paths[taxid], "gpcr": data["gpcr"], "controlRows": data["controls"]}, handle,
                      separators=(",", ":"), ensure_ascii=False)

    meta = {
        "families": [{"key": f, "display": PAGE_DISPLAY.get(f, f),
                      "lineage": next(l for l, fs in LINEAGE_OF_FAMILY.items() if f in fs)}
                     for f in FAMILIES],
        "controls": [{"key": f, "display": DISPLAY.get(f, f)} for f in CONTROLS],
        "lineages": list(LINEAGE_OF_FAMILY),
        "speciesLineages": lineage_order,
        "gpcrColumns": ["accession", "family", "evalue", "topologyScore", "status", "finalFamily",
                        "description", "gene", "overlapFamilies"],
        "controlColumns": ["accession", "control", "evalue", "description", "gene"],
        "species": index,
    }
    with open(os.path.join(args.out, "index.json"), "w", encoding="utf-8") as handle:
        json.dump(meta, handle, separators=(",", ":"), ensure_ascii=False)
    # the small list the site-wide species search reads: taxID, name, common name, lineage, GPCRs
    search = [[e["taxID"], e["name"], common.get(MERGED_TAXIDS.get(e["taxID"], e["taxID"]), ""), e["lineage"],
               e["total"]] for e in index]
    with open(os.path.join(args.out, "search.json"), "w", encoding="utf-8") as handle:
        json.dump(search, handle, separators=(",", ":"), ensure_ascii=False)
    print(f"{len(index)} genomes written to {args.out}")


if __name__ == "__main__":
    main()
