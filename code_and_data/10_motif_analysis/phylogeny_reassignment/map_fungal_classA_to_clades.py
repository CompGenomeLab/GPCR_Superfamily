"""Map every fungal Class A sequence of the Sep24 final membership to the GPR1 or Class A clade.

References are the tree representatives (Sep22 all set): the GPR1 clade is the GPR1 representatives
plus the Class A representatives on the GPR1 side in >= 6 of 8 best trees (tree_side_per_sequence.tsv);
the Class A clade is the other Class A-like representatives. The three Dikarya representatives on the
Class A side are left out as references. Each query is assigned to the clade of its best local
alignment (BLOSUM62, gap open 11, extend 1), self-hits excluded.

Reads  04_overlap_assignment/final_membership.tsv and acc_taxid_map.tsv,
       07_representative_selection/taxonomy_eukaryotes_filtered.json,
       02_clade_extraction/clade_sequences/classA.fasta.gz,
       10_motif_analysis/phylogeny_reassignment/inputs_Sep22/ (Sep22 representatives and their alignment)
       and tree_side_per_sequence.tsv (run tree_side_per_sequence.py first)
Writes 10_motif_analysis/phylogeny_reassignment/fungal_classA_members.tsv and fungal_classA_clade_map.tsv;
the Dikarya rows assigned to GPR1 are what phylogeny_reassignments.tsv lists.

Usage:
    python map_fungal_classA_to_clades.py      (needs Biopython)
"""
import csv
import json
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from gpcr_paths import open_text, stage  # noqa: E402

S4 = stage("04_overlap_assignment")
OUT = stage("10_motif_analysis", "phylogeny_reassignment")
IN = OUT / "inputs_Sep22"
A_LIKE = {"classA", "Olfactory", "Vomeronasal1", "classT", "Nematode"}
EXCEPT = {"A0A8H3VW39", "A0A9P8ZU35", "A0A4S8LVH9"}

# fungal Class A members, Sep24 final membership
with open_text(stage("07_representative_selection", "taxonomy_eukaryotes_filtered.json")) as fh:
    tax = {r["taxID"]: r for r in json.load(fh)}
with open_text(S4 / "final_membership.tsv") as fh:
    rows = list(csv.reader(fh, delimiter="\t"))
i = rows[0].index("classA")
classA = {r[i] for r in rows[1:] if i < len(r) and r[i]}
with open_text(S4 / "acc_taxid_map.tsv") as fh:
    tid = {r["accession"]: int(r["taxID"]) for r in csv.DictReader(fh, delimiter="\t")}
queries = {}
for a in sorted(classA):
    t = tax.get(tid.get(a))
    if t and t["Kingdom"] == "Fungi":
        queries[a] = {"phylum": t["Phylum"], "class": t["Class"], "order": t["Order"], "species": t["Species"], "taxID": t["taxID"]}
with open(OUT / "fungal_classA_members.tsv", "w", encoding="utf-8", newline="\n") as out:
    out.write("family\taccession\ttaxID\tphylum\tclass\torder\tspecies\n")
    for a, q in queries.items():
        out.write("\t".join(["classA", a, str(q["taxID"]), q["phylum"], q["class"], q["order"], q["species"]]) + "\n")

# references
with open_text(IN / "representatives_all_families_Sep22_topo0.8.tsv") as fh:
    reps = {r["accession"]: r for r in csv.DictReader(fh, delimiter="\t")}
with open_text(OUT / "tree_side_per_sequence.tsv") as fh:
    side = {r["acc"]: r for r in csv.DictReader(fh, delimiter="\t")}
ref = {}
with open_text(IN / "representatives_all_Sep22_topo0.8_ginsi.fasta") as fh:
    msa = list(SeqIO.parse(fh, "fasta"))
for rec in msa:
    a = rec.id.split("|")[1] if "|" in rec.id else rec.id
    f = reps.get(a, {}).get("gpcrfamily")
    if a in EXCEPT:
        continue
    if f == "GPR1" or (f == "classA" and a in side and int(side[a]["best_GPR1"]) >= 6):
        clade = "GPR1"
    elif f in A_LIKE:
        clade = "ClassA"
    else:
        continue
    ref[a] = (clade, f, reps[a]["kingdom"], reps[a]["species"], str(rec.seq).replace("-", "").upper())

seqs = {}
with open_text(stage("02_clade_extraction", "clade_sequences", "classA.fasta")) as fh:
    for rec in SeqIO.parse(fh, "fasta"):
        a = rec.id.split("|")[1] if "|" in rec.id else rec.id
        if a in queries:
            seqs[a] = str(rec.seq).upper().replace("*", "")

al = PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                     open_gap_score=-11, extend_gap_score=-1)
with open(OUT / "fungal_classA_clade_map.tsv", "w", encoding="utf-8", newline="\n") as out:
    out.write("accession\tphylum\tclass\torder\tspecies\tis_rep\ttree_side\tbest_GPR1\tbest_GPR1_hit\t"
              "best_ClassA\tbest_ClassA_hit\tassigned\tmargin\n")
    for a, q in queries.items():
        s = seqs.get(a)
        if not s:
            continue
        best = {"GPR1": (-1, ""), "ClassA": (-1, "")}
        for ra, (clade, f, k, spn, rs) in ref.items():
            if ra == a:
                continue
            sc = al.score(s, rs)
            if sc > best[clade][0]:
                best[clade] = (sc, f"{f}:{k}:{spn}")
        g, c = best["GPR1"][0], best["ClassA"][0]
        ts = f"GPR1 {side[a]['best_GPR1']}/8, ClassA {side[a]['best_ClassA']}/8" if a in side else ""
        out.write("\t".join([a, q["phylum"], q["class"], q["order"], q["species"], str(a in reps), ts,
                             f"{g:.0f}", best["GPR1"][1], f"{c:.0f}", best["ClassA"][1],
                             "GPR1" if g > c else "ClassA", f"{(g - c) / max(g, c):.3f}"]) + "\n")
print("wrote", OUT / "fungal_classA_members.tsv", "and", OUT / "fungal_classA_clade_map.tsv")
