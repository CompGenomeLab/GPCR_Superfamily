"""Supplementary table of the literature comparison: one row per receptor listed by a study, each step
of the comparison and of our pipeline in its own column.
Inputs: literature_sequences.tsv (compare_literature.py: genome, sequence source, ID matches and their
pipeline steps) and literature_all_outcomes.tsv (classify_all_literature.py: the MMseqs2 matching of
every receptor not matched through an ID). Pipeline steps of MMseqs2-matched proteins are read from the
same stage tables, search outputs and trees as compare_literature.py.
Writes 05_literature_comparison/Supplementary_Table_literature_receptors.tsv and .xlsx (README sheet + table)."""
import csv, os, sys

STEP6 = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, STEP6)
import compare_literature as C

LO = C.OUT_DIR
csv.field_size_limit(10**9)
read = lambda n: list(csv.DictReader(open(os.path.join(LO, n), encoding="utf-8"), delimiter="\t"))
seqs, outs = read("literature_sequences.tsv"), read("literature_all_outcomes.tsv")
assert len(seqs) == len(outs) and all((a["study"], a["their_family"], a["id"]) == (b["study"], b["their_family"], b["id"])
                                      for a, b in zip(seqs, outs))

mm = {o["our_accession"] for o in outs if o["route"] in ("species match", "genus match")}
stages = [(name, C.read_stage(path)) for name, path in C.STAGES]
_, hits, _, _ = C.scan_searches(set(), mm)
trees = C.read_tree_members(mm)

COLS = [
    ("study", "Study"), ("study_family", "Family as named by the study"),
    ("our_family", "Our corresponding family (none = no equivalent among our families; any = study's non-canonical candidates)"),
    ("receptor_id", "Receptor identifier as given by the study"),
    ("genome", "Genome (resolved to the NCBI name)"),
    ("genome_in_our_set", "yes / sister genome (same genus) / no"),
    ("sequence_source", "Where the receptor's sequence comes from"),
    ("step1_id_match", "Step 1. Matched through a UniProt accession or identifier the study gives: the identifier, or 'no'"),
    ("step2_mmseqs2_query", "Step 2. MMseqs2 query against our full proteome set (not run when step 1 matched; 'no sequence' when the study gives none we could use)"),
    ("step3_sequence_match", "Step 3. same species / same genus (>= 0.9 identical residues over the longer sequence) / no match / no significant hit"),
    ("matched_protein", "Our protein the receptor is traced to (step 1 or step 3)"),
    ("matched_species", "Species of the matched protein (step 3)"),
    ("identity_over_longer", "Identical residues over the length of the longer sequence (step 3)"),
    ("evalue", "MMseqs2 E-value of the match or of the nearest hit"),
    ("query_coverage", "MMseqs2 query coverage of the match or of the nearest hit"),
    ("step4_our_searches", "Step 4. Our searches (HMMER/MMseqs2) that hit the matched protein"),
    ("step5_clade", "Step 5. Family clade(s) of our search trees holding the protein ('-' = none)"),
    ("step6_7TM_filter", "Step 6. Families in which it passed the 7TM filter ('-' = removed / not reached)"),
    ("step7_topology_filter", "Step 7. Families in which it passed the topology filter (0.8)"),
    ("step8_final_family", "Step 8. Final family after overlap resolution"),
    ("outcome", "Where the receptor ends in the comparison"),
    ("nearest_family_indicator", "Unmatched receptors only: family of the best significant hit (E <= 1e-5, query coverage >= 0.5); an indicator, not membership"),
    ("explanation", "Each step the receptor passed, in words"),
]

dash = lambda v: v if v else "-"
rows = []
for s, o in zip(seqs, outs):
    fams = None if s["our_families"].startswith("none") else (["*"] if s["our_families"] == "any" else s["our_families"].split(","))
    gin = "no" if s["stopped_at"].startswith(("genome not in our set", "no species")) else (
        "sister genome" if "sister" in s["genome_how"] else "yes")
    row = dict(study=s["study"], study_family=s["their_family"], our_family=s["our_families"], receptor_id=s["id"],
               genome=s["genome"], genome_in_our_set=gin, sequence_source=s["sequence_source"],
               step1_id_match="no", step2_mmseqs2_query="", step3_sequence_match="", matched_protein="",
               matched_species="", identity_over_longer="", evalue="", query_coverage="", step4_our_searches="",
               step5_clade="", step6_7TM_filter="", step7_topology_filter="", step8_final_family="", outcome="",
               nearest_family_indicator="", explanation="")
    route = o["route"]
    if route == "ID match":
        row.update(step1_id_match=s["match"], step2_mmseqs2_query="not run (step 1 matched)", step3_sequence_match="not run",
                   matched_protein=s["our_accession"], step4_our_searches=dash(s["our_search_hits"]),
                   step5_clade=dash(s["our_clade"]), step6_7TM_filter=dash(s["our_7TM"]),
                   step7_topology_filter=dash(s["our_topology"]), step8_final_family=dash(s["our_final"]),
                   outcome=s["stopped_at"], explanation=s["explanation"])
    elif route in ("species match", "genus match"):
        acc = o["our_accession"]
        sf = {n: sorted(m.get(acc, set())) for n, m in stages}
        searched = sorted(hits.get(acc, ()))
        r = {"match": "MMseqs2, %s, %s identical over the longer sequence" % (route.replace(" match", ""), o["identity_over_longer"]),
             "named": True, "genome": s["genome"]}
        stopped, expl = C.explain(r, fams, acc, True, sf, trees, searched)
        if route == "genus match":
            expl = expl.replace("genome in our set; ", (
                "genome in our set but no same-species protein at >= 0.9 identity; " if gin == "yes" else "genome not in our set; ")
                + "matched in the same genus (%s); " % o["hit_species"], 1)
        row.update(step2_mmseqs2_query=o["hit_query_id"], step3_sequence_match="same " + route.split()[0],
                   matched_protein=acc, matched_species=o["hit_species"], identity_over_longer=o["identity_over_longer"],
                   evalue=o["evalue"], query_coverage=o["query_coverage"], step4_our_searches=dash(",".join(searched)),
                   step5_clade=dash(",".join(sf["clade"])), step6_7TM_filter=dash(",".join(sf["7TM"])),
                   step7_topology_filter=dash(",".join(sf["topology"])), step8_final_family=dash(",".join(sf["final"])),
                   outcome=stopped, explanation=expl)
    elif route == "indicator":
        row.update(step2_mmseqs2_query=o["query_ids"], step3_sequence_match="no match (best hit below 0.9 identity or in another genus)",
                   evalue=o["evalue"], query_coverage=o["query_coverage"],
                   outcome="not matched; nearest family is an indicator only",
                   nearest_family_indicator="%s (%s, %s)" % (o["outcome"].replace("nearest: ", ""), o["our_accession"], o["hit_species"]),
                   explanation="No protein of ours in the same species or genus at >= 0.9 identity over the longer sequence; "
                               "best significant hit %s (%s), E = %s, query coverage %s." % (
                                   o["our_accession"], o["hit_species"], o["evalue"], o["query_coverage"]))
    elif route == "no significant hit":
        row.update(step2_mmseqs2_query=o["query_ids"], step3_sequence_match="no significant hit",
                   outcome="no significant hit in our proteome set",
                   explanation="MMseqs2 found no hit at E <= 1e-5 with query coverage >= 0.5 in our proteome set.")
    else:
        row.update(step2_mmseqs2_query="no sequence", step3_sequence_match="not run", outcome="no sequence; identifier not in our data",
                   explanation="The study gives no usable sequence and its identifier matches no protein in our data.")
    rows.append(row)

keys = [k for k, _ in COLS]
with open(os.path.join(LO, "Supplementary_Table_literature_receptors.tsv"), "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, keys, delimiter="\t")
    w.writeheader()
    w.writerows(rows)

from openpyxl import Workbook
from openpyxl.styles import Font
wb = Workbook()
ws = wb.active
ws.title = "README"
for line in [
    ["Literature receptors traced through our data, one row per receptor listed by a study."],
    [],
    ["Order of the steps"],
    ["Step 1: a receptor is matched through a UniProt accession or an identifier the study gives (in our search outputs or the UniProt cross-references)."],
    ["Steps 2-3: every other receptor with a sequence is searched with MMseqs2 against all proteomes we searched. It is matched to our protein when a hit in the same species (or, failing that, the same genus) has >= 0.9 identical residues over the longer of the two sequences."],
    ["Steps 4-8: the matched protein is followed through our pipeline: family searches, family clade of the search tree, 7TM filter, topology filter (0.8), final family after overlap resolution."],
    ["Unmatched receptors carry the family of their best significant hit (E <= 1e-5, query coverage >= 0.5) as an indicator only; it is not counted as membership."],
    [],
    ["Column", "Meaning"],
] + [[k, d] for k, d in COLS]:
    ws.append(line)
for c in (ws["A1"], ws["A3"], ws["A9"], ws["B9"]):
    c.font = Font(bold=True)
ws.column_dimensions["A"].width = 26
ws.column_dimensions["B"].width = 120
t = wb.create_sheet("Receptors")
t.append(keys)
for c in t[1]:
    c.font = Font(bold=True)
for r in rows:
    t.append([r[k] for k in keys])
t.freeze_panes = "A2"
t.auto_filter.ref = t.dimensions
wb.save(os.path.join(LO, "Supplementary_Table_literature_receptors.xlsx"))

from collections import Counter
print("rows", len(rows))
for k, v in Counter(r["outcome"].split(" (")[0] for r in rows).most_common():
    print("%7d  %s" % (v, k))
