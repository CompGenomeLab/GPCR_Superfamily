"""Class A or GPR1 states of the fungal "Class A" sequences at the 28 superfamily motif positions.

Positions: the 28 of Fig. 3 (make_fig3.py: >= 12 of 18 families conserved at 60% on the treein
G-INS-i --ep 0.123 mapping JSONs). Each family's alignment column at a position is the website's
(top-5 voted consensus mapping), read in consensus_trimmed_msas_Oct11/ (Oct09 membership).

Reference states are recomputed without the sequences being tested: Class A without any of the
337 fungal Class A of Sep24 (fungal_classA_members.tsv), GPR1 without the 136 moved Dikarya
(phylogeny_reassignments.tsv). A family's state is its top residue's similarity group
(D/E, F/W/Y/H, K/R, N/Q, I/V, L/M), conserved when the group reaches 60% of sequences, gaps
included (the manuscript rule). Each fungal sequence is read in the family alignment it sits in
now (the 136 in GPR1, the rest in Class A) at the same superfamily columns.

A diagnostic position is one where both references are conserved and their states differ.

With --all, the positions are instead every superfamily column where the Class A or the GPR1
reference is conserved (60%) and both families have a mapped column; each sequence is also scored
on the family-specific sites (conserved in one family and not shared with the other in the same
state): the share of Class A-specific sites carrying the Class A state against the share of
GPR1-specific sites carrying the GPR1 state.

Writes fungal_motif_states[_all60]_positions.tsv and fungal_motif_states[_all60]_sequences.tsv.
"""
import csv
import math
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE.parent))
import superfamily_motif_positions as smp  # noqa: E402  (also puts gpcr_paths on the path)
from gpcr_paths import open_text, stage  # noqa: E402

S10 = stage("10_motif_analysis")
HERE = S10 / "phylogeny_reassignment"
MAPPING = S10 / "consensus_mapping" / "representatives_all_Oct08_ginsi_ep0.123_treein_conf0.33_top5.tsv"
TRIMMED = S10 / "family_alignments" / "consensus_trimmed_msas"
MIN_FAMILIES = 12
ALL = "--all" in sys.argv
TAG = "_all60" if ALL else ""
AA = set("ACDEFGHIKLMNPQRSTVWY")
DIKARYA = {"Ascomycota", "Basidiomycota"}


def read_fasta(path):
    h, s = None, []
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if h:
                    yield h, "".join(s)
                h, s = line[1:], []
            else:
                s.append(line)
    if h:
        yield h, "".join(s)


def group(aa):
    return smp.GROUPS.get(aa, aa)


# the 28 positions, as make_fig3.py selects them
fams = smp.load_all()
chosen = []
for c in range(len(fams["classA"])):
    present = [fams[f][c] for f in smp.COUNTED if fams[f][c]]
    if sum(p["cons"] >= smp.THRESHOLD for p in present) >= MIN_FAMILIES:
        chosen.append(c)

consensus = defaultdict(dict)
for r in csv.DictReader(open(MAPPING), delimiter="\t"):
    if r["consensus_family_column"]:
        consensus[r["family"]][int(r["sup_column"]) - 1] = int(r["consensus_family_column"]) - 1
if ALL:   # every column both families map; the 60% filter on the references follows below
    chosen = sorted(set(consensus["classA"]) & set(consensus["GPR1"]))


def label(c):
    for f in ("classA", "GPR1"):
        p = fams[f][c] if c < len(fams[f]) else None
        if p and p.get("label"):
            return p["label"] + ("" if f == "classA" else "[GPR1]")
    return f"col{c + 1}"


labels = {c: label(c) for c in chosen}

members = {r["accession"]: r for r in csv.DictReader(open(HERE / "fungal_classA_members.tsv"), delimiter="\t")}
moved = {r["accession"] for r in csv.DictReader(open(CODE / "phylogeny_reassignments.tsv"), delimiter="\t")}

# residues at the 28 positions: reference counts per family, and every fungal sequence
ref = {}
fungal = {}
for fam, exclude in (("classA", set(members)), ("GPR1", moved)):
    kept = [int(x) - 1 for x in (TRIMMED / f"{fam}_kept_columns.txt").read_text().split()]
    full_to_trim = {f: t for t, f in enumerate(kept)}
    cols = {c: full_to_trim.get(consensus[fam].get(c)) for c in chosen}
    counts, n = {c: Counter() for c in chosen}, 0
    for h, s in read_fasta(TRIMMED / f"{fam}_genes_filtered_db_FAMSA.consensus_trimmed.fasta"):
        a = h.split("|")[1]
        res = {c: (s[t].upper() if t is not None else "-") for c, t in cols.items()}
        if a in members:
            fungal[a] = (fam, res)
        if a in exclude:
            continue
        n += 1
        for c, x in res.items():
            if x in AA:
                counts[c][group(x)] += 1
    ref[fam] = {}
    for c in chosen:
        if counts[c]:
            g, k = counts[c].most_common(1)[0]
            ref[fam][c] = (g, 100 * k / n)
        else:
            ref[fam][c] = (None, 0.0)
    print(f"{fam}: {n} reference sequences; unmapped positions: "
          f"{[labels[c] for c in chosen if cols[c] is None] or 'none'}")

missing = sorted(set(members) - set(fungal))
print(f"fungal Class A of Sep24: {len(members)}; found in the Oct09 alignments: {len(fungal)}"
      + (f"; missing {missing}" if missing else ""))

cons = lambda fam, c: ref[fam][c][1] >= smp.THRESHOLD
if ALL:
    n_mapped = len(chosen)
    chosen = [c for c in chosen if cons("classA", c) or cons("GPR1", c)]
    print(f"columns mapped in both families: {n_mapped}; conserved in Class A or GPR1: {len(chosen)}")
diag = [c for c in chosen if cons("classA", c) and cons("GPR1", c) and ref["classA"][c][0] != ref["GPR1"][c][0]]
shared = [c for c in chosen if cons("classA", c) and cons("GPR1", c) and ref["classA"][c][0] == ref["GPR1"][c][0]]
specific = {f: [c for c in chosen if cons(f, c) and c not in shared] for f in ("classA", "GPR1")}


def subgroup(a):
    m = members[a]
    if m["phylum"] in DIKARYA:
        return "Dikarya, moved to GPR1" if a in moved else "Dikarya, kept in Class A"
    return "early-diverging fungi"


SUBGROUPS = ["Dikarya, moved to GPR1", "Dikarya, kept in Class A", "early-diverging fungi"]
by_sub = defaultdict(list)
for a in fungal:
    by_sub[subgroup(a)].append(a)

with open(HERE / f"fungal_motif_states{TAG}_positions.tsv", "w", newline="") as out:
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    head = ["position", "type", "classA_state", "classA_pct", "GPR1_state", "GPR1_pct"]
    for s in SUBGROUPS:
        head += [f"{s} (n={len(by_sub[s])}): % Class A state", f"{s}: % GPR1 state"]
    w.writerow(head)
    for c in chosen:
        ga, pa = ref["classA"][c]
        gg, pg = ref["GPR1"][c]
        typ = ("diagnostic" if c in diag else "shared" if c in shared else
               "Class A only" if cons("classA", c) else "GPR1 only" if cons("GPR1", c) else "neither conserved")
        row = [labels[c], typ, ga, f"{pa:.0f}", gg, f"{pg:.0f}"]
        for s in SUBGROUPS:
            accs = by_sub[s]
            res = [fungal[a][1][c] for a in accs]
            row += [f"{100 * sum(ga is not None and x in AA and group(x) == ga for x in res) / len(accs):.0f}",
                    f"{100 * sum(gg is not None and x in AA and group(x) == gg for x in res) / len(accs):.0f}"]
        w.writerow(row)

with open(HERE / f"fungal_motif_states{TAG}_sequences.tsv", "w", newline="") as out:
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    w.writerow(["accession", "subgroup", "phylum", "class", "order", "species", "alignment",
                "n_diagnostic", "classA_state", "GPR1_state", "other_or_gap", "call",
                "n_classA_specific", "classA_specific_matched", "n_GPR1_specific", "GPR1_specific_matched",
                "pct_classA_specific", "pct_GPR1_specific", "call_specific",
                "residues_at_diagnostic (" + " ".join(labels[c] for c in diag) + ")"])
    for a in sorted(fungal, key=lambda a: (SUBGROUPS.index(subgroup(a)), members[a]["phylum"], members[a]["class"], a)):
        fam, res = fungal[a]
        na = sum(res[c] in AA and group(res[c]) == ref["classA"][c][0] for c in diag)
        ng = sum(res[c] in AA and group(res[c]) == ref["GPR1"][c][0] for c in diag)
        call = "Class A" if na > ng else "GPR1" if ng > na else "tie"
        sa = sum(res[c] in AA and group(res[c]) == ref["classA"][c][0] for c in specific["classA"])
        sg = sum(res[c] in AA and group(res[c]) == ref["GPR1"][c][0] for c in specific["GPR1"])
        pa, pg = 100 * sa / len(specific["classA"]), 100 * sg / len(specific["GPR1"])
        call_s = "Class A" if pa > pg else "GPR1" if pg > pa else "tie"
        m = members[a]
        w.writerow([a, subgroup(a), m["phylum"], m["class"], m["order"], m["species"], fam,
                    len(diag), na, ng, len(diag) - na - ng, call,
                    len(specific["classA"]), sa, len(specific["GPR1"]), sg, f"{pa:.0f}", f"{pg:.0f}", call_s,
                    "".join(res[c] for c in diag)])

print(f"\n{len(chosen)} positions; diagnostic ({len(diag)}): {', '.join(labels[c] for c in diag)}")
print(f"shared ({len(shared)}): {', '.join(labels[c] for c in shared)}")
print(f"Class A-specific sites: {len(specific['classA'])}; GPR1-specific sites: {len(specific['GPR1'])}")
