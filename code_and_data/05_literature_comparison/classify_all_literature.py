"""Every literature receptor, one outcome, with the same criteria throughout.
  ID match       : matched through a UniProt accession or an identifier the study gives (compare_literature.py)
  otherwise, from the MMseqs2 hits of its sequence against our whole proteome set
  (unmatched_ LITU, remaining_ LITX and unsearched_ LITS query sets):
  species match  : a hit in the same species with >= 0.9 identical residues over the longer sequence
  genus match    : no species match, a hit in the same genus by the same identity rule
  indicator      : best hit with E <= 1e-5 and query coverage >= 0.5; its family is an indicator, not membership
  no significant hit / no sequence
Writes 05_literature_comparison/literature_all_outcomes.tsv and prints the summary."""
import csv, os, re, sys
from collections import Counter, defaultdict

STEP6 = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, STEP6)
sys.path.insert(0, os.path.dirname(STEP6))
from gpcr_paths import open_text, stage  # noqa: E402
LO = str(stage("05_literature_comparison"))
import literature_readers as L
csv.field_size_limit(10**9)
MIN_ID, MAX_E, MIN_QCOV = 0.9, 1e-5, 0.5


def wide(path):
    rows = list(csv.reader(open(path, encoding="utf-8"), delimiter="\t"))
    return {acc: fam for row in rows[1:] for fam, acc in zip(rows[0], row) if acc}

final = wide(str(stage("04_overlap_assignment", "final_membership.tsv")))
clade = wide(str(stage("02_clade_extraction", "clade_membership.tsv")))
status = lambda a: ("SLT_TRNS" if final[a] == "SLT_TRNS" else "final") if a in final else ("clade_filtered" if a in clade else "no_clade")
fam = lambda a: final.get(a) or clade.get(a) or ""
sp2 = lambda s: " ".join(s.replace(".", "").lower().split()[:2])
aa = lambda s: re.sub(r"[^A-Z]", "", s.upper())


def fasta(path):
    out, h = {}, None
    for line in open_text(path):
        if line.startswith(">"):
            h = line[1:].split("|")[0].strip(); out[h] = []
        else:
            out[h].append(line.strip())
    return {k: "".join(v) for k, v in out.items()}

# sequence -> query ids (all three query sets)
SETS = ["unmatched_literature_sequences", "remaining_literature_sequences", "unsearched_literature_sequences"]
seq2q = defaultdict(set)
for s in SETS:
    for q, seq in fasta(os.path.join(LO, s + ".fasta")).items():
        seq2q[seq].add(q)
# identifier-only receptors whose sequence was fetched (LITX): (study, family, id) -> query
fetched = {(r["study"], r["their_family"], r["id"]): r["query_id"] for r in
           csv.DictReader(open(os.path.join(LO, "remaining_literature_sequences.tsv"), encoding="utf-8"), delimiter="\t")}
rseq = {(r["study"], r["their_family"], r["seq_id"]): aa(r["sequence"]) for r in L.read_all()
        if r["kind"] == "sequence" and r["sequence"]}

# hits: per query the best significant hit and every hit at >= 0.9 identity over the longer sequence
best, close = {}, defaultdict(list)
searched = set()
for s in SETS:
    for line in open_text(os.path.join(LO, s + ".m8")):
        f = line.rstrip("\n").split("\t")
        q = f[0].split("|")[0]
        searched.add(q)
        fid, aln, ql, tl, bits, qcov, ev = float(f[4]), int(f[6]), int(f[7]), int(f[8]), float(f[11]), float(f[12]), float(f[14])
        sp = (re.search(r"OS=(.*?) OX=", f[2]) or [None, ""])[1]
        if ev <= MAX_E and qcov >= MIN_QCOV and (q not in best or bits > best[q][0]):
            best[q] = (bits, f[1], sp, ev, qcov, fid, q)
        if fid >= MIN_ID:
            idl = round(fid * aln) / max(ql, tl)
            if idl >= MIN_ID:
                close[q].append((bits, f[1], sp, idl, ev, qcov, q))
    print("read", s, flush=True)


def outcome(target, ours):
    st = status(target)
    if st == "final":
        return "final in claimed family" if fam(target) in ours else "final in another family"
    if st == "clade_filtered":
        return "in claimed family's clade, removed by filters" if fam(target) in ours else "in another family's clade, removed by filters"
    return "SLT_TRNS (background)" if st == "SLT_TRNS" else "not in our clades"


rows = []
for r in csv.DictReader(open(os.path.join(LO, "literature_sequences.tsv"), encoding="utf-8"), delimiter="\t"):
    key = (r["study"], r["their_family"], r["id"])
    ours = [x for x in r["our_families"].split(",") if x and not x.startswith("none") and x != "any"]
    if r["study"].startswith(("Behrens", "Itoigawa")) and r["their_family"] == "T2R":
        ours = ["classT"]
    base = [r["study"], r["their_family"], ",".join(ours), r["id"], r["genome"]]
    if r["match"].startswith(("accession", "identifier")):
        rows.append(base + ["ID match", r["stopped_at"], r["our_accession"], "", "", "", "", "", ""])
        continue
    qs = set(seq2q.get(rseq.get(key, ""), ())) | ({fetched[key]} if key in fetched else set())
    if not qs:
        rows.append(base + ["no sequence", "", "", "", "", "", "", "", ""]); continue
    hits = [h for q in qs for h in close.get(q, [])]
    g = sp2(r["genome"])
    same_sp = [h for h in hits if g and sp2(h[2]) == g]
    same_gen = [h for h in hits if g and sp2(h[2]).split(" ")[0] == g.split(" ")[0]]
    if same_sp or same_gen:
        h = max(same_sp or same_gen)
        rows.append(base + ["species match" if same_sp else "genus match", outcome(h[1], ours), h[1], h[2], "%.3f" % h[3],
                            "%.1e" % h[4], "%.2f" % h[5], h[6], ",".join(sorted(qs))])
        continue
    b = max((best[q] for q in qs if q in best), default=None)
    if b:
        st = status(b[1])
        o = ("nearest: claimed family" if fam(b[1]) in ours else "nearest: another family") if st in ("final", "clade_filtered") \
            else ("nearest: SLT_TRNS" if st == "SLT_TRNS" else "nearest: no family of ours")
        rows.append(base + ["indicator", o, b[1], b[2], "", "%.1e" % b[3], "%.2f" % b[4], b[6], ",".join(sorted(qs))])
    else:  # a query absent from the m8 had no significant hit
        rows.append(base + ["no significant hit", "", "", "", "", "", "", "", ",".join(sorted(qs))])

with open(os.path.join(LO, "literature_all_outcomes.tsv"), "w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["study", "their_family", "our_families", "id", "genome", "route", "outcome", "our_accession", "hit_species", "identity_over_longer",
                "evalue", "query_coverage", "hit_query_id", "query_ids"])
    w.writerows(rows)

print("receptors", len(rows))
for k, v in Counter(x[5] for x in rows).most_common():
    print("%-20s %6d" % (k, v))
print()
for route in ("ID match", "species match", "genus match", "indicator"):
    print(route)
    for k, v in Counter(x[6] for x in rows if x[5] == route).most_common():
        print("   %6d  %s" % (v, k))
print()
by = defaultdict(Counter)
for x in rows:
    by[x[0]][x[5]] += 1
for s, c in sorted(by.items(), key=lambda kv: -sum(kv[1].values())):
    print("%-38s %6d  %s" % (s[:38], sum(c.values()), dict(c)))
