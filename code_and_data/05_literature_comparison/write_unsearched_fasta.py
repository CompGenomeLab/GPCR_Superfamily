"""Study sequences not yet searched with MMseqs2 and not matched through an identifier the study gives
(UniProt accession, or an identifier in our search headers): each distinct sequence once.
Earlier query sets: unmatched_literature_sequences.fasta, remaining_literature_sequences.fasta.
Writes 05_literature_comparison/unsearched_literature_sequences.fasta/.tsv (query ids LITS#####)."""
import csv, os, re, sys
from collections import Counter

STEP6 = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, STEP6)
import literature_readers as L

sys.path.insert(0, os.path.dirname(STEP6))
from gpcr_paths import stage  # noqa: E402
OUT = str(stage("05_literature_comparison"))
csv.field_size_limit(10**9)
aa = lambda s: re.sub(r"[^A-Z]", "", s.upper())
searched = set()
for name in ("unmatched_literature_sequences.fasta", "remaining_literature_sequences.fasta"):
    seq = []
    for line in open(os.path.join(OUT, name), encoding="utf-8"):
        if line.startswith(">"):
            if seq: searched.add("".join(seq))
            seq = []
        else:
            seq.append(line.strip())
    if seq: searched.add("".join(seq))
match = {(r["study"], r["their_family"], r["id"]): r["match"] for r in
         csv.DictReader(open(os.path.join(OUT, "literature_sequences.tsv"), encoding="utf-8"), delimiter="\t")}
by_id = lambda m: m.startswith(("accession", "identifier"))
clean = lambda s: re.sub(r"[^A-Za-z0-9.=_-]+", "_", s).strip("_")
n, by, written, skipped_id = 0, Counter(), {}, 0
with open(os.path.join(OUT, "unsearched_literature_sequences.fasta"), "w", encoding="utf-8") as fa, \
        open(os.path.join(OUT, "unsearched_literature_sequences.tsv"), "w", encoding="utf-8", newline="") as tv:
    w = csv.writer(tv, delimiter="\t")
    w.writerow(["query_id", "study", "their_family", "id", "genome", "current_match", "length"])
    for r in L.read_all():
        if r["kind"] != "sequence" or not r["sequence"]:
            continue
        s, m = aa(r["sequence"]), match.get((r["study"], r["their_family"], r["seq_id"]), "")
        if by_id(m):
            skipped_id += 1
            continue
        if len(s) < 30 or s in searched:
            continue
        if s not in written:
            n += 1
            written[s] = "LITS%05d" % n
            fa.write(">%s|%s|claimed=%s|%s\n%s\n" % (written[s], clean(r["study"]), clean(r["their_family"]),
                                                     clean(r["genome"] or "no_species"), s))
            by[r["study"]] += 1
        w.writerow([written[s], r["study"], r["their_family"], r["seq_id"], r["genome"], m, len(s)])
print("matched through an identifier, left out:", skipped_id)
print("distinct sequences written:", n)
for k, v in by.most_common(): print("  %6d  %s" % (v, k))
