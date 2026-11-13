"""Every literature receptor sequence searched with MMseqs2 (the receptors not matched by ID), as one search.
Joins the three query sets classify_all_literature.py reads (unmatched_ LITU, remaining_ LITX, unsearched_ LITS)
and keeps each sequence once: a sequence listed by several receptors, or searched in two sets, keeps its first
query. Queries with no significant hit stay in the FASTA although they do not appear in the .m8.
Writes in 05_literature_comparison/:
  literature_mmseqs_queries_all.fasta   one query per unique sequence
  literature_mmseqs_all.m8              the hits of those queries, from the three runs
  literature_mmseqs_query_map.tsv       every receptor query -> the query its sequence was searched under"""
import os, shutil, sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gpcr_paths import open_text, stage  # noqa: E402
LO = str(stage("05_literature_comparison"))
SETS = ["unmatched_literature_sequences", "remaining_literature_sequences", "unsearched_literature_sequences"]


def fasta(path):
    out, h = {}, None
    for line in open(path, encoding="utf-8"):
        if line.startswith(">"):
            h = line[1:].strip(); out[h] = []
        else:
            out[h].append(line.strip())
    return {k: "".join(v) for k, v in out.items()}

kept, by_seq, qmap = [], {}, []
for s in SETS:
    for q, seq in fasta(os.path.join(LO, s + ".fasta")).items():
        if seq not in by_seq:
            by_seq[seq] = q
            kept.append((q, seq))
        qmap.append((q, by_seq[seq]))

with open(os.path.join(LO, "literature_mmseqs_queries_all.fasta"), "w", encoding="utf-8", newline="\n") as fo:
    for q, seq in kept:
        fo.write(f">{q}\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
with open(os.path.join(LO, "literature_mmseqs_query_map.tsv"), "w", encoding="utf-8", newline="\n") as fo:
    fo.write("receptor_query\tsearched_as\n")
    fo.writelines(f"{a}\t{b}\n" for a, b in qmap)

keep = {q for q, _ in kept}
hit_q, n = set(), 0
with open(os.path.join(LO, "literature_mmseqs_all.m8"), "w", encoding="utf-8", newline="\n") as out:
    for s in SETS:
        for line in open_text(os.path.join(LO, s + ".m8")):
            q = line.split("\t", 1)[0]
            if q in keep:
                out.write(line); hit_q.add(q); n += 1
print("receptor queries", len(qmap), "unique sequences", len(kept), "with hits", len(hit_q),
      "without hits", len(keep - hit_q), "m8 lines", n)
