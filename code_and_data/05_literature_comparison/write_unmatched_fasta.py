"""Literature receptors with a sequence that were not matched to a protein of ours (no accession,
sequence or identifier match), for an MMseqs2 search against our families. Header style follows
deeper_origin_tested_sequences.fasta."""
import csv, os, re, sys
from collections import Counter

STEP6 = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, STEP6)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gpcr_paths import stage  # noqa: E402
import literature_readers as L

OUT = str(stage("05_literature_comparison"))
csv.field_size_limit(10**9)

seqs = {}
for r in L.read_all():
    if r["kind"] == "sequence" and r.get("sequence"):
        seqs.setdefault((r["study"], r["their_family"], r["seq_id"]), re.sub(r"[^A-Z]", "", r["sequence"].upper()))

rows = list(csv.DictReader(open(os.path.join(OUT, "literature_sequences.tsv"), encoding="utf-8"), delimiter="\t"))
clean = lambda s: re.sub(r"[^A-Za-z0-9.=_-]+", "_", s).strip("_")
n, kept, no_seq = 0, Counter(), Counter()
with open(os.path.join(OUT, "unmatched_literature_sequences.fasta"), "w", encoding="utf-8") as fa, \
        open(os.path.join(OUT, "unmatched_literature_sequences.tsv"), "w", encoding="utf-8", newline="") as tv:
    w = csv.writer(tv, delimiter="\t")
    w.writerow(["query_id", "study", "their_family", "our_families", "id", "genome", "stopped_at", "length"])
    for x in rows:
        if x["our_accession"]:
            continue  # matched to a protein of ours; its fate is already known
        seq = seqs.get((x["study"], x["their_family"], x["id"]))
        if not seq:
            no_seq[(x["study"], x["stopped_at"])] += 1
            continue
        n += 1
        qid = "LITU%05d" % n
        fa.write(">%s|%s|claimed=%s|%s|stopped=%s\n%s\n" % (
            qid, clean(x["study"]), clean(x["their_family"]), clean(x["genome"] or "no_species"),
            clean(x["stopped_at"]), seq))
        w.writerow([qid, x["study"], x["their_family"], x["our_families"], x["id"], x["genome"], x["stopped_at"], len(seq)])
        kept[(x["study"], "counterpart" if x["our_families"] else "no counterpart", x["stopped_at"])] += 1

print("written:", n)
for k, v in sorted(kept.items()):
    print("  %s\t%s\t%s\t%d" % (k + (v,)))
print("unmatched but no sequence:", sum(no_seq.values()))
for k, v in sorted(no_seq.items()):
    print("  %s\t%s\t%d" % (k + (v,)))
