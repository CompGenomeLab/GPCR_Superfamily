"""The literature sequences tested for deeper origins, as one FASTA for DeepTMHMM and MMseqs2.

Takes the sequences literature_deeper_origins.py tests (its receptors and translated TAS2R loci),
drops repeats of the same genome, identifier and sequence, and writes each with a numbered header
LIT<n>|study|claimed=<their family>|genome|from=<source>. Run after literature_deeper_origins.py;
the output is the input of deeptmhmm_literature.sh and of the fourth task of literature_mmseqs2.sh.
The archive keeps the file that was used, whose headers also end in best_hit=<family>, the closest
family by an earlier local-alignment comparison that deeper_origin_best_hits.py replaces. Each sequence
keeps its LIT number from that file (matched by genome, identifier and sequence), so the numbers stay
those of the MMseqs2 hits, the 7TM gate and the trees even when study labels change.

Recovered from the session of 3 October 2026 in which the file was written (inline code, saved here
unchanged apart from the output path and the best_hit label).

Reads  05_literature_comparison/literature_comparison.tsv
Writes 05_literature_comparison/deeper_origin_tested_sequences.fasta

Usage:
    python write_deeper_origin_fasta.py [--out <fasta>]
"""
import argparse
import csv
import os
import re

import compare_literature as C
import literature_deeper_origins as D
import literature_readers as L

O = D.OUT
ARCHIVED = os.path.join(O, "deeper_origin_tested_sequences.fasta")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=ARCHIVED)
    args = ap.parse_args()

    records = L.read_all()
    names = {row["genome_as_reported"] for row in csv.DictReader(
        open(os.path.join(O, "literature_comparison.tsv"), encoding="utf-8"), delimiter="\t")}
    m = C.GenomeMatcher({r["genome"] for r in records if r["genome"]} | names, set())
    cands = D.candidates(m, D.our_origins(m))
    D.attach_sequences(cands)
    queries = D.sort_receptors(cands) + D.tas2r_queries(cands)

    def clean(v):
        return re.sub(r"[|\s]+", "_", v.strip())[:60]

    # the LIT numbers of the archived file (the one the trees used) are kept, matched by genome,
    # identifier and sequence; a sequence not in it is numbered after the last
    used, head = {}, None
    if os.path.exists(ARCHIVED):
        for line in open(ARCHIVED, encoding="utf-8"):
            if line.startswith(">"):
                head = line[1:].rstrip("\n")
            else:
                f = head.split(" ", 1)
                used[(f[0].split("|")[3], f[1] if len(f) > 1 else "", line.strip())] = int(f[0].split("|")[0][3:])
    last = max(used.values(), default=0)

    seen, out = set(), []
    for q in queries:
        key = (q["genome"], q["id"], q["sequence"])
        if key in seen:
            continue
        seen.add(key)
        seq = re.sub("[^A-Z]", "", q["sequence"].upper())
        n = used.get((clean(q["genome"]), q["id"][:80], seq))
        if n is None:
            last += 1
            n = last
        head = "LIT%04d|%s|claimed=%s|%s|from=%s" % (
            n, clean(q["study"]), q["their_family"].replace(" ", "_"), clean(q["genome"]), clean(q["source"]))
        out.append((n, ">%s %s\n%s\n" % (head, q["id"][:80], seq)))
    out = [rec for _, rec in sorted(out)]
    open(args.out, "w", encoding="utf-8").write("".join(out))
    print(args.out, len(out), "sequences")


if __name__ == "__main__":
    main()
