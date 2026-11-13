"""The closest family of ours for each literature sequence tested for a deeper origin: its MMseqs2 best hit
among our family members.

The 365 sequences of deeper_origin_tested_sequences.fasta were searched with MMseqs2 against the searched
proteomes (literature_mmseqs2.sh, fourth task). As for the other literature receptors
(classify_all_literature.py), hits count when E <= 1e-5 and query coverage >= 0.5. The closest family is the
family of the highest-scoring such hit that is a final member of one of our families; a sequence with no
such hit has no hit.

Reads  05_literature_comparison/deeper_origin_tested_sequences.fasta, deeper_origin_tested_sequences.m8
       04_overlap_assignment/final_membership.tsv
Writes 05_literature_comparison/literature_closest_family.tsv

Usage:
    python deeper_origin_best_hits.py
"""
import csv
import os
import sys

STEP = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(STEP))
from gpcr_paths import open_text, stage  # noqa: E402

MAX_E, MIN_QCOV = 1e-5, 0.5   # as in classify_all_literature.py
LO = stage("05_literature_comparison")


def wide(path):
    rows = list(csv.reader(open(path, encoding="utf-8"), delimiter="\t"))
    return {acc: fam for row in rows[1:] for fam, acc in zip(rows[0], row) if acc}


def main():
    final = wide(stage("04_overlap_assignment", "final_membership.tsv"))

    queries = []   # (header, sequence length)
    head, n = None, 0
    for line in open_text(LO / "deeper_origin_tested_sequences.fasta"):
        if line.startswith(">"):
            if head:
                queries.append((head, n))
            head, n = line[1:].strip(), 0
        else:
            n += len(line.strip())
    queries.append((head, n))

    best = {}
    for line in open_text(LO / "deeper_origin_tested_sequences.m8"):
        f = line.rstrip("\n").split("\t")
        q, bits, qcov, ev = f[0].split("|")[0], float(f[11]), float(f[12]), float(f[14])
        if f[1] in final and ev <= MAX_E and qcov >= MIN_QCOV and (q not in best or bits > best[q][0]):
            best[q] = (bits, f[1], float(f[4]), qcov, ev)

    cols = ["id", "study", "claimed", "genome", "sequence_from", "length", "closest_family", "best_hit",
            "bits", "identity", "query_coverage", "evalue"]
    with open(LO / "literature_closest_family.tsv", "w", newline="", encoding="utf-8") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(cols)
        for head, length in queries:
            parts = head.split(" ")[0].split("|")
            lit, study, claimed, genome, source = parts[0], parts[1], parts[2][len("claimed="):], parts[3], parts[4]
            row = [lit, study.replace("_", " "), claimed, genome.replace("_", " "), source[len("from="):].replace("_", " "), length]
            if lit not in best:
                w.writerow(row + ["no hit", "", "", "", "", ""])
                continue
            bits, acc, fid, qcov, ev = best[lit]
            w.writerow(row + [final[acc], acc, "%g" % bits, "%.3f" % fid, "%.3f" % qcov, "%.2e" % ev])
    print("wrote literature_closest_family.tsv:", len(queries), "sequences,", len(best), "with a hit to a family member")


if __name__ == "__main__":
    main()
