"""The 18 literature sequences added to the superfamily trees, trimmed to their 7TM core.

They test the literature presences that would place the origin of a family deeper than our data
(literature_deeper_origins.py). Eligible were sequences that pass the 7TM gate (literature_7tm_gate.py)
and come from genomes outside our dataset. Within each lineage (layer) outside the family's origin, the
candidates were ranked by the bit score of their best MMseqs2 hit to a final member of the claimed family
(deeper_origin_tested_sequences.m8), preferring those whose best hits are members of that family rather
than cAMP, and up to two were taken per lineage:

    GP157   Amoebozoa (Arcella x2), Apusozoa (1), CRuMs (Rigifila, 1), Discoba (Entosiphon, Euglena),
            SAR (Mataza, Nonionella), Haptista (Raphidiophrys x2), Metamonada (Anaeramoeba, Trimastix)
    Taste-2 the two cartilaginous fish loci with the highest scores (Scyliorhinus, Rhincodon)
    Vomeronasal 1  two lamprey receptors (Kowatschew & Korsching 2022)
    Acanthamoeba   the two sequences reported as Rhodopsin-like, closest to GPR1 or equally GPR1 and cAMP

The choice is recorded below (SELECTION), as made on 4 October 2026. Each sequence is cut to its 7TM core
(first to seventh helix of the single OWIMOWIMOWIMOWI pattern) with five flanking residues on either side,
as for our representatives (08_superfamily_alignment/extract_tm_region.py).

Reads  05_literature_comparison/deeper_origin_7TM_gate.tsv, predicted_topologies.3line
Writes 05_literature_comparison/literature_reps_7TM_flank5.fasta

Usage:
    python select_literature_representatives.py
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import open_text, stage  # noqa: E402

# (query id, candidate family label, lineage), in the order written
SELECTION = [("LIT0024", "GP157", "Amoebozoa"), ("LIT0025", "GP157", "Amoebozoa"), ("LIT0022", "GP157", "Apusozoa"),
             ("LIT0171", "GP157", "CRuMs"), ("LIT0063", "GP157", "Discoba"), ("LIT0081", "GP157", "Discoba"),
             ("LIT0114", "GP157", "Sar"), ("LIT0119", "GP157", "Sar"), ("LIT0161", "GP157", "Haptista"),
             ("LIT0162", "GP157", "Haptista"), ("LIT0015", "GP157", "Metamonada"), ("LIT0197", "GP157", "Metamonada"),
             ("LIT0361", "classT", "Chondrichthyes"), ("LIT0360", "classT", "Chondrichthyes")]
LAMPREY_V1R = ["LIT0234", "LIT0217"]
ACANTHAMOEBA = [("LIT0202", "L8GRC2"), ("LIT0203", "L8H6T9")]   # with the protein of ours each matches


def main():
    d = stage("05_literature_comparison")
    with open(d / "deeper_origin_7TM_gate.tsv", encoding="utf-8") as fh:
        gate = {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    with open_text(d / "predicted_topologies.3line") as fh:
        lines = [line.rstrip("\n") for line in fh]
    seq = {lines[i][1:].split("|")[0]: lines[i + 1] for i in range(0, len(lines), 3)}

    def core(q):
        g, s = gate[q], seq[q]
        st = max(1, int(g["core_start"]) - 5)
        en = min(len(s), int(g["core_stop"]) + 5)
        return g, st, en, s[st - 1:en].replace("*", "")

    with open(d / "literature_reps_7TM_flank5.fasta", "w") as out:
        for q, family, lineage in SELECTION:
            g, st, en, sub = core(q)
            out.write(f">LIT_{family}|{q}|{g['genome']}|claimed={g['claimed']}|layer={lineage}/{st}-{en}\n{sub}\n")
        for q in LAMPREY_V1R:
            g, st, en, sub = core(q)
            out.write(f">LIT_Vomeronasal1|{q}|{g['genome']}|claimed={g['claimed']}|study={g['study']}"
                      f"|layer=Cyclostomata/{st}-{en}\n{sub}\n")
        for q, ours in ACANTHAMOEBA:
            g, st, en, sub = core(q)
            out.write(f">LIT_GPR1like|{q}|{g['genome']}|claimed=Rhodopsin|matches_our={ours}"
                      f"|layer=Amoebozoa/{st}-{en}\n{sub}\n")
    print("wrote", len(SELECTION) + len(LAMPREY_V1R) + len(ACANTHAMOEBA), "sequences")


if __name__ == "__main__":
    main()
