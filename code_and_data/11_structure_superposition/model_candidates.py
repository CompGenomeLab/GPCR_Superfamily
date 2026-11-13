"""Rank family members as structure candidates for the motif figure.

For each family without an experimental structure, the motifs it must show are the final list
positions where the family carries the motif (conserved or consensus, in
final_motif_list_per_family.tsv). Every member is scored by how many of those it carries
(2x46 must be L; positions without a fixed class accept the family's own similarity group). The
models already in the Figure 4 session are scored the same way. Members with an experimental
structure (SIFTS) are flagged. Kingdom rows (cAMP_Plants, classC_Plants) rank only that kingdom's
members of the family, against the kingdom row's motifs.

Writes model_candidates_<family>.tsv, sorted best first.

Usage: python model_candidates.py <sifts uniprot_pdb.tsv.gz>
"""
import csv
import gzip
import json
import sys

from build_superposition import STRUCTURES
from final_motif_list import KINGDOM_ROWS, POSITIONS, family_accept
from motif_columns import HERE, MAPS, stage, accession, family_motif_columns, read_fasta, trimmed_alignment

FAMILIES = ["Vomeronasal2", "cAMP", "Mth", "FSLB", "GP143", "GP157", "GPR1", "STE3", "TM116",
            "Vomeronasal1", "Nematode", "STE2", "cAMP_Plants", "classC_Plants"]
KINGDOM = {"cAMP_Plants": "Viridiplantae", "classC_Plants": "Viridiplantae"}
ACC_TAXID = stage("10_motif_analysis", "phylogeny_reassignment", "acc_taxid_map_reassigned.tsv")
TAXONOMY = stage("07_representative_selection", "taxonomy_eukaryotes_filtered.json")
# AlphaFold models in the current session, by family
SESSION_MODELS = {}
for _, fam, _, _, acc, src in STRUCTURES:
    if src == "AF":
        SESSION_MODELS.setdefault(fam, []).append(acc)


def family_motifs(family):
    rows = csv.DictReader(open(HERE / "final_motif_list_per_family.tsv"), delimiter="\t")
    return [r["position"] for r in rows
            if r["family"] == family and r["status"] in ("conserved", "consensus")]


def kingdom_members(family, kingdom):
    kingdom_of = {r["taxID"]: r["Kingdom"] for r in json.load(open(TAXONOMY))}
    return {r["accession"] for r in csv.DictReader(open(ACC_TAXID), delimiter="\t")
            if r["family"] == family and kingdom_of.get(int(r["taxID"])) == kingdom}


def main():
    pdb = {}
    with gzip.open(sys.argv[1], "rt") as f:
        for line in f:
            if not line.startswith("#") and not line.startswith("SP_PRIMARY"):
                acc, ids = line.rstrip("\n").split("\t")
                pdb[acc] = ids
    cols = {label: col for label, col, *_ in POSITIONS}
    accepted = {label: acc for label, _, _, acc, _ in POSITIONS}
    for key in FAMILIES:
        fam = KINGDOM_ROWS.get(key, key)
        members = kingdom_members(fam, KINGDOM[key]) if key in KINGDOM else None
        site = json.load(open(MAPS / f"{key}.json"))["positions"]
        want = family_motifs(key)
        accept = {m: "L" if m == "2x46" else family_accept(accepted[m], site[cols[m] - 1])
                  for m in want}
        fcols = family_motif_columns(fam, None, None, {m: cols[m] for m in want})
        rows = []
        # residue letters only: counting from the header range drifts in trimmed alignments;
        # af_confidence.py maps real residue numbers against the model sequence
        for h, s in read_fasta(trimmed_alignment(fam)):
            acc = accession(h)
            if members is not None and acc not in members:
                continue
            cells, hits, missing = [], 0, []
            for m in want:
                c = fcols[m]
                aa = s[c].upper() if c is not None and s[c] != "-" else "-"
                cells.append(aa)
                if aa in accept[m]:
                    hits += 1
                else:
                    missing.append(f"{m}:{aa}")
            rows.append((hits, acc, h.split(" ", 1)[-1][:90], pdb.get(acc, ""),
                         "session" if acc in SESSION_MODELS.get(key, []) else "",
                         ",".join(missing), cells))
        rows.sort(key=lambda r: (-r[0], r[1]))
        with open(HERE / f"model_candidates_{key}.tsv", "w") as out:
            out.write("motifs_carried\tof\taccession\tdescription\tpdb\tin_session\tmissing\t"
                      + "\t".join(want) + "\n")
            for r in rows:
                out.write(f"{r[0]}\t{len(want)}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\t"
                          + "\t".join(r[6]) + "\n")
        full = sum(r[0] == len(want) for r in rows)
        print(f"\n{key}: {len(want)} motifs, {len(rows)} members, {full} carry all")
        for r in rows:
            if r[4]:
                print(f"  session model {r[1]}: {r[0]}/{len(want)}  missing {r[5] or '-'}")
        for r in [r for r in rows if r[3]][:5]:
            print(f"  PDB member {r[1]} ({r[3][:40]}): {r[0]}/{len(want)}  missing {r[5] or '-'}")


if __name__ == "__main__":
    main()
