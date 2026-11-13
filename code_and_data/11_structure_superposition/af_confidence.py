"""AlphaFold confidence for the best motif-carrying candidates of each family.

Takes the members that carry the most family motifs (model_candidates_<family>.tsv; Swiss-Prot
entries first, at most MAX_QUERY per family) plus the current session models, looks up each in
AlphaFold DB, and for the TOP by global pLDDT downloads the model to read pLDDT at every motif
residue and to confirm the model carries the residue the alignment shows there.

Residue numbers are found by aligning the member's alignment row to the model sequence, not by
counting from the header range: the trimmed family alignments drop columns, so counting
drifts after the first removed insertion.

Writes af_candidates.tsv.
"""
import csv
import json
import urllib.request

from Bio.Align import PairwiseAligner, substitution_matrices

from final_motif_list import KINGDOM_ROWS, POSITIONS
from model_candidates import FAMILIES
from motif_columns import HERE, accession, family_motif_columns, read_fasta, trimmed_alignment

MAX_QUERY, TOP = 40, 3
CACHE = HERE / "af_models"


def af_entry(acc):
    try:
        with urllib.request.urlopen(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}",
                                    timeout=30) as r:
            return json.load(r)[0]
    except Exception:
        return None


def plddt_by_residue(entry):
    """Residue number -> (one-letter residue, pLDDT) from the model's mmCIF atom records."""
    CACHE.mkdir(exist_ok=True)
    path = CACHE / entry["cifUrl"].rsplit("/", 1)[-1]
    if not path.exists():
        urllib.request.urlretrieve(entry["cifUrl"], path)
    three = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
             "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
             "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
    cols, out = [], {}
    for line in open(path):
        if line.startswith("_atom_site."):
            cols.append(line.strip().split(".", 1)[1])
        elif line.startswith("ATOM") and cols:
            f = dict(zip(cols, line.split()))
            if f["label_atom_id"] == "CA":
                out[int(f["label_seq_id"])] = (three.get(f["label_comp_id"], "X"),
                                               float(f["B_iso_or_equiv"]))
    return out


def column_to_residue(row, model):
    """Alignment column -> model residue number, by aligning the row's residues to the model."""
    cols = [i for i, a in enumerate(row) if a != "-"]
    query = "".join(row[i].upper() for i in cols)
    seq_nums = sorted(model)
    target = "".join(model[n][0] for n in seq_nums)
    aligner = PairwiseAligner(mode="global", open_gap_score=-10, extend_gap_score=-0.5,
                              end_gap_score=0)
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aln = aligner.align(target, query)[0]
    out = {}
    for (t0, t1), (q0, q1) in zip(*aln.aligned):
        for k in range(t1 - t0):
            out[cols[q0 + k]] = seq_nums[t0 + k]
    return out


def main():
    colnum = {label: col for label, col, *_ in POSITIONS}
    out = open(HERE / "af_candidates.tsv", "w")
    out.write("family\taccession\tin_session\tmotifs_carried\tof\tmissing\tglobal_pLDDT\t"
              "motif_mean_pLDDT\tmotif_min_pLDDT\tresidue_check\tdescription\tmodel_residues\n")
    for fam in FAMILIES:
        rows = list(csv.DictReader(open(HERE / f"model_candidates_{fam}.tsv"), delimiter="\t"))
        motifs = list(rows[0].keys())[7:]
        best = int(rows[0]["motifs_carried"])
        pool = [r for r in rows if int(r["motifs_carried"]) == best]
        pool.sort(key=lambda r: "Uncharacterized" in r["description"])
        pool = pool[:MAX_QUERY] + [r for r in rows if r["in_session"] and r not in pool[:MAX_QUERY]]
        scored = []
        for r in pool:
            e = af_entry(r["accession"])
            if e:
                scored.append((e.get("globalMetricValue") or 0, r, e))
        scored.sort(key=lambda x: -x[0])
        keep = scored[:TOP] + [s for s in scored[TOP:] if s[1]["in_session"]]
        base = KINGDOM_ROWS.get(fam, fam)  # kingdom rows read their family's alignment
        fcols = family_motif_columns(base, None, None, {m: colnum[m] for m in motifs})
        wanted = {r["accession"] for _, r, _ in keep}
        rows_by_acc = {accession(h): s for h, s in
                       read_fasta(trimmed_alignment(base))
                       if accession(h) in wanted}
        print(f"\n{fam}: {len(pool)} queried, {len(scored)} in AlphaFold DB")
        for g, r, e in keep:
            res = plddt_by_residue(e)
            to_res = column_to_residue(rows_by_acc[r["accession"]], res)
            vals, bad, cells = [], [], []
            for m in motifs:
                if r[m] == "-":
                    continue
                aa, num = r[m][0], to_res.get(fcols[m])
                got = res.get(num) if num else None
                if got is None or got[0] != aa:
                    bad.append(f"{m}:{aa}->{got[0] if got else 'none'}")
                else:
                    vals.append(got[1])
                    cells.append(f"{m}={aa}{num}")
            check = "ok" if not bad else "MISMATCH " + ",".join(bad)
            mean = sum(vals) / len(vals) if vals else 0
            low = min(vals) if vals else 0
            out.write(f"{fam}\t{r['accession']}\t{r['in_session']}\t{r['motifs_carried']}\t"
                      f"{r['of']}\t{r['missing']}\t{g:.1f}\t{mean:.1f}\t{low:.1f}\t{check}\t"
                      f"{r['description']}\t{' '.join(cells)}\n")
            print(f"  {r['accession']:12s} {'SESSION ' if r['in_session'] else '        '}"
                  f"{r['motifs_carried']}/{r['of']} global {g:.1f} motif mean {mean:.1f} "
                  f"min {low:.1f} {check} | {r['description'][:60]}")
    out.close()


if __name__ == "__main__":
    main()
