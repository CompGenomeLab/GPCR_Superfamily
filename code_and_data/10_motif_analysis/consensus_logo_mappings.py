"""Superfamily logo mapping JSONs built on the family column assigned to every superfamily
column by consensus_column_mapping.py (one-to-one, confidence at or above the threshold)
instead of the anchor's column.

Residue counts come from consensus_trimmed_msas/ (trim_family_alignments_consensus.py on the
cluster): every family alignment column where the GPCRdb numbering reference has a residue or
that the 0.33 mapping assigns in any alignment, with <family>_kept_columns.txt giving each kept
column's full alignment column. A consensus column that is not among them is left empty and
counted.

Each JSON is the website's (same keys, read from its folder) with the positions rebuilt: residue
counts, information content, the GPCRdb number and reference residue of acc2 (acc1 without a
reference) at the consensus column, and `confidence`, `votes`, `nReps`. Kingdom rows use their
family's mapping and count only the family's sequences of that kingdom
(family_column_counts_kingdom_groups_Oct09.tsv, as the website's kingdom rows).

python consensus_logo_mappings.py [threshold ...]   (default 0.33)
python consensus_logo_mappings.py --treein [--top5] [threshold ...]   the tree-guided alignments (TREEIN), with
    consensus_trimmed_msas_Oct11/ (every column any representative fills, consensus_trim_columns_allreps_Oct11.tsv)

Writes consensus_mapping/logo_conf<threshold>/<website mapping folder>/<family>.json and index.json.
"""
import csv
import itertools
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

from consensus_column_mapping import ALIGNMENTS, DATA, HERE, OUT, PUBLIC, SOURCE, TREEIN, read_fasta, residues

TREEIN_RUN = "--treein" in sys.argv
TAG = "_top5" if "--top5" in sys.argv else ""  # mapping voted by the top 5 representatives (consensus_mapping_top_reps.py)
THRESHOLDS = [float(x) for x in sys.argv[1:] if not x.startswith("--")] or [0.33]
if TREEIN_RUN:
    ALIGNMENTS = TREEIN
TRIMMED = DATA / "family_alignments" / "consensus_trimmed_msas"
# numbering reference (acc2, or acc1 without one) per family, as the alignments were trimmed
REFERENCE = {r["family"]: r["reference"] for r in csv.DictReader(
    open(DATA / "family_alignments" / "consensus_trim_columns_allreps.tsv"),
    delimiter="	")}
ANCHOR = {l.split("	")[1].strip(): l.split("	")[0].strip()
          for l in (SOURCE / "trim_info.tsv").read_text().splitlines()[1:] if l.strip()}
# Class C and cAMP kingdom rows, as the website's (Oct09 membership)
KINGDOM = {(f, a): g for a, f, g in (l.split("	") for l in
           (DATA / "family_alignments" / "family_column_counts_kingdom_groups.tsv").read_text().splitlines() if l)}
AA = "ACDEFGHIKLMNPQRSTVWY"
CONSERVATION = PUBLIC / "conservation_files"


def gpcrdb_numbers(header):
    parts = header.split("|")
    path = CONSERVATION / f"{parts[2].split('_')[0]}_conservation.txt" if len(parts) > 2 else None
    if not path or not path.exists():
        return {}
    out = {}
    for line in path.read_text().splitlines()[1:]:
        cols = line.split("\t")
        if len(cols) >= 6:
            out[int(cols[0])] = cols[5].strip()
    return out


def trimmed_family(family, reference):
    """Residue counts per trimmed column (whole family and per kingdom group), the full->trimmed
    column map, and the numbering reference's residues by trimmed column."""
    rows = list(read_fasta(TRIMMED / f"{family}_genes_filtered_db_FAMSA.consensus_trimmed.fasta"))
    kept = [int(x) - 1 for x in (TRIMMED / f"{family}_kept_columns.txt").read_text().split()]
    full_to_trim = {f: t for t, f in enumerate(kept)}
    h, s = next((h, s) for h, s in rows if h.split("|")[1] == reference)
    numbers = gpcrdb_numbers(h) if reference != ANCHOR[family] else {}
    _, res = residues(h, s)
    ref = {c: (n, aa, numbers.get(n)) for n, (c, aa) in res.items()}

    arr = np.frombuffer("".join(s for _, s in rows).upper().encode(), dtype=np.uint8).reshape(len(rows), -1)
    accs = [h.split("|")[1] for h, _ in rows]
    counts = {family: (len(rows), {aa: (arr == ord(aa)).sum(axis=0) for aa in AA})}
    for group in ("Animals", "Fungi", "Plants"):
        mask = np.array([KINGDOM.get((family, a)) == group for a in accs])
        if mask.any():
            sub = arr[mask]
            counts[f"{family}_{group}"] = (int(mask.sum()), {aa: (sub == ord(aa)).sum(axis=0) for aa in AA})
    return full_to_trim, ref, counts


def position(nseq, counts, t, ref, conf):
    rc = {aa: int(counts[aa][t]) for aa in AA if counts[aa][t]}
    nongap = sum(rc.values())
    p = {"residueCounts": rc, "totalSequences": nseq, "informationContent": 0, "letterHeights": {}}
    if nongap:
        freqs = {aa: c / nongap for aa, c in rc.items()}
        ic = max(0.0, math.log2(20) + sum(f * math.log2(f) for f in freqs.values()))
        p["informationContent"], p["letterHeights"] = ic, {aa: f * ic for aa, f in freqs.items()}
    if t in ref:
        n, aa, num = ref[t]
        p.update(refResidue=n, refAA=aa, gpcrdb=num)
    p.update(conf)
    return p


def main():
    trim_info = [l.split("\t") for l in (SOURCE / "trim_info.tsv").read_text().splitlines()[1:] if l.strip()]
    rows = [(r[1].strip(), r[0].strip(), (r[2].strip() if len(r) > 2 else "") or None) for r in trim_info]
    rows = [r for r in rows if r[0] != "SLT_TRNS"]
    base = lambda key: key.split("_")[0] if key.split("_")[-1] in ("Animals", "Fungi", "Plants") else key

    maps, counts = {}, {}
    for key, acc1, acc2 in rows:
        if key == base(key):
            full_to_trim, ref, c = trimmed_family(key, REFERENCE[key])
            maps[key] = (full_to_trim, ref)
            counts.update(c)
            print(f"{key}: {len(full_to_trim)} trimmed columns, {c[key][0]} sequences, "
                  + ", ".join(f"{k.split('_')[1]} {v[0]}" for k, v in c.items() if k != key))

    for threshold, (name, mapping_dir) in itertools.product(THRESHOLDS, ALIGNMENTS.items()):
        consensus = defaultdict(dict)
        for r in csv.DictReader(open(OUT / f"{name}_conf{threshold}{TAG}.tsv"), delimiter="\t"):
            if r["consensus_family_column"]:
                consensus[r["family"]][int(r["sup_column"]) - 1] = (int(r["consensus_family_column"]) - 1, {
                    "confidence": float(r["confidence"]), "votes": int(r["votes"]),
                    "nReps": int(r["n_residue_at_column"])})
        out_dir = OUT / f"logo_conf{threshold}{TAG}" / mapping_dir
        out_dir.mkdir(parents=True, exist_ok=True)
        index, missing = [], defaultdict(int)
        for key, acc1, acc2 in rows:
            src = PUBLIC / mapping_dir / f"{key}.json"
            if not src.exists():
                continue
            data = json.loads(src.read_text())
            full_to_trim, ref = maps[base(key)]
            nseq, cnt = counts[key]
            empty = {"residueCounts": {}, "totalSequences": nseq, "informationContent": 0, "letterHeights": {}}
            positions = []
            for c in range(data["length"]):
                f, conf = consensus[base(key)].get(c, (None, None))
                t = full_to_trim.get(f) if f is not None else None
                if f is not None and t is None:
                    missing[key] += 1
                positions.append(position(nseq, cnt, t, ref, conf) if t is not None else dict(empty))
            data["positions"], data["mapping"] = positions, "consensus"
            (out_dir / f"{key}.json").write_text(json.dumps(data))
            entry = {"familyKey": key, "length": data["length"], "acc1": acc1}
            if key != base(key):
                entry["kingdomOf"] = base(key)
            index.append(entry)
        (out_dir / "index.json").write_text(json.dumps(index))
        print(f"{name}, threshold {threshold}: consensus columns not in the trimmed alignment: "
              + (", ".join(f"{k} {v}" for k, v in missing.items() if k == base(k)) or "none"))


if __name__ == "__main__":
    main()
