#!/usr/bin/env python3
"""Build the sequence downloads of gpcrevolution.org, full length and unaligned, as searched (UniProtKB 2025_02),
with their UniProt headers. Together they hold every sequence the species pages list:

    <out>/<family>.fasta            the final members of each family (after the GPR1 reassignment,
                                    10_motif_analysis/phylogeny_reassignment/final_membership_reassigned.tsv)
    <out>/removed/<family>.fasta    the members of each family's clade removed by the filters: no single GPCR-like
                                    7TM core (03_topology_filter/core_membership.tsv) or topology score below 0.8
                                    (topology_membership.tsv); a sequence removed from two clades is in both files
    <out>/controls/<control>.fasta  the members of the control clades (G alpha, G beta, G gamma, TBP, RPB1)

GPCR sequences come from 02_clade_extraction/clade_sequences/, controls from their search's sequences
(01_homology_search/<control>/<control>_hmm_E1000_db.fasta.gz).

Each file is sorted by the genome's taxID (the OX= of the header, the taxID of the species pages) and then by
accession, so the sequences of one genome are one contiguous block. A file over MAX_BYTES is split into equal
parts at genome boundaries (<name>_1.fasta, <name>_2.fasta, ...). The site then indexes the files
(node scripts/index-sequence-fasta.mjs in the GPCRevolution repository: for each taxID the file, byte offset,
length and number of sequences), so it fetches one genome's sequences with an HTTP range request instead of the
whole file.

Usage:
    python build_family_fasta.py --out <GPCRevolution>/public/sequences
    (then, in GPCRevolution: node scripts/index-sequence-fasta.mjs)
"""
import argparse
import csv
import gzip
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from gpcr_paths import CONTROLS, FAMILIES, open_text, stage  # noqa: E402

OX = re.compile(r"\sOX=(\d+)")
MAX_BYTES = 95_000_000                                 # GitHub refuses files over 100 MB


def wide(path):
    """{column: set of accessions} of a membership table with one column per family."""
    with open_text(path) as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    return {name: {r[i].strip() for r in rows[1:] if i < len(r) and r[i].strip()} for i, name in enumerate(rows[0])}


def read_fasta(paths, wanted):
    """accession -> (header, sequence) for the wanted accessions (the first file holding one wins)."""
    found = {}
    for path in paths:
        acc = None
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith(">"):
                    acc = line[1:].split()[0].split("|")[1]
                    if acc in wanted and acc not in found:
                        found[acc] = [line[1:].rstrip("\n"), []]
                    else:
                        acc = None
                elif acc:
                    found[acc][1].append(line.strip())
    missing = wanted - set(found)
    if missing:
        raise SystemExit(f"{len(missing)} sequences not found, e.g. {sorted(missing)[:3]}")
    return {acc: (header, "".join(parts)) for acc, (header, parts) in found.items()}


def write_set(folder, name, accessions, seqs):
    """<folder>/<name>.fasta (or its numbered parts), the accessions' sequences grouped by genome."""
    folder.mkdir(parents=True, exist_ok=True)
    blocks = {}                                        # taxID -> its records, in file order
    for taxid, acc in sorted((int(OX.search(seqs[a][0]).group(1)), a) for a in accessions):
        header, seq = seqs[acc]
        blocks.setdefault(taxid, []).append(
            (">" + header + "\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n").encode())
    blocks = [b"".join(records) for records in blocks.values()]
    total = sum(len(b) for b in blocks)
    parts = max(1, -(-total // MAX_BYTES))             # number of files, at most MAX_BYTES each
    names = [f"{name}.fasta"] if parts == 1 else [f"{name}_{i + 1}.fasta" for i in range(parts)]
    # fill the files in turn with whole genomes, moving on once a file holds its share
    files, part, done = [[] for _ in names], 0, 0
    for block in blocks:
        if part < parts - 1 and done >= total * (part + 1) / parts:
            part += 1
        files[part].append(block)
        done += len(block)
    for old in list(folder.glob(f"{name}.fasta")) + list(folder.glob(f"{name}_*.fasta")):
        old.unlink()
    for file_name, chunks in zip(names, files):
        (folder / file_name).write_bytes(b"".join(chunks))
    sizes = " + ".join(f"{sum(len(b) for b in chunks) / 1e6:.1f}" for chunks in files)
    print(f"{folder.name + '/' + name:26s} {len(accessions):7d} sequences, {len(blocks):5d} genomes, {sizes} MB")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    out = Path(ap.parse_args().out)

    clade = wide(stage("02_clade_extraction", "clade_membership.tsv"))
    core = wide(stage("03_topology_filter", "core_membership.tsv"))
    topo = wide(stage("03_topology_filter", "topology_membership.tsv"))
    final = wide(stage("10_motif_analysis", "phylogeny_reassignment", "final_membership_reassigned.tsv"))
    # removed as on the species pages: no single core, or a core but below the topology cut-off
    removed = {f: (clade[f] - core[f]) | (core[f] - topo[f]) for f in FAMILIES}
    controls = wide(stage("02_clade_extraction", "clade_membership_controls.tsv"))

    gpcr = read_fasta(sorted(stage("02_clade_extraction", "clade_sequences").glob("*.fasta.gz")),
                      set().union(*final.values(), *removed.values()))
    ctrl = read_fasta([stage("01_homology_search", c, f"{c}_hmm_E1000_db.fasta.gz") for c in CONTROLS],
                      set().union(*controls.values()))

    for fam in FAMILIES:
        write_set(out, fam, final[fam], gpcr)
    for fam in FAMILIES:
        if removed[fam]:
            write_set(out / "removed", fam, removed[fam], gpcr)
    for c in CONTROLS:
        write_set(out / "controls", c, controls[c], ctrl)


if __name__ == "__main__":
    main()
