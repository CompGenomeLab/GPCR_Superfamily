"""Write the four representative sequence sets as the inputs of the superfamily alignments.

Each representative's sequence and DeepTMHMM topology are taken from the stage 03 predictions
and trimmed by extract_tm_region.py to its GPCR-like 7TM core (the first to the seventh helix of
the single OWIMOWIMOWIMOWI pattern) with 5 flanking residues on either side, so the trimmed
headers carry the kept range as '/start-stop'. The sets are the three stage 07 groups plus the
core group with the 18 literature sequences of stage 05, trimmed the same way.

One folder per set in 08_superfamily_alignment/fasta/, holding
    <set>.fasta                    full-length sequences
    <set>_TMregion.fasta           7TM core plus flank, unaligned (the alignment input)
    <set>_TMregion_boundaries.tsv  accession, start, stop of the kept region (1-based)
plus, beside the folders, one 3line with every sequence's topology under its FASTA header,
for projecting the topology onto any of the alignments.

Usage:
    python build_representative_fastas.py
"""
import argparse
import os
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from gpcr_paths import open_text, stage  # noqa: E402
GROUPS = ["all", "no_ste2", "no_ste2ste3classCV2"]
RANGE_RE = re.compile(r"/(\d+)-(\d+)$")


def read_3line(path):
    with open_text(path) as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip()]
    return [(lines[i][1:].split()[0], lines[i + 1], lines[i + 2]) for i in range(0, len(lines), 3)]


def read_fasta(path):
    records = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                records.append([line[1:].split()[0], ""])
            elif records:
                records[-1][1] += line.strip()
    return records


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--topology", default=str(stage("03_topology_filter", "deeptmhmm_predictions_oriented.3line")))
    p.add_argument("--literature", default=str(stage("05_literature_comparison", "literature_reps_7TM_flank5.fasta")))
    p.add_argument("--literature-topology", default=str(stage("05_literature_comparison", "predicted_topologies.3line")))
    p.add_argument("--output-dir", default=str(stage("08_superfamily_alignment", "fasta")))
    args = p.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    lists = {g: [l.strip() for l in open(stage("07_representative_selection", f"representatives_{g}.txt")) if l.strip()]
             for g in GROUPS}
    wanted = set(lists["all"])
    full = {}  # accession -> (header, sequence, topology)
    for h, s, t in read_3line(args.topology):
        acc = h.split("|")[1] if "|" in h else h
        if acc in wanted and acc not in full:
            full[acc] = (h, s, t)
    missing = [a for a in lists["all"] if a not in full]
    if missing:
        sys.exit(f"{len(missing)} representatives have no topology: {missing[:10]}")

    # Trim all representatives once; the sets are subsets of the same trimmed records.
    trimmed, bounds = {}, {}
    with tempfile.TemporaryDirectory() as tmp:
        fa, tl, out = (os.path.join(tmp, n) for n in ("in.fasta", "in.3line", "out.fasta"))
        with open(fa, "w") as f, open(tl, "w") as t:
            for a in lists["all"]:
                h, s, topo = full[a]
                f.write(f">{h}\n{s}\n")
                t.write(f">{h}\n{s}\n{topo}\n")
        subprocess.run([sys.executable, os.path.join(HERE, "extract_tm_region.py"),
                        "--fasta", fa, "--topology", tl, "--out", out, "--boundaries", out + ".tsv"], check=True)
        for h, s in read_fasta(out):
            trimmed[h.split("|")[1]] = (h, s)
        for line in open(out + ".tsv").read().splitlines()[1:]:
            k, start, stop = line.split("\t")
            bounds[k.split("|")[1] if "|" in k else k] = (k, start, stop)
    dropped = [a for a in lists["all"] if a not in trimmed]
    if dropped:
        sys.exit(f"extract_tm_region skipped {len(dropped)} representatives: {dropped}")

    # Literature: trimmed in step 6; full sequence and topology from its DeepTMHMM run, whose
    # headers share no identifier field with the FASTA headers, so they are renamed to them.
    lit_topo = {h.split("|")[0]: (s, t) for h, s, t in read_3line(args.literature_topology)}
    literature = []  # (full header, full sequence, topology, trimmed header, trimmed sequence, start, stop)
    for h, s in read_fasta(args.literature):
        m = RANGE_RE.search(h)
        seq, topo = lit_topo[h.split("|")[1]]
        literature.append((h[:m.start()], seq, topo, h, s, m.group(1), m.group(2)))

    def write_set(name, accessions, lit=()):
        d = os.path.join(args.output_dir, name)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, f"{name}.fasta"), "w") as f_full, \
                open(os.path.join(d, f"{name}_TMregion.fasta"), "w") as f_tm, \
                open(os.path.join(d, f"{name}_TMregion_boundaries.tsv"), "w") as f_b:
            f_b.write("accession\tstart\tstop\n")
            for a in accessions:
                h, s, _ = full[a]
                f_full.write(f">{h}\n{s}\n")
                th, ts = trimmed[a]
                f_tm.write(f">{th}\n{ts}\n")
                f_b.write("\t".join(bounds[a]) + "\n")
            for fh, fs, _, th, ts, start, stop in lit:
                f_full.write(f">{fh}\n{fs}\n")
                f_tm.write(f">{th}\n{ts}\n")
                f_b.write(f"{fh}\t{start}\t{stop}\n")
        print(f"{d}: {len(accessions) + len(lit)} sequences")

    for g in GROUPS:
        write_set(f"representatives_{g}", lists[g])
    write_set("representatives_and_literature", lists["no_ste2ste3classCV2"], literature)

    # One topology file for projecting every set's alignment, keyed by the FASTA headers.
    topology_out = os.path.join(args.output_dir, "representatives_and_literature.3line")
    with open(topology_out, "w") as tl:
        for a in lists["all"]:
            h, s, topo = full[a]
            tl.write(f">{h}\n{s}\n{topo}\n")
        for fh, fs, topo, *_ in literature:
            tl.write(f">{fh}\n{fs}\n{topo}\n")
    print(f"{topology_out}: {len(lists['all']) + len(literature)} topologies")


if __name__ == "__main__":
    main()
