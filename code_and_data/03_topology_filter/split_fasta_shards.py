#!/usr/bin/env python3
"""Split a large FASTA into shards with balanced total residue counts.

DeepTMHMM predicts topologies one sequence at a time, so a shard's runtime
tracks its total residues rather than its sequence count. Splitting the file
sequentially would give shards with very uneven runtimes whenever the length
distribution is uneven, so sequences are assigned longest-first to whichever
shard currently holds the fewest residues.

Usage:
    python3 split_fasta_shards.py <fasta> --shard-size 25000
    python3 split_fasta_shards.py <fasta> --num-shards 12

Shards are written as <fasta_stem>_shard<N>.fasta with N starting at 1,
into a "shards" directory beside the input unless --outdir is given.
"""

import argparse
import heapq
import os
import sys


def read_fasta(path):
    records = []
    header = None
    seq = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq)))
                header = line
                seq = []
            elif line:
                seq.append(line)
    if header is not None:
        records.append((header, "".join(seq)))
    return records


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fasta")
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--shard-size", type=int, help="max sequences per shard")
    group.add_argument("--num-shards", type=int, help="exact number of shards")
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()

    records = read_fasta(args.fasta)
    if not records:
        sys.exit("ERROR: no sequences read from " + args.fasta)

    n = len(records)
    if args.num_shards:
        n_shards = args.num_shards
    else:
        n_shards = -(-n // args.shard_size)
    if n_shards < 1:
        sys.exit("ERROR: shard count must be at least 1")
    if n_shards > n:
        sys.exit("ERROR: more shards (%d) than sequences (%d)" % (n_shards, n))

    stem = os.path.basename(args.fasta)
    for ext in (".fasta", ".fa", ".faa"):
        if stem.endswith(ext):
            stem = stem[: -len(ext)]
            break
    outdir = args.outdir or os.path.join(os.path.dirname(os.path.abspath(args.fasta)), "shards")
    os.makedirs(outdir, exist_ok=True)

    # Assign longest sequences first to the lightest shard, balancing residues.
    order = sorted(range(n), key=lambda i: len(records[i][1]), reverse=True)
    bins = [[0, s] for s in range(n_shards)]   # [residues, shard index]
    heapq.heapify(bins)
    assigned = [[] for _ in range(n_shards)]
    for i in order:
        residues, s = heapq.heappop(bins)
        assigned[s].append(i)
        heapq.heappush(bins, [residues + len(records[i][1]), s])

    total_seqs = 0
    total_res = 0
    print("%-40s %8s %12s %10s" % ("shard", "seqs", "residues", "mean_len"))
    for s in range(n_shards):
        idx = sorted(assigned[s])          # restore input order within the shard
        path = os.path.join(outdir, "%s_shard%d.fasta" % (stem, s + 1))
        res = 0
        with open(path, "w") as out:
            for i in idx:
                header, seq = records[i]
                out.write(header + "\n")
                for j in range(0, len(seq), 60):
                    out.write(seq[j : j + 60] + "\n")
                res += len(seq)
        total_seqs += len(idx)
        total_res += res
        print("%-40s %8d %12d %10.1f" % (os.path.basename(path), len(idx), res, res / len(idx)))

    print()
    print("input : %d sequences, %d residues" % (n, sum(len(s) for _, s in records)))
    print("output: %d sequences, %d residues across %d shards" % (total_seqs, total_res, n_shards))
    if total_seqs != n or total_res != sum(len(s) for _, s in records):
        sys.exit("ERROR: sequence or residue count mismatch after splitting")
    print("OK: all sequences accounted for")


if __name__ == "__main__":
    main()
