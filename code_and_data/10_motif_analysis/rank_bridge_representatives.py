#!/usr/bin/env python3
"""Rank each family's representatives by how well their topology agrees with both
their own family alignment and the superfamily alignments.

A sequence that anchors a family inside a superfamily comparison, such as the
reference sequence a family's logo is placed onto the superfamily columns through,
has to be a good member of both worlds. It must reproduce its family's helix bundle,
or it is a poor stand-in for the family, and it must reproduce the superfamily's, or
the columns it maps the family onto are the wrong ones. This script measures both.

Every alignment, the family alignments and the superfamily alignments alike, is
scored with `helix_features` from `03_topology_filter/per_helix_features.py`,
the same function behind `topology_helix_scores.tsv` and the topology viewer. The
consensus of each alignment is taken over all of its rows, at the same occupancy
threshold and minimum block width, so a family number here is the number that table
carries for the same alignment, and a superfamily number is measured on the same
scale. Only the representatives' rows are reported.

A representative is scored in its family alignment and in every superfamily alignment
that contains it, and ranked within its family by

1. its weakest helix across all of those alignments, highest first. A good anchor has
   to match every TM everywhere; a high average with one helix out of register is not
   one, which is why the minimum leads (the rule `find_bridge_sequences.py` introduced);
2. the mean of every per-helix Jaccard across those alignments;
3. its coverage of the superfamily core columns, averaged over the superfamily
   alignments. A column where the anchor has a gap is a column the family contributes
   nothing to, so among equally good anchors the one leaving fewer holes is preferred.

Ties after that break on accession, so a rerun on the same inputs gives the same order.

`--rank-by rank-sum` orders instead by family rank plus superfamily rank on mean
Jaccard, with the weakest helix breaking ties. Both ranks are written in either mode.

Family membership comes from `--groups`, not from which family alignment a sequence
appears in, because alignments built before overlap resolution can hold one accession
in two families.

Usage:
    python rank_bridge_representatives.py \\
        --superfamily ../step9_ml_hierarchy/MSA/representatives_all_Sep22_topo0.8_ginsi_ep0.123_topology_converted.fasta \\
                      ../step9_ml_hierarchy/MSA/representatives_no_ste2_Sep22_topo0.8_ginsi_ep0.123_topology_converted.fasta \\
        --names all no_ste2 \\
        --groups ../step8_representative_selection/representatives_all_families_Sep22_topo0.8.tsv \\
        --family-dir final_topology_msa_gt0.1 \\
        --family-suffix _filtered_genes_db_FAMSA_topology_converted_gt0.1.fasta \\
        --tag Sep24_topo0.8
"""

import argparse
import csv
import glob
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PIPELINE = os.path.dirname(HERE)
sys.path.insert(0, PIPELINE)
sys.path.insert(0, os.path.join(PIPELINE, "03_topology_filter"))
from gpcr_paths import stage  # noqa: E402
from per_helix_features import accession_of, helix_features, parse_fasta  # noqa: E402

DEFAULT_FAMILY_DIR = str(stage("10_motif_analysis", "family_alignments", "topology_alignments"))
DEFAULT_FAMILY_SUFFIX = "_filtered_genes_db_FAMSA_topology_converted_gt0.1.fasta"
RANGE_RE = re.compile(r"/(\d+)-(\d+)")


def load_groups(path):
    """accession -> family, from a two-column map or from a table with `accession` and
    `gpcrfamily` (or `family`) columns, such as step 6's representatives table."""
    with open(path, newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows:
        sys.exit(f"{path} is empty")
    header = [c.strip().lower() for c in rows[0]]
    if "accession" in header and ({"gpcrfamily", "family"} & set(header)):
        i_acc = header.index("accession")
        i_fam = header.index("gpcrfamily" if "gpcrfamily" in header else "family")
        body = rows[1:]
    else:
        i_acc, i_fam, body = 0, 1, rows
    groups = {}
    for row in body:
        if len(row) > max(i_acc, i_fam) and row[i_acc].strip():
            accession, family = row[i_acc].strip(), row[i_fam].strip()
            if groups.get(accession, family) != family:
                sys.exit(f"{accession} is given two families in {path}")
            groups[accession] = family
    return groups


def score_alignment(path, wanted, min_occupancy, min_block):
    """Score every row of one topology alignment against its own consensus and return
    {accession: (per-helix Jaccards, core coverage, header)} for the wanted accessions,
    along with the helix and core-column counts."""
    records = list(parse_fasta(path))
    if not records:
        sys.exit(f"No sequences in {path}")
    features, _, _, _, n_core = helix_features(records, min_occupancy, min_block)

    # Coverage needs the core columns, which helix_features does not return; the mask
    # is recomputed the same way score_all_family_topologies.py does it.
    width = len(records[0][1])
    matrix = np.frombuffer("".join(seq for _, seq in records).encode(),
                           dtype="S1").reshape(len(records), width)
    occupancy = (matrix != b"-").sum(axis=0) / len(records)
    core = np.where(occupancy >= min_occupancy)[0]

    out = {}
    for i, (header, _) in enumerate(records):
        accession = accession_of(header)
        if accession in wanted and accession not in out:
            coverage = float((matrix[i, core] != b"-").sum()) / n_core
            out[accession] = (features[i], coverage, header)
    return out, features.shape[1], n_core


def record_ranks(records, field, out_field):
    """Write each record's rank on `field`, highest first, as `out_field`. Tied values
    share the mean of the ranks they span, so a tie is never broken arbitrarily here."""
    ordered = sorted(records, key=lambda r: -r[field])
    i = 0
    while i < len(ordered):
        j = i
        while j + 1 < len(ordered) and ordered[j + 1][field] == ordered[i][field]:
            j += 1
        for record in ordered[i:j + 1]:
            record[out_field] = (i + j) / 2 + 1
        i = j + 1


def seq_range(header):
    match = RANGE_RE.search(header.split()[0])
    return f"{match.group(1)}-{match.group(2)}" if match else ""


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--superfamily", nargs="+", required=True,
                        help="superfamily topology alignments")
    parser.add_argument("--names", nargs="+", required=True,
                        help="a short label per superfamily alignment, used in column names")
    parser.add_argument("--groups", required=True,
                        help="accession -> family map of the representatives")
    parser.add_argument("--family-dir", default=DEFAULT_FAMILY_DIR,
                        help="folder of per-family topology alignments")
    parser.add_argument("--family-suffix", default=DEFAULT_FAMILY_SUFFIX,
                        help="filename after the family name in --family-dir")
    parser.add_argument("--min-occupancy", type=float, default=0.30,
                        help="core-column threshold, the one topology_helix_scores.tsv uses")
    parser.add_argument("--min-block", type=int, default=8,
                        help="shortest consensus helix run that counts as a block")
    parser.add_argument("--rank-by", choices=["weakest", "rank-sum"], default="weakest",
                        help="weakest: weakest helix across every alignment, then the mean. "
                             "rank-sum: family rank plus superfamily rank on mean Jaccard, "
                             "ties broken by the weakest helix")
    parser.add_argument("--top", type=int, default=5, help="representatives to list per family")
    parser.add_argument("--output-dir", default=str(stage("10_motif_analysis", "anchors")))
    parser.add_argument("--tag", required=True, help="filename suffix, e.g. Sep24_topo0.8")
    args = parser.parse_args()
    if len(args.names) != len(args.superfamily):
        sys.exit("--names must give one label per superfamily alignment")

    groups = load_groups(args.groups)
    wanted = set(groups)
    print(f"{len(wanted)} representatives across {len(set(groups.values()))} families")

    # accession -> {alignment label: (per-helix Jaccards, coverage)}
    scored = {accession: {} for accession in wanted}
    ranges = {}

    for name, path in zip(args.names, args.superfamily):
        print(f"[superfamily {name}] {os.path.basename(path)}", flush=True)
        found, n_helix, n_core = score_alignment(path, wanted, args.min_occupancy, args.min_block)
        print(f"  {n_helix} consensus helices over {n_core} core columns; "
              f"{len(found)} of {len(wanted)} representatives present")
        if n_helix != 7:
            print(f"  WARNING: superfamily consensus has {n_helix} helices, not 7")
        for accession, (js, coverage, header) in found.items():
            scored[accession][name] = (js, coverage)
            ranges.setdefault(accession, seq_range(header))

    families = sorted(set(groups.values()), key=str.lower)
    for family in families:
        path = os.path.join(args.family_dir, family + args.family_suffix)
        members = {a for a, f in groups.items() if f == family}
        if not os.path.exists(path):
            print(f"[{family}] no family alignment at {path}; its representatives are not ranked")
            continue
        print(f"[{family}] family alignment", flush=True)
        found, n_helix, n_core = score_alignment(path, members, args.min_occupancy, args.min_block)
        print(f"  {n_helix} consensus helices over {n_core} core columns; "
              f"{len(found)} of {len(members)} representatives present")
        if n_helix != 7:
            print(f"  WARNING: {family} family consensus has {n_helix} helices, not 7")
        for accession, (js, coverage, _) in found.items():
            scored[accession]["family"] = (js, coverage)

    labels = ["family"] + args.names
    records = []
    missing_family = missing_superfamily = 0
    for accession, by_alignment in scored.items():
        if "family" not in by_alignment:
            missing_family += 1
            continue
        superfamily = [n for n in args.names if n in by_alignment]
        if not superfamily:
            missing_superfamily += 1
            continue
        every = np.concatenate([by_alignment[label][0] for label in by_alignment])
        record = {
            "family": groups[accession],
            "accession": accession,
            "superfamily_range": ranges.get(accession, ""),
            "n_alignments": len(by_alignment),
            "weakest_tm": float(every.min()),
            "mean_tm": float(every.mean()),
            "superfamily_coverage": float(np.mean([by_alignment[n][1] for n in superfamily])),
        }
        for label in labels:
            js, coverage = by_alignment.get(label, (None, None))
            record[f"{label}_mean_jaccard"] = float(js.mean()) if js is not None else None
            record[f"{label}_weakest_tm"] = float(js.min()) if js is not None else None
            record[f"{label}_coverage"] = coverage
            record[f"{label}_tms"] = ",".join(f"{j:.2f}" for j in js) if js is not None else ""
        records.append(record)

    if missing_family:
        print(f"{missing_family} representatives were not found in their family alignment")
    if missing_superfamily:
        print(f"{missing_superfamily} representatives were in no superfamily alignment")

    # Rank on each side separately, within the family, so both orderings are always in
    # the output whichever one the final order is taken from.
    by_family = {}
    for record in records:
        by_family.setdefault(record["family"], []).append(record)
    for members in by_family.values():
        record_ranks(members, "family_mean_jaccard", "family_rank")
        per_alignment = []
        for name in args.names:
            present = [r for r in members if r[f"{name}_mean_jaccard"] is not None]
            record_ranks(present, f"{name}_mean_jaccard", f"_{name}_rank")
            per_alignment.append(name)
        for record in members:
            ranks = [record[f"_{n}_rank"] for n in per_alignment if f"_{n}_rank" in record]
            record["superfamily_rank"] = sum(ranks) / len(ranks)
            record["rank_sum"] = record["family_rank"] + record["superfamily_rank"]

    if args.rank_by == "rank-sum":
        key = lambda r: (r["family"].lower(), r["rank_sum"], -r["weakest_tm"],
                         -r["superfamily_coverage"], r["accession"])
    else:
        key = lambda r: (r["family"].lower(), -r["weakest_tm"], -r["mean_tm"],
                         -r["superfamily_coverage"], r["accession"])
    records.sort(key=key)
    rank = {}
    for record in records:
        rank[record["family"]] = rank.get(record["family"], 0) + 1
        record["rank_in_family"] = rank[record["family"]]

    columns = ["family", "rank_in_family", "accession", "superfamily_range", "n_alignments",
               "family_rank", "superfamily_rank", "rank_sum",
               "weakest_tm", "mean_tm", "superfamily_coverage"]
    for label in labels:
        columns += [f"{label}_mean_jaccard", f"{label}_weakest_tm", f"{label}_coverage",
                    f"{label}_tms"]

    def fmt(value):
        if value is None:
            return ""
        return f"{value:.4f}" if isinstance(value, float) else str(value)

    os.makedirs(args.output_dir, exist_ok=True)
    suffix = "_ranksum" if args.rank_by == "rank-sum" else ""
    stem = os.path.join(args.output_dir, f"bridge_representatives_{args.tag}{suffix}")
    for path, subset in ((stem + ".tsv", records),
                         (stem + "_top.tsv", [r for r in records if r["rank_in_family"] <= args.top])):
        with open(path, "w", newline="") as out:
            out.write("\t".join(columns) + "\n")
            for record in subset:
                out.write("\t".join(fmt(record[c]) for c in columns) + "\n")
        print(f"Wrote {path} ({len(subset)} rows)")

    print(f"\n{'family':<14}{'best':<13}{'famRk':>6}{'sfRk':>6}{'weakest':>8}{'mean':>7}"
          f"{'famMJ':>7}" + "".join(f"{n[:7]:>8}" for n in args.names) + f"{'sfCov':>7}")
    for record in records:
        if record["rank_in_family"] != 1:
            continue
        sf = "".join(f"{fmt(record[f'{n}_mean_jaccard']) or '-':>8}" for n in args.names)
        print(f"{record['family']:<14}{record['accession']:<13}{record['family_rank']:>6g}"
              f"{record['superfamily_rank']:>6g}{record['weakest_tm']:>8.3f}"
              f"{record['mean_tm']:>7.3f}{record['family_mean_jaccard']:>7.3f}{sf}"
              f"{record['superfamily_coverage']:>7.3f}")


if __name__ == "__main__":
    main()
