#!/usr/bin/env python3
"""Trim sequences down to their transmembrane core using DeepTMHMM topology.

By default the window is the GPCR-like core. The membrane state is split by the
direction the helix crosses, M for inside to outside and W for outside to
inside, matching project_topology_onto_alignment.py, and a GPCR then reads as
OWIMOWIMOWIMOWI once repeats are collapsed. The seven helices of that pattern
are located and the window runs from the first of them to the last, so the helix
count is never tested directly and a sequence is judged on its architecture
instead. A sequence whose pattern is absent, or present in more than one phase
and therefore not uniquely placed, is skipped and counted.

Anchoring on the pattern means an extra helix in front of the core, which is
usually a signal peptide predicted as a helix, is left outside the window rather
than dragging the start of it backwards. The older behaviour, --first-last,
would trim such a sequence from the spurious helix and carry roughly one helix
of extra sequence into the alignment; it is kept for the case where the input is
known to be clean and nothing should be dropped.

Either way a flank of a few residues is added on each side so the helix ends are
not cut flush, and everything outside the window is dropped, which is what makes
a global alignment of the remaining cores sensible. The helix counts of the
sequences actually written are reported alongside those of the input, so the
composition of the output is visible without testing a count to build it.

The kept region is recorded twice, in whichever form the next step prefers: a
three-column table of accession, start and stop, and, unless --plain-headers is
given, the same range appended to each identifier as '/start-stop'. Both are
1-based inclusive and describe the span actually written, flank included.
"""

import argparse
import re
import sys

GAP_CHARS = "-."
SIDES = ("I", "O")
# A GPCR once repeats are collapsed: N terminus outside, seven helices alternating
# W (outside to inside) and M (inside to outside), C terminus inside.
GPCR_PATTERN = "OWIMOWIMOWIMOWI"
# Where TM1 and TM7 sit inside that pattern, used to map a match back to residues.
CORE_FIRST_HELIX = 1
CORE_LAST_HELIX = 13


def parse_fasta(path):
    """Yield (header, sequence) pairs, header without the leading '>'."""
    header = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].strip()
                chunks = []
            elif header is not None:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def parse_3line(path):
    """Yield (header, sequence, topology) triples from a DeepTMHMM 3line file."""
    header = None
    lines = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, lines[0], lines[1]
                header = line[1:].strip()
                lines = []
            else:
                lines.append(line)
                if len(lines) > 2:
                    raise ValueError(f"More than two lines under header {header!r}")
    if header is not None:
        if len(lines) != 2:
            raise ValueError(f"Header {header!r} is not followed by two lines")
        yield header, lines[0], lines[1]


def id_keys(header):
    """Return the identifiers a header can be matched by, most specific first.

    A DeepTMHMM header looks like 'tr|A0AAD4YJB4|A0AAD4YJB4_OVIAM | TM' and a
    FASTA header like 'tr|A0AAD4YJB4|A0AAD4YJB4_OVIAM', so the first
    whitespace-delimited token matches directly. The bare accession and the
    entry name are added as fallbacks for files that were renamed.
    """
    token = header.split()[0] if header.split() else header
    keys = [token]
    fields = token.split("|")
    if len(fields) >= 3:
        keys.append(fields[1])
        keys.append(fields[2])
    return keys


def build_lookup(topology_path):
    lookup = {}
    for header, sequence, topology in parse_3line(topology_path):
        if len(sequence) != len(topology):
            raise ValueError(
                f"{header}: sequence is {len(sequence)} residues but topology is "
                f"{len(topology)} labels"
            )
        for key in id_keys(header):
            lookup.setdefault(key, (sequence, topology))
    return lookup


def tm_segments(topology, label="M"):
    """Return [(start, end)] half-open, 0-based spans of consecutive TM labels."""
    segments = []
    start = None
    for index, char in enumerate(topology):
        if char == label:
            if start is None:
                start = index
        elif start is not None:
            segments.append((start, index))
            start = None
    if start is not None:
        segments.append((start, len(topology)))
    return segments


def orient_helices(topology):
    """Split the membrane state M into two, by the direction the helix crosses.

    Verbatim in behaviour with project_topology_onto_alignment.py: 'IMO' becomes M
    and 'OMI' becomes W. A run at either terminus is decided from whichever side it
    has, and a re-entrant helix flanked by the same state is left as M.
    """
    out = list(topology)
    for match in re.finditer(r"M+", topology):
        start, end = match.start(), match.end()
        before = topology[start - 1] if start > 0 else ""
        after = topology[end] if end < len(topology) else ""
        if before in SIDES and after in SIDES:
            if before == "O" and after == "I":
                direction = "W"
            elif before == "I" and after == "O":
                direction = "M"
            else:
                continue
        elif before in SIDES:
            direction = "W" if before == "O" else "M"
        elif after in SIDES:
            direction = "W" if after == "I" else "M"
        else:
            continue
        out[start:end] = direction * (end - start)
    return "".join(out)


def find_cores(topology):
    """Return [(tm1_start, tm7_end)] for every GPCR-like core in the label string.

    Occurrences are counted overlapping, because a sequence with more helices than
    one core matches the pattern in two phases and neither phase is the right
    answer; such a sequence has no uniquely located core and is refused upstream.
    """
    runs = [(match.group()[0], match.start(), match.end())
            for match in re.finditer(r"(.)\1*", orient_helices(topology))]
    collapsed = "".join(label for label, _, _ in runs)
    cores = []
    index = collapsed.find(GPCR_PATTERN)
    while index != -1:
        cores.append((runs[index + CORE_FIRST_HELIX][1], runs[index + CORE_LAST_HELIX][2]))
        index = collapsed.find(GPCR_PATTERN, index + 1)
    return cores


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True, help="unaligned FASTA to trim")
    parser.add_argument("--topology", required=True, help="DeepTMHMM .3line file")
    parser.add_argument("--out", required=True, help="output FASTA of trimmed sequences")
    parser.add_argument(
        "--boundaries",
        help="output TSV of boundaries (default: <out>.boundaries.tsv)",
    )
    parser.add_argument(
        "--flank",
        type=int,
        default=5,
        help="residues kept outside the first and last TM helix (default 5)",
    )
    parser.add_argument(
        "--first-last",
        action="store_true",
        help="trim to the first and last TM helix instead of to the GPCR-like core, "
             "keeping every sequence that has any helix at all",
    )
    parser.add_argument(
        "--require-tm",
        type=int,
        help="skip sequences that do not have exactly this many TM helices",
    )
    parser.add_argument(
        "--plain-headers",
        action="store_true",
        help="keep headers unchanged instead of appending the kept region",
    )
    parser.add_argument(
        "--allow-missing",
        action="store_true",
        help="skip sequences that have no topology instead of failing",
    )
    parser.add_argument(
        "--wrap",
        type=int,
        default=0,
        help="wrap output lines at this width (0 = one line per sequence)",
    )
    args = parser.parse_args()

    lookup = build_lookup(args.topology)
    boundaries_path = args.boundaries or args.out + ".boundaries.tsv"

    missing = []
    length_mismatch = []
    residue_mismatch = []
    no_tm = []
    wrong_count = []
    no_core = []
    multi_core = []
    tm_counts = {}
    kept_tm_counts = {}
    written = 0

    with open(args.out, "w") as out, open(boundaries_path, "w") as table:
        table.write("accession\tstart\tstop\n")

        for header, raw in parse_fasta(args.fasta):
            sequence = "".join(c for c in raw if c not in GAP_CHARS)
            entry = None
            for key in id_keys(header):
                if key in lookup:
                    entry = lookup[key]
                    break
            if entry is None:
                missing.append(header)
                continue

            reference, topology = entry
            if len(sequence) != len(topology):
                length_mismatch.append((header, len(sequence), len(topology)))
                continue
            if sequence != reference:
                differences = sum(1 for a, b in zip(sequence, reference) if a != b)
                residue_mismatch.append((header, differences))

            segments = tm_segments(topology)
            tm_counts[len(segments)] = tm_counts.get(len(segments), 0) + 1
            if not segments:
                no_tm.append(header)
                continue
            if args.require_tm is not None and len(segments) != args.require_tm:
                wrong_count.append((header, len(segments)))
                continue

            if args.first_last:
                first_start = segments[0][0]
                last_end = segments[-1][1]
            else:
                cores = find_cores(topology)
                if not cores:
                    no_core.append((header, len(segments)))
                    continue
                if len(cores) > 1:
                    multi_core.append((header, len(cores)))
                    continue
                first_start, last_end = cores[0]
            kept_tm_counts[len(segments)] = kept_tm_counts.get(len(segments), 0) + 1

            kept_start = max(0, first_start - args.flank)
            kept_end = min(len(sequence), last_end + args.flank)
            trimmed = sequence[kept_start:kept_end]

            identifier = header.split()[0] if header.split() else header
            # 1-based inclusive, and the span actually written, flank included.
            start, stop = kept_start + 1, kept_end
            table.write(f"{identifier}\t{start}\t{stop}\n")

            if args.plain_headers:
                out.write(f">{header}\n")
            else:
                # The range goes on the identifier itself, as '/start-stop', which
                # is what project_topology_onto_alignment.py reads back. Put it in
                # a trailing field instead and the parser would not see it.
                rest = header[len(identifier):]
                out.write(f">{identifier}/{start}-{stop}{rest}\n")
            if args.wrap > 0:
                for start in range(0, len(trimmed), args.wrap):
                    out.write(trimmed[start : start + args.wrap] + "\n")
            else:
                out.write(trimmed + "\n")
            written += 1

    print(f"Wrote {written} trimmed sequences to {args.out}")
    print(f"Wrote boundaries to {boundaries_path}")
    print(
        "TM helix counts, all input: "
        + ", ".join(
            f"{count} helices: {number}" for count, number in sorted(tm_counts.items())
        )
    )
    print(
        "TM helix counts, written:   "
        + ", ".join(
            f"{count} helices: {number}" for count, number in sorted(kept_tm_counts.items())
        )
    )

    for header, differences in residue_mismatch:
        print(
            f"WARNING {header}: lengths agree but {differences} residues differ "
            "from the DeepTMHMM sequence",
            file=sys.stderr,
        )
    if not args.first_last:
        print(f"Skipped {len(no_core)} without a GPCR-like core, "
              f"{len(multi_core)} with more than one")
    for header, count in no_core:
        print(f"SKIPPED {header}: no GPCR-like core, {count} TM helices", file=sys.stderr)
    for header, count in multi_core:
        print(f"SKIPPED {header}: {count} GPCR-like cores, none uniquely placed", file=sys.stderr)
    for header, count in wrong_count:
        print(f"SKIPPED {header}: {count} TM helices", file=sys.stderr)
    for header in no_tm:
        print(f"SKIPPED {header}: no TM helix predicted", file=sys.stderr)
    for header, sequence_length, topology_length in length_mismatch:
        print(
            f"ERROR {header}: {sequence_length} residues but "
            f"{topology_length} topology labels",
            file=sys.stderr,
        )
    for header in missing:
        print(f"ERROR {header}: no topology found", file=sys.stderr)

    if length_mismatch or (missing and not args.allow_missing):
        sys.exit(1)


if __name__ == "__main__":
    main()
