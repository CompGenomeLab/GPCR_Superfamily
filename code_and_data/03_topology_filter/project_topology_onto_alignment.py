#!/usr/bin/env python3
"""Replace the residues of an alignment with their DeepTMHMM topology labels.

Every aligned residue becomes the label DeepTMHMM assigned to it and every gap
stays a gap, so the output is an alignment of topology strings that can be read
column by column.

The membrane state is split by direction: M for a helix crossing from inside to
outside and W for one crossing outside to inside, the two letters being mirror
images of each other. Two proteins with the same helix count but opposite
orientation therefore no longer look alike, which the single M state hid.
"""

import argparse
import re
import sys

GAP_CHARS = "-."
SIDES = ("I", "O")


def orient_helices(topology):
    """Split the membrane state M into two, by the direction the helix crosses.

    A helix running inside to outside keeps M, whose strokes rise; one running
    outside to inside becomes W, the same shape inverted. Direction is read from
    the states flanking each run, so 'IMO' becomes M and 'OMI' becomes W. Every
    run in the DeepTMHMM output is one of those two, but a run at either terminus
    is decided from whichever side it has, and anything else (a re-entrant helix
    flanked by the same state) is left as M rather than guessed at.
    """
    out = list(topology)
    for match in re.finditer(r"M+", topology):
        start, end = match.start(), match.end()
        before = topology[start - 1] if start > 0 else ""
        after = topology[end] if end < len(topology) else ""
        # Only I and O say which side the helix starts and ends on; a signal
        # peptide or a terminus says nothing, so it is skipped rather than read.
        if before in SIDES and after in SIDES:
            if before == "O" and after == "I":
                direction = "W"
            elif before == "I" and after == "O":
                direction = "M"
            else:
                continue                      # re-entrant, same state both sides
        elif before in SIDES:
            direction = "W" if before == "O" else "M"
        elif after in SIDES:
            direction = "W" if after == "I" else "M"
        else:
            continue
        out[start:end] = direction * (end - start)
    return "".join(out)


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


RANGE_RE = re.compile(r"[/:](\d+)-(\d+)$")


def split_range(token):
    """'tr|ACC|NAME/61-290' -> ('tr|ACC|NAME', 61, 290); no range -> (token, None, None).

    extract_tm_region.py stamps the residue range it kept onto the header, so a
    trimmed sequence carries its own coordinates and cannot drift out of step
    with a separate table the way a side file can.
    """
    match = RANGE_RE.search(token)
    if not match:
        return token, None, None
    return token[:match.start()], int(match.group(1)), int(match.group(2))


def header_range(header):
    token = header.split()[0] if header.split() else header
    _, start, end = split_range(token)
    return None if start is None else (start, end)


def id_keys(header):
    """Return the identifiers a header can be matched by, most specific first.

    A DeepTMHMM header looks like 'tr|A0AAD4YJB4|A0AAD4YJB4_OVIAM | TM' and an
    alignment header like 'tr|A0AAD4YJB4|A0AAD4YJB4_OVIAM', so the first
    whitespace-delimited token matches directly. Any trimming range appended to
    the token is stripped first, so a trimmed sequence still finds its full
    topology. The bare accession and the entry name are added as fallbacks for
    alignments that were renamed.
    """
    token = header.split()[0] if header.split() else header
    token, _, _ = split_range(token)
    keys = [token]
    fields = token.split("|")
    if len(fields) >= 3:
        keys.append(fields[1])
        keys.append(fields[2])
    return keys


def load_boundaries(path):
    """Read 'accession<TAB>start<TAB>stop' into a lookup keyed by every id form.

    Coordinates are 1-based inclusive and describe the region that was actually
    written, flank included. A header line is skipped when its second field is
    not a number.
    """
    bounds = {}
    with open(path) as handle:
        for line in handle:
            parts = line.rstrip("\n\r").split("\t")
            if len(parts) < 3:
                continue
            try:
                start, stop = int(parts[1]), int(parts[2])
            except ValueError:
                continue                        # header line
            for key in id_keys(parts[0]):
                bounds.setdefault(key, (start, stop))
    if not bounds:
        sys.exit(f"No 'accession<TAB>start<TAB>stop' rows found in {path}")
    return bounds


def region_of(header, bounds):
    """The trimmed region for an alignment row: header range first, file second."""
    from_header = header_range(header)
    if from_header is not None:
        return from_header
    if bounds:
        for key in id_keys(header):
            if key in bounds:
                return bounds[key]
    return None


def slice_topology(topology, start, stop):
    """Cut out the 1-based inclusive region the alignment was trimmed to."""
    return topology[start - 1:stop]


def build_lookup(topology_path, orient=True):
    lookup = {}
    for header, sequence, topology in parse_3line(topology_path):
        if len(sequence) != len(topology):
            raise ValueError(
                f"{header}: sequence is {len(sequence)} residues but topology is "
                f"{len(topology)} labels"
            )
        if orient:
            topology = orient_helices(topology)
        for key in id_keys(header):
            lookup.setdefault(key, (sequence, topology))
    return lookup


def project(aligned, sequence, topology):
    """Map the topology onto the aligned sequence, keeping gaps in place."""
    out = []
    position = 0
    for char in aligned:
        if char in GAP_CHARS:
            out.append("-")
            continue
        if position >= len(topology):
            return None, position
        out.append(topology[position])
        position += 1
    if position != len(topology):
        return None, position
    return "".join(out), position


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment", required=True, help="aligned FASTA")
    parser.add_argument("--topology", required=True, help="DeepTMHMM .3line file")
    parser.add_argument("--out", required=True, help="output FASTA of topology strings")
    parser.add_argument(
        "--boundaries",
        help="three-column TSV, accession<TAB>start<TAB>stop, 1-based inclusive, "
             "giving the region each sequence was trimmed to. Use it when the "
             "alignment is of TM-only sequences.",
    )
    parser.add_argument(
        "--no-orient",
        action="store_true",
        help="keep DeepTMHMM's single M state instead of splitting it into "
             "M (inside to outside) and W (outside to inside)",
    )
    parser.add_argument(
        "--allow-missing",
        action="store_true",
        help="skip alignment sequences that have no topology instead of failing",
    )
    parser.add_argument(
        "--wrap",
        type=int,
        default=0,
        help="wrap output lines at this width (0 = one line per sequence)",
    )
    args = parser.parse_args()

    lookup = build_lookup(args.topology, orient=not args.no_orient)

    # Helix direction is resolved on the full topology inside build_lookup, before
    # any trimming. Doing it after would misread helices whose flanking inside or
    # outside state was cut away, which is exactly what a small --flank leaves.
    bounds = None
    first_header = next((h for h, _ in parse_fasta(args.alignment)), "")
    trimmed = bool(args.boundaries) or header_range(first_header) is not None
    if trimmed:
        bounds = load_boundaries(args.boundaries) if args.boundaries else {}
        if not args.boundaries:
            print("Using the residue ranges stamped on the alignment headers")

    missing = []
    no_boundary = []
    length_mismatch = []
    residue_mismatch = []
    written = 0

    with open(args.out, "w") as out:
        for header, aligned in parse_fasta(args.alignment):
            entry = None
            for key in id_keys(header):
                if key in lookup:
                    entry = lookup[key]
                    break
            if entry is None:
                missing.append(header)
                continue

            sequence, topology = entry
            if trimmed:
                region = region_of(header, bounds)
                if region is None:
                    no_boundary.append(header)
                    continue
                start, stop = region
                sequence = slice_topology(sequence, start, stop)
                topology = slice_topology(topology, start, stop)
            ungapped = "".join(c for c in aligned if c not in GAP_CHARS)
            projected, consumed = project(aligned, sequence, topology)
            if projected is None:
                length_mismatch.append((header, len(ungapped), len(topology)))
                continue
            if ungapped != sequence:
                differences = sum(1 for a, b in zip(ungapped, sequence) if a != b)
                residue_mismatch.append((header, differences))

            out.write(f">{header}\n")
            if args.wrap > 0:
                for start in range(0, len(projected), args.wrap):
                    out.write(projected[start : start + args.wrap] + "\n")
            else:
                out.write(projected + "\n")
            written += 1

    print(f"Wrote {written} projected sequences to {args.out}")

    for header, differences in residue_mismatch:
        print(
            f"WARNING {header}: lengths agree but {differences} residues differ "
            "from the DeepTMHMM sequence",
            file=sys.stderr,
        )
    for header, aligned_length, topology_length in length_mismatch:
        print(
            f"ERROR {header}: {aligned_length} aligned residues but "
            f"{topology_length} topology labels",
            file=sys.stderr,
        )
    for header in missing:
        print(f"ERROR {header}: no topology found", file=sys.stderr)
    for header in no_boundary:
        print(f"ERROR {header}: no TM region in the boundaries file", file=sys.stderr)

    if length_mismatch or no_boundary or (missing and not args.allow_missing):
        sys.exit(1)


if __name__ == "__main__":
    main()
