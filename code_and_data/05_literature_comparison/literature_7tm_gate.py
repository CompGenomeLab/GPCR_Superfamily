"""Apply the topology gate of stage 03 to the literature sequences that could deepen an origin.

The sequences tested by literature_deeper_origins.py were run through DeepTMHMM
(predicted_topologies.3line). Each prediction is read with its helices split by direction (W, outside to
inside; M, inside to outside) and its consecutive identical states collapsed, exactly as for our own
sequences, and a sequence passes when the GPCR pattern OWIMOWIMOWIMOWI occurs exactly once. The 7TM core
coordinates are kept for trimming (select_literature_representatives.py).

Reads  05_literature_comparison/predicted_topologies.3line
Writes 05_literature_comparison/deeper_origin_7TM_gate.tsv

Usage:
    python literature_7tm_gate.py
"""
import csv
import re
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import open_text, stage  # noqa: E402

SIDES = ("I", "O")
GPCR_PATTERN = "OWIMOWIMOWIMOWI"
CORE_FIRST_HELIX, CORE_LAST_HELIX = 1, 13   # TM1 and TM7 within the pattern


def orient_helices(topology):
    """M -> M (inside to outside) or W (outside to inside), as in 03_topology_filter."""
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


def find_cores(oriented):
    runs = [(m.group()[0], m.start(), m.end()) for m in re.finditer(r"(.)\1*", oriented)]
    collapsed = "".join(label for label, _, _ in runs)
    cores, index = [], collapsed.find(GPCR_PATTERN)
    while index != -1:
        cores.append((runs[index + CORE_FIRST_HELIX][1], runs[index + CORE_LAST_HELIX][2]))
        index = collapsed.find(GPCR_PATTERN, index + 1)
    return collapsed, cores


def classify(collapsed, n_cores):
    if n_cores == 1:
        return "GPCR_like_core"
    if n_cores > 1:
        return "multi_core"
    if collapsed.endswith("IMOWIMOWIMOWI") and collapsed.lstrip("S") == "IMOWIMOWIMOWI":
        return "missing_TM1"
    if collapsed.lstrip("S") == "OWIMOWIMOWIMO":
        return "missing_TM7"
    if "IMOWIMOWIMOWIMO" in collapsed:
        return "inverted"
    return "other"


def main():
    out_dir = stage("05_literature_comparison")
    with open_text(out_dir / "predicted_topologies.3line") as fh:
        lines = [line.rstrip("\n") for line in fh]
    rows = []
    for i in range(0, len(lines), 3):
        header = lines[i][1:]
        query = header.split(" | ")[0]
        prediction = header.split(" | ")[1] if " | " in header else ""
        seq, labels = lines[i + 1], lines[i + 2]
        oriented = orient_helices(labels)
        collapsed, cores = find_cores(oriented)
        f = query.split("|")
        rows.append(dict(query=query, id=f[0], study=f[1], claimed=f[2][8:], genome=f[3],
                         prediction=prediction, length=len(seq),
                         n_tm=sum(1 for _ in re.finditer(r"[MW]+", oriented)), core_count=len(cores),
                         topology_class=classify(collapsed, len(cores)),
                         core_start=cores[0][0] + 1 if len(cores) == 1 else "",
                         core_stop=cores[0][1] if len(cores) == 1 else "",
                         collapsed=collapsed, passes=len(cores) == 1))
    with open(out_dir / "deeper_origin_7TM_gate.tsv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(len(rows), "sequences; pass", sum(r["passes"] for r in rows))
    print(Counter(r["topology_class"] for r in rows))


if __name__ == "__main__":
    main()
