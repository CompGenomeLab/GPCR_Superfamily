"""Map each superfamily alignment column to a family alignment column by the votes of all
the family's representatives, instead of through the single anchor (acc1).

Every representative sits in its family's TM-only FAMSA alignment
(family_reps_TMonly_FAMSA/, every column kept) and in the superfamily alignment; residue
numbers come from the /start-end range in each header, as in the website's
precompute-mappings.mjs. For a superfamily column c and family F, each representative with
a residue at c votes for the family column holding that residue. The consensus column f is
pair (c, f) has

    confidence = votes for f / representatives with a residue at c or at f

so a representative that has a residue at c but elsewhere in the family alignment, or at f
but a gap at c, counts against it. The mapping is one-to-one: all pairs at or above the
threshold are taken in order of confidence (then votes), each only if neither its superfamily
nor its family column is already used. A column that loses its most voted family column to a
better-supported one can therefore take its next one if that passes the threshold. The most
voted column (top_family_column) and the anchor's own column are written beside the assigned
one. Rows are written only where some representative has a residue at c.

python consensus_column_mapping.py [threshold ...]   (default 0.33)

Writes consensus_mapping/<alignment>_conf<threshold>_Oct10.tsv for each of the 8 Oct08
alignments and prints a per-family summary.
"""
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import open_text, stage  # noqa: E402
DATA = stage("10_motif_analysis")
PUBLIC = DATA / "website_inputs"
SOURCE = PUBLIC / "superfamily_logo_source_files"
REPS_DIR = DATA / "family_alignments" / "family_reps_TMonly_FAMSA"
OUT = DATA / "consensus_mapping"
import superfamily_motif_positions as smp  # noqa: E402

# superfamily alignment -> website mapping folder built from it (same column count)
ALIGNMENTS = {
    "representatives_all_Oct08_ginsi_ep0.123": "superfamily_logo_mappings",
    "representatives_all_Oct08_ginsi": "superfamily_logo_mappings_lower_gap",
    "representatives_no_ste2_Oct08_ginsi_ep0.123": "superfamily_logo_mappings_no_ste2",
    "representatives_no_ste2_Oct08_ginsi": "superfamily_logo_mappings_no_ste2_lower_gap",
    "representatives_no_ste2ste3classCV2_Oct08_ginsi_ep0.123": "superfamily_logo_mappings_no_ste2ste3classCV2",
    "representatives_no_ste2ste3classCV2_Oct08_ginsi": "superfamily_logo_mappings_no_ste2ste3classCV2_lower_gap",
    "representatives_and_literature_Oct08_ginsi_ep0.123": "superfamily_logo_mappings_and_literature",
    "representatives_and_literature_Oct08_ginsi": "superfamily_logo_mappings_and_literature_lower_gap",
}
# tree-guided G-INS-i (step9 guide_trees), built by build_treein_logo_inputs.sh
TREEIN = {
    "representatives_all_Oct08_ginsi_ep0.123_treein": "superfamily_logo_mappings_treein",
    "representatives_all_Oct08_ginsi_treein": "superfamily_logo_mappings_lower_gap_treein",
}


def read_fasta(path):
    header, seq = None, []
    for line in open_text(path):
        line = line.rstrip()
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq)
            header, seq = line[1:], []
        elif line:
            seq.append(line)
    if header is not None:
        yield header, "".join(seq)


def residues(header, seq):
    """{residue number: (column, residue)} from the /start-end range of the header."""
    acc = header.split("|")[1]
    start = int(header.split("/")[1].split("-")[0]) if "/" in header else 1
    out, n = {}, start - 1
    for col, aa in enumerate(seq):
        if aa != "-":
            n += 1
            out[n] = (col, aa.upper())
    return acc, out


def load_families():
    fams = {}
    for path in sorted(REPS_DIR.glob("*_family_reps_TMonly_FAMSA.fasta*")):
        family = path.name.replace(".gz", "").replace("_family_reps_TMonly_FAMSA.fasta", "")
        fams[family] = dict(residues(h, s) for h, s in read_fasta(path))
    return fams


def anchors():
    rows = [l.split("\t") for l in (SOURCE / "trim_info.tsv").read_text().splitlines()[1:] if l.strip()]
    return {r[1].strip(): r[0].strip() for r in rows}


def label(labelled, c):
    """GPCRdb number at superfamily column c, as smp.label, from the families present."""
    for f in smp.LABEL_ORDER:
        p = labelled.get(f, [None] * (c + 1))[c]
        if p and p["label"] and "x" in str(p["label"]):
            return str(p["label"])
    return f"col{c + 1}"


def map_alignment(name, mapping_dir, fams, acc1, threshold, tag=""):
    sup = dict(residues(h, s) for h, s in read_fasta(SOURCE / f"{name}_plus_SLT_TRNS_topology.fasta"))
    smp.MAPPINGS = PUBLIC / mapping_dir
    labelled = {f: smp.family_columns(f) for f in smp.LABEL_ORDER
                if (smp.MAPPINGS / f"{f}.json").exists()}
    rows, summary = [], {}
    for family, reps in fams.items():
        present = [a for a in reps if a in sup]
        mismatch = sum(1 for a in present for n, (_, aa) in sup[a].items()
                       if n in reps[a] and reps[a][n][1] != aa)
        if mismatch:
            print(f"  {family}: {mismatch} residues differ between family and superfamily rows")
        # per superfamily column: rep -> family column (None when the residue is outside the
        # family alignment); per family column: reps with a residue there
        at_sup = defaultdict(dict)
        at_fam = defaultdict(set)
        for a in present:
            for n, (c, _) in sup[a].items():
                at_sup[c][a] = reps[a][n][0] if n in reps[a] else None
            for n, (f, _) in reps[a].items():
                at_fam[f].add(a)
        anchor = acc1.get(family)
        # every (superfamily column, family column) pair voted for, with its confidence
        cand, top = [], {}
        for c in sorted(at_sup):
            votes = Counter(f for f in at_sup[c].values() if f is not None)
            for f, v in votes.items():
                cand.append((v / len(set(at_sup[c]) | at_fam[f]), v, c, f))
            if votes:
                f, v = votes.most_common(1)[0]
                top[c] = (f, v / len(set(at_sup[c]) | at_fam[f]))
        # one-to-one: pairs at or above the threshold taken in order of confidence, then votes,
        # each kept only if neither column is already used
        assigned, used = {}, set()
        for conf, v, c, f in sorted(cand, key=lambda x: (-x[0], -x[1], x[2], x[3])):
            if conf >= threshold and c not in assigned and f not in used:
                assigned[c] = (f, v, conf)
                used.add(f)
        stats = Counter()
        for c in sorted(at_sup):
            f, v, conf = assigned.get(c, (None, 0, 0.0))
            t_f, t_conf = top.get(c, (None, 0.0))
            a_col = at_sup[c].get(anchor, "-") if anchor in present else "-"
            agree = "" if a_col in ("-", None) else str(a_col == f)
            rows.append((c + 1, label(labelled, c), family, len(present), len(at_sup[c]), v,
                         "" if f is None else f + 1, f"{conf:.3f}",
                         "" if t_f is None else t_f + 1, f"{t_conf:.3f}",
                         "" if a_col in ("-", None) else a_col + 1, agree))
            if a_col not in ("-", None):
                stats["anchor_cols"] += 1
                stats["agree"] += a_col == f
                stats["mapped"] += f is not None
                stats["from_top"] += f is not None and f == t_f
        stats["assigned"] = len(assigned)
        stats["fallback"] = sum(1 for c, (f, _, _) in assigned.items() if f != top[c][0])
        summary[family] = (len(present), stats)
    OUT.mkdir(exist_ok=True)
    with open(OUT / f"{name}_conf{threshold}{tag}.tsv", "w") as fh:
        fh.write("sup_column\tlabel\tfamily\tn_reps\tn_residue_at_column\tvotes\tconsensus_family_column"
                 "\tconfidence\ttop_family_column\ttop_confidence\tanchor_family_column\tanchor_agrees\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    return summary


def main():
    thresholds = [float(x) for x in sys.argv[1:]] or [0.33]
    fams, acc1 = load_families(), anchors()
    for threshold in thresholds:
        for name, mapping_dir in ALIGNMENTS.items():
            print(f"{name}, threshold {threshold}")
            summary = map_alignment(name, mapping_dir, fams, acc1, threshold)
            print("  family        reps  mapped  fallback  anchor_cols  anchor_cols_mapped  anchor=consensus")
            for family, (n, s) in summary.items():
                k = s["anchor_cols"] or 1
                print(f"  {family:13} {n:4}  {s['assigned']:6}  {s['fallback']:8}  {s['anchor_cols']:11}  "
                      f"{s['mapped']:18}  {s['agree'] / k:16.3f}")


if __name__ == "__main__":
    main()
