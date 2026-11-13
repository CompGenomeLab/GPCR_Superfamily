#!/usr/bin/env python3
"""Build a self-contained HTML viewer for a topology-projected alignment.

Reads the output of project_topology_onto_alignment.py, computes per-column
occupancy and consensus topology, and writes one HTML file that draws the
alignment as a topology map, hides low-occupancy columns on a slider, and
scores every sequence by how well it matches the column consensus.
"""

import argparse
import base64
import collections
import gzip
import json
import os
import re
import sys


RANGE_RE = re.compile(r"[/:]\d+-\d+$")


def strip_range(token):
    """Drop a trimming range such as '/61-290' that extract_tm_region.py appends."""
    return RANGE_RE.sub("", token)

GAP = "-"
LABEL_ORDER = ["M", "W", "I", "O", "S", "B", "P"]
TOPOLOGY_LABELS = set(LABEL_ORDER)
AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")


def detect_alphabet(labels):
    """'topology' when every label is a topology state, else 'protein'.

    The two alphabets share the letters B, M, O, P, S and W, so a page must be
    told which one it is holding rather than guessing per column.
    """
    return "topology" if set(labels) <= TOPOLOGY_LABELS else "protein"


def leaf_accession(name, known):
    """Accession of a tree leaf, whatever decoration the label carries.

    Leaves are written either 'tr|A0A123|A0A123_HUMAN/3-277' or, once
    label_tree_families.py has run, 'classA|tr|A0A123|A0A123_HUMAN/3-277', so
    the field that is actually the accession is found by asking the alignment.
    """
    for field in name.split("|"):
        field = strip_range(field.strip())
        if field in known:
            return field
    return short_id(name)


def tree_structure(tree, accession_of):
    """The rooted tree as flat arrays the page can draw without a tree library.

    x is the distance from the root, y is the row the node sits on (its leaf
    rank for a leaf, the mean of its children for an internal node), so a
    rectangular phylogram is two line segments per node and nothing else.
    """
    nodes = list(tree.traverse("postorder"))     # children before parents
    index = {node: i for i, node in enumerate(nodes)}

    depth = {}
    for node in tree.traverse("preorder"):
        depth[node] = depth[node.up] + node.dist if node.up else 0.0

    rank_of = {index[leaf]: rank for rank, leaf in enumerate(tree.iter_leaves())}
    y = [0.0] * len(nodes)
    for i, node in enumerate(nodes):
        if node.is_leaf():
            y[i] = float(rank_of[i])
        else:
            kids = [index[child] for child in node.children]
            y[i] = sum(y[k] for k in kids) / len(kids)

    return {
        "x": [round(depth[node], 6) for node in nodes],
        "y": [round(v, 4) for v in y],
        "parent": [index[node.up] if node.up else -1 for node in nodes],
        "support": [round(float(node.support), 4) if not node.is_leaf() else -1.0
                    for node in nodes],
        "leafAcc": [accession_of[node.name] if node.is_leaf() else None for node in nodes],
        "maxDepth": round(max(depth.values()), 6),
    }


def load_tree_ranks(newick, known, root_spec=None, group_map=None):
    """accession -> position in the tree's leaf order, after optional rerooting.

    root_spec is 'midpoint', a comma-separated list of group names from
    --groups, or a file of accessions; the tree is rooted on the edge that
    separates that set from everything else, which is the edge a reader means
    by 'root it between these two clades'.
    """
    try:
        from ete3 import Tree
    except ImportError:
        sys.exit("--tree needs ete3 installed (pip install ete3)")
    # format=0 reads the numbers on internal nodes as support. format=1 would
    # read them as node *names* and silently leave every support at a default
    # of 1.0, so a tree carrying TBE values must not be read that way; format=1
    # is kept only as the fallback for a tree that has no support values.
    try:
        tree = Tree(newick, format=0)
    except Exception:
        tree = Tree(newick, format=1)
    accession_of = {leaf.name: leaf_accession(leaf.name, known) for leaf in tree}

    target = None
    if root_spec and root_spec != "midpoint":
        if os.path.exists(root_spec):
            with open(root_spec) as handle:
                target = {line.strip() for line in handle if line.strip()}
        elif group_map:
            wanted = {part.strip() for part in root_spec.split(",") if part.strip()}
            target = {a for a, g in group_map.items() if g in wanted}
            if not target:
                sys.exit(f"--tree-root named no group present in --groups: {root_spec}")
        else:
            sys.exit("--tree-root by group name needs --groups")

    if root_spec == "midpoint":
        tree.set_outgroup(tree.get_midpoint_outgroup())
    elif target:
        names = {n for n, a in accession_of.items() if a in target}
        if not names:
            sys.exit("--tree-root matched no leaf of the tree")
        # Rooting is a property of an edge, not of a side of it, and the file's
        # own root is arbitrary, so the wanted set may sit on either side. Both
        # are scored and the edge whose split best matches is the one used.
        everything = set(accession_of)
        complement = everything - names
        best = None
        for node in tree.traverse():
            if node.is_root():
                continue
            leaves = set(node.get_leaf_names())
            misplaced = min(len(leaves ^ names), len(leaves ^ complement))
            if best is None or misplaced < best[0]:
                best = (misplaced, node)
        misplaced, node = best
        if misplaced > 0.2 * len(names):
            sys.exit(f"--tree-root set is not a clade of the tree: the closest edge "
                     f"misplaces {misplaced} leaves")
        tree.set_outgroup(node)
        print(f"  rooted on the edge splitting {len(node)} leaves from "
              f"{len(everything) - len(node)}; {misplaced} leaves off the requested split")

    order = [accession_of[leaf.name] for leaf in tree]
    ranks = {accession: rank for rank, accession in enumerate(order)}
    return ranks, tree_structure(tree, accession_of)


def parse_fasta(path):
    header = None
    chunks = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
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


def short_id(header):
    """'tr|A0A8D2HSC9|A0A8D2HSC9_UROPR Calcitonin receptor ...' -> accession."""
    token = header.split()[0] if header.split() else header
    token = strip_range(token)
    fields = token.split("|")
    return fields[1] if len(fields) >= 3 else token


def gene_of(header):
    """The GN= field UniProt puts in the description, else an empty string."""
    if " GN=" not in header:
        return ""
    return header.split(" GN=", 1)[1].split()[0]


def load_family_scores(path):
    """(family, accession) -> the sequence's score in its OWN family's alignment.

    This is the `mean_jaccard` column of `score_all_family_topologies.py`, computed
    once per family against that family's own consensus. It answers a different
    question from the score the page computes live, which is against the consensus
    of whatever alignment is on screen, and the two can disagree sharply: an
    olfactory receptor can sit at 1.000 inside Olfactory and at 0.83 against the
    superfamily, because the family's architecture is genuinely divergent. Both are
    carried so a row can be read either way.

    Keyed by (family, accession), not by accession alone. 775 accessions are scored
    in more than one family, because a sequence contested between two families is
    scored in each one's alignment, and the two numbers can be far apart: A0A1Y2D3Q0
    is 0.963 in cAMP and 0.674 in class A. 33 of the representatives are such
    accessions, so collapsing on the accession would silently show a score from the
    wrong family for them.
    """
    scores = {}
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        try:
            i_fam = header.index("family")
            i_acc, i_val = header.index("accession"), header.index("mean_jaccard")
        except ValueError:
            sys.exit(f"{path}: needs 'family', 'accession' and 'mean_jaccard' columns; got {header}")
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(i_fam, i_acc, i_val):
                continue
            try:
                scores[(parts[i_fam], parts[i_acc])] = round(float(parts[i_val]), 4)
            except ValueError:
                continue
    if not scores:
        sys.exit(f"No usable rows in {path}")
    return scores


def load_group_map(path):
    """accession -> group, from a TSV or from a directory of per-group lists."""
    mapping = {}
    if os.path.isdir(path):
        for entry in sorted(os.listdir(path)):
            if not entry.endswith(".txt"):
                continue
            group = entry.split("_accessions")[0]
            with open(os.path.join(path, entry)) as handle:
                for line in handle:
                    accession = line.strip()
                    if accession:
                        mapping[accession] = group
    else:
        with open(path) as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2 and parts[0]:
                    mapping[parts[0]] = parts[1]
    if not mapping:
        sys.exit(f"No accession-to-group assignments found in {path}")
    return mapping


def load_descriptions(paths, wanted):
    """accession -> UniProt description, read from search outputs.

    Takes hmmsearch output (the '>> ' line of each hit) and m8 files (the full
    header in column 3), or directories of them. A description belongs to the
    accession, not the search, so every file can be read into one map. Only
    accessions in `wanted` are kept."""
    files = []
    for path in paths:
        if os.path.isdir(path):
            files += sorted(os.path.join(root, f) for root, _, names in os.walk(path)
                            for f in names if f.endswith(("_hmmsearch_E1000.txt", ".m8", ".m8.gz",
                                                          "_hmmsearch_E1000.txt.gz")))
        else:
            files.append(path)
    found = {}
    for path in files:
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith(">> "):
                    header = line[3:].strip()
                elif path.endswith((".m8", ".m8.gz")) and line.count("	") >= 2:
                    header = line.split("	", 3)[2].strip()
                else:
                    continue
                parts = header.split(None, 1)
                if len(parts) < 2:
                    continue
                accession = short_id(parts[0])
                if accession in wanted and accession not in found:
                    found[accession] = parts[1]
    return found


def _group_columns(ids, path, kind):
    """accession -> group index, for either --groups or --lineages.

    Both are listed largest first so the biggest blocks lead, and anything
    unassigned collects into a trailing "unassigned" block.
    """
    mapping = load_group_map(path)
    assigned = [mapping.get(accession, "unassigned") for accession in ids]
    sizes = collections.Counter(assigned)
    unassigned = sizes.pop("unassigned", 0)
    names = [name for name, _ in sizes.most_common()]
    if unassigned:
        names.append("unassigned")
    index_of = {name: i for i, name in enumerate(names)}
    of = [index_of[name] for name in assigned]
    print(f"{len(names)} {kind}s; {len(ids) - unassigned} of {len(ids)} sequences assigned")
    if unassigned:
        print(f"  {unassigned} sequences had no {kind} and are shown as 'unassigned'")
    return names, of


def build_payload(path, max_gap_payload, groups, descriptions=None,
                  alphabet=None, tree_ranks=None, tree_struct=None, family_scores=None,
                  lineages=None):
    """One alignment -> the compressed payload the page decodes on demand."""
    records = list(parse_fasta(path))
    if not records:
        sys.exit(f"No sequences found in {path}")

    # Alignments whose headers were cut to the ID get their description back
    # from the search outputs, so the tooltip and gene names still read.
    if descriptions:
        filled = []
        for header, seq in records:
            if len(header.split(None, 1)) < 2 and short_id(header) in descriptions:
                header = f"{header} {descriptions[short_id(header)]}"
            filled.append((header, seq))
        records = filled
        named = sum(1 for header, _ in records if len(header.split(None, 1)) > 1)
        print(f"{named} of {len(records)} sequences have a description")

    widths = {len(seq) for _, seq in records}
    if len(widths) != 1:
        sys.exit(f"Sequences are not all the same width: {sorted(widths)[:5]}...")
    width = widths.pop()
    n_seq = len(records)

    # Per-column label counts. Consensus is a property of the column alone, so it
    # is computed once here rather than on every slider move in the browser.
    counts = [collections.Counter() for _ in range(width)]
    for _, seq in records:
        for index, char in enumerate(seq):
            counts[index][char] += 1

    labels_seen = sorted(
        {c for column in counts for c in column if c != GAP},
        key=lambda c: LABEL_ORDER.index(c) if c in LABEL_ORDER else 99,
    )

    kept = []
    for index, column in enumerate(counts):
        gap_fraction = column[GAP] / n_seq
        if gap_fraction > max_gap_payload:
            continue
        non_gap = [(c, k) for c, k in column.items() if c != GAP]
        consensus, consensus_count = max(non_gap, key=lambda item: item[1])
        kept.append({
            "i": index,
            "gap": round(gap_fraction, 5),
            "cons": consensus,
            "consFrac": round(consensus_count / max(1, n_seq - column[GAP]), 5),
            "counts": {c: column[c] for c in labels_seen if column[c]},
        })

    if not kept:
        sys.exit("No columns survived --max-gap-payload; raise the cutoff.")

    kept_indices = [column["i"] for column in kept]
    matrix = "".join("".join(seq[i] for i in kept_indices) for _, seq in records)

    ids = [short_id(header) for header, _ in records]

    # Row grouping, by GPCR family and, independently, by taxonomic lineage.
    group_names, group_of = [], []
    if groups:
        group_names, group_of = _group_columns(ids, groups, "group")

    lineage_names, lineage_of = [], []
    if lineages:
        lineage_names, lineage_of = _group_columns(ids, lineages, "lineage")

    # Row order taken from a tree, as a rank per row; absent when no tree was
    # given, or for the rows a tree of a different sequence set does not hold.
    tree_rank = []
    if tree_ranks:
        missing = [a for a in ids if a not in tree_ranks]
        if missing:
            print(f"  {len(missing)} of {n_seq} sequences are not leaves of the tree; "
                  f"they sort to the end of tree order")
        tree_rank = [tree_ranks.get(a, len(tree_ranks) + i) for i, a in enumerate(ids)]

    # The drawn tree needs each leaf to name the row it belongs to, so the
    # accessions it was built with are resolved against this alignment's rows.
    drawn_tree = None
    if tree_struct and tree_rank:
        row_of = {accession: i for i, accession in enumerate(ids)}
        drawn_tree = {k: v for k, v in tree_struct.items() if k != "leafAcc"}
        drawn_tree["leaf"] = [row_of.get(a, -1) if a else -1
                              for a in tree_struct["leafAcc"]]

    # null where the table has no row for an accession, so the tooltip can say so
    # rather than printing a misleading zero.
    family_score_list = None
    if family_scores:
        family_score_list = []
        for k, accession in enumerate(ids):
            fam = group_names[group_of[k]] if group_names and group_of else None
            family_score_list.append(family_scores.get((fam, accession)) if fam else None)
        found = sum(1 for v in family_score_list if v is not None)
        print(f"  own-family score matched for {found} of {len(ids)} sequences")

    payload = {
        "nSeq": n_seq,
        "nCol": len(kept),
        "alignmentWidth": width,
        "alphabet": alphabet or detect_alphabet(labels_seen),
        "treeRank": tree_rank,
        "tree": drawn_tree,
        "labels": labels_seen,
        "columns": kept,
        "ids": ids,
        "descriptions": [header for header, _ in records],
        "genes": [gene_of(header) for header, _ in records],
        "groupNames": group_names,
        "groupOf": group_of,
        "lineageNames": lineage_names,
        "lineageOf": lineage_of,
        "familyScore": family_score_list,
        "matrix": matrix,
    }
    blob = base64.b64encode(
        gzip.compress(json.dumps(payload, separators=(",", ":")).encode(), 6)
    ).decode()
    print(f"{n_seq} sequences x {width} alignment columns, "
          f"{len(kept)} carried into the file (gap <= {max_gap_payload:.0%}), "
          f"alphabet {payload['alphabet']}")
    return blob


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology-alignment", required=True, nargs="+",
                        help="one or more projected topology FASTAs; with several, the page "
                             "gets a selector to switch between them")
    parser.add_argument("--names", nargs="+", default=None,
                        help="a display name per alignment (default: the file name)")
    parser.add_argument("--out", required=True, help="output HTML file")
    parser.add_argument("--title", default=None, help="title shown in the viewer")
    parser.add_argument("--default-max-gap", type=float, default=0.10,
                        help="columns with a gap fraction above this start hidden (default 0.10)")
    parser.add_argument("--max-gap-payload", type=float, default=0.90,
                        help="columns above this gap fraction are left out of the file entirely "
                             "(default 0.90), which bounds the slider's range")
    parser.add_argument("--groups", default=None,
                        help="group the rows: either a two-column accession<TAB>group file, or a "
                             "directory of per-group accession lists named <group>_accessions*.txt")
    parser.add_argument("--lineages", default=None,
                        help="group the rows by taxonomic lineage instead of family, in the same "
                             "accession<TAB>lineage two-column format as --groups. Independent of "
                             "--groups: the page gets 'lineage'/'lineage_median' row-order options "
                             "alongside the family 'group'/'group_median' ones, so lineage blocks "
                             "can be viewed regardless of family.")
    parser.add_argument("--default-group-sort", default="group",
                        choices=["group", "group_median", "group_mean", "lineage", "lineage_median",
                                 "score_desc", "score", "tree"],
                        help="row order the page opens with when --groups is given (default: group); "
                             "'tree' needs --tree. 'group_median'/'lineage_median' rank the families "
                             "or lineages by their median score; 'group_mean' is the old name for "
                             "'group_median' and behaves identically.")
    parser.add_argument("--alphabet", nargs="+", default=None,
                        choices=["topology", "protein"],
                        help="one per alignment, overriding the guess made from the letters "
                             "present (topology alignments and protein alignments share the "
                             "letters B, M, O, P, S and W, so they are coloured differently)")
    parser.add_argument("--tree", default=None,
                        help="Newick tree whose leaf order becomes a row-order option; leaves "
                             "may carry a family prefix and a trimming range")
    parser.add_argument("--trees", nargs="+", default=None,
                        help="one tree per alignment instead of a single shared --tree, for a "
                             "page whose alignments were inferred into different trees; '-' "
                             "leaves an alignment without one")
    parser.add_argument("--tree-root", default=None,
                        help="reroot the tree before taking its leaf order: 'midpoint', a "
                             "comma-separated list of --groups group names, or a file of "
                             "accessions. The root is put on the edge separating that set "
                             "from the rest.")
    parser.add_argument("--family-scores", default=None,
                        help="topology_helix_scores.tsv from score_all_family_topologies.py. "
                             "Adds each sequence's score in its own family's alignment to the "
                             "tooltip, beside the score computed live for the alignment on "
                             "screen. It is shown only, never used for ordering.")
    parser.add_argument("--descriptions", nargs="+", default=None,
                        help="hmmsearch outputs or m8 files (or directories of them) to take "
                             "UniProt descriptions from, for headers that carry only the ID")
    args = parser.parse_args()

    paths = args.topology_alignment
    names = args.names or [os.path.basename(path) for path in paths]
    if len(names) != len(paths):
        sys.exit("--names must give one name per alignment")
    alphabets = args.alphabet or [None] * len(paths)
    if len(alphabets) != len(paths):
        sys.exit("--alphabet must give one value per alignment")

    family_scores = None
    if args.family_scores:
        family_scores = load_family_scores(args.family_scores)
        print(f"Within-family scores for {len(family_scores)} accessions")

    descriptions = None
    if args.descriptions:
        wanted = {short_id(header) for path in paths for header, _ in parse_fasta(path)}
        print(f"Reading descriptions for {len(wanted)} accessions...", flush=True)
        descriptions = load_descriptions(args.descriptions, wanted)
        print(f"  found {len(descriptions)}")

    if args.tree and args.trees:
        sys.exit("give either --tree or --trees, not both")
    if args.trees and len(args.trees) != len(paths):
        sys.exit("--trees must give one tree (or '-') per alignment")

    # Tree order is matched up by accession, so a protein alignment and its
    # topology projection sort the same way. --tree is shared by every
    # alignment; --trees gives each its own, read once per distinct file.
    known = {short_id(header) for path in paths for header, _ in parse_fasta(path)} \
        if (args.tree or args.trees) else None
    group_map = load_group_map(args.groups) if args.groups and known else None
    loaded = {}

    def read_tree(tree_path):
        if tree_path not in loaded:
            print(f"Reading tree {os.path.basename(tree_path)}...", flush=True)
            ranks, struct = load_tree_ranks(tree_path, known, args.tree_root, group_map)
            print(f"  {len(ranks)} leaves in tree order, "
                  f"{len(struct['x'])} nodes, "
                  f"deepest tip {struct['maxDepth']:.3f} substitutions/site")
            loaded[tree_path] = (ranks, struct)
        return loaded[tree_path]

    tree_paths = args.trees or [args.tree] * len(paths)

    # Each family stays compressed in the page and is decoded only when chosen.
    payloads = {}
    for name, path, alphabet, tree_path in zip(names, paths, alphabets, tree_paths):
        print(f"[{name}]", flush=True)
        tree_ranks, tree_struct = (read_tree(tree_path) if tree_path and tree_path != "-"
                                   else (None, None))
        payloads[name] = build_payload(path, args.max_gap_payload, args.groups,
                                       descriptions, alphabet, tree_ranks, tree_struct,
                                       family_scores, args.lineages)

    title = args.title or (names[0] if len(names) == 1 else "topology alignments")
    html = HTML_TEMPLATE.replace("__TITLE__", title)
    html = html.replace("__DEFAULT_MAX_GAP__", repr(args.default_max_gap))
    html = html.replace("__MAX_GAP_PAYLOAD__", repr(args.max_gap_payload))
    # 'group_mean' is the old name for the median ranking; the page only carries
    # the new option value, so the old spelling is folded in here.
    default_sort = ("group_median" if args.default_group_sort == "group_mean"
                    else args.default_group_sort)
    html = html.replace("__DEFAULT_GROUP_SORT__", json.dumps(default_sort))
    html = html.replace("__PAYLOADS__", json.dumps(payloads))

    with open(args.out, "w", encoding="utf-8") as out:
        out.write(html)
    print(f"Wrote {args.out} ({os.path.getsize(args.out) / 1e6:.1f} MB, "
          f"{len(payloads)} alignment{'s' if len(payloads) > 1 else ''})")


HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Topology Map</title>
<style>
:root {
  color-scheme: light;
  --surface-0: #f4f3f0;
  --surface-1: #fcfcfb;
  --surface-2: #eceae5;
  --border: #d8d5ce;
  --text-primary: #0b0b0b;
  --text-secondary: #52514e;
  --text-muted: #86847d;
  --accent: #256abf;
  --hist: #7fb3e8;
  /* Membrane is green by request. The four states were searched as a set for the
     palest steps that still clear CVD separation on every pair in both modes. */
  --state-M: #63be86;
  --state-W: #1c7f4c;
  --state-I: #5f5d59;
  --state-O: #b5b2ac;
  --state-S: #c1688d;
  --state-B: #b39ddb;
  --state-P: #ef9a9a;
  /* Taylor residue colours, muted for print and for long reading. */
  --aa-A: #C7D787;
  --aa-C: #D7D787;
  --aa-D: #D78787;
  --aa-E: #D787A7;
  --aa-F: #87D7A7;
  --aa-G: #D7B787;
  --aa-H: #87A7D7;
  --aa-I: #A7D787;
  --aa-K: #A787D7;
  --aa-L: #97D787;
  --aa-M: #87D787;
  --aa-N: #C787D7;
  --aa-P: #D7C787;
  --aa-Q: #D787C7;
  --aa-R: #8787D7;
  --aa-S: #D79787;
  --aa-T: #D7A787;
  --aa-V: #B7D787;
  --aa-W: #87C7D7;
  --aa-Y: #87D7C7;
  --aa-other: #cfcdc6;
  --gap: #eeece7;
  --outlier: #d9534f;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    color-scheme: dark;
    --surface-0: #121211;
    --surface-1: #1a1a19;
    --surface-2: #232321;
    --border: #3a3a36;
    --text-primary: #ffffff;
    --text-secondary: #c3c2b7;
    --text-muted: #8d8b82;
    --gap: #262624;
    --hist: #6da7ec;
    --state-M: #6fd087;
    --state-W: #0a9645;
    --state-I: #d2cfc8;
    --state-O: #6f6d68;
    --state-S: #e885bf;
    --state-B: #9c86c9;
    --state-P: #e07b7b;
  /* Taylor residue colours, muted for print and for long reading. */
  --aa-A: #B5C86A;
  --aa-C: #C8C86A;
  --aa-D: #C86A6A;
  --aa-E: #C86A90;
  --aa-F: #6AC890;
  --aa-G: #C8A26A;
  --aa-H: #6A90C8;
  --aa-I: #90C86A;
  --aa-K: #906AC8;
  --aa-L: #7DC86A;
  --aa-M: #6AC86A;
  --aa-N: #B56AC8;
  --aa-P: #C8B56A;
  --aa-Q: #C86AB5;
  --aa-R: #6A6AC8;
  --aa-S: #C87D6A;
  --aa-T: #C8906A;
  --aa-V: #A2C86A;
  --aa-W: #6AB5C8;
  --aa-Y: #6AC8B5;
  --aa-other: #4a4a46;
    --outlier: #e08381;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --surface-0: #121211;
  --surface-1: #1a1a19;
  --surface-2: #232321;
  --border: #3a3a36;
  --text-primary: #ffffff;
  --text-secondary: #c3c2b7;
  --text-muted: #8d8b82;
  --gap: #262624;
  --hist: #6da7ec;
  --state-M: #6fd087;
  --state-W: #0a9645;
  --state-I: #d2cfc8;
  --state-O: #6f6d68;
  --state-S: #e885bf;
  --state-B: #9c86c9;
  --state-P: #e07b7b;
  /* Taylor residue colours, muted for print and for long reading. */
  --aa-A: #B5C86A;
  --aa-C: #C8C86A;
  --aa-D: #C86A6A;
  --aa-E: #C86A90;
  --aa-F: #6AC890;
  --aa-G: #C8A26A;
  --aa-H: #6A90C8;
  --aa-I: #90C86A;
  --aa-K: #906AC8;
  --aa-L: #7DC86A;
  --aa-M: #6AC86A;
  --aa-N: #B56AC8;
  --aa-P: #C8B56A;
  --aa-Q: #C86AB5;
  --aa-R: #6A6AC8;
  --aa-S: #C87D6A;
  --aa-T: #C8906A;
  --aa-V: #A2C86A;
  --aa-W: #6AB5C8;
  --aa-Y: #6AC8B5;
  --aa-other: #4a4a46;
  --outlier: #e08381;
}
* { box-sizing: border-box; }
body {
  margin: 0;
  padding: 24px 16px 64px;
  background: var(--surface-0);
  color: var(--text-primary);
  font: 14px/1.5 ui-sans-serif, system-ui, -apple-system, "Segoe UI", sans-serif;
}
.wrap { max-width: 1400px; margin: 0 auto; }
h1 { font-size: 20px; margin: 0 0 4px; letter-spacing: -0.01em; }
.sub { color: var(--text-secondary); font-size: 13px; margin: 0 0 20px; word-break: break-all; }
.panel {
  background: var(--surface-1);
  border: 1px solid var(--border);
  border-radius: 10px;
  padding: 16px;
  margin-bottom: 16px;
}
.tiles { display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 12px; margin-bottom: 16px; }
.tile { background: var(--surface-1); border: 1px solid var(--border); border-radius: 10px; padding: 12px 14px; }
.tile .k { font-size: 11px; text-transform: uppercase; letter-spacing: 0.06em; color: var(--text-muted); }
.tile .v { font-size: 22px; font-weight: 600; margin-top: 2px; font-variant-numeric: tabular-nums; }
.tile .n { font-size: 11px; color: var(--text-secondary); margin-top: 2px; }
.controls { display: flex; flex-wrap: wrap; gap: 20px; align-items: flex-end; }
[hidden] { display: none !important; }
.ctl { display: flex; flex-direction: column; gap: 6px; }
.ctl label { font-size: 11px; text-transform: uppercase; letter-spacing: 0.06em; color: var(--text-muted); }
.ctl input[type=range] { width: 220px; accent-color: var(--accent); }
.ctl select, .ctl input[type=number] {
  background: var(--surface-2); color: var(--text-primary);
  border: 1px solid var(--border); border-radius: 6px; padding: 6px 8px; font: inherit;
}
.ctl .val { font-variant-numeric: tabular-nums; font-size: 13px; color: var(--text-secondary); }
button {
  background: var(--surface-2); color: var(--text-primary); border: 1px solid var(--border);
  border-radius: 6px; padding: 7px 12px; font: inherit; cursor: pointer;
}
button:hover { border-color: var(--text-muted); }
button.primary { background: var(--accent); border-color: var(--accent); color: #fff; }
/* The three filter views are one control, so they are drawn as one segmented
   strip rather than as three loose buttons that happen to sit together. */
.seg { display: inline-flex; border: 1px solid var(--border); border-radius: 6px; overflow: hidden; }
.seg button { border: 0; border-radius: 0; padding: 7px 11px; background: var(--surface-2); }
.seg button + button { border-left: 1px solid var(--border); }
.seg button.on { background: var(--accent); color: #fff; }
.seg button.on .n { color: rgba(255,255,255,0.8); }
.seg .n { font-variant-numeric: tabular-nums; font-size: 11px; color: var(--text-muted); margin-left: 5px; }
.legend { display: flex; flex-wrap: wrap; gap: 14px; margin: 14px 0 0; font-size: 12px; color: var(--text-secondary); }
.legend span { display: inline-flex; align-items: center; gap: 6px; }
.sw { width: 12px; height: 12px; border-radius: 3px; display: inline-block; }
/* flex:1 with min-width:0 is what stops the box collapsing to the canvas's
   intrinsic pixel width, which is only as wide as the visible column count. */
.canvas-box {
  position: relative; overflow: hidden; flex: 1 1 auto; min-width: 0;
  border: 1px solid var(--border); border-radius: 8px; background: var(--surface-2);
}
.heat-wrap { display: flex; align-items: flex-start; gap: 8px; width: 100%; }
.side-box {
  position: relative; flex: 0 0 360px; width: 360px;
  border: 1px solid var(--border); border-radius: 8px; background: var(--surface-1);
}
@media (max-width: 1100px) { .side-box { flex: 0 0 180px; width: 180px; } }
.tree-box {
  position: relative; flex: 0 0 300px; width: 300px;
  border: 1px solid var(--border); border-radius: 8px; background: var(--surface-1);
}
.tree-box[hidden] { display: none; }
@media (max-width: 1100px) { .tree-box { flex: 0 0 190px; width: 190px; } }
.gutter { position: relative; width: 132px; flex: 0 0 132px; }
.gutter .g {
  position: absolute; left: 0; right: 0; font-size: 11px; line-height: 1.1;
  color: var(--text-secondary); white-space: nowrap; overflow: hidden; text-overflow: ellipsis;
  border-top: 1px solid var(--border); padding-top: 1px;
}
.gutter .g b { color: var(--text-primary); font-weight: 600; }
.gutter .g i { font-style: normal; color: var(--text-muted); }
.seps { position: absolute; inset: 0; pointer-events: none; }
.selbar { display: flex; flex-wrap: wrap; gap: 12px; align-items: center; margin-top: 12px;
  padding-top: 12px; border-top: 1px solid var(--border); font-size: 13px; }
.selbar .count { font-weight: 600; font-variant-numeric: tabular-nums; }
.selbar .status { color: var(--text-muted); font-size: 12px; }
.selbar .status.err { color: var(--outlier); }
.seps div { position: absolute; left: 0; right: 0; height: 1px; background: var(--surface-0); opacity: 0.85; }
canvas { display: block; image-rendering: pixelated; }
#consensus { display: block; width: 100%; height: 34px; }
#occupancy { display: block; width: 100%; height: 44px; }
#heat { width: 100%; }
.track-label { font-size: 11px; text-transform: uppercase; letter-spacing: 0.06em; color: var(--text-muted); margin: 14px 0 4px; }
table { width: 100%; border-collapse: collapse; font-size: 13px; }
th, td { text-align: left; padding: 7px 10px; border-bottom: 1px solid var(--border); }
th { font-size: 11px; text-transform: uppercase; letter-spacing: 0.06em; color: var(--text-muted); font-weight: 600; }
td.num { text-align: right; font-variant-numeric: tabular-nums; }
tr.flag td.num.score { color: var(--outlier); font-weight: 600; }
.mono { font-family: ui-monospace, "Cascadia Mono", Menlo, monospace; font-size: 12px; }
.tip {
  position: fixed; pointer-events: none; z-index: 10; display: none;
  background: var(--surface-1); border: 1px solid var(--border); border-radius: 8px;
  padding: 8px 10px; font-size: 12px; box-shadow: 0 4px 16px rgba(0,0,0,0.18); max-width: 380px;
}
.tip .t { font-weight: 600; margin-bottom: 3px; }
.tip .r { color: var(--text-secondary); }
.note { font-size: 12px; color: var(--text-secondary); margin-top: 10px; }
.top-row { display: flex; justify-content: space-between; align-items: flex-start; gap: 16px; flex-wrap: wrap; }
#busy { color: var(--text-muted); font-size: 12px; }
@media (max-width: 640px) { .ctl input[type=range] { width: 160px; } }
</style>
</head>
<body>
<div class="wrap">
  <div class="top-row">
    <div>
      <h1>Alignment map &mdash; <span id="famtitle">__TITLE__</span></h1>
      <p class="sub" id="subtitle">loading&hellip;</p>
    </div>
    <button id="theme">Toggle theme</button>
  </div>

  <div class="tiles" id="tiles"></div>

  <div class="panel">
    <div class="controls">
      <div class="ctl" id="famctl">
        <label for="family">Alignment</label>
        <select id="family"></select>
      </div>
      <div class="ctl">
        <label for="maxgap">Hide columns with more than</label>
        <input type="range" id="maxgap" min="0" max="100" step="1">
        <span class="val" id="maxgapval"></span>
      </div>
      <div class="ctl">
        <label for="sort">Row order</label>
        <select id="sort">
          <option value="score_desc">Helix Jaccard similarity, best first</option>
          <option value="score">Helix Jaccard similarity, worst first</option>
          <option value="group">Group, worst first within group</option>
          <option value="group_median">Group, groups by median score, best first</option>
          <option value="lineage">Lineage, worst first within lineage</option>
          <option value="lineage_median">Lineage, groups by median score, best first</option>
          <option value="tree">Tree order</option>
          <option value="file">File order</option>
          <option value="coverage">Coverage, lowest first</option>
        </select>
      </div>
      <div class="ctl" id="treectl">
        <label for="showtree">Tree panel</label>
        <select id="showtree">
          <option value="on">Show at left</option>
          <option value="off">Hide</option>
        </select>
      </div>
      <div class="ctl">
        <label for="helixpick">Side panel</label>
        <select id="helixpick"></select>
      </div>
      <div class="ctl">
        <label for="flagcut" id="flagcutlabel">Flag score below (%)</label>
        <input type="number" id="flagcut" min="0" max="100" step="1">
      </div>
      <div class="ctl">
        <label>Filter view</label>
        <div class="seg" id="rowview">
          <button data-view="all" class="on">All <span class="n"></span></button>
          <button data-view="kept">Retained <span class="n"></span></button>
          <button data-view="removed">Filtered out <span class="n"></span></button>
        </div>
      </div>
      <div class="ctl">
        <label for="mode">Mouse mode</label>
        <select id="mode">
          <option value="pan">Pan / zoom</option>
          <option value="single">Select sequences (click)</option>
          <option value="range">Select range (drag rows)</option>
        </select>
      </div>
      <div class="ctl">
        <label>&nbsp;</label>
        <button id="resetzoom">Reset zoom</button>
      </div>
      <div class="ctl">
        <label>&nbsp;</label>
        <button id="export" class="primary">Download flagged IDs</button>
      </div>
      <div class="ctl">
        <label>&nbsp;</label>
        <span id="busy"></span>
      </div>
    </div>
    <div class="legend" id="legend"></div>
    <div class="selbar">
      <span>Excluded from downstream: <span class="count" id="selcount">0</span></span>
      <button id="selclear">Clear this family</button>
      <button id="seldownload">Download list</button>
      <span class="status" id="selstatus"></span>
    </div>
  </div>

  <div class="panel">
    <div class="track-label">Column occupancy (bar height = fraction of sequences present)</div>
    <canvas id="occupancy"></canvas>
    <div class="track-label">Column consensus topology</div>
    <canvas id="consensus"></canvas>
    <div class="track-label">Alignment &mdash; one row per sequence, one pixel per visible column</div>
    <div class="heat-wrap">
      <div class="tree-box" id="treebox" hidden><canvas id="tree"></canvas></div>
      <div class="gutter" id="gutter"></div>
      <div class="canvas-box"><canvas id="heat"></canvas><div class="seps" id="seps"></div></div>
      <div class="side-box"><canvas id="side"></canvas></div>
    </div>
    <p class="note" id="heatnote"></p>
    <p class="note" id="treenote" hidden></p>
  </div>

  <div class="panel" id="grouppanel" hidden>
    <div class="track-label">Per-group summary</div>
    <table id="groups"><thead><tr>
      <th>Group</th><th class="num">Sequences</th><th class="num">Mean score</th><th class="num">Median score</th>
      <th class="num">Lowest</th><th class="num">Median coverage</th><th class="num">Below cutoff</th>
    </tr></thead><tbody></tbody></table>
  </div>

  <div class="panel">
    <div class="track-label">Helix Jaccard similarity distribution</div>
    <canvas id="hist" style="width:100%;height:120px"></canvas>
    <p class="note">Every sequence is scored against each consensus helix by Jaccard similarity,
      helix in both over helix in either, so agreeing that a column is loop earns nothing and a
      helix shifted out of register is charged twice, once for the consensus helix it fails to
      cover and once for the loop column it occupies instead. The score shown is the mean of those
      per-helix Jaccards, on the same scale from 0 to 1, and the side panel can show any single
      helix instead. Nothing is standardised or rescaled, so a value is read directly as how much
      of the helix bundle this sequence reproduces. Sequences in the left tail are the candidates
      to inspect.</p>
  </div>

  <div class="panel">
    <div class="track-label">Lowest-scoring sequences <span id="flagcount" style="text-transform:none;letter-spacing:0"></span></div>
    <table id="outliers"><thead><tr>
      <th>#</th><th>Accession</th><th id="grouphead" hidden>Group</th>
      <th class="num">Score</th><th class="num">Coverage</th>
      <th class="num">Residues</th><th>Topology composition</th>
    </tr></thead><tbody></tbody></table>
    <p class="note">Coverage is the fraction of the 30% occupancy columns where this sequence has a residue
      rather than a gap. A low score with high coverage means a genuinely divergent architecture;
      a low coverage means a fragment that simply does not span the family core.</p>
  </div>
</div>

<div class="tip" id="tip"></div>

<script>
const DEFAULT_MAX_GAP = __DEFAULT_MAX_GAP__;
const MAX_GAP_PAYLOAD = __MAX_GAP_PAYLOAD__;
const DEFAULT_GROUP_SORT = __DEFAULT_GROUP_SORT__;
const PAYLOADS = __PAYLOADS__;   // name -> gzipped, base64 payload

const STATE_COLORS = { M: "--state-M", W: "--state-W", I: "--state-I",
                       O: "--state-O", S: "--state-S", B: "--state-B", P: "--state-P" };
const STATE_NAMES = { M: "Membrane, inside to outside", W: "Membrane, outside to inside",
                      I: "Inside (cytoplasmic)", O: "Outside (extracellular)",
                      S: "Signal peptide", B: "Beta barrel", P: "Periplasm" };

/* Taylor's residue colours, muted. The hues carry the chemistry, so residues
   that behave alike sit close together and a column of conservative
   substitutions still reads as one band. */
const AA_COLORS = Object.fromEntries(
  "ACDEFGHIKLMNPQRSTVWY".split("").map(a => [a, `--aa-${a}`]));
const AA_NAMES = { A: "Alanine", C: "Cysteine", D: "Aspartate", E: "Glutamate",
  F: "Phenylalanine", G: "Glycine", H: "Histidine", I: "Isoleucine", K: "Lysine",
  L: "Leucine", M: "Methionine", N: "Asparagine", P: "Proline", Q: "Glutamine",
  R: "Arginine", S: "Serine", T: "Threonine", V: "Valine", W: "Tryptophan",
  Y: "Tyrosine", X: "Any residue", B: "Asn or Asp", Z: "Gln or Glu", U: "Selenocysteine" };

/* Topology and protein alignments share the letters B, M, O, P, S and W, so
   the palette is chosen by what the payload says it is, never per letter. */
function isProtein() { return !!D && D.alphabet === "protein"; }
function activePalette() { return isProtein() ? AA_COLORS : STATE_COLORS; }
function colorVar(label) {
  return activePalette()[label] || (isProtein() ? "--aa-other" : "--gap");
}
function labelName(label) {
  return (isProtein() ? AA_NAMES : STATE_NAMES)[label] || label;
}

let D = null;             // decoded payload of the family on screen
let family = null;        // its name
let selected = new Set(); // accessions of this family marked for exclusion
let serverOK = false;     // true when served by serve_topology_viewer.py
let colGap = null;        // gap fraction per payload column
let visible = null;       // indices into payload columns that are shown
let scores = null, coverage = null, residues = null, comps = null;
let order = null;         // row draw order
let rgb = {};             // label -> [r,g,b], refreshed on theme change
let hasGroups = false;    // true when the builder was given --groups
let hasLineages = false;  // true when the builder was given --lineages
let hasTree = false;      // true when the builder was given --tree
let scoreMax = 1;         // every score on the page is a Jaccard, so always 1
/* Row order must mean the same thing in both panels, so the protein alignment is
   ordered by the topology Jaccard of the same accession rather than by its own
   identity to consensus. Sequence similarity and helix architecture disagree often
   enough that ordering each panel by its own number would put the same sequence in
   two different places and make the panels uncomparable, which is the whole reason
   they are in one page. These hold the topology payload and the scores taken from
   it, computed once over its 30% occupancy columns. */
let decodedCache = new Map();   // payload name -> decoded payload
let topoScoreById = new Map();
let topoScoreKey = null;        // which topology payload and gap the scores are from
let orderScores = null;   // the per-row values computeOrder and groupMedians rank on
let blocks = [], helixWindows = [], helixScores = [];  // consensus helices

/* Which sequences the map is showing: every one, only those the cutoff keeps, or
   only those it removes. The three are one alignment seen three ways rather than
   three alignments, so the view changes which rows `order` carries and nothing
   else -- the columns, the consensus, the scores and the cutoff are identical in
   all three, and a sequence sits at the same score in whichever view shows it.
   'all' additionally tints the removed rows, so the before picture and the two
   after pictures can be read against each other without switching. */
let rowView = "all";

/* A sequence is removed when it scores below the cutoff box. This is the one
   definition on the page: the tint, the two filtered views, the counts on the
   buttons and the flag list all ask this function, so they cannot disagree. */
function isRemoved(s) { return scores[s] < cutValue(); }

/* Rows currently drawn. Every drawing path counts rows through this rather than
   through D.nSeq, because in a filtered view the map is shorter than the family. */
function nRows() { return order ? order.length : 0; }

function cssRGB(varName) {
  const v = getComputedStyle(document.documentElement).getPropertyValue(varName).trim();
  const m = v.match(/^#([0-9a-f]{6})$/i);
  if (m) { const n = parseInt(m[1], 16); return [n >> 16 & 255, n >> 8 & 255, n & 255]; }
  const p = v.match(/\d+/g); return p ? p.slice(0, 3).map(Number) : [128, 128, 128];
}
function refreshColors() {
  rgb = {};
  for (const [label, varName] of Object.entries(activePalette())) rgb[label] = cssRGB(varName);
  rgb["-"] = cssRGB("--gap");
  if (D) for (const label of D.labels) if (!rgb[label]) rgb[label] = cssRGB(colorVar(label));
}

async function decode(b64) {
  const text = atob(b64), bytes = new Uint8Array(text.length);
  for (let i = 0; i < text.length; i++) bytes[i] = text.charCodeAt(i);
  const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream("gzip"));
  return JSON.parse(await new Response(stream).text());
}

function fmtInt(n) { return n.toLocaleString(); }
function pct(x) { return (100 * x).toFixed(1) + "%"; }

/* Every score on the page is a Jaccard, so it is shown as a percentage of the
   consensus helix bundle it reproduces. Values stay fractions internally -- the
   cutoff comparisons, the ramps and the ordering all work in 0 to 1 -- and only the
   printing changes, so there is one place to get the rounding right. */
function fmtScore(v) { return (100 * v).toFixed(1) + "%"; }

/* The cutoff box holds a percentage because that is what it displays; everything
   that compares against it wants the fraction. */
function cutValue() {
  const v = parseFloat(document.getElementById("flagcut").value);
  return isFinite(v) ? v / 100 : 0;
}

function computeVisible(maxGap) {
  visible = [];
  for (let c = 0; c < D.nCol; c++) if (colGap[c] <= maxGap + 1e-9) visible.push(c);
}

/* Score every sequence against the column consensus over the visible columns.
   Walks the packed matrix once; ~6M character reads for a family this size. */
/* A protein alignment has no helices to be in or out of register with, so a
   row is scored by how often it carries its column's consensus residue, over
   the columns it actually occupies. That lands on 0 to 1 like a Jaccard, so
   the cutoff box, the histogram and the flag list keep working unchanged. */
function computeIdentityScores() {
  const n = D.nSeq, w = D.nCol, m = D.matrix;
  scores = new Float32Array(n); coverage = new Float32Array(n); residues = new Int32Array(n);
  comps = [];
  const cons = visible.map(c => D.columns[c].cons);
  const vis = visible, vn = vis.length;
  for (let s = 0; s < n; s++) {
    const base = s * w;
    let present = 0, same = 0;
    const comp = {};
    for (let k = 0; k < vn; k++) {
      const ch = m[base + vis[k]];
      if (ch === "-") continue;
      present++;
      comp[ch] = (comp[ch] || 0) + 1;
      if (ch === cons[k]) same++;
    }
    residues[s] = present;
    coverage[s] = vn ? present / vn : 0;
    scores[s] = present ? same / present : 0;
    comps.push(comp);
  }
  blocks = []; helixWindows = []; helixScores = [];
  scoreMax = 1;
  refreshHelixPicker();
  syncCutoffScale();
}

/* Runs of consensus helix at least eight columns wide, as index pairs into
   whichever column list the mask was built from. */
function helixBlocksFor(consIsHelix) {
  const vn = consIsHelix.length, blocks = [];
  for (let k = 0, start = -1; k <= vn; k++) {
    const on = k < vn && consIsHelix[k];
    if (on && start < 0) start = k;
    else if (!on && start >= 0) { if (k - start >= 8) blocks.push([start, k]); start = -1; }
  }
  return blocks;
}

/* Every score on this page is computed over the columns that clear 30% occupancy
   and nothing else. The occupancy slider is a display control: it changes which
   columns are drawn, never what a sequence scores, so a number here means the same
   thing as the same number in topology_helix_scores.tsv and on the tree pages. */
const SCORE_MAX_GAP = 0.70;
function scoreColumnsOf(P) {
  const cols = [];
  for (let c = 0; c < P.nCol; c++) if (P.columns[c].gap <= SCORE_MAX_GAP + 1e-9) cols.push(c);
  return cols;
}

/* The helix scoring, written against a payload and a column list rather than
   against the globals, so the protein panel can borrow the topology payload's
   scores without a second copy of the arithmetic living somewhere else. */
function helixScoresFor(P, vis) {
  const n = P.nSeq, w = P.nCol, m = P.matrix;
  const scores = new Float32Array(n), coverage = new Float32Array(n);
  const residues = new Int32Array(n), comps = [];
  const cons = vis.map(c => P.columns[c].cons);
  const consIsHelix = cons.map(c => c === "M" || c === "W");
  const vn = vis.length;

  /* Consensus helix blocks, and a window per helix running to the midpoint of
     the loop on either side, so a helix that has slipped out of register is
     charged inside its own window rather than its neighbour's. */
  const blocks = helixBlocksFor(consIsHelix);
  const edges = [0];
  for (let b = 0; b + 1 < blocks.length; b++)
    edges.push(Math.floor((blocks[b][1] + blocks[b + 1][0]) / 2));
  edges.push(vn);
  const helixWindows = blocks.map((_, b) => [edges[b], edges[b + 1]]);
  const helixScores = helixWindows.map(() => new Float32Array(n));
  for (let s = 0; s < n; s++) {
    const base = s * w;
    let present = 0, both = 0, helixOnly = 0, consensusOnly = 0;
    const comp = {};
    for (let k = 0; k < vn; k++) {
      const ch = m[base + vis[k]];
      if (ch === "-") continue;
      present++;
      comp[ch] = (comp[ch] || 0) + 1;
      // M and W are different states, inside-to-outside and outside-to-inside, so a
      // helix running the wrong way is off register just as a helix in a loop is.
      const isHelix = ch === "M" || ch === "W";
      if (isHelix && ch === cons[k]) both++;
      else if (isHelix) helixOnly++;              // loop column, or the wrong direction
      else if (consIsHelix[k]) consensusOnly++;   // family helix left uncovered
    }
    residues[s] = present;
    coverage[s] = vn ? present / vn : 0;
    comps.push(comp);

    let jaccardSum = 0;
    for (let b = 0; b < helixWindows.length; b++) {
      const [a, e] = helixWindows[b];
      let hb = 0, ho = 0, hc = 0;
      for (let k = a; k < e; k++) {
        const ch = m[base + vis[k]];
        if (ch === "-") continue;
        const isHelix = ch === "M" || ch === "W";
        if (isHelix && ch === cons[k]) hb++;
        else if (isHelix) ho++;
        else if (consIsHelix[k]) hc++;
      }
      const u = hb + ho + hc;
      const j = u ? hb / u : 0;
      helixScores[b][s] = j;
      jaccardSum += j;
    }

    /* The score is the plain mean of the sequence's Jaccard against each consensus
       helix, so it is a similarity on the same 0 to 1 scale as the per-helix numbers
       beside it. The Euclidean distance of those Jaccards from the origin was tried
       and dropped: squaring lets the helices that match well dominate, so a sequence
       with two intact helices and five collapsed scores about 0.51 under it against
       0.27 under the mean, and the mean is the reading that places such a sequence
       where it belongs. */
    if (helixWindows.length) {
      scores[s] = jaccardSum / helixWindows.length;
    } else {
      const union = both + helixOnly + consensusOnly;   // no helix blocks resolved
      scores[s] = union ? both / union : 0;
    }
  }
  return { scores, coverage, residues, comps, blocks, helixWindows, helixScores };
}

function computeScores() {
  if (isProtein()) { computeIdentityScores(); return; }
  const r = helixScoresFor(D, scoreColumnsOf(D));
  scores = r.scores;
  coverage = r.coverage; residues = r.residues; comps = r.comps;
  helixWindows = r.helixWindows; helixScores = r.helixScores;
  /* The consensus track is drawn against the visible columns, so its helix
     brackets are found in that frame. Only the drawing uses them; every number
     comes from the fixed frame above. */
  blocks = helixBlocksFor(visible.map(c => {
    const ch = D.columns[c].cons; return ch === "M" || ch === "W";
  }));
  scoreMax = 1;   // a mean of Jaccards, like every other number on the page
  refreshHelixPicker();
  syncCutoffScale();
}

/* Every score on the page is a Jaccard on 0 to 1, so the box needs no rescaling
   when the resolved helix count changes. The function is kept so the one caller
   does not have to know that. */
function syncCutoffScale() {
  const box = document.getElementById("flagcut");
  const previous = parseFloat(box.dataset.scale || "1");
  box.max = (100 * scoreMax).toFixed(0);
  box.step = "1";
  // 80% is the cutoff the filter is being applied at, so the page opens on the
  // split that is actually in force rather than on a nearby illustrative one.
  if (!box.dataset.scale) box.value = (80 * scoreMax).toFixed(0);
  box.dataset.scale = String(scoreMax);
}

/* The side panel can show any one helix, or the whole-bundle score. The mean is
   the default because it is the number the filter is argued from; TM3 remains a
   click away and is where the class C intruders break. */
function refreshHelixPicker() {
  const pick = document.getElementById("helixpick");
  const previous = pick.value;
  pick.innerHTML = helixWindows.map((_, b) => `<option value="${b}">TM${b + 1}</option>`)
    .join("") + '<option value="all">All helices (mean Jaccard)</option>';
  if (previous && [...pick.options].some(o => o.value === previous)) pick.value = previous;
  else pick.value = "all";
  // Nothing to choose between when no consensus helices were resolved.
  pick.closest(".ctl").hidden = helixWindows.length === 0;
}

function sidePanelValues() {
  const pick = document.getElementById("helixpick").value;
  if (pick === "all" || !helixScores.length)
    return { values: scores, label: "mean Jaccard", max: 1 };
  const b = Math.min(helixScores.length - 1, parseInt(pick, 10));
  return { values: helixScores[b], label: `TM${b + 1}`, max: 1 };
}

/* The knee of the sorted mean-Jaccard curve.

   Sorting the score turns it into a quantile function, so the derivative of that
   curve is one over the density and the row where it peaks is where the curve
   falls away from the body of the distribution. Measured across families the knee
   sits in the bottom few percent, matching what a smooth unimodal distribution
   produces, so it marks the start of the tail rather than a boundary between two
   populations. It is drawn to show where the tail begins, not to be cut on.

   Always computed from the mean Jaccard, whichever helix the panel is showing,
   because the whole-bundle score is what the rows are ranked by. */
const kneeCache = new Map();
function kneeScore() {
  if (kneeCache.has(family)) return kneeCache.get(family);
  let result = null;
  const n = scores.length;
  if (n >= 50) {
    const v = Float64Array.from(scores).sort();          // ascending
    const w = Math.max(11, Math.floor(0.02 * n) | 1);
    const m = n - w + 1;
    if (m > 4) {
      const s = new Float64Array(m);
      let sum = 0;
      for (let i = 0; i < w; i++) sum += v[i];
      s[0] = sum / w;
      for (let i = 1; i < m; i++) { sum += v[i + w - 1] - v[i - 1]; s[i] = sum / w; }
      const dx = 1 / (m - 1);
      const half = Math.max(2, Math.floor(m / 2));
      let best = -Infinity, bi = 1;
      for (let i = 1; i < half; i++) {
        const d = (s[i + 1] - s[i - 1]) / (2 * dx);
        if (d > best) { best = d; bi = i; }
      }
      result = s[bi];
    }
  }
  kneeCache.set(family, result);
  return result;
}

/* The row the knee falls at, or null when the rows are not ranked by score and a
   single horizontal line would therefore not separate anything. */
function kneeRow() {
  const k = kneeScore();
  if (k === null || !order) return null;
  let previous = Infinity, crossing = -1;
  for (let row = 0; row < order.length; row++) {
    const v = scores[order[row]];
    if (v > previous + 1e-9) return null;               // not sorted by score
    previous = v;
    if (crossing < 0 && v < k) crossing = row;
  }
  return crossing > 0 ? crossing : null;
}

/* The row the cutoff falls at, on the same terms as kneeRow: null unless the rows
   are ranked by score, since a horizontal line across rows in tree or group order
   would cut through sequences on both sides of the threshold and separate nothing.
   This is the boundary between what the filter keeps and what it removes, so it is
   the one line the side panel draws. */
function cutoffRow() {
  const cut = cutValue();
  if (!order || !scores) return null;
  let previous = Infinity, crossing = -1;
  for (let row = 0; row < order.length; row++) {
    const v = scores[order[row]];
    if (v > previous + 1e-9) return null;               // not sorted by score
    previous = v;
    if (crossing < 0 && v < cut) crossing = row;
  }
  return crossing > 0 ? crossing : null;
}

async function payloadNamed(name) {
  if (name === family && D) return D;
  if (!decodedCache.has(name)) decodedCache.set(name, await decode(PAYLOADS[name]));
  return decodedCache.get(name);
}

/* Which topology alignment this protein alignment should borrow its order from.
   A page can hold several dates at once, and each date's protein alignment must
   be ranked by its own date's topology, not by whichever topology payload
   happens to come first. The pairing is decided by how many accessions the two
   actually share, which needs no naming convention to hold. */
async function topologyPartner() {
  let best = null, bestShared = 0;
  const mine = new Set(D.ids);
  for (const name of Object.keys(PAYLOADS)) {
    const candidate = await payloadNamed(name);
    if (candidate.alphabet !== "topology") continue;
    let shared = 0;
    for (const id of candidate.ids) if (mine.has(id)) shared++;
    if (shared > bestShared) { best = { name, payload: candidate }; bestShared = shared; }
  }
  return bestShared ? best : null;
}

/* The partner topology alignment scored on its own 30% occupancy columns and keyed
   by accession. Returns quietly when the page holds no topology alignment sharing
   sequences with this one. */
async function ensureTopoScores() {
  const partner = await topologyPartner();
  if (!partner) { topoScoreById = new Map(); topoScoreKey = null; return; }
  const key = partner.name;
  if (topoScoreKey === key && topoScoreById.size) return;
  const P = partner.payload;
  const r = helixScoresFor(P, scoreColumnsOf(P));
  topoScoreById = new Map();
  for (let i = 0; i < P.nSeq; i++) topoScoreById.set(P.ids[i], r.scores[i]);
  topoScoreKey = key;
}

/* What the rows are ranked by. The topology panel ranks on its own score; the
   protein panel borrows the topology score for the same accession, and falls back
   to its own only for an accession the topology alignment does not contain. */
/* The sort options name whichever score is actually doing the ranking, so the
   protein panel does not appear to be sorted by identity when it is not. */
function refreshSortLabels() {
  const borrowed = isProtein() && topoScoreById.size > 0;
  const name = borrowed ? "Topology Jaccard similarity"
             : isProtein() ? "Identity to consensus" : "Helix Jaccard similarity";
  const sel = document.getElementById("sort");
  for (const opt of sel.options) {
    if (opt.value === "score_desc") opt.textContent = `${name}, best first`;
    else if (opt.value === "score") opt.textContent = `${name}, worst first`;
    else if (opt.value === "group_median") opt.textContent = `Group, groups by median ${borrowed ? "topology " : ""}score, best first`;
    else if (opt.value === "lineage_median") opt.textContent = `Lineage, groups by median ${borrowed ? "topology " : ""}score, best first`;
  }
}

function computeOrderScores() {
  if (!isProtein() || !topoScoreById.size) { orderScores = scores; return; }
  const v = new Float32Array(D.nSeq);
  for (let i = 0; i < D.nSeq; i++) {
    const t = topoScoreById.get(D.ids[i]);
    v[i] = t === undefined ? scores[i] : t;
  }
  orderScores = v;
}

function computeOrder(mode) {
  const idx = Array.from({ length: D.nSeq }, (_, i) => i);
  const rank = orderScores || scores;
  if (mode === "score") idx.sort((a, b) => rank[a] - rank[b] || a - b);
  else if (mode === "score_desc") idx.sort((a, b) => rank[b] - rank[a] || a - b);
  else if (mode === "coverage") idx.sort((a, b) => coverage[a] - coverage[b] || a - b);
  else if (mode === "tree" && hasTree) idx.sort((a, b) => D.treeRank[a] - D.treeRank[b] || a - b);
  else if (mode === "group" && hasGroups)
    idx.sort((a, b) => D.groupOf[a] - D.groupOf[b] || rank[a] - rank[b] || a - b);
  else if ((mode === "group_median" || mode === "group_mean") && hasGroups) {
    // Groups ranked by their mean score, best group on top, best sequence first within each.
    const medians = groupMedians("group"), grank = new Array(medians.length);
    medians.map((m, g) => [m, g]).sort((a, b) => b[0] - a[0] || a[1] - b[1])
      .forEach(([, g], r) => { grank[g] = r; });
    idx.sort((a, b) => grank[D.groupOf[a]] - grank[D.groupOf[b]] || rank[b] - rank[a] || a - b);
  }
  else if (mode === "lineage" && hasLineages)
    idx.sort((a, b) => D.lineageOf[a] - D.lineageOf[b] || rank[a] - rank[b] || a - b);
  else if (mode === "lineage_median" && hasLineages) {
    // Same ranking as group_median, but by lineage instead of family, and
    // regardless of which family each sequence belongs to.
    const medians = groupMedians("lineage"), grank = new Array(medians.length);
    medians.map((m, g) => [m, g]).sort((a, b) => b[0] - a[0] || a[1] - b[1])
      .forEach(([, g], r) => { grank[g] = r; });
    idx.sort((a, b) => grank[D.lineageOf[a]] - grank[D.lineageOf[b]] || rank[b] - rank[a] || a - b);
  }
  /* Filtering after the sort, not before, so a row keeps the position the chosen
     sort gives it and the filtered views are the full view with rows taken out
     rather than a differently ordered alignment. */
  order = rowView === "all" ? idx
        : rowView === "kept" ? idx.filter(s => !isRemoved(s))
        : idx.filter(s => isRemoved(s));
  refreshViewCounts();
}

/* The counts live on the buttons themselves, so the split is readable without
   switching views. They are recomputed whenever the cutoff or the family moves. */
function refreshViewCounts() {
  if (!scores) return;
  let removed = 0;
  for (let s = 0; s < D.nSeq; s++) if (isRemoved(s)) removed++;
  const n = { all: D.nSeq, kept: D.nSeq - removed, removed };
  for (const button of document.querySelectorAll("#rowview button"))
    button.querySelector(".n").textContent = fmtInt(n[button.dataset.view]);
}

/* A group's standing is its median score, not its mean. A handful of badly
   aligning members drags a mean down far enough to reorder the groups around
   it, which makes the ranking a statement about the outliers rather than about
   the group, and the outliers are already visible in their own right. Works on
   either axis: 'group' is GPCR family, 'lineage' is taxonomic lineage. */
function groupMedians(axis) {
  const names = axis === "lineage" ? D.lineageNames : D.groupNames;
  const of = axis === "lineage" ? D.lineageOf : D.groupOf;
  const rank = orderScores || scores;
  const per = names.map(() => []);
  for (let i = 0; i < D.nSeq; i++) per[of[i]].push(rank[i]);
  return per.map(v => v.length ? medianOf(v) : 0);
}

/* Where each block starts and ends in the current row order, on whichever axis
   the sort mode names. Only meaningful when the rows are grouped that way, so
   the gutter is hidden otherwise. */
function groupBlocks() {
  const mode = document.getElementById("sort").value;
  const byLineage = mode.startsWith("lineage");
  if (byLineage ? !hasLineages : !(hasGroups && mode.startsWith("group"))) return [];
  const names = byLineage ? D.lineageNames : D.groupNames;
  const of = byLineage ? D.lineageOf : D.groupOf;
  const blocks = [];
  let start = 0;
  for (let row = 1; row <= nRows(); row++) {
    if (row === nRows() || of[order[row]] !== of[order[start]]) {
      blocks.push({ name: names[of[order[start]]], start, end: row });
      start = row;
    }
  }
  return blocks;
}

function drawGutter() {
  const gutter = document.getElementById("gutter"), seps = document.getElementById("seps");
  const blocks = groupBlocks();
  if (!blocks.length) {
    gutter.innerHTML = ""; seps.innerHTML = "";
    gutter.style.display = "none";
    return;
  }
  gutter.style.display = "block";
  const height = VIEW_HEIGHT;
  const top0 = view ? view.y : 0, rows = view ? view.h : nRows();
  // A label only fits if its block is tall enough; smaller blocks keep the
  // separator line but drop the text rather than overprinting a neighbour.
  gutter.innerHTML = blocks.map(b => {
    const top = ((b.start - top0) / rows) * height;
    const h = ((b.end - b.start) / rows) * height;
    if (top + h < 0 || top > height) return "";
    const members = order.slice(b.start, b.end).map(i => scores[i]);
    const stat = medianOf(members);
    const statName = "median";
    const text = h >= 11
      ? `<b>${b.name}</b> <i>${b.end - b.start} · ${fmtScore(stat)}</i>` : "";
    return `<div class="g" style="top:${top}px;height:${Math.max(1, h)}px"
      title="${b.name} — ${b.end - b.start} sequences, ${statName} score ${fmtScore(stat)}">${text}</div>`;
  }).join("");
  seps.innerHTML = blocks.slice(1).map(b =>
    `<div style="top:${((b.start - top0) / rows) * height}px"></div>`).join("");
}

function medianOf(values) {
  if (!values.length) return 0;
  const s = [...values].sort((a, b) => a - b);
  return s[Math.floor(s.length / 2)];
}

function drawGroupTable() {
  if (!hasGroups) return;
  const cut = cutValue();
  const rows = D.groupNames.map((name, g) => {
    const members = [];
    for (let i = 0; i < D.nSeq; i++) if (D.groupOf[i] === g) members.push(i);
    const sc = members.map(i => scores[i]);
    return {
      name, n: members.length,
      mean: sc.reduce((a, v) => a + v, 0) / sc.length, median: medianOf(sc), lowest: Math.min(...sc),
      cov: medianOf(members.map(i => coverage[i])),
      below: sc.filter(v => v < cut).length,
    };
  }).sort((a, b) => b.median - a.median);
  document.getElementById("groups").querySelector("tbody").innerHTML = rows.map(r => `
    <tr class="${r.below ? "flag" : ""}">
      <td>${r.name}</td><td class="num">${fmtInt(r.n)}</td>
      <td class="num">${fmtScore(r.mean)}</td>
      <td class="num">${fmtScore(r.median)}</td>
      <td class="num score">${fmtScore(r.lowest)}</td>
      <td class="num">${pct(r.cov)}</td>
      <td class="num">${r.below || "—"}</td></tr>`).join("");
}

/* The map is painted once at its natural size, one pixel per cell, into an
   offscreen canvas. What the page shows is a rectangle of that buffer blitted
   to the visible canvas, so zooming and panning only move the rectangle and
   never repaint the data. `view` is that rectangle, in natural cell units.
   Browsers cap a canvas at a few tens of thousands of pixels a side, which
   class A and Olfactory exceed by rows alone, so the buffer is a stack of
   horizontal tiles; `base` keeps the full width and height for the viewport. */
let base = null;
let view = null;
const VIEW_HEIGHT = 620;
const TILE_ROWS = 8192;
/* How far a removed row is pulled towards the flag colour in the combined view.
   Enough that the band of removed sequences is obvious when the whole family is
   zoomed out, little enough that M and W are still told apart inside it. */
const TINT = 0.45;

function drawHeat(keepView) {
  const w = Math.max(1, visible.length), h = nRows();
  const m = D.matrix, stride = D.nCol;
  const tintRemoved = rowView === "all";
  const cut = cssRGB("--outlier");
  /* A view can legitimately be empty -- a family where the cutoff removes nothing
     has nothing to show under "Filtered out" -- so that is said plainly rather
     than left as a blank canvas or allowed to divide by a zero row count. */
  if (!h) {
    base = null; view = null;
    for (const id of ["heat", "side"]) {
      const cv = document.getElementById(id);
      cv.getContext("2d").clearRect(0, 0, cv.width, cv.height);
    }
    document.getElementById("heatnote").textContent =
      rowView === "removed"
        ? `No sequence of ${family} scores below ${fmtScore(cutValue())}; nothing is filtered out.`
        : `Every sequence of ${family} scores below ${fmtScore(cutValue())}; nothing is retained.`;
    document.getElementById("treebox").hidden = true;
    return;
  }
  base = { width: w, height: h, tiles: [] };
  for (let y0 = 0; y0 < h; y0 += TILE_ROWS) {
    const rows = Math.min(TILE_ROWS, h - y0);
    const tile = document.createElement("canvas");
    tile.width = w; tile.height = rows;
    const tctx = tile.getContext("2d");
    const img = tctx.createImageData(w, rows);
    const data = img.data;
    for (let r = 0; r < rows; r++) {
      const s = order[y0 + r], from = s * stride, out = r * w * 4;
      /* In the combined view a removed row is washed towards the flag colour, so
         it stays legible as topology while reading as struck out. The tint is a
         blend rather than a replacement because the point of that view is to see
         what the removed sequences look like, not merely where they are. */
      const wash = tintRemoved && isRemoved(s);
      for (let k = 0; k < visible.length; k++) {
        const c = rgb[m[from + visible[k]]] || rgb["-"];
        const o = out + k * 4;
        if (wash) {
          data[o] = c[0] + (cut[0] - c[0]) * TINT;
          data[o + 1] = c[1] + (cut[1] - c[1]) * TINT;
          data[o + 2] = c[2] + (cut[2] - c[2]) * TINT;
        } else {
          data[o] = c[0]; data[o + 1] = c[1]; data[o + 2] = c[2];
        }
        data[o + 3] = 255;
      }
    }
    tctx.putImageData(img, 0, 0);
    base.tiles.push({ y0, rows, canvas: tile });
  }
  if (!keepView || !view) view = { x: 0, y: 0, w: base.width, h: base.height };
  clampView();
  renderView();
}

function clampView() {
  if (!view || !base) return;
  view.w = Math.min(base.width, Math.max(2, view.w));
  view.h = Math.min(base.height, Math.max(2, view.h));
  view.x = Math.min(base.width - view.w, Math.max(0, view.x));
  view.y = Math.min(base.height - view.h, Math.max(0, view.h ? view.y : 0));
  view.y = Math.min(base.height - view.h, Math.max(0, view.y));
}

function renderView() {
  const cv = document.getElementById("heat");
  if (!base || !view) return;
  syncTreeBox();
  const dpr = window.devicePixelRatio || 1;
  const cssW = cv.clientWidth || 900;
  cv.style.height = VIEW_HEIGHT + "px";
  cv.width = Math.round(cssW * dpr); cv.height = Math.round(VIEW_HEIGHT * dpr);
  const ctx = cv.getContext("2d");
  ctx.imageSmoothingEnabled = false;
  ctx.clearRect(0, 0, cv.width, cv.height);
  const yScale = cv.height / view.h;
  for (const t of base.tiles) {
    const top = Math.max(view.y, t.y0), bottom = Math.min(view.y + view.h, t.y0 + t.rows);
    if (bottom <= top) continue;
    ctx.drawImage(t.canvas, view.x, top - t.y0, view.w, bottom - top,
                  0, (top - view.y) * yScale, cv.width, (bottom - top) * yScale);
  }

  /* Once a cell is a few pixels across it is drawn as a cell: a separating grid
     so the squares read individually, and at larger sizes the topology letter
     itself, so a single sequence can be read residue by residue. */
  const pxCol = cssW / view.w, pxRow = VIEW_HEIGHT / view.h;
  if (pxCol >= 6 && pxRow >= 6) {
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    const first = Math.max(0, Math.floor(view.x)), lastC = Math.min(visible.length, Math.ceil(view.x + view.w));
    const firstR = Math.max(0, Math.floor(view.y)), lastR = Math.min(nRows(), Math.ceil(view.y + view.h));
    const grid = cssRGB("--surface-0");
    ctx.strokeStyle = `rgba(${grid.join(",")},0.85)`;
    ctx.lineWidth = 1;
    ctx.beginPath();
    for (let c = first; c <= lastC; c++) {
      const x = Math.round((c - view.x) * pxCol) + 0.5;
      ctx.moveTo(x, 0); ctx.lineTo(x, VIEW_HEIGHT);
    }
    for (let r = firstR; r <= lastR; r++) {
      const y = Math.round((r - view.y) * pxRow) + 0.5;
      ctx.moveTo(0, y); ctx.lineTo(cssW, y);
    }
    ctx.stroke();

    if (pxCol >= 11 && pxRow >= 11) {
      const size = Math.min(pxCol, pxRow) * 0.62;
      ctx.font = `600 ${size.toFixed(1)}px ui-monospace, "Cascadia Mono", Menlo, monospace`;
      ctx.textAlign = "center"; ctx.textBaseline = "middle";
      ctx.fillStyle = cssRGB("--text-primary").length
        ? `rgba(${cssRGB("--text-primary").join(",")},0.92)` : "#000";
      for (let r = firstR; r < lastR; r++) {
        const s = order[r], from = s * D.nCol;
        const y = (r - view.y) * pxRow + pxRow / 2;
        for (let c = first; c < lastC; c++) {
          const ch = D.matrix[from + visible[c]];
          if (ch === "-") continue;
          ctx.fillText(ch, (c - view.x) * pxCol + pxCol / 2, y);
        }
      }
      ctx.textAlign = "start"; ctx.textBaseline = "alphabetic";
    }
    ctx.setTransform(1, 0, 0, 1, 0, 0);
  }

  // Selected rows are tinted, one band per run of consecutive selected rows.
  if (selected.size) {
    const out = cssRGB("--outlier"), yScale = cv.height / view.h;
    ctx.fillStyle = `rgba(${out.join(",")},0.55)`;
    const firstR = Math.max(0, Math.floor(view.y)), lastR = Math.min(nRows(), Math.ceil(view.y + view.h));
    for (let r = firstR; r < lastR; r++) {
      if (!selected.has(D.ids[order[r]])) continue;
      let e = r + 1;
      while (e < lastR && selected.has(D.ids[order[e]])) e++;
      const y = (r - view.y) * yScale;
      ctx.fillRect(0, y, cv.width, Math.max(1, (e - r) * yScale));
      r = e;
    }
  }

  // The gutter is shown or hidden first, since it moves the map and the tracks align to the map.
  drawGutter();
  renderTree();
  renderSide();
  drawTracks();

  const rowsPerPixel = view.h / VIEW_HEIGHT;
  document.getElementById("heatnote").textContent =
    `Showing sequences ${fmtInt(Math.floor(view.y) + 1)}–${fmtInt(Math.ceil(view.y + view.h))} ` +
    `of ${fmtInt(nRows())} and columns ${fmtInt(Math.floor(view.x) + 1)}–` +
    `${fmtInt(Math.ceil(view.x + view.w))} of ${fmtInt(visible.length)}, ` +
    `${rowsPerPixel < 1 ? (1 / rowsPerPixel).toFixed(1) + " px per sequence"
                        : rowsPerPixel.toFixed(1) + " sequences per pixel"}. ` +
    "Scroll to zoom, drag to pan.";
}

/* The tree is drawn only while the rows are actually in tree order, since a
   dendrogram beside rows sorted by anything else would draw edges between
   sequences that are not neighbours and quietly misstate the topology. */
function rowsAreTreeOrdered() {
  if (!D || !D.tree || !hasTree || !order) return false;
  for (let r = 0; r < nRows(); r++) if (D.treeRank[order[r]] !== r) return false;
  return true;
}

let treeKids = null;   // children per node, rebuilt when the payload changes
function buildTreeKids() {
  treeKids = null;
  if (!D || !D.tree) return;
  treeKids = Array.from({ length: D.tree.x.length }, () => []);
  for (let i = 0; i < D.tree.x.length; i++)
    if (D.tree.parent[i] >= 0) treeKids[D.tree.parent[i]].push(i);
}

function treeVisible() {
  return !!(D && D.tree) && document.getElementById("showtree").value === "on";
}

/* Shown or hidden before the map is measured, because the box takes its width
   out of the same flex row and the map would otherwise be a frame stale. */
function syncTreeBox() {
  const wanted = treeVisible();
  document.getElementById("treebox").hidden = !wanted;
  document.getElementById("treectl").hidden = !(D && D.tree);
  document.getElementById("treenote").hidden = !wanted;
}

/* Rectangular phylogram sharing the map's vertical viewport, so a leaf sits at
   exactly the height of its row. Leaf labels appear only once a row is tall
   enough to hold one; below that the panel is the tree alone. */
function renderTree() {
  if (!treeVisible() || !view || !order) return;
  const box = document.getElementById("treebox");
  const cv = document.getElementById("tree");
  const dpr = window.devicePixelRatio || 1;
  const cssW = Math.max(60, box.clientWidth - 2);
  cv.style.width = cssW + "px"; cv.style.height = VIEW_HEIGHT + "px";
  cv.width = Math.round(cssW * dpr); cv.height = Math.round(VIEW_HEIGHT * dpr);
  const ctx = cv.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, cssW, VIEW_HEIGHT);

  const note = document.getElementById("treenote");
  const ink = cssRGB("--text-secondary"), muted = cssRGB("--text-muted");
  const strong = cssRGB("--text-primary");

  if (!rowsAreTreeOrdered()) {
    ctx.fillStyle = "rgb(" + muted.join(",") + ")";
    ctx.font = "12px ui-sans-serif, system-ui, sans-serif";
    ctx.textBaseline = "middle";
    ctx.fillText("Set Row order to “Tree order”", 10, VIEW_HEIGHT / 2 - 9);
    ctx.fillText("to draw the tree.", 10, VIEW_HEIGHT / 2 + 9);
    note.textContent = "The tree is drawn only in tree order; any other row order would "
      + "put non-neighbouring sequences side by side.";
    return;
  }

  const T = D.tree, n = T.x.length;
  const pxPerRow = VIEW_HEIGHT / view.h;
  const labelled = pxPerRow >= 9;
  const labelW = labelled ? Math.min(180, cssW * 0.62) : 0;
  const plotW = Math.max(24, cssW - labelW - 8);
  const xmax = T.maxDepth || 1;
  const sx = v => 4 + (v / xmax) * (plotW - 8);
  const sy = row => (row - view.y) * pxPerRow + pxPerRow / 2;

  ctx.strokeStyle = "rgb(" + ink.join(",") + ")";
  ctx.lineWidth = pxPerRow >= 4 ? 1.1 : 0.8;
  ctx.beginPath();
  for (let i = 0; i < n; i++) {
    const p = T.parent[i];
    if (p < 0) continue;
    const y = sy(T.y[i]);
    if (y < -2 || y > VIEW_HEIGHT + 2) continue;
    ctx.moveTo(sx(T.x[p]), y);
    ctx.lineTo(sx(T.x[i]), y);
  }
  // Verticals are drawn from each node's own child span, so a clade taller than
  // the viewport still connects across it instead of vanishing.
  for (let i = 0; i < n; i++) {
    const kids = treeKids[i];
    if (!kids || kids.length < 2) continue;
    let lo = Infinity, hi = -Infinity;
    for (const c of kids) { if (T.y[c] < lo) lo = T.y[c]; if (T.y[c] > hi) hi = T.y[c]; }
    const yl = sy(lo), yh = sy(hi);
    if (yh < -2 || yl > VIEW_HEIGHT + 2) continue;
    const x = sx(T.x[i]);
    ctx.moveTo(x, Math.max(-2, yl));
    ctx.lineTo(x, Math.min(VIEW_HEIGHT + 2, yh));
  }
  ctx.stroke();

  if (labelled) {
    const first = Math.max(0, Math.floor(view.y) - 1);
    const last = Math.min(nRows(), Math.ceil(view.y + view.h) + 1);
    ctx.textBaseline = "middle";
    for (let i = 0; i < n; i++) {
      const s = T.leaf[i];
      if (s < 0) continue;
      const row = T.y[i];
      if (row < first || row > last) continue;
      const y = sy(row);
      if (y < -2 || y > VIEW_HEIGHT + 2) continue;
      const fam = hasGroups ? D.groupNames[D.groupOf[s]] : "";
      let x = plotW + 4;
      if (fam) {
        ctx.font = "600 10px ui-sans-serif, system-ui, sans-serif";
        ctx.fillStyle = "rgb(" + strong.join(",") + ")";
        ctx.fillText(fam, x, y);
        x += ctx.measureText(fam).width + 5;
      }
      ctx.font = "10px ui-monospace, SFMono-Regular, Menlo, monospace";
      ctx.fillStyle = "rgb(" + muted.join(",") + ")";
      ctx.fillText(D.ids[s], x, y);
    }
  }

  note.textContent =
    "Tree panel: rectangular phylogram, x axis in substitutions per site, deepest tip "
    + T.maxDepth.toFixed(2) + ". Rooted where the build was told to root it. "
    + (labelled ? "Family and accession are shown beside each leaf."
                : "Zoom in past nine pixels per row for family and accession labels.");
}

/* The side panel shares the map's vertical viewport exactly, so a row sits at
   the same height in both and they zoom and pan together. Bars are drawn at
   display resolution rather than blitted, so the gene names stay crisp. */
function renderSide() {
  const cv = document.getElementById("side");
  if (!view || !order) return;
  const dpr = window.devicePixelRatio || 1;
  const cssW = cv.parentElement.clientWidth || 360;
  cv.style.width = cssW + "px"; cv.style.height = VIEW_HEIGHT + "px";
  cv.width = Math.round(cssW * dpr); cv.height = Math.round(VIEW_HEIGHT * dpr);
  const ctx = cv.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, cssW, VIEW_HEIGHT);

  const { values, label, max } = sidePanelValues();
  const barW = Math.min(110, Math.max(60, cssW * 0.3));
  const pxPerRow = VIEW_HEIGHT / view.h;
  const ink = cssRGB("--text-secondary"), muted = cssRGB("--text-muted");
  const good = cssRGB("--hist"), bad = cssRGB("--outlier");
  const cut = cutValue();

  // axis: 0 on the left, 1 at the bar's full width
  ctx.strokeStyle = `rgb(${muted.join(",")})`; ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(0.5, 0); ctx.lineTo(0.5, VIEW_HEIGHT); ctx.stroke();

  const first = Math.max(0, Math.floor(view.y));
  const last = Math.min(nRows(), Math.ceil(view.y + view.h));
  const labelled = pxPerRow >= 9;
  ctx.font = "11px ui-sans-serif, system-ui, sans-serif";
  ctx.textBaseline = "middle";

  // One line running down the panel: x is the score, y is the row, so the trace
  // sits beside the sequences it describes and moves with them.
  ctx.strokeStyle = `rgb(${good.join(",")})`;
  ctx.lineWidth = 1.5;
  ctx.lineJoin = "round";
  ctx.beginPath();
  let started = false;
  if (pxPerRow >= 1) {
    for (let row = first; row < last; row++) {
      const v = values[order[row]];
      const y = (row - view.y) * pxPerRow + pxPerRow / 2;
      const x = 1 + barW * (v / max);
      if (started) ctx.lineTo(x, y); else { ctx.moveTo(x, y); started = true; }
    }
  } else {
    // More rows than pixels: one point per screen row, the mean of what it covers.
    for (let py = 0; py < VIEW_HEIGHT; py++) {
      const a = Math.floor(view.y + py * view.h / VIEW_HEIGHT);
      const b = Math.max(a + 1, Math.floor(view.y + (py + 1) * view.h / VIEW_HEIGHT));
      let sum = 0, n = 0;
      for (let row = a; row < b && row < nRows(); row++) { sum += values[order[row]]; n++; }
      if (!n) continue;
      const x = 1 + barW * (sum / n / max);
      if (started) ctx.lineTo(x, py + 0.5); else { ctx.moveTo(x, py + 0.5); started = true; }
    }
  }
  ctx.stroke();

  // The cutoff, drawn across the panel at the row it falls on, so the boundary
  // between what is kept and what is removed lines up with the same sequences in
  // the alignment map beside it. It is the only dashed line here: the knee and the
  // vertical score line were both dropped, because three references to the same
  // axis made the panel harder to read than the one that decides the membership.
  const cutRow = cutoffRow();
  if (cutRow !== null) {
    const cy = (cutRow - view.y) * pxPerRow;
    if (cy >= -1 && cy <= VIEW_HEIGHT + 1) {
      ctx.strokeStyle = `rgb(${ink.join(",")})`;
      ctx.lineWidth = 1;
      ctx.setLineDash([6, 4]);
      ctx.beginPath();
      ctx.moveTo(0, cy + 0.5); ctx.lineTo(cssW, cy + 0.5);
      ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillStyle = `rgb(${ink.join(",")})`;
      ctx.textAlign = "right";
      ctx.fillText(fmtScore(cut), cssW - 4, cy - 7);
      ctx.textAlign = "left";
    }
  }

  // Rows below the threshold get a marker, so they stay findable once the line
  // is dense, and the labels beside them are written in the same colour.
  if (pxPerRow >= 3) {
    for (let row = first; row < last; row++) {
      const s = order[row], v = values[s];
      if (v >= (max === 1 ? cut / scoreMax : cut)) continue;
      const y = (row - view.y) * pxPerRow + pxPerRow / 2;
      ctx.fillStyle = `rgb(${bad.join(",")})`;
      ctx.beginPath();
      ctx.arc(1 + barW * (v / max), y, Math.min(3, Math.max(1.5, pxPerRow / 3)), 0, 6.284);
      ctx.fill();
    }
  }

  if (labelled) {
    for (let row = first; row < last; row++) {
      const s = order[row], v = values[s];
      const y = (row - view.y) * pxPerRow + pxPerRow / 2;
      ctx.fillStyle = `rgb(${(v < (max === 1 ? cut / scoreMax : cut) ? bad : ink).join(",")})`;
      const gene = D.genes && D.genes[s] ? " " + D.genes[s] : "";
      const name = (D.descriptions[s].split(/\s+/).slice(1).join(" ") || "")
        .replace(/ (OS|OX|GN|PE|SV)=.*$/, "");
      ctx.fillText(`${fmtScore(v)}  ${D.ids[s]}${gene}  ${name}`.slice(0, 64),
                   barW + 8, y);
    }
  }

  ctx.fillStyle = `rgb(${muted.join(",")})`;
  ctx.textBaseline = "top";
  ctx.fillText(`${label}  0`, 2, 2);
  ctx.fillText(fmtScore(max), barW - 16, 2);
  if (!labelled) ctx.fillText("zoom in for gene names", barW + 8, 2);
}

/* Screen point to cell, honouring the current viewport. */
function cellAt(clientX, clientY) {
  const cv = document.getElementById("heat");
  const r = cv.getBoundingClientRect();
  const fx = (clientX - r.left) / r.width, fy = (clientY - r.top) / r.height;
  if (fx < 0 || fx > 1 || fy < 0 || fy > 1) return null;
  return { col: Math.floor(view.x + fx * view.w), row: Math.floor(view.y + fy * view.h) };
}

function setupZoom() {
  const cv = document.getElementById("heat");
  cv.addEventListener("wheel", event => {
    if (!view || !base) return;
    event.preventDefault();               // zoom the map, do not scroll the page
    const r = cv.getBoundingClientRect();
    const fx = (event.clientX - r.left) / r.width, fy = (event.clientY - r.top) / r.height;
    const anchorX = view.x + fx * view.w, anchorY = view.y + fy * view.h;
    const factor = event.deltaY < 0 ? 1 / 1.2 : 1.2;
    /* The matrix is far taller than it is wide, so zooming both axes by the same
       factor would keep cells as thin horizontal lines however far you went in.
       Instead the row count is pulled toward the value that makes a cell square,
       closing most of the excess on every step, so square cells arrive while
       there are still tens of columns on screen rather than three. Zooming out
       reverses the pull and returns to the full matrix just as quickly. */
    view.w *= factor;
    const squareRows = VIEW_HEIGHT * view.w / (cv.clientWidth || 900);
    const PULL = 0.55;
    const scaled = view.h * factor;
    view.h = factor < 1 ? squareRows + (scaled - squareRows) * PULL
                        : squareRows + (scaled - squareRows) / PULL;
    clampView();
    // Keep the cell under the cursor under the cursor.
    view.x = anchorX - fx * view.w; view.y = anchorY - fy * view.h;
    clampView();
    renderView();
  }, { passive: false });

  /* Pan mode drags the view. Select mode toggles the row clicked. Range mode
     selects every row between press and release; a drag that starts on a row
     already selected clears the range instead, so a mistake is undone the
     same way it was made. */
  const mode = () => document.getElementById("mode").value;
  const rowAt = event => {
    const r = cv.getBoundingClientRect();
    const at = cellAt(r.left + 1, Math.min(Math.max(event.clientY, r.top), r.bottom - 1));
    return at ? Math.min(nRows() - 1, Math.max(0, at.row)) : null;
  };
  let dragging = null, ranging = null;
  cv.addEventListener("mousedown", event => {
    if (mode() === "pan") {
      dragging = { x: event.clientX, y: event.clientY, vx: view.x, vy: view.y };
      cv.style.cursor = "grabbing";
      return;
    }
    const row = rowAt(event);
    if (row === null) return;
    const id = D.ids[order[row]];
    if (mode() === "single") {
      if (selected.has(id)) selected.delete(id); else selected.add(id);
      selectionChanged();
    } else {
      ranging = { from: row, to: row, add: !selected.has(id), before: new Set(selected) };
      applyRange();
    }
    event.preventDefault();
  });
  function applyRange() {
    selected = new Set(ranging.before);
    const a = Math.min(ranging.from, ranging.to), b = Math.max(ranging.from, ranging.to);
    for (let r = a; r <= b; r++) {
      const id = D.ids[order[r]];
      if (ranging.add) selected.add(id); else selected.delete(id);
    }
    document.getElementById("selstatus").textContent =
      `${ranging.add ? "selecting" : "clearing"} rows ${fmtInt(a + 1)}–${fmtInt(b + 1)}`;
    renderView();
  }
  addEventListener("mousemove", event => {
    if (ranging) {
      const row = rowAt(event);
      if (row !== null && row !== ranging.to) { ranging.to = row; applyRange(); }
      return;
    }
    if (!dragging) return;
    const r = cv.getBoundingClientRect();
    view.x = dragging.vx - (event.clientX - dragging.x) / r.width * view.w;
    view.y = dragging.vy - (event.clientY - dragging.y) / r.height * view.h;
    clampView(); renderView();
  });
  addEventListener("mouseup", () => {
    if (ranging) { ranging = null; selectionChanged(); }
    dragging = null;
    cv.style.cursor = mode() === "pan" ? "grab" : "crosshair";
  });
  document.getElementById("mode").addEventListener("change", () => {
    cv.style.cursor = mode() === "pan" ? "grab" : "crosshair";
  });
  cv.style.cursor = "grab";

  document.getElementById("resetzoom").addEventListener("click", () => {
    view = { x: 0, y: 0, w: base.width, h: base.height };
    renderView();
  });
}

/* The two tracks sit exactly over the alignment canvas and show the same
   columns it does, so they follow its zoom and pan column for column. */
function drawTracks() {
  const occ = document.getElementById("occupancy");
  const cons = document.getElementById("consensus");
  const heat = document.getElementById("heat");
  const dpr = window.devicePixelRatio || 1;
  for (const cv of [occ, cons]) {
    cv.style.marginLeft = "0px";
    cv.style.width = heat.clientWidth + "px";
    cv.style.marginLeft = (heat.getBoundingClientRect().left - cv.getBoundingClientRect().left) + "px";
  }
  const boxW = heat.clientWidth || 900;
  const x0 = view ? view.x : 0, span = view ? view.w : Math.max(1, visible.length);
  const bw = boxW / span;
  const first = Math.max(0, Math.floor(x0)), last = Math.min(visible.length, Math.ceil(x0 + span));

  occ.width = boxW * dpr; occ.height = 44 * dpr;
  const octx = occ.getContext("2d"); octx.scale(dpr, dpr);
  octx.clearRect(0, 0, boxW, 44);
  cons.width = boxW * dpr; cons.height = 34 * dpr;
  const cctx = cons.getContext("2d"); cctx.scale(dpr, dpr);
  cctx.clearRect(0, 0, boxW, 34);
  for (let k = first; k < last; k++) {
    const col = D.columns[visible[k]];
    const c = rgb[col.cons] || rgb["-"];
    const x = (k - x0) * bw;
    octx.fillStyle = cctx.fillStyle = `rgb(${c[0]},${c[1]},${c[2]})`;
    const bh = Math.max(1, (1 - col.gap) * 40);
    octx.fillRect(x, 44 - bh, Math.max(0.5, bw - 0.4), bh);
    cctx.fillRect(x, 0, Math.max(0.5, bw + 0.2), 18);
  }
  drawSegmentNames(cctx, boxW, x0, bw, first, last);
}

/* The consensus read as a receptor rather than as a run of letters: the helix
   blocks the score already uses become TM1..TMn, and each stretch between them
   is named from its own consensus side, ICL when the family sits inside and ECL
   when it sits outside, numbered in the order they appear. Nothing is assumed
   about the canonical alternation, so a family whose bundle is not the textbook
   seven-helix arrangement is labelled by what its own consensus actually says.
   Anything before the first helix or after the last is a terminus. */
function consensusSegments() {
  if (!blocks || !blocks.length) return [];
  const sideOf = (a, b) => {
    let inside = 0, outside = 0;
    for (let k = a; k < b; k++) {
      const c = D.columns[visible[k]].cons;
      if (c === "I") inside++; else if (c === "O") outside++;
    }
    return inside === 0 && outside === 0 ? null : (inside >= outside ? "I" : "O");
  };
  const segs = [];
  let icl = 0, ecl = 0;
  if (blocks[0][0] > 0) segs.push([0, blocks[0][0], sideOf(0, blocks[0][0]) === "I" ? "N-term in" : "N-term"]);
  blocks.forEach(([a, b], i) => {
    segs.push([a, b, `TM${i + 1}`]);
    const next = i + 1 < blocks.length ? blocks[i + 1][0] : visible.length;
    if (next > b) {
      const side = sideOf(b, next);
      const last = i + 1 === blocks.length;
      if (side === null) segs.push([b, next, ""]);
      else if (last) segs.push([b, next, side === "I" ? "C-term in" : "C-term out"]);
      else segs.push([b, next, side === "I" ? `ICL${++icl}` : `ECL${++ecl}`]);
    }
  });
  return segs;
}

function drawSegmentNames(ctx, boxW, x0, bw, first, last) {
  const segs = consensusSegments();
  if (!segs.length) return;
  ctx.save();
  ctx.font = "600 10px ui-monospace, SFMono-Regular, Menlo, monospace";
  ctx.textBaseline = "middle";
  ctx.textAlign = "center";
  const ink = cssRGB("--text-secondary"), rule = cssRGB("--border");
  for (const [a, b, name] of segs) {
    if (!name || b <= first || a >= last) continue;
    const xa = (a - x0) * bw, xb = (b - x0) * bw;
    // A tick at each boundary, so a segment stays readable once it is too
    // narrow for its own name.
    ctx.strokeStyle = `rgb(${rule[0]},${rule[1]},${rule[2]})`;
    ctx.beginPath(); ctx.moveTo(xa, 18); ctx.lineTo(xa, 22); ctx.stroke();
    const w = ctx.measureText(name).width;
    if (xb - xa < w + 6) continue;
    ctx.fillStyle = `rgb(${ink[0]},${ink[1]},${ink[2]})`;
    ctx.fillText(name, Math.max(w / 2 + 1, Math.min(boxW - w / 2 - 1, (xa + xb) / 2)), 27);
  }
  ctx.restore();
}

function drawHist() {
  const cv = document.getElementById("hist");
  const dpr = window.devicePixelRatio || 1;
  const w = cv.clientWidth || 900, h = 120;
  cv.width = w * dpr; cv.height = h * dpr;
  const ctx = cv.getContext("2d"); ctx.scale(dpr, dpr);
  ctx.clearRect(0, 0, w, h);
  const bins = new Array(50).fill(0);
  for (let i = 0; i < scores.length; i++)
    bins[Math.min(49, Math.floor(scores[i] / scoreMax * 50))]++;
  const max = Math.max(...bins);
  const bar = cssRGB("--hist"), flag = cssRGB("--outlier");
  const cut = cutValue();
  const bw = w / bins.length;
  for (let i = 0; i < bins.length; i++) {
    const bh = max ? (bins[i] / max) * (h - 22) : 0;
    const c = ((i + 1) / 50 * scoreMax <= cut) ? flag : bar;
    ctx.fillStyle = `rgb(${c[0]},${c[1]},${c[2]})`;
    const r = Math.min(4, bw / 2, bh);
    ctx.beginPath();
    ctx.roundRect(i * bw + 1, h - 18 - bh, Math.max(1, bw - 2), bh, [r, r, 0, 0]);
    ctx.fill();
  }
  ctx.fillStyle = `rgb(${cssRGB("--text-muted").join(",")})`;
  ctx.font = "11px ui-sans-serif, system-ui, sans-serif";
  ctx.fillText("0", 2, h - 4);
  ctx.fillText("helix Jaccard similarity", w / 2 - 62, h - 4);
  ctx.textAlign = "right"; ctx.fillText(fmtScore(scoreMax), w - 2, h - 4);
  ctx.textAlign = "left";
}

function drawTable() {
  const cut = cutValue();
  const idx = Array.from({ length: D.nSeq }, (_, i) => i)
    .sort((a, b) => scores[a] - scores[b] || a - b);
  const flagged = idx.filter(i => scores[i] < cut);
  document.getElementById("flagcount").textContent =
    `— ${fmtInt(flagged.length)} of ${fmtInt(D.nSeq)} below ${fmtScore(cut)}`;
  const rows = idx.slice(0, 100);
  const body = document.getElementById("outliers").querySelector("tbody");
  body.innerHTML = rows.map((i, rank) => {
    const comp = comps[i];
    const total = Object.values(comp).reduce((a, b) => a + b, 0) || 1;
    const chips = D.labels.filter(l => comp[l]).map(l =>
      `<span style="display:inline-flex;align-items:center;gap:4px;margin-right:8px">
         <span class="sw" style="background:var(${colorVar(l)})"></span>${l} ${pct(comp[l] / total)}</span>`).join("");
    return `<tr class="${scores[i] < cut ? "flag" : ""}">
      <td class="num">${rank + 1}</td>
      <td class="mono" title="${D.descriptions[i].replace(/"/g, "&quot;")}">${D.ids[i]}</td>
      ${hasGroups ? `<td>${D.groupNames[D.groupOf[i]]}</td>` : ""}
      <td class="num score">${fmtScore(scores[i])}</td>
      <td class="num">${pct(coverage[i])}</td>
      <td class="num">${fmtInt(residues[i])}</td>
      <td>${chips}</td></tr>`;
  }).join("");
  window.__flagged = flagged;
}

function drawTiles() {
  /* The tiles stay family-wide on purpose. They describe the alignment, which is
     the same alignment in all three views, so they must not move when the view
     does -- a median that changed with the button would read as a property of the
     subset rather than of the family the cutoff is being argued about. */
  const mean = scores.reduce((a, b) => a + b, 0) / scores.length;
  const sorted = Array.from(scores).sort((a, b) => a - b);
  const median = sorted[Math.floor(sorted.length / 2)];
  const p01 = sorted[Math.floor(sorted.length * 0.01)];
  document.getElementById("tiles").innerHTML = `
    <div class="tile"><div class="k">Sequences</div><div class="v">${fmtInt(D.nSeq)}</div></div>
    <div class="tile"><div class="k">Visible columns</div><div class="v">${fmtInt(visible.length)}</div>
      <div class="n">of ${fmtInt(D.alignmentWidth)} in the alignment</div></div>
    <div class="tile"><div class="k">Median score</div><div class="v">${fmtScore(median)}</div>
      <div class="n">mean ${mean.toFixed(3)}</div></div>
    <div class="tile"><div class="k">1st percentile</div><div class="v">${p01.toFixed(3)}</div>
      <div class="n">bottom 1% fall below this</div></div>`;
}

function redrawAll() {
  drawTracks(); drawHeat(); drawHist(); drawTable(); drawGroupTable(); drawTiles();
}

function recompute() {
  const busy = document.getElementById("busy");
  busy.textContent = "computing…";
  setTimeout(async () => {
    const maxGap = parseInt(document.getElementById("maxgap").value, 10) / 100;
    document.getElementById("maxgapval").textContent =
      `${(100 * maxGap).toFixed(0)}% gaps — ${fmtInt(visibleCount(maxGap))} columns shown`;
    computeVisible(maxGap);
    computeScores();
    // Only a protein panel borrows, and finding the partner decodes payloads, so
    // the per-family viewer with its nineteen topology alignments never pays for it.
    if (isProtein()) await ensureTopoScores();
    computeOrderScores();
    refreshSortLabels();
    computeOrder(document.getElementById("sort").value);
    redrawAll();
    busy.textContent = "";
  }, 0);
}
function visibleCount(maxGap) {
  let n = 0; for (let c = 0; c < D.nCol; c++) if (colGap[c] <= maxGap + 1e-9) n++; return n;
}

/* What the live score is of, named so the tooltip cannot be misread. For a topology
   alignment it is the helix Jaccard against that alignment's consensus, over its 30%
   occupancy columns whatever the slider shows. */
function scoreNoun() {
  return isProtein() ? "identity to this alignment" : "this alignment, helix Jaccard";
}

/* The precomputed score from the sequence's own family alignment, which is a
   different measurement from the one above and is display-only. */
function famScoreLine(s) {
  if (!D.familyScore) return "own-family score not loaded";
  const v = D.familyScore[s];
  if (v === null || v === undefined) return "own-family score unavailable for this accession";
  const delta = v - scores[s];
  const sign = delta >= 0 ? "+" : "\u2212";
  /* Both sides are now computed over the 30% occupancy columns, so the difference is
     between the two consensuses -- this alignment's and the sequence's own family's --
     and not between two slider settings. */
  return `own family, helix Jaccard <b>${fmtScore(v)}</b> · ${sign}${(100 * Math.abs(delta)).toFixed(1)} points vs this alignment`;
}

function setupHover() {
  const cv = document.getElementById("heat"), tip = document.getElementById("tip");
  cv.addEventListener("mousemove", e => {
    const at = cellAt(e.clientX, e.clientY);
    if (!at) { tip.style.display = "none"; return; }
    const row = at.row, col = at.col;
    if (row < 0 || row >= nRows() || col < 0 || col >= visible.length) { tip.style.display = "none"; return; }
    const s = order[row], column = D.columns[visible[col]];
    const ch = D.matrix[s * D.nCol + visible[col]];
    tip.innerHTML = `<div class="t">${D.ids[s]}${hasGroups ? " &middot; " + D.groupNames[D.groupOf[s]] : ""}${selected.has(D.ids[s]) ? " &middot; excluded" : ""}</div>
      <div class="r">${D.descriptions[s].split(" ").slice(1).join(" ").slice(0, 90) || "&mdash;"}</div>
      <div class="r">${scoreNoun()} <b>${fmtScore(scores[s])}</b> · coverage ${pct(coverage[s])} · ${fmtInt(residues[s])} residues</div>
      <div class="r">${famScoreLine(s)}</div>
      <div class="r">alignment column ${fmtInt(column.i + 1)} · this sequence ${ch === "-" ? "gap" : labelName(ch)}</div>
      <div class="r">column consensus ${labelName(column.cons)} (${pct(column.consFrac)}) · occupancy ${pct(1 - column.gap)}</div>`;
    tip.style.display = "block";
    const x = Math.min(e.clientX + 14, window.innerWidth - tip.offsetWidth - 10);
    const y = Math.min(e.clientY + 14, window.innerHeight - tip.offsetHeight - 10);
    tip.style.left = x + "px"; tip.style.top = y + "px";
  });
  cv.addEventListener("mouseleave", () => { tip.style.display = "none"; });
}

function setupLegend() {
  document.getElementById("legend").innerHTML =
    D.labels.map(l => `<span><span class="sw" style="background:var(${colorVar(l)})"></span>
      ${l} &mdash; ${labelName(l)}</span>`).join("") +
    `<span><span class="sw" style="background:var(--gap)"></span>gap</span>`;
}

/* The exclusion list is saved per family by serve_topology_viewer.py on every
   change, so closing the page loses nothing. Opened as a plain file there is no
   server, and the Download button is the way to keep the list. */
let saveTimer = null, pendingSave = null;
function selectionChanged() {
  document.getElementById("selcount").textContent = fmtInt(selected.size);
  renderView();
  const status = document.getElementById("selstatus");
  if (!serverOK) {
    status.className = "status err";
    status.textContent = "not saved: open through serve_topology_viewer.py, or use Download";
    return;
  }
  status.className = "status"; status.textContent = "saving…";
  clearTimeout(saveTimer);
  pendingSave = { name: family, ids: [...selected].sort() };
  saveTimer = setTimeout(flushSave, 400);
}

async function flushSave() {
  clearTimeout(saveTimer);
  if (!pendingSave) return;
  const { name, ids } = pendingSave;
  pendingSave = null;
  const status = document.getElementById("selstatus");
  try {
    const res = await fetch(`selections/${encodeURIComponent(name)}`, {
      method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(ids) });
    const info = await res.json();
    if (!res.ok) throw new Error(info.error || res.status);
    if (name === family) status.textContent = `saved ${fmtInt(ids.length)} to ${info.path}`;
  } catch (err) {
    status.className = "status err"; status.textContent = `save failed for ${name}: ${err.message}`;
  }
}

async function loadSelection(name) {
  await flushSave();   // a save still pending for the previous family goes out first
  selected = new Set();
  const status = document.getElementById("selstatus");
  serverOK = false;
  if (location.protocol.startsWith("http")) {
    try {
      const res = await fetch(`selections/${encodeURIComponent(name)}`);
      if (res.ok) {
        const info = await res.json();
        const present = new Set(D.ids);
        selected = new Set(info.ids.filter(id => present.has(id)));
        serverOK = true;
        status.className = "status";
        status.textContent = info.ids.length ? `loaded ${fmtInt(selected.size)} from ${info.path}`
                                             : `will save to ${info.path}`;
      }
    } catch (err) { /* no server; fall through */ }
  }
  if (!serverOK) {
    status.className = "status err";
    status.textContent = "no server: selections are not saved automatically";
  }
  document.getElementById("selcount").textContent = fmtInt(selected.size);
}

/* Swap in another family. Everything drawn is derived from D, so decoding the
   new payload and recomputing is enough; the zoom resets to the full matrix. */
async function loadFamily(name) {
  const busy = document.getElementById("busy");
  busy.textContent = "loading " + name + "…";
  D = await decode(PAYLOADS[name]);
  family = name;
  await loadSelection(name);
  colGap = D.columns.map(c => c.gap);
  hasGroups = Array.isArray(D.groupNames) && D.groupNames.length > 0;
  hasLineages = Array.isArray(D.lineageNames) && D.lineageNames.length > 0;
  hasTree = Array.isArray(D.treeRank) && D.treeRank.length === D.nSeq;
  buildTreeKids();
  const showTree = document.getElementById("showtree");
  document.getElementById("treectl").hidden = !D.tree;
  if (!D.tree) showTree.value = "off";
  else if (!showTree.dataset.touched) showTree.value = "on";
  // The two alphabets colour differently, so the cached RGB has to be rebuilt
  // whenever the alignment on screen changes, not only when the theme does.
  refreshColors();
  const sortSelect = document.getElementById("sort");
  document.getElementById("grouppanel").hidden = !hasGroups;
  document.getElementById("grouphead").hidden = !hasGroups;
  document.getElementById("gutter").style.display = (hasGroups || hasLineages) ? "" : "none";
  for (const value of ["group", "group_median"]) {
    const option = sortSelect.querySelector(`option[value="${value}"]`);
    option.hidden = option.disabled = !hasGroups;
  }
  for (const value of ["lineage", "lineage_median"]) {
    const option = sortSelect.querySelector(`option[value="${value}"]`);
    option.hidden = option.disabled = !hasLineages;
  }
  const treeOption = sortSelect.querySelector('option[value="tree"]');
  treeOption.hidden = treeOption.disabled = !hasTree;
  // The score means a different thing in each alphabet, so the menu says which.
  const scoreLabel = isProtein() ? "Identity to consensus" : "Helix Jaccard similarity";
  sortSelect.querySelector('option[value="score_desc"]').textContent = scoreLabel + ", best first";
  sortSelect.querySelector('option[value="score"]').textContent = scoreLabel + ", worst first";
  document.getElementById("flagcutlabel").textContent =
    isProtein() ? "Flag identity below" : "Flag score below";
  if (!sortSelect.dataset.touched)
    sortSelect.value = (DEFAULT_GROUP_SORT === "tree" && hasTree) ? "tree"
                     : DEFAULT_GROUP_SORT.startsWith("lineage") && hasLineages ? DEFAULT_GROUP_SORT
                     : hasGroups ? DEFAULT_GROUP_SORT : hasLineages ? "lineage_median" : "score_desc";
  else if (!hasGroups && sortSelect.value.startsWith("group")) sortSelect.value = "score_desc";
  else if (!hasLineages && sortSelect.value.startsWith("lineage")) sortSelect.value = "score_desc";
  else if (!hasTree && sortSelect.value === "tree") sortSelect.value = "score_desc";
  if (Object.keys(PAYLOADS).length > 1) document.getElementById("famtitle").textContent = name;
  document.getElementById("subtitle").textContent =
    `${fmtInt(D.nSeq)} sequences · alignment ${fmtInt(D.alignmentWidth)} columns wide · ` +
    `${fmtInt(D.nCol)} columns with at most ${(100 * MAX_GAP_PAYLOAD).toFixed(0)}% gaps carried into this file · ` +
    `${isProtein() ? "residues, Taylor colours" : "predicted topology"}` +
    `${hasTree ? " · tree order available" : ""}`;
  setupLegend();
  recompute();
}

(async function init() {
  refreshColors();
  const names = Object.keys(PAYLOADS);
  const pick = document.getElementById("family");
  pick.innerHTML = names.map(n => `<option value="${n}">${n}</option>`).join("");
  if (names.length < 2) document.getElementById("famctl").hidden = true;
  pick.addEventListener("change", () => loadFamily(pick.value));
  document.getElementById("selclear").addEventListener("click", () => {
    if (!selected.size || !confirm(`Clear all ${selected.size} excluded sequences of ${family}?`)) return;
    selected = new Set(); selectionChanged();
  });
  document.getElementById("seldownload").addEventListener("click", () => {
    const url = URL.createObjectURL(new Blob([[...selected].sort().join("\n") + "\n"], { type: "text/plain" }));
    const a = document.createElement("a");
    a.href = url; a.download = `${family}_excluded.txt`; a.click();
    URL.revokeObjectURL(url);
  });
  const slider = document.getElementById("maxgap");
  slider.max = Math.round(MAX_GAP_PAYLOAD * 100);
  slider.value = Math.round(DEFAULT_MAX_GAP * 100);
  setupHover(); setupZoom();
  slider.addEventListener("input", () => {
    document.getElementById("maxgapval").textContent =
      `${slider.value}% gaps — ${fmtInt(visibleCount(slider.value / 100))} columns shown`;
  });
  slider.addEventListener("change", recompute);
  document.getElementById("sort").addEventListener("change", () => {
    document.getElementById("sort").dataset.touched = "1";
    computeOrder(document.getElementById("sort").value); drawHeat(true);
  });
  /* The cutoff defines what "removed" means, so moving it has to rebuild the row
     set as well as the panels; otherwise a filtered view would keep showing the
     split from the previous cutoff. */
  document.getElementById("flagcut").addEventListener("change", () => {
    computeOrder(document.getElementById("sort").value);
    drawHeat(rowView === "all"); drawHist(); drawTable(); drawGroupTable(); renderSide();
  });
  document.getElementById("rowview").addEventListener("click", event => {
    const button = event.target.closest("button[data-view]");
    if (!button || button.dataset.view === rowView) return;
    rowView = button.dataset.view;
    for (const b of document.querySelectorAll("#rowview button"))
      b.classList.toggle("on", b.dataset.view === rowView);
    // The three views have different row counts, so the vertical zoom cannot be
    // carried across; the horizontal one is about columns and is left alone.
    computeOrder(document.getElementById("sort").value);
    /* Only the row-bearing panels are redrawn. The tiles, the histogram, the
       occupancy curve and the consensus track are properties of the alignment
       rather than of the rows on screen, so they are left exactly as they are
       and read the same in all three views. */
    drawHeat(false); drawGutter(); renderSide();
  });
  document.getElementById("helixpick").addEventListener("change", renderSide);
  document.getElementById("showtree").addEventListener("change", () => {
    document.getElementById("showtree").dataset.touched = "1";
    renderView();
  });
  document.getElementById("theme").addEventListener("click", () => {
    const dark = matchMedia("(prefers-color-scheme: dark)").matches;
    const cur = document.documentElement.getAttribute("data-theme") || (dark ? "dark" : "light");
    document.documentElement.setAttribute("data-theme", cur === "dark" ? "light" : "dark");
    refreshColors(); redrawAll();
  });
  document.getElementById("export").addEventListener("click", () => {
    const ids = (window.__flagged || []).map(i => D.descriptions[i].split(" ")[0]).join("\n");
    const url = URL.createObjectURL(new Blob([ids + "\n"], { type: "text/plain" }));
    const a = document.createElement("a");
    a.href = url; a.download = `${family}_flagged_sequences.txt`; a.click();
    URL.revokeObjectURL(url);
  });
  addEventListener("resize", () => { drawTracks(); drawHist(); renderView(); });
  // The side panel scrolls the map too, so the two never drift apart.
  document.getElementById("side").addEventListener("wheel", event => {
    document.getElementById("heat").dispatchEvent(new WheelEvent("wheel", {
      deltaY: event.deltaY, clientX: document.getElementById("heat")
        .getBoundingClientRect().left + 10, clientY: event.clientY, bubbles: false }));
    event.preventDefault();
  }, { passive: false });
  await loadFamily(names[0]);
})();
</script>
</body>
</html>
"""


if __name__ == "__main__":
    main()
