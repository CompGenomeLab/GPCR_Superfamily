"""
Build the data of the family-assignment page of gpcrevolution.org (transparency/family-assignment).

For every family in ``gene_boundaries`` (read out of ../identify_gene_clade.py, so the two never
drift apart) this reads the ordered search outputs of the archive -- the score-ordered tree, the
tree-ordered MSA and the score-and-index line plot of stage 01 -- and writes one compact JSON bundle
per family (<out>/<family>.json) and the family list (<out>/index.json). The page draws, along one
x-axis (leaf order in the score-ordered tree), the binned MSA, the search score and index, the tree
and the extracted clade, re-tinted for the members left after each step:

    initial            the clade as extracted (stage 02)
    after_7tm          single GPCR-like 7TM core (03_topology_filter/core_membership.tsv)
    after_topology     topology score >= 0.8 (03_topology_filter/topology_membership.tsv)
    after_assignment   overlaps resolved (04_overlap_assignment/final_membership.tsv)
    after_reassignment Dikarya Class A moved to GPR1
                       (10_motif_analysis/phylogeny_reassignment/final_membership_reassigned.tsv)

Adapted from the standalone figure builder of the working pipeline; only the inputs (archive),
the stages and the output (JSON for the website, no HTML) changed.

Usage:
    python build_family_assignment_data.py --out <GPCRevolution>/public/family-assignment [family ...]
"""

import ast
import base64
import csv
import glob
import json
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from gpcr_paths import open_text, stage  # noqa: E402

TARGET_BINS = 1000     # leaves are aggregated into at most this many x-axis bins

# How each family is named everywhere a reader sees it -- the page title, the
# statistics bar, the clade legend, the family switcher, the index table and the
# cross-annotation gutter. The keys stay as they are, because they are the
# `gene_boundaries` names the source filenames are built from, so renaming a
# family here never touches a path.
DISPLAY_NAMES = {
    "classA": "Class A",
    "classB1": "Class B1",
    "classB2": "Class B2",
    "classC": "Class C",
    "classF": "Class F",
    "classT": "Taste-2",
    "Galpha": "Gα",
    "Gbeta": "Gβ",
    "Ggamma": "Gγ",
    "FSLB_DICDI": "FSL",
    "GPR1_YEAST": "GPR1",
    "GP143": "GP143",
    "GP157": "GP157",
    "TM116": "TM116",
    "Vomeronasal1": "Vomeronasal 1",
    "Vomeronasal2": "Vomeronasal 2",
    "Mth": "Mth",
    "GIT3_SCHPO": "Git3",
}

# Families still present in gene_boundaries but not part of the final analysis.
# They get no figure of their own and are never used to cross-annotate another
# family's tree.
EXCLUDED_FAMILIES = {"LRLA_DICDI"}

# The G protein subunits and the non-receptor controls. They still get their own
# pages, but they take no part in the overlap analysis: they never mark another
# family's tree, and another family's clade is never judged against theirs.
# GIT3_SCHPO is deliberately NOT here: its proposed clade is shown as a track on
# every other family's tree, to see where Git3 lands across the superfamily. Put
# it back in this set to take it off the other pages again.
OVERLAP_EXCLUDED = {"Galpha", "Gbeta", "Ggamma", "TBP", "RBP1"}

# Cross-annotation: how another family's extracted clade is reported when its
# sequences land inside this family's tree. Several searches share one tree, so
# these blocks are usually clean; across different searches they can be patchy.
# Some families were sampled into their own tree as representatives, so the clade
# read off that tree understates the family's real membership. Where a complete
# accession list is available it is unioned with the clade for cross-annotation
# purposes, which is what lets the family be marked in full in other families'
# trees. The family's own highlighted clade still comes from its tree boundaries.
COMPLETE_MEMBER_FILES = {
    "classA": "02_clade_extraction/clade_membership.tsv",
    "Olfactory": "02_clade_extraction/clade_membership.tsv",
}

# Sequence sets that are not extracted families but are worth tracking across the
# trees anyway, each defined as a leaf range of one family's tree so that it
# regenerates rather than depending on a hand-made list. They annotate every tree,
# including the one they were taken from, where they simply mark themselves.
EXTRA_CROSS_SETS = {}

# Families that are not in `gene_boundaries` at all, given with the search and
# the boundary they are read from. Entries here are merged with `setdefault`, so
# a family that later gains a real `gene_boundaries` entry takes it from there
# instead. TMEM116 was defined here while it was figure-only and now lives in
# `gene_boundaries`, so nothing remains.
EXTRA_FAMILIES = {}
PENDING_FAMILIES = {}
VARIANT_FAMILIES = set()
VARIANT_BOUNDARIES = {}

CROSS_MAX_GAP_BINS = 2       # merge runs separated by at most this many bins
CROSS_MIN_PURITY = 0.5       # drop blocks where under half the span belongs to it
# There is no minimum block size: a run of any length is drawn as long as it is
# solid. Purity is about this tree, not about the family the block came from, which
# would hide a real stretch just because that family's clade is large elsewhere.

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = None          # set from --out in main()
LEAF_ORDER_DIR = None    # a work folder beside the output, not published
CLADE_SCRIPT = os.path.join(os.path.dirname(SCRIPT_DIR), "identify_gene_clade.py")
SEARCH_DIR = stage("01_homology_search")

# -----------------------------------------------------------------------------
# 7TM / rank-based filtering overlay -- which of the initially extracted leaves
# survive DeepTMHMM's 7TM filter, and which additionally survive the rank-based
# cross-family dedup. Both tables are the wide family-by-accession TSVs written
# by DeepTMHMM/analyze_family_overlaps.py and build_deeptmhmm_tables.py.
#
# Column names in those tables are the family's own name except for the three
# seed searches, whose gene_boundaries key carries the source organism.
# -----------------------------------------------------------------------------
FILTER_TABLE_ALIASES = {"GPR1_YEAST": "GPR1", "FSLB_DICDI": "FSLB", "LRLA_DICDI": "LRLA"}

# Mirrors NESTED_MODEL_EXCLUSIONS in step1_clade_extraction/write_baselines.py, in the
# wide table's own column names (so "GPR1", not "GPR1_YEAST"). GPR1 and Vomeronasal 1 are
# more specific models than the class A search they also match, and Vomeronasal 2 is more
# specific than the class C search; their shared sequences are resolved in favour of the
# specific family, not because the 7TM filter removed them from the broad one. That
# resolution belongs at the "Overlaps resolved" stage, so the "After 7TM filter" highlight
# must not read it as a 7TM removal -- see seven_tm_membership_for_highlight().
NESTED_MODEL_EXCLUSIONS = {
    "classA": ["GPR1", "Vomeronasal1"],
    "classC": ["Vomeronasal2"],
}


FILTER_7TM_TSV = str(stage("03_topology_filter", "core_membership.tsv"))
FILTER_TOPOLOGY_TSV = str(stage("03_topology_filter", "topology_membership.tsv"))
FILTER_ASSIGNMENT_TSV = str(stage("04_overlap_assignment", "final_membership.tsv"))
FILTER_REASSIGNMENT_TSV = str(stage("10_motif_analysis", "phylogeny_reassignment",
                                    "final_membership_reassigned.tsv"))
CLADE_MEMBERSHIP_TSV = str(stage("02_clade_extraction", "clade_membership.tsv"))


_WIDE_TABLE_CACHE = {}


def load_wide_accession_table(path):
    """Read a family-per-column accession TSV into {family: set(accessions)}.

    Cached on the path, because every family asks for the same two 7 MB tables and
    re-parsing them per family was most of the run's file reading."""
    if not path or not os.path.isfile(path):
        return {}
    if path in _WIDE_TABLE_CACHE:
        return _WIDE_TABLE_CACHE[path]
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows:
        return {}
    header = rows[0]
    table = {name: set() for name in header}
    for row in rows[1:]:
        for name, value in zip(header, row):
            value = value.strip()
            if value:
                table[name].add(value)
    _WIDE_TABLE_CACHE[path] = table
    return table


_MEMBERSHIP_CACHE = {}


def seven_tm_membership():
    """Each family's membership after the 7TM topology filter.

    This is read straight from the wide table, which is now the whole story. The
    earlier pipeline deduplicated the families against each other after filtering,
    so the post-7TM membership had to be reconstructed by adding every contested
    accession back to each family that claimed it. Nothing is deduplicated any
    more. An accession claimed by more than one family stays in all of them, the
    overlap analysis reports it, and which family keeps it is decided downstream
    rather than applied here.
    """
    if "7tm" not in _MEMBERSHIP_CACHE:
        table = load_wide_accession_table(FILTER_7TM_TSV)
        if table:
            print("  post-7TM membership from %s (%d families)"
                  % (os.path.basename(FILTER_7TM_TSV), len(table)))
        else:
            print("  WARNING: no 7TM table at %s -- the filtered view will be empty"
                  % FILTER_7TM_TSV)
        _MEMBERSHIP_CACHE["7tm"] = table
    return _MEMBERSHIP_CACHE["7tm"]


def seven_tm_membership_for_highlight():
    """The "After 7TM filter" membership as the filter alone would leave it.

    `FILTER_7TM_TSV` already has the nested-model exclusion applied -- GPR1 and
    Vomeronasal 1 removed from class A, Vomeronasal 2 removed from class C -- because
    that rule is applied once, upstream, to the same wide tables everything else reads.
    Read straight through, that makes class A and class C look like the 7TM filter
    itself stripped those sequences, when the filter left them in and it was the
    nested-model resolution (an overlap decision, not a topology one) that moved them
    to the specific family. This adds them back for the two broad families they were
    removed from, so the "After 7TM filter" highlight reflects only the 7TM filter; the
    resolution still shows up, correctly, as a further narrowing at the "Overlaps
    resolved" stage.
    """
    if "7tm_for_highlight" not in _MEMBERSHIP_CACHE:
        base = seven_tm_membership()
        table = {family: set(accessions) for family, accessions in base.items()}
        for broad, specific_list in NESTED_MODEL_EXCLUSIONS.items():
            if broad not in table:
                continue
            for specific in specific_list:
                table[broad] |= base.get(specific, set())
        _MEMBERSHIP_CACHE["7tm_for_highlight"] = table
    return _MEMBERSHIP_CACHE["7tm_for_highlight"]


def topology_membership():
    """Each family's membership after the topology-viewer exclusions."""
    if "topology" not in _MEMBERSHIP_CACHE:
        table = load_wide_accession_table(FILTER_TOPOLOGY_TSV)
        if table:
            print("  post-topology membership from %s (%d families)"
                  % (os.path.basename(FILTER_TOPOLOGY_TSV), len(table)))
        else:
            print("  WARNING: no topology table at %s -- the topology view will be empty"
                  % FILTER_TOPOLOGY_TSV)
        _MEMBERSHIP_CACHE["topology"] = table
    return _MEMBERSHIP_CACHE["topology"]


def topology_membership_for_highlight():
    """The "After topology filter" membership, with the nested-model sequences added
    back to class A and class C for the reason given in
    seven_tm_membership_for_highlight()."""
    if "topology_for_highlight" not in _MEMBERSHIP_CACHE:
        base = topology_membership()
        table = {family: set(accessions) for family, accessions in base.items()}
        for broad, specific_list in NESTED_MODEL_EXCLUSIONS.items():
            if broad not in table:
                continue
            for specific in specific_list:
                table[broad] |= base.get(specific, set())
        _MEMBERSHIP_CACHE["topology_for_highlight"] = table
    return _MEMBERSHIP_CACHE["topology_for_highlight"]


def assignment_membership():
    """Each family's membership once the overlaps were resolved (stage 04)."""
    return load_wide_accession_table(FILTER_ASSIGNMENT_TSV)


def reassignment_membership():
    """The final membership, after the Dikarya Class A sequences moved to GPR1."""
    return load_wide_accession_table(FILTER_REASSIGNMENT_TSV)


# -----------------------------------------------------------------------------
# gene_boundaries, lifted from identify_gene_clade.py without importing it
# (that module pulls in Bio and ete3, which this script does not need).
# -----------------------------------------------------------------------------
def load_gene_boundaries():
    with open(CLADE_SCRIPT, encoding="utf-8") as f:
        tree = ast.parse(f.read())
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == "gene_boundaries":
                    return ast.literal_eval(node.value)
    raise RuntimeError("gene_boundaries not found in " + CLADE_SCRIPT)


def resolve_sources(gene_name, gene_data):
    """Mirror _resolve_paths() in identify_gene_clade.py, and add the matching
    line plot. Returns (tree_path, fasta_path, lineplot_path_or_None)."""
    hmmregion = gene_data["hmmregion"]
    is_m8 = str(hmmregion).upper() == "NA"
    search_name = gene_data.get("search_name", gene_name)
    hmmregion_txt = "hmmregions_" if hmmregion == 1 else ""

    if is_m8:
        stem = search_name + "_db_FAMSA_gt0.1"
        ordered_suffix = "m8ordered"
        lineplot_kind = "bitsANDindex"
    else:
        stem = search_name + "_hmm_E1000_db_" + hmmregion_txt + "FAMSA_gt0.1"
        ordered_suffix = "hmmordered"
        lineplot_kind = "hmmscoreANDindex"

    folder = os.path.join(str(SEARCH_DIR), (gene_data.get("search_folders") or [search_name])[0])
    tree_path = os.path.join(folder, stem + "_" + ordered_suffix + ".tree")
    fasta_path = os.path.join(folder, stem + "_treeordered.fasta.gz")

    # Line plots are named for the search they came from. Only a plot built from
    # this same search may be used: the series are keyed by sequence ID, so a plot
    # from a different search would silently contribute another search's values.
    lineplot_name = stem + "_" + ordered_suffix + "_score_lineplot.html"
    lineplot_path = os.path.join(folder, lineplot_name)

    if not os.path.exists(tree_path):
        raise FileNotFoundError(gene_name + ": missing tree " + tree_path)
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(gene_name + ": missing MSA " + fasta_path)
    if not os.path.exists(lineplot_path):
        return tree_path, fasta_path, None, lineplot_name
    return tree_path, fasta_path, lineplot_path, lineplot_name


# -----------------------------------------------------------------------------
# Newick parsing (iterative -- avoids Python recursion limits on ~100,000 leaves)
# -----------------------------------------------------------------------------
def parse_newick(text):
    s = text.strip()
    if s.endswith(";"):
        s = s[:-1]
    n = len(s)
    i = 0
    stack = [[]]

    def read_label_and_blen(i):
        start = i
        while i < n and s[i] not in "(),;:":
            i += 1
        name = s[start:i] if i > start else None
        blen = 0.0
        if i < n and s[i] == ":":
            i += 1
            start2 = i
            while i < n and s[i] not in "(),;":
                i += 1
            try:
                blen = float(s[start2:i])
            except ValueError:
                blen = 0.0
        return name, blen, i

    while i < n:
        c = s[i]
        if c == "(":
            stack.append([])
            i += 1
        elif c == ")":
            children = stack.pop()
            i += 1
            name, blen, i = read_label_and_blen(i)
            stack[-1].append({"name": name, "blen": blen, "children": children})
        elif c == ",":
            i += 1
        else:
            name, blen, i = read_label_and_blen(i)
            stack[-1].append({"name": name, "blen": blen, "children": []})

    return stack[0][0]


def analyze_tree(root):
    """Iterative Euler tour: assigns depth, leaf_start, leaf_end to every node
    and returns the left-to-right leaf name order."""
    leaf_order = []
    root["depth"] = 0
    stack = [[root, 0]]
    while stack:
        frame = stack[-1]
        node, ci = frame
        if ci == 0:
            if not node["children"]:
                node["leaf_start"] = len(leaf_order)
                node["leaf_end"] = len(leaf_order)
                leaf_order.append(node["name"])
                stack.pop()
                continue
            node["leaf_start"] = len(leaf_order)
        if ci < len(node["children"]):
            child = node["children"][ci]
            child["depth"] = node["depth"] + 1
            frame[1] += 1
            stack.append([child, 0])
        else:
            node["leaf_end"] = len(leaf_order) - 1
            stack.pop()
    return leaf_order


def iter_nodes(root):
    stack = [root]
    while stack:
        node = stack.pop()
        yield node
        stack.extend(node["children"])


# -----------------------------------------------------------------------------
# Clade boundaries -> leaf-index ranges
#
# A node's leaves are contiguous in the tree's own left-to-right order, so the
# most recent common ancestor of the two boundary sequences is exactly one
# [leaf_start, leaf_end] interval. This reproduces the ete3 get_common_ancestor
# and get_leaf_names calls in identify_gene_clade.py without the ete3
# dependency, and yields the highlight ranges directly. Entries ending in -1
# are subtracted.
# -----------------------------------------------------------------------------
def mrca_range(internal_nodes, leaf_index, left_id, right_id):
    if left_id not in leaf_index or right_id not in leaf_index:
        raise KeyError("boundary sequence not in tree: " + left_id + " / " + right_id)
    lo = min(leaf_index[left_id], leaf_index[right_id])
    hi = max(leaf_index[left_id], leaf_index[right_id])
    best = None
    best_span = None
    for node in internal_nodes:
        if node["leaf_start"] <= lo and node["leaf_end"] >= hi:
            span = node["leaf_end"] - node["leaf_start"]
            if best_span is None or span < best_span:
                best = (node["leaf_start"], node["leaf_end"])
                best_span = span
    if best is None:
        raise RuntimeError("no common ancestor found")
    return best


def merge_ranges(ranges):
    ranges = sorted(r for r in ranges if r[1] >= r[0])
    merged = []
    for start, end in ranges:
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def subtract_ranges(included, excluded):
    """Both are lists of inclusive [start, end] leaf-index ranges."""
    result = list(included)
    for ex_start, ex_end in excluded:
        nxt = []
        for start, end in result:
            if ex_end < start or ex_start > end:
                nxt.append((start, end))
                continue
            if start < ex_start:
                nxt.append((start, ex_start - 1))
            if end > ex_end:
                nxt.append((ex_end + 1, end))
        result = nxt
    return merge_ranges(result)


def highlight_ranges_for(root, leaf_index, boundaries):
    internal_nodes = [n for n in iter_nodes(root) if n["children"]]
    included = []
    excluded = []
    for b in boundaries:
        if len(b) < 2:
            continue
        rng = mrca_range(internal_nodes, leaf_index, b[0], b[1])
        if len(b) == 3 and b[2] == -1:
            excluded.append(rng)
        else:
            included.append(rng)
    return subtract_ranges(merge_ranges(included), merge_ranges(excluded))


def crown_depths_for(root, ranges):
    """The depth of the ancestral node each highlight segment hangs from.

    A segment is a contiguous run of leaves, so the shallowest node lying wholly
    inside it is that segment's crown. Subtracting an excluded span can split one
    clade into several segments, and each piece then hangs from its own ancestor
    rather than from the crown of the whole clade, which is why this is computed
    per segment. A segment holding a single leaf has no internal node inside it
    and is given the tip depth.
    """
    depths = [None] * len(ranges)
    for node in iter_nodes(root):
        if not node["children"]:
            continue
        for i, (start, end) in enumerate(ranges):
            if node["leaf_start"] >= start and node["leaf_end"] <= end:
                if depths[i] is None or node["depth"] < depths[i]:
                    depths[i] = node["depth"]
    return depths


# -----------------------------------------------------------------------------
# Cladogram segments
# -----------------------------------------------------------------------------
def build_dendrogram_segments(root, bin_size):
    """Rectangular cladogram: leaf tips aligned to a common depth, internal
    structure pruned below bin_size leaves (matching the render resolution of
    the other tracks)."""
    all_nodes = list(iter_nodes(root))
    max_leaf_depth = max(n["depth"] for n in all_nodes if not n["children"])
    kept = set()
    for n in all_nodes:
        if n["children"] and (n["leaf_end"] - n["leaf_start"] + 1) >= bin_size:
            kept.add(id(n))

    # Anything drawn as a tip -- a real leaf, or an internal node collapsed for
    # holding less than one bin -- sits at the centre of its leaf span. Every
    # internal node that survives sits at the mean of its own children instead,
    # which is how the phyletic-distribution tree in GPCRevolution places its
    # nodes. On an unbalanced split the span centre falls inside the larger
    # subtree, so the apex reads as belonging to that side; the mean of the
    # children puts it between the two, and the vertical stub rising from a
    # parent then meets the middle of the horizontal bar below it.
    def span_mid(node):
        return (node["leaf_start"] + node["leaf_end"] + 1) / 2.0 / bin_size

    x_by_id = {}
    for node in reversed(all_nodes):   # preorder reversed: children before parents
        if not node["children"] or id(node) not in kept:
            x_by_id[id(node)] = span_mid(node)
        else:
            cx = [x_by_id[id(c)] for c in node["children"]]
            x_by_id[id(node)] = sum(cx) / len(cx)

    def x_of(node):
        return x_by_id[id(node)]

    h_segments = []  # [x0, x1, y]
    v_segments = []  # [x, y0, y1]
    for node in all_nodes:
        if id(node) not in kept:
            continue
        child_x = [x_of(c) for c in node["children"]]
        y = node["depth"]
        h_segments.append([round(min(child_x), 3), round(max(child_x), 3), y])
        for c, cx in zip(node["children"], child_x):
            if not c["children"]:
                y_end = max_leaf_depth
            elif id(c) in kept:
                y_end = c["depth"]
            else:
                y_end = max_leaf_depth
            v_segments.append([round(cx, 3), y, y_end])

    return {"max_depth": max_leaf_depth, "h": h_segments, "v": v_segments}


# -----------------------------------------------------------------------------
# Plotly HTML data extraction
# -----------------------------------------------------------------------------
def extract_plotly_traces(html_path):
    with open(html_path, encoding="utf-8") as f:
        content = f.read()
    idx = content.rfind("Plotly.newPlot")
    if idx < 0:
        raise ValueError("no Plotly.newPlot call in " + html_path)
    sub = content[idx:]
    start = sub.find("[{")
    depth = 0
    i = start
    while True:
        c = sub[i]
        if c == "[":
            depth += 1
        elif c == "]":
            depth -= 1
            if depth == 0:
                end = i + 1
                break
        i += 1
    return json.loads(sub[start:end])


def decode_plotly_array(value):
    """Plotly serializes numeric trace arrays either as a plain JSON list or as
    a dict carrying a dtype and base64 bdata. Handle both."""
    if isinstance(value, list):
        return value
    if isinstance(value, dict) and "bdata" in value:
        raw = base64.b64decode(value["bdata"])
        dtype_map = {
            "i1": np.int8, "u1": np.uint8,
            "i2": np.int16, "u2": np.uint16,
            "i4": np.int32, "u4": np.uint32,
            "i8": np.int64, "u8": np.uint64,
            "f4": np.float32, "f8": np.float64,
        }
        return np.frombuffer(raw, dtype=dtype_map[value["dtype"]]).tolist()
    raise ValueError("unrecognized plotly array encoding: " + str(type(value)))


# -----------------------------------------------------------------------------
# MSA raster
# -----------------------------------------------------------------------------
CAT_LUT = np.full(256, 5, dtype=np.uint8)          # 5 = other
for _ch in "-.":
    CAT_LUT[ord(_ch)] = 0                          # gap
for _ch in "AVLIPFMW":
    CAT_LUT[ord(_ch)] = 1                          # hydrophobic
for _ch in "GSTCYNQ":
    CAT_LUT[ord(_ch)] = 2                          # polar / neutral
for _ch in "DE":
    CAT_LUT[ord(_ch)] = 3                          # acidic
for _ch in "KRH":
    CAT_LUT[ord(_ch)] = 4                          # basic
HEXDIGITS = np.frombuffer(b"0123456789abcdef", dtype=np.uint8)


def read_fasta(path):
    seqs = {}
    with open_text(path) as f:
        header = None
        buf = []
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(buf)
                header = line[1:].strip()
                buf = []
            else:
                buf.append(line)
        if header is not None:
            seqs[header] = "".join(buf)
    return seqs


def build_msa_raster(fasta_path, leaf_order, bin_size, n_bins):
    seqs = read_fasta(fasta_path)
    aln_len = len(next(iter(seqs.values())))
    n_leaves = len(leaf_order)

    char_matrix = np.empty((n_leaves, aln_len), dtype="S1")
    missing = 0
    for i, name in enumerate(leaf_order):
        seq = seqs.get(name)
        if seq is None or len(seq) != aln_len:
            missing += 1
            char_matrix[i] = b"-"
        else:
            char_matrix[i] = np.frombuffer(seq.encode("ascii", "replace"), dtype="S1")

    codes = CAT_LUT[char_matrix.view(np.uint8)]
    del char_matrix, seqs
    pad = n_bins * bin_size - n_leaves
    if pad:
        codes = np.pad(codes, ((0, pad), (0, 0)), constant_values=0)
    codes = codes.reshape(n_bins, bin_size, aln_len)

    counts = np.empty((6, n_bins, aln_len), dtype=np.int32)
    for k in range(6):
        counts[k] = (codes == k).sum(axis=1, dtype=np.int32)
    del codes

    majority = counts.argmax(axis=0).astype(np.uint8)          # (n_bins, aln_len)
    frac = counts.max(axis=0) / float(bin_size)
    del counts

    cat_bytes = np.ascontiguousarray((majority + 48).T)        # (aln_len, n_bins)
    frac_idx = np.rint(frac * 15).astype(np.int32).clip(0, 15)
    frac_bytes = np.ascontiguousarray(HEXDIGITS[frac_idx].T)

    cat_rows = [row.tobytes().decode("ascii") for row in cat_bytes]
    frac_rows = [row.tobytes().decode("ascii") for row in frac_bytes]
    return aln_len, cat_rows, frac_rows, missing


# -----------------------------------------------------------------------------
# Per-family build
# -----------------------------------------------------------------------------
# Category colors follow the Taylor (1997) amino-acid color scheme, matched by
# physicochemical group -- hydrophobic to Taylor's L/V/I green, polar/neutral to
# its N/Q cyan, acidic to its D/E red, basic to its K/R blue, and other/ambiguous
# to its C/M yellow -- desaturated toward a muted, dusty tone rather than a bright
# or pale one, so a dense bin x column heatmap doesn't visually overpower the rest
# of the figure.
PALETTE = ["#f2f2f0", "#9AAF9A", "#A7D0D0", "#C79E9E", "#A7B3D2", "#D4D4A9"]
SCORE_TRACES = ("Score", "Score_smoothed", "Normalized_Index", "Normalized_Index_smoothed")


def binned_mean(value_by_id, order, bin_size, n_bins):
    arr = np.array([value_by_id.get(name, np.nan) for name in order], dtype=float)
    pad = n_bins * bin_size - len(arr)
    if pad:
        arr = np.pad(arr, (0, pad), constant_values=np.nan)
    with np.errstate(invalid="ignore"):
        return np.nanmean(arr.reshape(n_bins, bin_size), axis=1)


def clean(values, digits):
    return [None if np.isnan(v) else round(float(v), digits) for v in values]


def per_leaf(value_by_id, order, digits):
    """The unaggregated series, one value per tree leaf in tree order. Binning
    is itself an averaging step, so a bin mean of the unsmoothed series is not
    raw; this is what the score track draws underneath the smoothed line."""
    out = []
    for name in order:
        v = value_by_id.get(name)
        out.append(None if v is None else round(float(v), 3))
    return out


def _ranges_from_positions(positions):
    """Sorted leaf indices -> contiguous (start, end) runs, breaking at every gap.

    Unlike the cross-annotation grouping this never merges across a gap, because a
    filtered-out leaf inside a run is a real hole in what survived, not noise to
    smooth over.
    """
    if not positions:
        return []
    positions = sorted(positions)
    out = []
    start = prev = positions[0]
    for p in positions[1:]:
        if p == prev + 1:
            prev = p
        else:
            out.append((start, prev))
            start = prev = p
    out.append((start, prev))
    return out


def build_highlight_sets(root, leaf_order, ranges, initial_depths, max_depth, bin_size,
                         family_key, n_clade):
    """The initial-extraction highlight, plus the same clade re-tinted for which
    leaves survive the 7TM filter and which additionally survive the rank-based
    cross-family dedup.

    Both filtering stages only ever remove leaves that the initial extraction
    already included, so each later set is a subset of the initial clade's leaf
    positions -- the highlight for a later stage can only shrink, never move
    outside the boundary-picked clade.

    A filtered set therefore carries no depths of its own. Filtering drops leaves
    scattered through the clade, which cuts it into thousands of runs, and giving
    each run its own crown would both dominate the run time (one pass over every
    internal node per run) and redraw the tint at a new height for every gap. The
    tint keeps the initial extraction's height and loses only the vertical strips
    of the leaves that were removed, so what changes between the three views is
    which leaves are still in, not how deep the clade is drawn.
    """
    def to_bins(index_ranges):
        return [[start / bin_size, (end + 1) / bin_size] for start, end in index_ranges]

    sets = {
        "initial": {
            "label": "Initial extraction", "n": n_clade,
            "ranges": to_bins(ranges), "depths": initial_depths, "no_data": False,
        }
    }

    clade_positions = [i for start, end in ranges for i in range(start, end + 1)]
    stages = [
        ("after_7tm", "After 7TM filter", seven_tm_membership_for_highlight),
        ("after_topology", "After topology filter", topology_membership_for_highlight),
        ("after_assignment", "Overlaps resolved", assignment_membership),
        ("after_reassignment", "Final (Dikarya Class A to GPR1)", reassignment_membership),
    ]
    table_key = FILTER_TABLE_ALIASES.get(family_key, family_key)
    for mode, label, membership in stages:
        table = membership() if clade_positions else {}
        if not clade_positions or table_key not in table:
            sets[mode] = {"label": label, "n": 0, "ranges": [], "depths": [], "no_data": True}
            continue
        survivors = table[table_key]
        kept = [i for i in clade_positions if accession(leaf_order[i]) in survivors]
        kept_ranges = _ranges_from_positions(kept)
        sets[mode] = {
            "label": label, "n": len(kept),
            "ranges": to_bins(kept_ranges), "depths": [], "no_data": False,
        }
    return sets


def refresh_highlight_sets_only(gene_name, gene_data):
    """Recompute just the highlight rectangles -- the initial-extraction tint and
    its two filtered re-tints -- and patch them into an already-built bundle.

    The MSA raster and the score/index series never depend on which highlight
    mode is selected, so re-rasterizing an 85 MB alignment to move a rectangle is
    wasted work. This reads only the small ordered-tree file, not the alignment
    or the line plot, and reuses the tree's own max_depth from the existing
    bundle rather than rebuilding the full dendrogram. Requires that family to
    already have a bundle on disk from a normal build_family() run.
    """
    output_name = gene_data.get("output_name", gene_name)
    bundle_path = os.path.join(DATA_DIR, output_name + ".json")
    if not os.path.exists(bundle_path):
        raise FileNotFoundError(
            "no existing bundle for %s -- run a normal build first" % output_name)
    with open(bundle_path, encoding="utf-8") as f:
        bundle = json.load(f)

    tree_path, _fasta_path, _lineplot_path, _lineplot_name = resolve_sources(gene_name, gene_data)
    print("\n=== " + output_name + " (highlights only) ===")
    with open(tree_path, encoding="utf-8") as f:
        root = parse_newick(f.read())
    leaf_order = analyze_tree(root)
    leaf_index = {name: i for i, name in enumerate(leaf_order)}
    n_leaves = len(leaf_order)
    bin_size = max(1, -(-n_leaves // TARGET_BINS))

    pending = not gene_data.get("boundary")
    ranges = []
    if not pending:
        try:
            ranges = highlight_ranges_for(root, leaf_index, gene_data["boundary"])
        except KeyError:
            pending = True
    n_clade = sum(end - start + 1 for start, end in ranges)

    max_depth = bundle["tree"]["max_depth"]
    crown_depths = [d if d is not None else max_depth for d in crown_depths_for(root, ranges)]
    highlight_sets = build_highlight_sets(
        root, leaf_order, ranges, crown_depths, max_depth, bin_size, output_name, n_clade)
    print("        %d leaves in %d range(s)" % (n_clade, len(ranges)))
    for mode, hs in highlight_sets.items():
        print("        %-16s n=%-8d no_data=%s" % (mode, hs["n"], hs["no_data"]))

    bundle["pending"] = pending
    bundle["n_clade_leaves"] = n_clade
    bundle["highlight_ranges"] = highlight_sets["initial"]["ranges"]
    bundle["highlight_depths"] = highlight_sets["initial"]["depths"]
    bundle["highlight_sets"] = highlight_sets
    with open(bundle_path, "w", encoding="utf-8") as f:
        json.dump(bundle, f, separators=(",", ":"))

    # Same shape as build_family()'s return, taking whatever this pass did not
    # recompute (aln_len, has_score, lineplot_expected) from the untouched bundle.
    return {
        "gene": output_name,
        "file": "figures/figure_" + output_name + ".html",
        "n_leaves": n_leaves,
        "n_clade_leaves": n_clade,
        "n_ranges": len(ranges),
        "aln_len": bundle["aln_len"],
        "has_score": bundle["score"] is not None,
        "pending": pending,
        "variant_of": gene_data.get("variant_of"),
        "lineplot_expected": bundle["sources"]["lineplot_expected"],
        "bin_size": bin_size,
        "n_bins": -(-n_leaves // bin_size),
        "leaf_ranges": [[a, b] for a, b in ranges],
    }


def build_family(gene_name, gene_data):
    output_name = gene_data.get("output_name", gene_name)
    tree_path, fasta_path, lineplot_path, lineplot_name = resolve_sources(gene_name, gene_data)
    print("\n=== " + output_name + " ===")

    print("  [1/5] parsing ordered tree ...")
    with open(tree_path, encoding="utf-8") as f:
        root = parse_newick(f.read())
    leaf_order = analyze_tree(root)
    n_leaves = len(leaf_order)
    bin_size = max(1, -(-n_leaves // TARGET_BINS))
    n_bins = -(-n_leaves // bin_size)
    leaf_index = {name: i for i, name in enumerate(leaf_order)}
    print("        %d leaves -> %d bins of %d" % (n_leaves, n_bins, bin_size))

    pending = not gene_data.get("boundary")
    if pending:
        print("  [2/5] no clade boundary defined yet -- nothing highlighted")
        ranges = []
        n_clade = 0
    else:
        print("  [2/5] locating the extracted clade ...")
        try:
            ranges = highlight_ranges_for(root, leaf_index, gene_data["boundary"])
        except KeyError as exc:
            print("        boundary sequence absent from this tree (%s)" % exc)
            print("        nothing highlighted")
            ranges = []
            pending = True
        n_clade = sum(end - start + 1 for start, end in ranges)
        if ranges:
            print("        %d leaves in %d range(s): %s" % (
                n_clade, len(ranges), ", ".join("[%d, %d]" % r for r in ranges)))

    print("  [3/5] building the cladogram ...")
    dendro = build_dendrogram_segments(root, bin_size)
    crown_depths = [d if d is not None else dendro["max_depth"]
                    for d in crown_depths_for(root, ranges)]
    print("        %d horizontal / %d vertical segments" % (
        len(dendro["h"]), len(dendro["v"])))

    highlight_sets = build_highlight_sets(
        root, leaf_order, ranges, crown_depths, dendro["max_depth"], bin_size,
        output_name, n_clade)

    print("  [4/5] rasterizing the alignment ...")
    aln_len, cat_rows, frac_rows, missing = build_msa_raster(
        fasta_path, leaf_order, bin_size, n_bins)
    if missing:
        print("        WARNING: %d leaves had no MSA row (filled with gaps)" % missing)
    print("        raster: %d columns x %d bins" % (aln_len, n_bins))

    score = None
    if lineplot_path:
        print("  [5/5] extracting score / normalized-index series ...")
        by_name = {}
        required = ("Score", "Score_smoothed", "Normalized_Index", "Normalized_Index_smoothed")
        for tr in extract_plotly_traces(lineplot_path):
            if tr.get("name") in SCORE_TRACES:
                xs = decode_plotly_array(tr["x"])
                ys = decode_plotly_array(tr["y"])
                by_name[tr["name"]] = dict(zip(xs, ys))
        if all(k in by_name for k in required):
            score = {
                "score_smoothed_mean": clean(binned_mean(
                    by_name["Score_smoothed"], leaf_order, bin_size, n_bins), 4),
                "index_smoothed_mean": clean(binned_mean(
                    by_name["Normalized_Index_smoothed"], leaf_order, bin_size, n_bins), 4),
                # The raw series are per leaf, not per bin, and are drawn faintly
                # underneath so the smoothed line reads against real scatter.
                "score_raw_leaf": per_leaf(by_name["Score"], leaf_order, 3),
                "index_raw_leaf": per_leaf(by_name["Normalized_Index"], leaf_order, 3),
            }
        else:
            print("        WARNING: expected traces absent -- score track omitted")
    else:
        print("  [5/5] score track omitted -- no line plot at")
        print("        " + lineplot_name)

    bundle = {
        "gene": output_name,
        "search_name": gene_data.get("search_name", gene_name),
        "hmmregion": gene_data["hmmregion"],
        "n_leaves": n_leaves,
        "n_clade_leaves": n_clade,
        "pending": pending,
        "variant_of": gene_data.get("variant_of"),
        "bin_size": bin_size,
        "n_bins": n_bins,
        "aln_len": aln_len,
        "highlight_ranges": [[start / bin_size, (end + 1) / bin_size]
                             for start, end in ranges],
        "highlight_depths": crown_depths,
        "highlight_sets": highlight_sets,
        "msa": {
            "categories": ["gap", "hydrophobic", "polar", "acidic", "basic", "other"],
            "colors": PALETTE,
            "cat_rows": cat_rows,
            "frac_rows": frac_rows,
        },
        "score": score,
        "tree": dendro,
        "sources": {
            "tree": os.path.basename(tree_path),
            "msa_fasta": os.path.basename(fasta_path),
            "lineplot": os.path.basename(lineplot_path) if lineplot_path else None,
            "lineplot_expected": lineplot_name,
        },
    }

    # Kept so the cross-annotation pass can map other families' extracted clades
    # onto this tree without parsing every Newick file a second time.
    with open(os.path.join(LEAF_ORDER_DIR, output_name + ".txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(leaf_order))

    bundle_path = os.path.join(DATA_DIR, output_name + ".json")
    with open(bundle_path, "w", encoding="utf-8") as f:
        json.dump(bundle, f, separators=(",", ":"))
    print("        wrote %s (%.2f MB)" % (bundle_path, os.path.getsize(bundle_path) / 1e6))

    return {
        "gene": output_name,
        "file": "figures/figure_" + output_name + ".html",
        "n_leaves": n_leaves,
        "n_clade_leaves": n_clade,
        "n_ranges": len(ranges),
        "aln_len": aln_len,
        "has_score": score is not None,
        "pending": pending,
        "variant_of": gene_data.get("variant_of"),
        "lineplot_expected": lineplot_name,
        "bin_size": bin_size,
        "n_bins": n_bins,
        "leaf_ranges": [[a, b] for a, b in ranges],
    }


# -----------------------------------------------------------------------------
# Cross-annotation
#
# Several families are extracted from the same search, so one tree often carries
# more than one named clade. Labelling the other families' clades where they fall
# in this tree is what explains why a similar-looking stretch was left out: it was
# taken by a different family. Blocks are found by mapping the other family's
# extracted sequence IDs into this tree's leaf order, so this also works across
# searches wherever the two hit lists overlap.
# -----------------------------------------------------------------------------
def load_leaf_order(gene):
    with open(os.path.join(LEAF_ORDER_DIR, gene + ".txt"), encoding="utf-8") as f:
        return f.read().splitlines()


def group_runs(positions, max_gap):
    """Sorted positions -> [start, end, count] runs, merging gaps up to max_gap."""
    runs = []
    start = end = positions[0]
    count = 1
    for p in positions[1:]:
        if p - end <= max_gap:
            end = p
            count += 1
        else:
            runs.append([start, end, count])
            start = end = p
            count = 1
    runs.append([start, end, count])
    return runs


def accession(leaf_name):
    """'tr|A0A7L0ZC94|A0A7L0ZC94_RAVEN' -> 'A0A7L0ZC94'. Matching is done on the
    accession because the complete member lists are bare accessions, and because
    accessions are unique within every tree."""
    parts = leaf_name.split("|")
    return parts[1] if len(parts) > 2 else leaf_name


def member_accessions(entry, order):
    ids = set(accession(n) for a, b in entry["leaf_ranges"] for n in order[a:b + 1])
    from_clade = len(ids)
    if entry["gene"] in COMPLETE_MEMBER_FILES:
        ids.update(load_wide_accession_table(CLADE_MEMBERSHIP_TSV).get(entry["gene"], set()))
        print("  %s: clade %d + clade membership table -> %d members"
              % (entry["gene"], from_clade, len(ids)))
    return ids


# The three memberships a family can annotate another tree with: the clade as
# first extracted, and what is left of it after each filtering stage.
#
# The blocks themselves are found ONCE, on the widest membership a family has, and
# the later stages only ever subtract from those blocks. Recomputing the runs per
# stage looked right but was not: group_runs merges positions up to a gap and the
# purity gate then drops dilute runs, so removing members could break one wide
# dilute run into several short dense ones and make a rectangle appear at a
# stricter stage. On the GPR1 tree, class A over leaves 4467-4651 held 44 of 185
# leaves (purity 0.24, dropped); after the dedup the survivors formed [4642, 4648],
# 6 of 7 leaves (purity 0.86, kept), so a block materialised out of a filter that
# had only taken sequences away.
CROSS_STAGE_NAMES = ("initial", "after_7tm", "after_topology", "after_assignment",
                     "after_reassignment")


def _stage_positions(ids, idx):
    return sorted(idx[x] for x in ids if x in idx)


def cross_blocks_by_stage(entry, order, sources):
    """Where each source family's sequences fall in this family's tree, per stage.

    Returns {stage: [block, ...]}. Every later stage's blocks sit inside an
    initial block: the span is found once and then trimmed to whichever of its
    members survive, and a block whose members are all gone disappears.
    """
    idx = {accession(name): i for i, name in enumerate(order)}
    bin_size = entry["bin_size"]
    max_gap = max(1, bin_size * CROSS_MAX_GAP_BINS)
    gene = entry["gene"]

    out = {stage: [] for stage in CROSS_STAGE_NAMES}
    for other_name, ids_by_stage, min_purity, only_on in sources:
        if other_name == gene or (only_on and gene not in only_on):
            continue
        pos_by_stage = {stage: _stage_positions(ids_by_stage.get(stage, ()), idx)
                        for stage in CROSS_STAGE_NAMES}
        if not pos_by_stage["initial"]:
            continue

        for start, end, n in group_runs(pos_by_stage["initial"], max_gap):
            # A run too dilute to be honestly painted as a bar is not discarded:
            # it is marked sparse and carries its member positions, and the page
            # draws one line per sequence instead of a rectangle. Nothing a family
            # contributes to another tree is invisible, and nothing dilute is drawn
            # as though it were a contiguous clade.
            initial_sparse = n / float(end - start + 1) < min_purity
            for stage in CROSS_STAGE_NAMES:
                kept = [p for p in pos_by_stage[stage] if start <= p <= end]
                if not kept:
                    continue            # nothing of this family left here
                first, last = kept[0], kept[-1]
                width = last - first + 1
                # The gate is re-tested on what survives, because a block can be
                # emptied out without being emptied: on the GPR1 tree class A held
                # 1,059 of the 1,391 leaves over [3053, 4443], and the dedup left
                # 29 of them spread across 955 leaves. Drawn from first to last
                # survivor that is still a full-width bar, painted across the
                # middle of the GPR1 clade, reading as a large overlap where there
                # are 29 scattered sequences sitting in holes of the clade. This
                # can only ever drop a block: the spans still come from the
                # initial run, so nothing new can appear at a stricter stage.
                # A run that was dilute to begin with stays sparse at every later
                # stage even if its survivors happen to be dense, because the span
                # still comes from the initial run and a bar appearing only after a
                # filter removed sequences would read as an overlap that grew.
                sparse = initial_sparse or len(kept) / float(width) < min_purity
                total = len(ids_by_stage.get(stage, ()))
                block = {
                    "gene": other_name,
                    "bin_start": first / bin_size,
                    "bin_end": (last + 1) / bin_size,
                    "leaf_start": first,
                    "leaf_end": last,
                    "n": len(kept),
                    "sparse": sparse,
                    "purity": round(len(kept) / float(width), 3),
                    # Share of that family's membership at this stage, which for a
                    # family with a complete accession list is the whole family,
                    # not just the part sampled into its own tree.
                    "coverage": round(len(kept) / float(total), 3) if total else 0.0,
                    "complete_set": other_name in COMPLETE_MEMBER_FILES,
                    "extra_set": other_name in EXTRA_CROSS_SETS,
                }
                if sparse:
                    block["positions"] = kept
                out[stage].append(block)

    for stage in CROSS_STAGE_NAMES:
        out[stage].sort(key=lambda b: (b["bin_start"], -b["n"]))
    return out


def compute_cross_annotations(entries):
    clade_ids = {}
    leaf_orders = {}
    for e in entries:
        order = load_leaf_order(e["gene"])
        leaf_orders[e["gene"]] = order
        clade_ids[e["gene"]] = member_accessions(e, order)

    annotating = [e["gene"] for e in entries
                  if not e.get("variant_of") and e["gene"] not in OVERLAP_EXCLUDED]
    # The highlight-corrected table, not the raw one: class A and class C's "after 7TM"
    # cross-annotation must show GPR1/Vomeronasal1/Vomeronasal2 overlap the same way the
    # families' own highlight does, or the cross track and the highlight would disagree
    # about whether the 7TM filter removed those sequences. See
    # seven_tm_membership_for_highlight().
    seven_tm = seven_tm_membership_for_highlight()
    topology = topology_membership_for_highlight()
    assigned = assignment_membership()
    reassigned = reassignment_membership()

    sources = []
    for gene in annotating:
        key = FILTER_TABLE_ALIASES.get(gene, gene)
        after_7tm = seven_tm.get(key, set())
        # The initial membership is the widest the family has anywhere, being its
        # extracted clade, its complete accession list where it has one, and the
        # filtered set. The filtered set should already sit inside the clade, but
        # class A carries accessions that appear in the 7TM table and in neither,
        # and a stage must never hold a sequence the initial stage lacks.
        after_assignment = assigned.get(key, set())
        initial = clade_ids[gene] | after_7tm
        sources.append((gene, {
            "initial": initial,
            "after_7tm": after_7tm,
            "after_topology": topology.get(key, set()),
            "after_assignment": after_assignment,
            "after_reassignment": reassigned.get(key, set()),
        }, CROSS_MIN_PURITY, None))

    # Sets carved out of one family's leaf order carry no filtering status of
    # their own, so they annotate every stage with the same sequences.
    for name, spec in EXTRA_CROSS_SETS.items():
        order = leaf_orders.get(spec["from_tree"])
        if order is None:
            print("  %s: WARNING source tree %s not built" % (name, spec["from_tree"]))
            continue
        ids = set(accession(n) for n in order[spec["leaf_start"]:spec["leaf_end"] + 1])
        sources.append((name, {stage: ids for stage in CROSS_STAGE_NAMES},
                        spec.get("min_purity", CROSS_MIN_PURITY), spec.get("only_on")))
        print("  %s: %d sequences from %s leaves %d-%d"
              % (name, len(ids), spec["from_tree"], spec["leaf_start"], spec["leaf_end"]))

    summary = {}
    for e in entries:
        gene = e["gene"]
        sets = cross_blocks_by_stage(e, leaf_orders[gene], sources)
        summary[gene] = sets["initial"]

        bundle_path = os.path.join(DATA_DIR, gene + ".json")
        with open(bundle_path, encoding="utf-8") as f:
            bundle = json.load(f)
        # `cross` stays the initial set so nothing that reads the old key changes.
        bundle["cross"] = sets["initial"]
        bundle["cross_sets"] = sets
        # Families sampled into their own tree carry a complete accession list, so
        # the page can say how much larger the family is than the tree drawn on it.
        # This is the same set that annotates every other family's tree.
        if gene in COMPLETE_MEMBER_FILES:
            bundle["complete_set"] = {
                "file": COMPLETE_MEMBER_FILES[gene],
                "n": len(clade_ids[gene]),
            }
            # These two families are sampled into their own trees, so a count of
            # leaves is not a count of the family and a share of the tree is not
            # the share a filter actually removed. Class A keeps 73.7% of its
            # complete set through the 7TM filter against 58.0% of its tree, a
            # sixteen-point difference that put it among the most heavily filtered
            # families when it is not. Every stage therefore also carries the count
            # over the complete membership, which is what the page reports.
            key = FILTER_TABLE_ALIASES.get(gene, gene)
            members = clade_ids[gene]
            by_stage = {
                "initial": members,
                "after_7tm": members & seven_tm.get(key, set()),
                "after_topology": members & topology.get(key, set()),
                "after_assignment": members & assigned.get(key, set()),
                "after_reassignment": members & reassigned.get(key, set()),
            }
            for stage, ids in by_stage.items():
                entry = bundle["highlight_sets"].get(stage)
                if entry and not entry.get("no_data"):
                    entry["n_complete"] = len(ids)
            bundle["complete_set"]["n_leaves_note"] = True
        else:
            bundle.pop("complete_set", None)
        with open(bundle_path, "w", encoding="utf-8") as f:
            json.dump(bundle, f, separators=(",", ":"))

    return summary


# -----------------------------------------------------------------------------
# Human sequences
#
# Where the human members of a family fall in its tree is the quickest way for a
# reader to orient themselves in an axis that is otherwise 10,000 to 200,000
# anonymous accessions. Every human leaf is recorded here, with its UniProt entry
# name, and the figure draws a tick for each one. Which of them get a written
# name is decided in the page itself, from the pixel width the track actually has,
# so a family with hundreds of human sequences shows a spread-out subset rather
# than an unreadable pile.
# -----------------------------------------------------------------------------
HUMAN_SUFFIX = "_HUMAN"

# GPR1 is fungi-only, so it has no human orthologs to orient a reader with. It gets
# a small named set of literature reference receptors in the same slot instead: the
# yeast GPR1 seed itself (Q12361) and the one true GprD ortholog in its tree
# (B8NHG3, Aspergillus flavus). This is deliberately not every TrEMBL entry that
# homology auto-annotation happens to label "GPR1" in some other fungus -- those are
# propagated names, not distinct characterized reference receptors.
REFERENCE_OVERRIDES = {
    "GPR1_YEAST": {
        "accessions": {"Q12361": 1, "B8NHG3": 1},
        "noun": "reference",
        "group_title": "Reference sequences (GPR1 / GprD)",
        "track_label": "Reference (GPR1 / GprD)",
    },
}


# Guide sequences: named receptors whose position in a tree is evidence for
# where the clade boundary should be cut. They are not tied to one family,
# because the point of a guide is to be findable in every tree that contains it
# -- GPR1 itself is in the cAMP tree, and Git3 is in the class A tree, and both
# of those placements bear on where the GPR1 boundary belongs.
#
# They ride in the human-marker slot, so the page's existing "Show ... sequences"
# button reveals them, and they carry a priority the page reads as "always write
# this name" and "never drop at a later filter stage". A guide that thins out at
# some window widths, or vanishes when the filter stage changes, cannot guide
# anything. To retire one, delete its line here and rebuild.
GUIDE_MARKERS = {
    "Q12361": "GPR1_YEAST",
    "O94744": "GIT3_SCHPO",
}
GUIDE_PRIORITY = 2


def guide_markers(leaf_order):
    """[[leaf index, entry name, GUIDE_PRIORITY], ...] for the guides present."""
    out = []
    for i, name in enumerate(leaf_order):
        parts = name.split("|")
        acc = parts[1] if len(parts) > 2 else name
        if acc in GUIDE_MARKERS:
            out.append([i, GUIDE_MARKERS[acc], GUIDE_PRIORITY])
    return out


def human_markers(leaf_order):
    """[[leaf index, entry name, 1 if reviewed], ...] for every human leaf."""
    out = []
    for i, name in enumerate(leaf_order):
        parts = name.split("|")
        entry = parts[2] if len(parts) > 2 else name
        if not entry.endswith(HUMAN_SUFFIX):
            continue
        # TrEMBL human entries are named after their own accession, so they carry
        # no gene name to print. They are still human sequences and still get a
        # tick; the page just prefers a reviewed neighbour when it picks which
        # markers to label.
        reviewed = 1 if parts[0] == "sp" else 0
        out.append([i, entry, reviewed])
    return out


def reference_markers(leaf_order, accessions):
    """[[leaf index, entry name, 1 if priority-labelled], ...] for a fixed accession
    set, used in place of human_markers for the families in REFERENCE_OVERRIDES."""
    out = []
    for i, name in enumerate(leaf_order):
        parts = name.split("|")
        accession = parts[1] if len(parts) > 2 else name
        if accession not in accessions:
            continue
        entry = parts[2] if len(parts) > 2 else name
        out.append([i, entry, accessions[accession]])
    return out


def marker_sets_by_stage(gene, leaf_order, markers):
    """{stage: [marker, ...]} for one family's markers.

    A marker is dropped from a stage when the filter for that stage removed its
    sequence from this family, so switching the page to "After 7TM filter" takes
    the tick and its written name away together rather than leaving a marker
    standing for a sequence the filter discarded.
    """
    key = FILTER_TABLE_ALIASES.get(gene, gene)
    kept = {
        "initial": None,                                  # None means "keep all"
        "after_7tm": seven_tm_membership_for_highlight().get(key),
        "after_topology": topology_membership_for_highlight().get(key),
        "after_assignment": assignment_membership().get(key),
        "after_reassignment": reassignment_membership().get(key),
    }
    out = {}
    for stage, ids in kept.items():
        if ids is None:
            out[stage] = markers
        else:
            # A guide is pinned across the stages. It is there to say where a
            # named receptor sits, which is a fact about the tree and not about
            # whether this family kept the sequence at this filter stage.
            out[stage] = [m for m in markers
                          if len(m) > 2 and m[2] >= GUIDE_PRIORITY
                          or accession(leaf_order[m[0]]) in ids]
    return out


def annotate_human(entries):
    """Add the human-sequence (or, for the families in REFERENCE_OVERRIDES, named
    reference-receptor) markers to every family's stored bundle."""
    counts = {}
    for e in entries:
        leaf_order = load_leaf_order(e["gene"])
        override = REFERENCE_OVERRIDES.get(e["gene"])
        markers = (reference_markers(leaf_order, override["accessions"]) if override
                   else human_markers(leaf_order))
        # A guide wins any collision, since GPR1_YEAST is already a reference
        # marker on its own page and would otherwise be listed twice.
        guides = guide_markers(leaf_order)
        if guides:
            pinned = {m[0] for m in guides}
            markers = sorted([m for m in markers if m[0] not in pinned] + guides,
                             key=lambda m: m[0])
        counts[e["gene"]] = len(markers)
        bundle_path = os.path.join(DATA_DIR, e["gene"] + ".json")
        with open(bundle_path, encoding="utf-8") as f:
            bundle = json.load(f)
        bundle["human"] = markers
        bundle["human_sets"] = marker_sets_by_stage(e["gene"], leaf_order, markers)
        bundle.pop("temp_strip", None)   # clears the retired lineage strip
        if override:
            bundle["human_noun"] = override["noun"]
            bundle["human_group_title"] = override["group_title"]
            bundle["human_track_label"] = override["track_label"]
        else:
            # Clear any labels a previous override left behind. Without this the
            # markers revert to human sequences while the page keeps describing them
            # as whatever the override called them, which is how a removed track can
            # go on being named on a page long after it stopped existing.
            for key in ("human_noun", "human_group_title", "human_track_label"):
                bundle.pop(key, None)
        with open(bundle_path, "w", encoding="utf-8") as f:
            json.dump(bundle, f, separators=(",", ":"))
    return counts


def merge_entries(new_entries):
    """Carry forward families built by earlier runs, so that building a single
    family still leaves every page's family switcher complete."""
    store = os.path.join(DATA_DIR, "index.json")
    by_gene = {}
    if os.path.exists(store):
        with open(store, encoding="utf-8") as f:
            for entry in json.load(f):
                by_gene[entry["gene"]] = entry
    for entry in new_entries:
        by_gene[entry["gene"]] = entry
    # Drop excluded families, and any whose bundle has since been deleted.
    entries = [e for e in by_gene.values()
               if e["gene"] not in EXCLUDED_FAMILIES
               and os.path.exists(os.path.join(DATA_DIR, e["gene"] + ".json"))]
    entries.sort(key=lambda e: e["gene"].lower())
    with open(store, "w", encoding="utf-8") as f:
        json.dump(entries, f, indent=1)
    return entries


def variant_families(gene_boundaries):
    """The other hmmregion run of each search, as a family in its own right.

    A variant is a separate analysis, not another view of the same one: it has its
    own tree, its own tree-ordered alignment and its own line plot, and its tree
    orders the sequences differently. So it gets its own page built entirely from
    its own files, rather than having its alignment mapped onto the other tree.
    A variant never inherits the boundary of the run it is paired with: that
    boundary was read off the other tree, and this tree orders the sequences by a
    different signal, so the clade has to be read off this page in its own right.
    A variant therefore starts pending, highlighting nothing, until a boundary
    read off its own page is entered in VARIANT_BOUNDARIES.
    """
    out = {}
    for gene_name, gene_data in gene_boundaries.items():
        if gene_name not in VARIANT_FAMILIES or gene_data.get("variant_of"):
            continue
        hmmregion = gene_data.get("hmmregion")
        if str(hmmregion).upper() == "NA":
            continue                      # m8 searches have no hmmregion variants
        base = gene_data.get("output_name", gene_name)
        alt = dict(gene_data)
        alt["hmmregion"] = 0 if hmmregion == 1 else 1
        # The page is keyed by the variant's own name, so the search it came from
        # has to be carried explicitly or the filenames are built from the key.
        alt["search_name"] = gene_data.get("search_name", gene_name)
        alt["output_name"] = base + ("_fullseq" if hmmregion == 1 else "_hmmregions")
        alt["variant_of"] = base
        alt["boundary"] = VARIANT_BOUNDARIES.get(alt["output_name"], [])
        try:
            resolve_sources(gene_name, alt)
        except FileNotFoundError:
            continue                      # that variant was never built
        out[alt["output_name"]] = alt
    return out


def main():
    global DATA_DIR, LEAF_ORDER_DIR
    argv = sys.argv[1:]
    if "--out" not in argv:
        raise SystemExit("usage: build_family_assignment_data.py --out <folder> [family ...]")
    i = argv.index("--out")
    DATA_DIR = os.path.abspath(argv[i + 1])
    argv = argv[:i] + argv[i + 2:]
    LEAF_ORDER_DIR = tempfile.mkdtemp(prefix="family_assignment_leaf_order_")
    os.makedirs(DATA_DIR, exist_ok=True)
    os.makedirs(LEAF_ORDER_DIR, exist_ok=True)
    gene_boundaries = load_gene_boundaries()

    requested = argv or list(gene_boundaries)
    unknown = [g for g in requested if g not in gene_boundaries]
    if unknown:
        raise SystemExit("unknown family names: " + ", ".join(unknown))

    entries, failed = [], []
    for gene_name in requested:
        try:
            entry = build_family(gene_name, gene_boundaries[gene_name])
            entry["display"] = DISPLAY_NAMES.get(entry["gene"], entry["gene"])
            entry["control"] = entry["gene"] in OVERLAP_EXCLUDED
            entry.pop("file", None)
            entries.append(entry)
        except Exception as exc:                       # keep going over the whole set
            print("  SKIPPED %s: %s: %s" % (gene_name, type(exc).__name__, exc))
            failed.append((gene_name, "%s: %s" % (type(exc).__name__, exc)))

    entries = merge_entries(entries)

    print("\ncross-annotating trees with the other families' clades ...")
    compute_cross_annotations(entries)
    print("\nmarking the human sequences in each tree ...")
    human_counts = annotate_human(entries)
    for entry in entries:
        print("  %-16s %d human leaves" % (entry["gene"], human_counts[entry["gene"]]))
    print("wrote %d family bundles and index.json to %s" % (len(entries), DATA_DIR))
    if failed:
        print("\nFamilies that could not be built:")
        for gene_name, msg in failed:
            print("  %s: %s" % (gene_name, msg))
        sys.exit(1)


if __name__ == "__main__":
    main()
