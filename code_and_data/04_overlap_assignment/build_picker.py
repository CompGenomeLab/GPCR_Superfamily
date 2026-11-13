"""Build the family-assignment picker from the annotated tree.

Stage 04, after annotate_tree_with_families.py. This turns the annotated tree into a single self-contained HTML
page, `picker_template.html` with the tree data injected, which is the tool the
assignment is actually made in.

Three things happen here.

The tree is re-rooted on a family that the tree already resolves cleanly, so that
every other family is read against a fixed, trustworthy reference rather than against
an arbitrary midpoint. The anchor is the largest clade that contains all of that
family's representatives and no leaf belonging to any other family, which is a group
the tree itself defines rather than one chosen by eye. Methuselah is the default
because its representatives and the sequences it shares with class B2 form one such
block.

Each leaf is given the median-normalized score it holds in every family claiming it.
That is the search bit score divided by the median bit score of that family's own
post-filter members, so a typical member of any family is 1.0 and higher is a better
match. It is the tie-breaker the assignment falls back on where the tree gives no
clustering signal, and it is computed the same way as in
analyze_family_overlaps.py.

Layout is precomputed, both as tree depth and as cumulative branch length, so the page
can switch between a cladogram and a phylogram without recomputing anything.

Usage (from this folder, after `source ../../config/paths.sh`):
    python build_picker.py         --tree $GPCR_ZENODO_DIR/04_overlap_assignment/assignment_tree/overlaps_and_reps.with_families.nwk         --after-tsv $GPCR_ZENODO_DIR/03_topology_filter/topology_membership.tsv         --diagnostics $GPCR_ZENODO_DIR/04_overlap_assignment/overlap_diagnostics.tsv         --search-dir $GPCR_ZENODO_DIR/01_homology_search         --selections mrca_selections.tsv         --output overlap_mrca_picker.html
    python serve_picker.py      # saves every decision to mrca_selections.tsv
"""

import argparse
import collections
import csv
import json
import os
import sys
from statistics import median

from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import exists, open_text  # noqa: E402

import newick

SHARED_SEARCH = {"classB1": "classB", "classB2": "classB", "Vomeronasal1": "classT"}
M8_SEARCH = {"GPR1": "GPR1_YEAST", "FSLB": "FSLB_DICDI"}

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TEMPLATE = os.path.join(SCRIPT_DIR, "picker_template.html")


def accession_of(token):
    fields = token.split("|")
    return fields[1].strip() if len(fields) >= 3 else token.strip()


def tag_of(label):
    return label.split("|")[-1]


def read_wide_table(path):
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    header = [name.strip() for name in rows[0]]
    table = {name: set() for name in header}
    for row in rows[1:]:
        for name, value in zip(header, row):
            value = value.strip()
            if value:
                table[name].add(value)
    return header, table


# --------------------------------------------------------------------------- #
# search scores
# --------------------------------------------------------------------------- #
def hmm_scores(path):
    scores = {}
    in_scores = header_seen = reading = False
    with open_text(path) as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("Domain annotation for each sequence:"):
                break
            if not in_scores:
                if line.startswith("Scores for complete sequences"):
                    in_scores = True
                continue
            if not header_seen:
                if line.lstrip().startswith("E-value"):
                    header_seen = True
                continue
            if not reading:
                if set(line.strip()) <= {"-", " "} and "-" in line:
                    reading = True
                continue
            parts = line.split()
            if len(parts) < 9:
                continue
            accession = accession_of(parts[8])
            if accession in scores:
                continue
            try:
                scores[accession] = float(parts[1])
            except ValueError:
                pass
    return scores


def m8_scores(path):
    scores = {}
    with open_text(path) as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            accession = accession_of(parts[1])
            if accession in scores:
                continue
            try:
                scores[accession] = float(parts[11])
            except ValueError:
                pass
    return scores


def scores_for(family, search_dir):
    if family in M8_SEARCH:
        folder = os.path.join(search_dir, M8_SEARCH[family])
        found = sorted(f for f in os.listdir(folder) if f.endswith((".m8", ".m8.gz"))) \
            if os.path.isdir(folder) else []
        if not found:
            sys.exit("no .m8 file for %s in %s" % (family, folder))
        return m8_scores(os.path.join(folder, found[0]))
    search = SHARED_SEARCH.get(family, family)
    path = os.path.join(search_dir, search, "%s_hmmsearch_E1000.txt" % search)
    if not exists(path):
        sys.exit("missing search output for %s: %s" % (family, path))
    return hmm_scores(path)


# --------------------------------------------------------------------------- #
# re-rooting
# --------------------------------------------------------------------------- #
def anchor_node(root, nodes, family):
    """Largest clade holding every representative of `family` and nothing foreign."""
    reps = {n.label for n in nodes
            if not n.children and tag_of(n.label) == family + "*"}
    if not reps:
        sys.exit("no %s* representatives in the tree" % family)
    candidates = [nd for nd in nodes
                  if nd is not root and reps <= nd.leaves
                  and all(family in tag_of(leaf) for leaf in nd.leaves)]
    if not candidates:
        sys.exit("no clean %s clade to root on" % family)
    return max(candidates, key=lambda nd: len(nd.leaves))


def reroot(nodes, adjacency, target, parent, target_length, target_support,
           parent_support):
    """Rebuild the tree rooted on a new node inserted in the target/parent edge."""
    out = []

    def add(parent_id, label, length, support):
        index = len(out)
        out.append({"id": index, "p": parent_id, "c": [], "lab": label,
                    "len": length, "sup": support})
        if parent_id is not None:
            out[parent_id]["c"].append(index)
        return index

    def walk(u, came_from, parent_id, length, support):
        node = nodes[u]
        me = add(parent_id, node.label if not node.children else None, length, support)
        for v, weight, sup in adjacency[u]:
            if v == came_from:
                continue
            walk(v, u, me, weight, sup)
        return me

    half = (target_length or 0.0) / 2.0
    new_root = add(None, None, 0.0, None)
    walk(target, parent, new_root, half, target_support)
    walk(parent, target, new_root, half, parent_support)
    return out, new_root


def main():
    parser = argparse.ArgumentParser(description="Build the assignment picker page")
    parser.add_argument("--tree", required=True, help="annotated newick from annotate_tree_with_families.py")
    parser.add_argument("--after-tsv", required=True)
    parser.add_argument("--diagnostics", default=None,
                        help="optional. Without it the claiming families "
                             "are read off --after-tsv directly.")
    parser.add_argument("--search-dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--root-family", default="Mth",
                        help="family whose clean clade the tree is rooted on")
    parser.add_argument("--selections", default=None,
                        help="an existing mrca_selections.tsv to seed the page with, "
                             "so work already done is not lost on a rebuild")
    parser.add_argument("--selections-other-tree", action="store_true",
                        help="the selections were made on a different tree: carry every "
                             "clade pick over as per-accession decisions, never as a node")
    args = parser.parse_args()

    sys.setrecursionlimit(20000)
    root = newick.leafset(newick.parse(open(args.tree, encoding="utf-8").read()))
    nodes = list(newick.walk(root))

    index_of = {node: i for i, node in enumerate(nodes)}
    adjacency = collections.defaultdict(list)
    for node in nodes:
        for child in node.children:
            weight = child.length if child.length is not None else 0.0
            adjacency[index_of[node]].append((index_of[child], weight, child.support))
            adjacency[index_of[child]].append((index_of[node], weight, node.support))

    anchor = anchor_node(root, nodes, args.root_family)
    composition = collections.Counter(tag_of(leaf) for leaf in anchor.leaves)
    print("rooted on the %s group: %d leaves %s"
          % (args.root_family, len(anchor.leaves), dict(composition)))

    out, new_root = reroot(nodes, adjacency, index_of[anchor],
                           index_of[anchor.parent], anchor.length, anchor.support,
                           anchor.parent.support)

    leaves = [0]

    def layout(i, depth, cumulative):
        node = out[i]
        node["d"] = depth
        node["x"] = cumulative
        if not node["c"]:
            node["y"] = leaves[0]
            leaves[0] += 1
        else:
            for child in node["c"]:
                layout(child, depth + 1, cumulative + (out[child]["len"] or 0.0))
            node["y"] = (out[node["c"][0]]["y"] + out[node["c"][-1]]["y"]) / 2.0

    layout(new_root, 0, 0.0)

    for node in out:
        if node["lab"]:
            parts = node["lab"].split("|")
            node["acc"] = parts[1] if len(parts) > 1 else node["lab"]
            node["nm"] = parts[2] if len(parts) > 2 else ""
            node["t"] = parts[-1]

    # median-normalized score of every contested accession in each claiming family
    families, members = read_wide_table(args.after_tsv)
    ratios = {}
    for family in families:
        scores = scores_for(family, args.search_dir)
        mine = {a: scores[a] for a in members[family] if a in scores}
        if not mine:
            print("  WARNING: no scores for %s" % family)
            ratios[family] = {}
            continue
        middle = median(mine.values())
        ratios[family] = {a: v / middle for a, v in mine.items()} if middle > 0 else {}
        print("  %-14s %6d members, median %8.1f" % (family, len(mine), middle))

    # Which families claim each contested accession. Read from the diagnostics when one
    # is given, otherwise off the wide table, which holds the same information and
    # cannot fall out of step with the table the tree was built from.
    claims = {}
    if args.diagnostics:
        with open(args.diagnostics, newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                claims[row["accession"].strip()] = [
                    f.strip() for f in row["families"].split(",") if f.strip()]
    else:
        for family in families:
            for accession in members[family]:
                claims.setdefault(accession, []).append(family)
        claims = {a: f for a, f in claims.items() if len(f) > 1}

    per_accession = {}
    for accession, claiming in claims.items():
        values = [(f, ratios[f][accession]) for f in claiming
                  if accession in ratios.get(f, {})]
        if values:
            per_accession[accession] = sorted(values, key=lambda v: -v[1])

    for node in out:
        if node.get("acc") in per_accession:
            node["m"] = [[f, round(v, 3)] for f, v in per_accession[node["acc"]]]

    family_names = sorted({node["t"][:-1] for node in out
                           if node.get("t", "").endswith("*")})
    data = {"nodes": out, "root": new_root, "nleaf": leaves[0],
            "families": family_names,
            "maxdepth": max(n["d"] for n in out),
            "maxx": max(n["x"] for n in out)}

    seed = {"picks": {}, "sweep": {}, "override": {}}
    if args.selections and os.path.isfile(args.selections):
        leaf_cache = {}

        def leaves_under(node_id):
            """Accessions of every leaf below node_id, in the tree just built."""
            if node_id in leaf_cache:
                return leaf_cache[node_id]
            node = out[node_id]
            if not node["c"]:
                result = {node["acc"]} if node.get("acc") else set()
            else:
                result = set()
                for child in node["c"]:
                    result |= leaves_under(child)
            leaf_cache[node_id] = result
            return result

        # A clade pick is only a node index, and node indices are assigned fresh by
        # this build, so a pick from a selections file built off a different tree
        # (this one included, once the tree is rebuilt with new or removed leaves)
        # can silently point at the wrong clade instead of the one it named. Each
        # clade row also carries the accession it covered, so that is used to check
        # the old node id still covers the same ground in the new tree; a pick that
        # fails the check is not dropped, it is kept as a per-accession override
        # instead, so nothing already decided is lost or silently misassigned.
        clade_rows = collections.defaultdict(list)
        picks = collections.defaultdict(list)
        fallback = {}
        manual = {}
        ambiguous = {}
        with open(args.selections, newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                kind = row.get("kind")
                if kind == "clade":
                    clade_rows[(row["family"], int(row["node"]))].append(
                        row["accession"])
                elif kind == "sweep":
                    seed["sweep"][row["accession"]] = {
                        "f": row["family"], "node": int(row["node"])}
                elif kind == "manual":
                    # Changed during review, which outranks the clade and sweep passes.
                    manual[row["accession"]] = row["family"]
                elif kind == "excluded":
                    # Set aside on purpose. Without this the next clade covering it
                    # would silently claim it again on the following build.
                    ambiguous[row["accession"]] = None
                # An `ambiguous` row is the other kind of unassigned, the one that
                # simply has not been decided yet, and it seeds nothing. Reading it
                # as an exclusion, which is what the page used to write for both,
                # pinned every undecided sequence shut: a sweep would write the
                # assignment and the exclusion would delete it again on the same
                # redraw, so the label never changed.

        stale = 0
        for (family, node_id), accessions in clade_rows.items():
            # From another tree a node id means nothing, and a large node that happens
            # to contain the old accessions would pass the check below while claiming
            # leaves nobody picked, so every clade pick is carried per accession.
            valid = (not args.selections_other_tree
                     and 0 <= node_id < len(out) and set(accessions) <= leaves_under(node_id))
            if valid:
                if node_id not in picks[family]:
                    picks[family].append(node_id)
            else:
                stale += 1
                for accession in accessions:
                    fallback[accession] = family
        seed["picks"] = dict(picks)

        # Manual review decisions always win, however the old clade pick fared.
        seed["override"] = {**fallback, **manual, **ambiguous}

        # A decision for a family that no longer claims the accession (the claims can
        # change between builds, e.g. after more sequences are filtered out) would
        # assign it somewhere it cannot go, so it is dropped and listed instead.
        invalid = []
        for store in (seed["override"], seed["sweep"]):
            for accession in list(store):
                value = store[accession]
                family = value["f"] if isinstance(value, dict) else value
                if accession not in claims:
                    del store[accession]          # no longer contested at all
                elif family is not None and family not in claims[accession]:
                    invalid.append((accession, family, claims[accession]))
                    del store[accession]
        for accession, family, claiming in invalid:
            print("  dropped: %s was assigned to %s, now claimed only by %s"
                  % (accession, family, ", ".join(claiming)))
        print("seeded %d families, %d sweep, %d review decisions"
              % (len(seed["picks"]), len(seed["sweep"]), len(seed["override"])))
        if stale:
            print("  %d clade pick(s) no longer matched the new tree; their "
                  "%d accession(s) were carried over as individual overrides instead"
                  % (stale, len(fallback)))

    template = open(TEMPLATE, encoding="utf-8").read()
    page = (template.replace("__TREE_DATA__", json.dumps(data, separators=(",", ":")))
                    .replace("__SEED__", json.dumps(seed, separators=(",", ":"))))
    with open(args.output, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(page)
    print("\n%d nodes, %d leaves, %d families -> %s"
          % (len(out), leaves[0], len(family_names), args.output))


if __name__ == "__main__":
    main()
