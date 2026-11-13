# -*- coding: utf-8 -*-
"""
Created on Sun Sep 21 10:48:17 2025

@author: selcuk.1
"""

import argparse
import os
from pathlib import Path

INPUT_FASTA = None
INPUT_NEWICK = None
OUT_PREFIX = None



PIPE_TO_UNDERSCORE = True   # replace '|' with '_' in both FASTA and Newick labels
KEEP_CASE          = True   # keep case in labels (False = lowercase everything)
ZEROIZE_TINY_SCI   = True   # convert "1e-10" style lengths to 0 (like the Ruby script)
PRUNE_TO_INTERSECTION = True  # if FASTA and tree donâ€™t match, prune tree to overlap

SCALE = 1.0   # multiply branch lengths by this factor (e.g., 1/len_alignment if lengths are in substitutions)
WARN_LIMIT = 10.0  # warn/stop if any scaled branch length > this (MAFFT expects substitutions/site)
# ---------------------------------------------------------------------------

import re, os, sys
from collections import OrderedDict

# ---------------------- helpers: FASTA & names ------------------------------
def norm_name(s, pipe_to_underscore=False, keep_case=True):
    token = s.split(" ", 1)[0]  # everything before first space
    if pipe_to_underscore:
        token = token.replace("|", "_")
    token = token.strip().strip("'").strip('"')
    if not keep_case:
        token = token.lower()
    return token

def read_fasta_headers(path, pipe_to_underscore=False, keep_case=True):
    """
    Returns (headers_dict, order_list)
      headers_dict: normalized_name -> (raw_header_line, sequence_string)
      order_list:   [normalized_name, ...] in file order
    """
    headers = OrderedDict()
    name_order = []
    with open(path, "r", encoding="utf-8") as fh:
        header = None
        seq_chunks = []
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    name = norm_name(header[1:], pipe_to_underscore, keep_case)
                    if name in headers:
                        # make unique if duplicates appear
                        i = 2
                        uniq = f"{name}_{i}"
                        while uniq in headers:
                            i += 1
                            uniq = f"{name}_{i}"
                        name = uniq
                    headers[name] = (header, "".join(seq_chunks))
                    name_order.append(name)
                header = line
                seq_chunks = []
            else:
                if header is None:
                    raise ValueError("FASTA does not start with '>' header.")
                seq_chunks.append(line.strip())
        if header is not None:
            name = norm_name(header[1:], pipe_to_underscore, keep_case)
            if name in headers:
                i = 2
                uniq = f"{name}_{i}"
                while uniq in headers:
                    i += 1
                    uniq = f"{name}_{i}"
            else:
                uniq = name
            headers[uniq] = (header, "".join(seq_chunks))
            name_order.append(uniq)
    return headers, name_order

def write_fasta(pairs, out_path):
    with open(out_path, "w", encoding="utf-8") as out:
        for name, (raw_header, seq) in pairs:
            out.write(f">{name}\n")
            out.write(seq + "\n")

# ------------------------- helpers: Newick parse ----------------------------
def read_newick_until_semicolon(path):
    with open(path, "r", encoding="utf-8") as fh:
        buf = []
        for line in fh:
            buf.append(line.strip())
            joined = "".join(buf)
            if ";" in joined:
                return joined.split(";", 1)[0] + ";"
    raise ValueError("No semicolon found in Newick file.")

def sanitize_newick(s, pipe_to_underscore=False, zeroize_tiny_sci=False):
    s = re.sub(r"\[.*?\]", "", s)     # strip bracketed comments
    s = s.replace(" ", "")
    if pipe_to_underscore:
        s = s.replace("|", "_")
    if zeroize_tiny_sci:
        s = re.sub(r"[0-9\.]*[eE]-[0-9]*", "0", s)  # 1e-06 -> 0
    return s

def tokenize_newick(s):
    tokens = []
    i = 0
    buf = ""
    while i < len(s):
        ch = s[i]
        if ch in "(),:;":
            if buf:
                tokens.append(("TEXT", buf))
                buf = ""
            tokens.append((ch, ch))
            i += 1
        else:
            buf += ch
            i += 1
    if buf:
        tokens.append(("TEXT", buf))
    return tokens

def parse_tree(tokens, idx=0):
    """
    Returns (node, next_idx). A node is:
      {"children":[...], "label":str|None, "length":str|None}
    length kept as string to preserve exact numeric text.
    """
    if tokens[idx][0] == "(":
        idx += 1
        children = []
        while True:
            child, idx = parse_tree(tokens, idx)
            children.append(child)
            if tokens[idx][0] == ",":
                idx += 1
                continue
            elif tokens[idx][0] == ")":
                idx += 1
                break
            else:
                raise ValueError("Unexpected token in list: " + str(tokens[idx]))
        # optional internal label
        label = None
        if idx < len(tokens) and tokens[idx][0] == "TEXT":
            label = tokens[idx][1]; idx += 1
        length = None
        if idx < len(tokens) and tokens[idx][0] == ":":
            idx += 1
            if tokens[idx][0] != "TEXT":
                raise ValueError("Expected length after ':' for internal node")
            length = tokens[idx][1]; idx += 1
        return {"children":children, "label":label, "length":length}, idx
    else:
        # leaf
        if tokens[idx][0] != "TEXT":
            raise ValueError("Expected label for leaf")
        label = tokens[idx][1]; idx += 1
        length = None
        if idx < len(tokens) and tokens[idx][0] == ":":
            idx += 1
            if tokens[idx][0] != "TEXT":
                raise ValueError("Expected length after ':' for leaf")
            length = tokens[idx][1]; idx += 1
        return {"children":[], "label":label, "length":length}, idx

def serialize_tree(node):
    if node["children"]:
        inner = ",".join(serialize_tree(ch) for ch in node["children"] if ch is not None)
        s = f"({inner})"
        if node["label"] is not None:
            s += node["label"]
        if node["length"] is not None:
            s += f":{node['length']}"
        return s
    else:
        s = node["label"]
        if node["length"] is not None:
            s += f":{node['length']}"
        return s

def leaf_order_left_to_right(node, out=None):
    if out is None: out = []
    if not node["children"]:
        out.append(node["label"])
        return out
    for ch in node["children"]:
        leaf_order_left_to_right(ch, out)
    return out

def collect_leaf_set(node, out=None):
    if out is None: out = set()
    if not node["children"]:
        out.add(node["label"])
        return out
    for ch in node["children"]:
        collect_leaf_set(ch, out)
    return out

def prune_tree_to_keep(node, keep):
    if not node["children"]:
        return node if node["label"] in keep else None
    pruned_children = []
    for ch in node["children"]:
        kept = prune_tree_to_keep(ch, keep)
        if kept is not None:
            pruned_children.append(kept)
    if not pruned_children:
        return None
    if len(pruned_children) == 1:
        # collapse unary
        return pruned_children[0]
    return {"children": pruned_children, "label": node["label"], "length": node["length"]}

# ---------------------- emit MAFFT steps (Ruby-compatible) ------------------
def resolve_polytomies_to_binary(tree_core):
    """ textual hack similar to Ruby code: for numbered-leaf strings """
    while True:
        hit = re.search(r",([0-9]+):(-?[0-9\.eE]+),([0-9]+):(-?[0-9\.eE]+)", tree_core)
        if not hit:
            break
        i1, l1, i2, l2 = hit.group(1), hit.group(2), hit.group(3), hit.group(4)
        pos = hit.start()
        height = 0
        k = pos
        while k >= 0:
            ch = tree_core[k]
            if ch == ')': height += 1
            elif ch == '(':
                if height == 0: break
                height -= 1
            k -= 1
        tree_core = tree_core[:k] + "(" + tree_core[k+1:]
        tree_core = re.sub(r",([0-9]+):(-?[0-9\.eE]+),([0-9]+):(-?[0-9\.eE]+)", ",XXX", tree_core, count=1)
        tree_core = tree_core.replace("XXX", f"{i1}:{l1}):0,{i2}:{l2}", 1)
    return tree_core

def emit_mafft_steps_and_collapse(tree_core, scale=1.0, warn_limit=10.0):
    """
    Return list of strings "i j len_i len_j" and ensure N-1 lines (if possible).
    """
    steps = []
    core = tree_core[:-1] if tree_core.endswith(";") else tree_core
    pat = re.compile(r"\(([0-9]+):(-?[0-9\.eE]+),([0-9]+):(-?[0-9\.eE]+)\)")
    # pre-resolve aggressively
    prev = None
    while prev != core:
        prev = core
        core = resolve_polytomies_to_binary(core)
    guard = 0
    while "(" in core and guard < 10_000_000:
        m = pat.search(core)
        if not m:
            prev2 = core
            core = resolve_polytomies_to_binary(core)
            m = pat.search(core)
            if not m:
                break
        i, l1, j, l2 = int(m.group(1)), float(m.group(2)), int(m.group(3)), float(m.group(4))
        l1 *= scale; l2 *= scale
        if l1 > warn_limit or l2 > warn_limit:
            raise ValueError(
                "Please check the scale of branch length! MAFFT expects lengths in substitutions/site.\n"
                "If your tree uses 'substitutions', set SCALE = 1/(alignment_length)."
            )
        if j < i:
            i, j = j, i; l1, l2 = l2, l1
        steps.append(f"{i:5d} {j:5d} {l1:10.5f} {l2:10.5f}")
        core = core[:m.start()] + str(i) + core[m.end():]
        guard += 1

    # Add a final merge when the expression ends with a dangling "i:len,j:len)".
    m2 = re.match(r"^\s*([0-9]+):([-0-9\.eE]+),([0-9]+):([-0-9\.eE]+)\)\s*$", core)
    if m2:
        i = int(m2.group(1)); l1 = float(m2.group(2))
        j = int(m2.group(3)); l2 = float(m2.group(4))
        l1 *= scale; l2 *= scale
        if j < i:
            i, j = j, i; l1, l2 = l2, l1
        steps.append(f"{i:5d} {j:5d} {l1:10.5f} {l2:10.5f}")
    return steps

def count_leaves(node):
    if not node["children"]: return 1
    return sum(count_leaves(ch) for ch in node["children"])

def replace_leaf_labels_with_numbers(node, name_to_num):
    if not node["children"]:
        lab = node["label"]
        if lab not in name_to_num:
            raise ValueError(f"Leaf label '{lab}' not in numbering map")
        node["label"] = str(name_to_num[lab])
        return
    for ch in node["children"]:
        replace_leaf_labels_with_numbers(ch, name_to_num)

# ------------------------------- main flow ----------------------------------
def run_pipeline():
    # Read FASTA
    headers_dict, fasta_order = read_fasta_headers(
        INPUT_FASTA, pipe_to_underscore=PIPE_TO_UNDERSCORE, keep_case=KEEP_CASE
    )

    # Read and sanitize Newick
    raw_newick = read_newick_until_semicolon(INPUT_NEWICK)
    newick = sanitize_newick(
        raw_newick,
        pipe_to_underscore=PIPE_TO_UNDERSCORE,
        zeroize_tiny_sci=ZEROIZE_TINY_SCI
    )

    # Parse newick
    tokens = tokenize_newick(newick)
    tree, next_idx = parse_tree(tokens, 0)
    # next token should be ';'
    # Build normalized leaf set/order (labels already had '|' handled in sanitize)
    # Apply the same "before first space" rule and case behavior to each leaf label:
    def renorm_label_inplace(node):
        if not node["children"]:
            node["label"] = norm_name(node["label"], pipe_to_underscore=False, keep_case=KEEP_CASE)
            return
        for ch in node["children"]:
            renorm_label_inplace(ch)
    renorm_label_inplace(tree)

    tree_leaf_set = collect_leaf_set(tree)
    fasta_leaf_set = set(headers_dict.keys())

    # Prune (optional) and choose leaf order for numbering
    if PRUNE_TO_INTERSECTION:
        keep = tree_leaf_set & fasta_leaf_set
        missing_from_fasta = sorted(list(tree_leaf_set - keep))
        missing_from_tree  = sorted(list(fasta_leaf_set - keep))
        if missing_from_fasta:
            print(f"[info] Pruning: {len(missing_from_fasta)} leaves in tree are absent from FASTA (saved lists below).")
        if missing_from_tree:
            print(f"[info] Pruning: {len(missing_from_tree)} FASTA headers absent from tree.")
        pruned_tree = prune_tree_to_keep(tree, keep)
        if pruned_tree is None:
            raise ValueError("After pruning, no leaves remain. Check your inputs.")
        tree = pruned_tree
        # Save mismatch reports (next to OUT_PREFIX)
        base = os.path.splitext(OUT_PREFIX)[0]
        with open(base + ".tree_only_names.txt", "w", encoding="utf-8") as f:
            for x in missing_from_fasta: f.write(x + "\n")
        with open(base + ".fasta_only_names.txt", "w", encoding="utf-8") as f:
            for x in missing_from_tree: f.write(x + "\n")
    else:
        # Require exact match
        if tree_leaf_set != fasta_leaf_set:
            extra_tree = sorted(list(tree_leaf_set - fasta_leaf_set))
            extra_fa   = sorted(list(fasta_leaf_set - tree_leaf_set))
            raise ValueError(
                "FASTA/tree mismatch. Disable PRUNE_TO_INTERSECTION or fix inputs.\n"
                f"  In tree only:  {len(extra_tree)}\n"
                f"  In FASTA only: {len(extra_fa)}"
            )

    # Newick leftâ†’right leaf order (for reordering FASTA)
    leaf_order = leaf_order_left_to_right(tree, out=[])

    # Reorder FASTA to this order
    missing_for_order = [x for x in leaf_order if x not in headers_dict]
    if missing_for_order:
        raise ValueError("Internal mismatch after pruning; missing: " + ", ".join(missing_for_order[:8]))
    reordered = [(name, headers_dict[name]) for name in leaf_order]
    fasta_out = OUT_PREFIX + ".reordered.fasta"
    os.makedirs(os.path.dirname(fasta_out), exist_ok=True)
    write_fasta(reordered, fasta_out)

    # Numbering map: 1..N in reordered FASTA order
    numbered = {name: i+1 for i, (name, _) in enumerate(reordered)}

    # Numbered Newick
    #   (replace leaf labels in the parsed tree, then serialize)
    replace_leaf_labels_with_numbers(tree, numbered)
    numbered_newick_core = serialize_tree(tree)
    numbered_newick = numbered_newick_core + ";"
    newick_out = OUT_PREFIX + ".numbered.newick"
    with open(newick_out, "w", encoding="utf-8") as outn:
        outn.write(numbered_newick + "\n")

    # Ensure bifurcating & emit MAFFT steps
    steps = emit_mafft_steps_and_collapse(numbered_newick_core + ";", scale=SCALE, warn_limit=WARN_LIMIT)

    # Sanity check: steps must be N-1
    N = len(reordered)
    if len(steps) != N - 1:
        print(f"[warn] Steps lines = {len(steps)}; expected N-1 = {N-1}. (A root shape might have hidden a final merge.)")

    steps_out = OUT_PREFIX + ".mafft_steps.txt"
    with open(steps_out, "w", encoding="utf-8") as f:
        for line in steps:
            f.write(line + "\n")

    print("Done.")
    print(f"scale = {SCALE}")
    print(f"FASTA reordered -> {fasta_out}")
    print(f"Numbered Newick -> {newick_out}")
    print(f"MAFFT steps     -> {steps_out}")
    print(f"N (sequences)   -> {N}  |  Steps lines -> {len(steps)} (should be N-1)")

# -------------------------- optional: run in Spyder -------------------------
if __name__ == "__main__":
    phylogeny_dir = Path(__file__).resolve().parent
    run_name = "All_clust0.7_cov1_Oct25_minsize2"

    parser = argparse.ArgumentParser(
        description="Convert a preliminary RAxML-NG tree into MAFFT --treein merge steps."
    )
    parser.add_argument(
        "--fasta",
        default=phylogeny_dir / "data" / "alignments" / f"{run_name}_linsi_ep0.123_trimends.fasta",
        help="Manually trimmed 7TM FASTA used for the preliminary tree.",
    )
    parser.add_argument(
        "--tree",
        default=phylogeny_dir / "data" / "trees" / f"{run_name}_linsi_ep0.123_trimends_kpicgappy.raxml.bestTree",
        help="Preliminary RAxML-NG Newick tree.",
    )
    parser.add_argument(
        "--out-prefix",
        default=phylogeny_dir / "data" / "alignments" / f"{run_name}_linsi_ep0.123_trimends_treein",
        help="Prefix for the reordered FASTA, numbered tree, merge steps, and mismatch reports.",
    )
    args = parser.parse_args()

    INPUT_FASTA = str(Path(args.fasta).resolve())
    INPUT_NEWICK = str(Path(args.tree).resolve())
    OUT_PREFIX = str(Path(args.out_prefix).resolve())
    run_pipeline()

    # If you also want to run MAFFT from inside Spyder (optional), uncomment:
    # import subprocess, shlex
    # cmd = f"mafft --treein {OUT_PREFIX}.mafft_steps.txt {OUT_PREFIX}.reordered.fasta"
    # result = subprocess.run(shlex.split(cmd), capture_output=True, text=True)
    # with open(OUT_PREFIX + ".aligned.fasta", "w", encoding="utf-8") as out_aln:
    #     out_aln.write(result.stdout)
    # print("MAFFT finished ->", OUT_PREFIX + ".aligned.fasta")
