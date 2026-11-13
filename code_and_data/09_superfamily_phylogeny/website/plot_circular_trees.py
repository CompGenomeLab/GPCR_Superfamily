"""Circular SVG and interactive HTML of each best ML tree with annotation rings, the first step of the data of
gpcrevolution.org/superfamily-trees.

Rings, from the tree outward: GPCR family, lineage (highest available taxonomic rank from the
step6 representative table), accession labels, and family names bracketing runs of leaves.
A tree whose booster output <tree>.booster.tbe_stats.txt sits in this folder gets a third ring
between lineage and labels: booster's per-taxon transfer index (tIndex), a heatmap from 0 to the
tree's maximum.
Each tree is rooted on the branch whose side best matches the class A-like, GPR1 and fungal cAMP
leaves (F1, tolerating a few leaves placed elsewhere). The consensus trees carry no branch lengths, so they are drawn as cladograms
with the leaves aligned on the outer circle; the best ML trees keep their branch lengths.
Branches of an MRE tree carry Bipartitions (the % of ML search trees holding the split, read from
its .consensusTreeMRE) and, for a _bsTBE tree, TBE; branches of a best ML tree carry TBE only.
Reads each *.families.nwk in 09_superfamily_phylogeny/trees/ and writes <tree>.circular.svg and <tree>.circular.html
to website/build/circular/ (local, not part of the archive), the input of build_tree_viewer_data.py.
The HTML highlights the root-to-leaf path and shows a tooltip when hovering a leaf, its ring
cells or its label; hovering a branch highlights its clade and shows TBE support and composition.
Clicking a branch marks it; the marked-branch panel copies their IDs (B<n>) as text.
Each branch ID is resolved to its leaves in <tree>.branch_ids.tsv, written alongside.
"""
import csv
import html
import json
import re
import math
from collections import Counter
from pathlib import Path
from ete3 import Tree

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from gpcr_paths import stage  # noqa: E402

HERE = stage("09_superfamily_phylogeny")
TREE_DIR = HERE / "trees"
# drawings are an intermediate of the website data, kept in a local build folder, not in the archive
BUILD = Path(__file__).resolve().parent / "build" / "circular"

# The outgroup is the class A-like + GPR1 + fungal cAMP clade, so the root sits between it and
# the rest of the cAMP-related families. It is found tolerantly: the side of any branch that best
# matches those leaves (F1), so a few unstable leaves placed elsewhere do not drag the root away.
OUTGROUP_FAMILIES = {"classA", "Olfactory", "Vomeronasal1", "classT", "Nematode", "GPR1"}


def in_outgroup(name):
    fam, _, acc = (name.split("|") + ["", ""])[:3]
    return fam in OUTGROUP_FAMILIES or (fam == "cAMP" and lineage(acc) == "Fungi")


def root_on_outgroup(t):
    """Root t on the branch whose side best matches the outgroup leaves (in_outgroup), by F1.

    Literature leaves are left out of the count. Returns the outgroup side after rooting.
    """
    ours = [l for l in t.get_leaves() if not l.name.startswith("LIT_")]
    targets = {id(l) for l in ours if in_outgroup(l.name)}
    if not targets:
        raise SystemExit("no class A, GPR1 or fungal cAMP leaves to root on")
    total, nt = len(ours), len(targets)
    below = {}
    for n in t.traverse("postorder"):
        if n.is_leaf():
            below[n] = (0, 0) if n.name.startswith("LIT_") else (1, int(id(n) in targets))
        else:
            below[n] = tuple(map(sum, zip(*(below[c] for c in n.children))))
    best = None
    for n in t.traverse():
        if n.is_root():
            continue
        k, hit = below[n]
        for size, inter in ((k, hit), (total - k, nt - hit)):
            f1 = 2 * inter / (size + nt) if size else 0
            if best is None or f1 > best[0]:
                best = (f1, n)
    t.set_outgroup(best[1])
    return max(t.children, key=lambda c: sum(id(l) in targets for l in c))


def split_key(names, everyone):
    """One key per bipartition, whichever side it is given by: the side without the first leaf."""
    names = frozenset(names)
    return names if min(everyone) not in names else frozenset(everyone - names)


def load_bipartitions(nwk_path):
    """Split -> % of ML search trees holding it, from the set's .consensusTreeMRE.

    A _bsTBE tree is that consensus with TBE mapped onto it, so every one of its splits is
    found there; a plain consensus tree gets its own labels back.
    """
    if "_MRE" not in nwk_path.name:          # a best ML tree has no consensus to read them from
        return {}
    path = TREE_DIR / (nwk_path.name.split("_MRE")[0] + "_MRE.raxml.consensusTreeMRE")
    t = Tree(path.read_text().strip())
    everyone = frozenset(t.get_leaf_names())
    return {split_key(n.get_leaf_names(), everyone): n.support
            for n in t.traverse() if not n.is_leaf() and not n.is_root()}


SAR = "SAR"
TAXONOMY_JSON = stage("07_representative_selection", "taxonomy_eukaryotes_filtered.json")

RANKS = ["kingdom", "phylum", "subphylum", "superclass", "class", "superorder", "order", "suborder",
         "infraorder", "parvorder", "superfamily", "family_tax", "subfamily", "genus", "species"]

# The lineage ring is read from the representative table (stage 07).
TAX = {}
TABLE = stage("07_representative_selection", "representatives_all_families.tsv")


def load_tax(tree_name):
    with open(TABLE, newline="") as fh:
        return "representatives", {r["accession"]: r for r in csv.DictReader(fh, delimiter="	")}


def load_supergroups():
    """taxID -> the taxonomy JSON's supergroup, which is what defines SAR here."""
    if not TAXONOMY_JSON.exists():
        raise SystemExit(f"{TAXONOMY_JSON} not found; it defines which lineages are SAR")
    with open(TAXONOMY_JSON, encoding="utf-8") as fh:
        return {str(r["taxID"]): str(r.get("Clade", "NA")) for r in json.load(fh)}


SUPERGROUP = load_supergroups()


def lineage(acc):
    r = TAX.get(acc)
    if r is None:
        return "Unknown"
    if SUPERGROUP.get(str(r.get("taxid", "")), "NA").lower() == "sar":
        return SAR
    return next((r[k] for k in RANKS if k in r and r[k] not in UNKNOWN), "Unknown")


FAMILY_COLORS = {
    "classA": "#1f77b4", "Olfactory": "#aec7e8", "Nematode": "#17becf", "classT": "#9edae5",
    "Vomeronasal1": "#393b79", "GPR1": "#6b6ecf", "classC": "#d62728", "Vomeronasal2": "#ff9896",
    "STE2": "#8c564b", "STE3": "#c49c94", "classB1": "#2ca02c", "classB2": "#98df8a", "Mth": "#637939",
    "classF": "#ff7f0e", "FSLB": "#ffbb78", "cAMP": "#9467bd", "GP143": "#e377c2", "TM116": "#f7b6d2",
    "GP157": "#bcbd22",
    # Literature sequences from lineages our genomes do not sample, labelled with the family the
    # study claimed. They take that family's colour, and are told apart from our own leaves by a
    # black outline on their ring cells and a star at their tip.
    "LIT_GP157": "#bcbd22", "LIT_classT": "#9edae5", "LIT_Vomeronasal1": "#393b79",
    "LIT_GPR1like": "#6b6ecf",
}
LIT_MARK = "#000000"
# SAR is drawn as one lineage. Its members are spread across a dozen phyla and classes
# that each contribute a handful of leaves -- Sep21 has thirteen of them over 45 leaves --
# and splitting the ring that finely says less about the tree than grouping them does.
# Membership comes from the taxonomy JSON's supergroup column rather than a list kept
# here, so a SAR lineage this run happens not to sample is still grouped in the next.
LINEAGE_COLORS = {
    "Metazoa": "#4477aa", "Fungi": "#cc6677", "Viridiplantae": "#228833",
    SAR: "#ddcc77",
    "Haptophyta": "#44aa99", "Evosea": "#88ccee", "Discosea": "#0077bb",
    "Heterolobosea": "#332288", "Filasterea": "#117733", "Anaeramoebidae": "#555555",
    # Lineages the Sep21 selection reached that the earlier runs never sampled. Without
    # their own entries they would all fall back to one grey and read as a single group.
    "Choanoflagellata": "#66ccee", "Rhodophyta": "#994455", "Apusomonadida": "#aa7744",
    "Fornicata": "#777733",
    # Lineages that only literature sequences bring in, named by their header's layer.
    "Amoebozoa": "#6699cc", "Apusozoa": "#cc9966", "CRuMs": "#ee7733", "Discoba": "#aa4499",
    "Haptista": "#66bbaa", "Metamonada": "#bbbb55",
}
# A literature leaf's lineage is its header's layer, put on the ring's scheme where it maps onto it.
LIT_LINEAGE = {"Sar": SAR, "Chondrichthyes": "Metazoa", "Cyclostomata": "Metazoa"}


def literature_fields(name):
    """'LIT_GP157|LIT_GP157|LIT0024|Arcella_intermedia|claimed=GPR157|layer=Amoebozoa/13-294'
    -> {'study': ..., 'claimed': ..., 'layer': ...} plus 'species'."""
    parts = name.split("/")[0].split("|")
    f = dict(x.split("=", 1) for x in parts[4:] if "=" in x)
    f["species"] = parts[3].replace("_", " ")
    return f
# Grey, and the same grey every time, for a lineage no palette entry covers. A new run
# can reach lineages the earlier ones never sampled, and an unknown colour should not
# stop the figure being drawn.
LINEAGE_FALLBACK = "#9a9a9a"

# The two selection scripts spell an unknown rank differently: the mmseqs one wrote
# "NA", the topology one writes "__unk__". Both mean the same thing here.
UNKNOWN = {"NA", "__unk__", "", "nan", "None"}

SIZE, TREE_R = 2600, 760
FAM_R0, FAM_R1 = 766, 798
LIN_R0, LIN_R1 = 802, 834
ACC_R, NAME_R = 842, 976
# tIndex ring; when drawn, the labels and family names move out by TIX_SHIFT
TIX_R0, TIX_R1, TIX_SHIFT = 838, 862, 28
# single-hue sequential ramp, white (0) to dark red (the tree's maximum tIndex)
TIX_RAMP = ["#fff5f0", "#fcbba1", "#fb6a4a", "#cb181d", "#67000d"]


def load_tindex(tree_name):
    """Taxon -> booster tIndex from booster_stats/<tree>.tbe_stats.txt (or the older
    <tree>.booster.tbe_stats.txt in this folder), or None if there is none."""
    stem = tree_name.replace(".raxml.families.nwk", "")
    path = HERE / "booster_stats" / f"{stem}.tbe_stats.txt"
    if not path.exists():
        path = HERE / f"{stem}.booster.tbe_stats.txt"
    if not path.exists():
        return None
    out, on = {}, False
    for line in path.read_text().splitlines():
        if line.startswith("Taxon\ttIndex"):
            on = True
        elif on and line.startswith("Edge\t"):
            break
        elif on and line.strip():
            name, val = line.rsplit("\t", 1)
            out[name] = float(val)
    return out


def tix_color(v, vmax):
    x = min(max(v / vmax, 0.0), 1.0) * (len(TIX_RAMP) - 1) if vmax > 0 else 0.0
    k = min(int(x), len(TIX_RAMP) - 2)
    f = x - k
    c0, c1 = TIX_RAMP[k], TIX_RAMP[k + 1]
    return "#" + "".join(f"{round(int(c0[i:i + 2], 16) * (1 - f) + int(c1[i:i + 2], 16) * f):02x}"
                         for i in (1, 3, 5))


def polar(r, a):
    return SIZE / 2 + r * math.cos(a), SIZE / 2 + r * math.sin(a)


def arc_path(r, a0, a1):
    x0, y0 = polar(r, a0); x1, y1 = polar(r, a1)
    large = 1 if abs(a1 - a0) > math.pi else 0
    sweep = 1 if a1 >= a0 else 0
    return f"M{x0:.2f},{y0:.2f} A{r},{r} 0 {large} {sweep} {x1:.2f},{y1:.2f}"


def wedge(r0, r1, a0, a1, color, attrs=""):
    x0, y0 = polar(r0, a0); x1, y1 = polar(r1, a0); x2, y2 = polar(r1, a1); x3, y3 = polar(r0, a1)
    return (f'<path {attrs} d="M{x0:.2f},{y0:.2f} L{x1:.2f},{y1:.2f} A{r1},{r1} 0 0 1 {x2:.2f},{y2:.2f} '
            f'L{x3:.2f},{y3:.2f} A{r0},{r0} 0 0 0 {x0:.2f},{y0:.2f} Z" fill="{color}"/>')


def star_path(cx, cy, r):
    """A five-pointed star centred on (cx, cy), outer radius r; marks a literature leaf's tip."""
    pts = []
    for k in range(10):
        a = -math.pi / 2 + k * math.pi / 5
        rr = r if k % 2 == 0 else r * 0.45
        pts.append(f"{cx + rr * math.cos(a):.2f},{cy + rr * math.sin(a):.2f}")
    return "M" + " L".join(pts) + " Z"


def rotated_text(r, a, text, size, weight="normal", attrs=""):
    x, y = polar(r, a)
    deg = math.degrees(a)
    flip = 90 < deg % 360 < 270
    rot = deg + 180 if flip else deg
    anchor = "end" if flip else "start"
    return (f'<text {attrs} x="{x:.2f}" y="{y:.2f}" font-size="{size}" font-weight="{weight}" '
            f'text-anchor="{anchor}" dominant-baseline="central" '
            f'transform="rotate({rot:.2f} {x:.2f} {y:.2f})">{text}</text>')


def build(nwk_path):
    """Return (svg markup, per-node data) for one tree."""
    global OUTGROUP_COMP
    t = Tree(nwk_path.read_text().strip())
    outgroup = root_on_outgroup(t)
    OUTGROUP_COMP = Counter(l.name.split("|")[0] for l in outgroup).most_common()
    t.ladderize()
    # mirror the outgroup (class A-like + GPR1 + fungal cAMP) side, all the way down, so that
    # side reads from Nematode at the top; display order only
    for n in outgroup.traverse():
        n.children.reverse()
    leaves = t.get_leaves()
    n = len(leaves)
    acc_r, name_r = (ACC_R + TIX_SHIFT, NAME_R + TIX_SHIFT) if TIX else (ACC_R, NAME_R)
    if TIX:
        missing = [l.name for l in leaves if l.name.split("|", 1)[1] not in TIX]
        if missing:
            raise SystemExit(f"{nwk_path.name}: {len(missing)} leaves not in the booster stats, "
                             f"e.g. {missing[0]}")
        tix_max = max(TIX.values())
    step = 2 * math.pi / n
    start = -math.pi / 2

    nodes = list(t.traverse("preorder"))
    nid = {node: i for i, node in enumerate(nodes)}
    if PHYLO:
        # a best ML tree keeps its branch lengths
        depth = {t: 0.0}
        for node in nodes[1:]:
            depth[node] = depth[node.up] + node.dist
    else:
        # Cladogram: a node sits as many steps in from the rim as its longest path down to a
        # leaf, so every leaf lands on the outer circle.
        height = {}
        for node in t.traverse("postorder"):
            height[node] = 0 if node.is_leaf() else 1 + max(height[c] for c in node.children)
        depth = {node: float(height[t] - height[node]) for node in nodes}
        for node in nodes[1:]:
            node.dist = depth[node] - depth[node.up]
    scale = TREE_R / max(depth[l] for l in leaves)
    ang = {l: start + (i + 0.5) * step for i, l in enumerate(leaves)}
    for node in t.traverse("postorder"):
        if not node.is_leaf():
            ang[node] = (ang[node.children[0]] + ang[node.children[-1]]) / 2

    data = []
    for node in nodes:
        d = {"id": f"B{nid[node]}", "p": nid[node.up] if node.up else -1, "bl": round(node.dist, 4),
             "rt": round(depth[node], 4)}
        if node.is_leaf() and node.name.startswith("LIT_"):
            fam, _, acc = node.name.split("|")[:3]
            f = literature_fields(node.name)
            lin = LIT_LINEAGE.get(f.get("layer"), f.get("layer", "Unknown"))
            d.update(leaf=1, lit=1, name=node.name, fam=fam, acc=acc, lin=lin,
                     tax=[[k, f[k].replace("_", " ")] for k in ("species", "claimed", "study", "layer", "matches_our")
                          if k in f],
                     cs="—", csl="literature sequence")
        elif node.is_leaf():
            fam, _, acc = node.name.split("|")[:3]
            r = TAX.get(acc, {})
            if "cluster_size" in r:
                stat, label = r["cluster_size"], "cluster size"
            elif "score" in r:
                stat, label = r["score"], "topology score"
            else:
                stat, label = "\u2014", "selection statistic"
            d.update(leaf=1, name=node.name, fam=fam, acc=acc, lin=lineage(acc),
                     tax=[[k, r[k]] for k in RANKS if k in r and r[k] not in UNKNOWN],
                     cs=stat, csl=label)
        if node.is_leaf() and TIX:
            d["ti"] = TIX[node.name.split("|", 1)[1]]
        if not node.is_leaf():
            below = [l.name.split("|", 1)[1] for l in node]
            d.update(sup=round(node.support, 3), n=len(node),
                     bp=BIPARTITIONS.get(split_key(below, EVERYONE)) if node.up else None,
                     comp=Counter(l.name.split("|")[0] for l in node).most_common())
        data.append(d)

    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{SIZE}" height="{SIZE}" viewBox="0 0 {SIZE} {SIZE}" '
           f'font-family="Arial, Helvetica, sans-serif">', f'<rect width="{SIZE}" height="{SIZE}" fill="white"/>']

    out.append('<g class="guides" stroke="#d9d9d9" stroke-width="0.4">')
    for l in leaves:
        x0, y0 = polar(depth[l] * scale, ang[l]); x1, y1 = polar(FAM_R0, ang[l])
        out.append(f'<line data-n="{nid[l]}" x1="{x0:.2f}" y1="{y0:.2f}" x2="{x1:.2f}" y2="{y1:.2f}"/>')
    out.append('</g>')

    # one path per node: arc along the parent's radius, then the radial branch out to the node
    edges = []
    for node in nodes[1:]:
        rp, r = depth[node.up] * scale, depth[node] * scale
        x1, y1 = polar(r, ang[node])
        edges.append((nid[node], f'{arc_path(rp, ang[node.up], ang[node])} L{x1:.2f},{y1:.2f}'))
    out.append('<g class="edges" stroke="#222" stroke-width="0.9" fill="none">'
               + "".join(f'<path data-n="{i}" d="{d}"/>' for i, d in edges) + '</g>')
    # wider invisible copies make branches easier to hover in the HTML
    out.append('<g class="hit" stroke="transparent" stroke-width="7" fill="none">'
               + "".join(f'<path data-n="{i}" d="{d}"/>' for i, d in edges) + '</g>')

    fams, lins, ring, text, stars = [], [], [], [], []
    for i, l in enumerate(leaves):
        d = data[nid[l]]
        fams.append(d["fam"]); lins.append(d["lin"])
        a0, a1 = start + i * step, start + (i + 1) * step
        attr = f'data-n="{nid[l]}"'
        if TIX:
            ring.append(wedge(TIX_R0, TIX_R1, a0, a1, tix_color(d["ti"], tix_max), attr))
        if d.get("lit"):
            lit_attr = attr + f' stroke="{LIT_MARK}" stroke-width="1.2"'
            ring.append(wedge(FAM_R0, FAM_R1, a0, a1, FAMILY_COLORS[d["fam"]], lit_attr))
            ring.append(wedge(LIN_R0, LIN_R1, a0, a1,
                              LINEAGE_COLORS.get(d["lin"], LINEAGE_FALLBACK), lit_attr))
            label = (f'★ {d["acc"]} · {html.escape(d["tax"][0][1])} · '
                     f'{html.escape(d["tax"][1][1])} (literature)')
            text.append(rotated_text(acc_r, ang[l], label, 8, weight="bold", attrs=attr + f' fill="{LIT_MARK}"'))
            x, y = polar(depth[l] * scale, ang[l])
            stars.append(f'<path d="{star_path(x, y, 9)}" fill="{FAMILY_COLORS[d["fam"]]}" '
                         f'stroke="{LIT_MARK}" stroke-width="1.2"/>')
            continue
        ring.append(wedge(FAM_R0, FAM_R1, a0, a1, FAMILY_COLORS[d["fam"]], attr))
        ring.append(wedge(LIN_R0, LIN_R1, a0, a1,
                          LINEAGE_COLORS.get(d["lin"], LINEAGE_FALLBACK), attr))
        text.append(rotated_text(acc_r, ang[l], f'{d["acc"]} · {html.escape(d["lin"])}', 7, attrs=attr))
    out.append('<g class="ring" stroke="white" stroke-width="0.15">' + "".join(ring) + '</g>')
    out.append('<g class="labels" fill="#333">' + "".join(text) + '</g>')
    if stars:
        out.append('<g class="litstars" pointer-events="none">' + "".join(stars) + '</g>')

    out.append('<g class="brackets">')
    i = 0
    while i < n:
        j = i
        while j + 1 < n and fams[j + 1] == fams[i]:
            j += 1
        # literature leaves carry only their own leaf label
        if j - i + 1 >= 3 and not fams[i].startswith("LIT_"):
            a0, a1 = start + i * step + step * 0.2, start + (j + 1) * step - step * 0.2
            out.append(f'<path d="{arc_path(name_r - 8, a0, a1)}" stroke="{FAMILY_COLORS[fams[i]]}" '
                       f'stroke-width="4" fill="none"/>')
            out.append(rotated_text(name_r, (a0 + a1) / 2, f"{fams[i]} ({j - i + 1})", 16, weight="bold"))
        i = j + 1
    out.append('</g>')

    def legend(x, y, title, colors, present, key=""):
        g = [f'<text x="{x}" y="{y}" font-size="22" font-weight="bold">{title}</text>']
        for k, name in enumerate([c for c in colors if c in present]):
            yy = y + 16 + k * 26
            attr = f' class="lg" data-key="{key}" data-val="{name}"' if key else ""
            lit = name.startswith("LIT_")
            edge = f' stroke="{LIT_MARK}" stroke-width="2"' if lit else ""
            g.append(f'<g{attr}><rect x="{x}" y="{yy}" width="20" height="20" fill="{colors[name]}"{edge}/>'
                     f'<text x="{x + 28}" y="{yy + 15}" font-size="18">{name} ({present.count(name)})</text></g>')
        return "".join(g)

    out.append(legend(30, 50, "Family (ring 1)", FAMILY_COLORS, fams, key="fam"))
    out.append(legend(SIZE - 330, 50, "Lineage (ring 2)", LINEAGE_COLORS, lins, key="lin"))
    if TIX:
        x, y, w = SIZE - 430, SIZE - 150, 360
        stops = "".join(f'<stop offset="{k / (len(TIX_RAMP) - 1):.2f}" stop-color="{c}"/>'
                        for k, c in enumerate(TIX_RAMP))
        out.append(f'<defs><linearGradient id="tixg">{stops}</linearGradient></defs>'
                   f'<text x="{x}" y="{y}" font-size="22" font-weight="bold">Transfer index (ring 3)</text>'
                   f'<text x="{x}" y="{y + 22}" font-size="15" fill="#555">booster tIndex: % of well-supported '
                   f'branches the taxon moves in</text>'
                   f'<rect x="{x}" y="{y + 34}" width="{w}" height="20" fill="url(#tixg)" stroke="#999" '
                   f'stroke-width="0.5"/>'
                   f'<text x="{x}" y="{y + 74}" font-size="16">0</text>'
                   f'<text x="{x + w}" y="{y + 74}" font-size="16" text-anchor="end">{tix_max:.2f}</text>')
    out.append(f'<text x="30" y="{SIZE - 60}" font-size="18">{nwk_path.name}</text>')
    out.append(f'<text x="30" y="{SIZE - 34}" font-size="16" fill="#555">{n} leaves; {NOTE}; '
               f'lineage = highest available rank in the step6 taxonomy table</text>')
    out.append('</svg>')
    ids = [(f"B{nid[node]}", node) for node in nodes[1:]]
    return "\n".join(out), data, ids


HTML = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>__TITLE__</title>
<style>
  html, body { margin: 0; height: 100%; background: #fff; color: #222; font-family: Arial, Helvetica, sans-serif; }
  #bar { position: fixed; top: 0; left: 0; right: 0; padding: 8px 12px; background: #f4f4f4;
         border-bottom: 1px solid #ccc; display: flex; gap: 10px; align-items: center; z-index: 2; font-size: 13px; }
  #bar input { width: 260px; padding: 4px 6px; }
  #wrap { position: absolute; top: 44px; left: 0; right: 0; bottom: 0; overflow: hidden; cursor: grab; }
  #wrap svg { width: 100%; height: 100%; }
  .edges path.sub { stroke: #0a7cff; stroke-width: 2; }
  .guides line.hl { stroke: #e6007e; stroke-width: 1.2; }
  .ring path.hl { stroke: #000; stroke-width: 1.5; }
  .labels text.hl { fill: #e6007e; font-weight: bold; font-size: 11px; }
  .ring path, .labels text, .edges path { cursor: pointer; }
  .hit path { pointer-events: stroke; cursor: pointer; }
  .edges path.mark { stroke: #ff8c00; stroke-width: 4; }
  .edges path.hl { stroke: #e6007e; stroke-width: 3.2; }
  .dim { opacity: 0.1; }
  #marks { position: fixed; left: 10px; bottom: 10px; z-index: 2; background: rgba(255,255,255,.96); border: 1px solid #999;
           border-radius: 4px; padding: 8px 10px; font-size: 12px; max-width: 460px; max-height: 40vh; overflow: auto; display: none; }
  #marks ol { margin: 4px 0; padding-left: 18px; } #marks li { margin: 2px 0; }
  .lg { cursor: pointer; } .lg.off { opacity: 0.3; } .lg.on text { font-weight: bold; }
  #tip { position: fixed; pointer-events: none; background: rgba(255,255,255,.97); border: 1px solid #999;
         border-radius: 4px; padding: 8px 10px; font-size: 12px; max-width: 380px; box-shadow: 0 2px 8px rgba(0,0,0,.2);
         display: none; z-index: 3; line-height: 1.4; }
  #tip b { font-size: 13px; } #tip table { border-collapse: collapse; } #tip td { padding: 0 6px 0 0; vertical-align: top; }
  #tip td:first-child { color: #666; }
</style></head><body>
<div id="bar">
  <strong>__TITLE__</strong>
  <input id="q" placeholder="Search accession, family, lineage or species">
  <span id="hits"></span>
  <span style="color:#666">Scroll to zoom, drag to pan, double-click to reset. Click families or lineages in the legends to show only those.</span>
  <button id="clearSel" style="display:none">Show all</button>
</div>
<div id="wrap">__SVG__</div>
<div id="tip"></div>
<div id="marks"><b>Marked branches (<span id="mcount">0</span>)</b>
  <button id="mcopy">Copy</button> <button id="mclear">Clear</button><ol id="mlist"></ol></div>
<script>
const D = __DATA__;
const svg = document.querySelector('#wrap svg'), tip = document.getElementById('tip');
const byNode = {};
svg.querySelectorAll('[data-n]').forEach(el => (byNode[el.dataset.n] ??= []).push(el));
const children = D.map(() => []);
D.forEach((d, i) => { if (d.p >= 0) children[d.p].push(i); });
const esc = s => String(s).replace(/[&<>"]/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));
let lit = [];
function clear() { lit.forEach(([el, c]) => el.classList.remove(c)); lit = []; }
function mark(i, c) { (byNode[i] || []).forEach(el => { el.classList.add(c); lit.push([el, c]); }); }
function pathTo(i) { for (let k = i; k > 0; k = D[k].p) mark(k, 'hl'); }
function subtree(i) { const st = [i]; while (st.length) { const k = st.pop(); if (k !== i) mark(k, D[k].leaf ? 'hl' : 'sub'); children[k].forEach(c => st.push(c)); } mark(i, 'hl'); }
function leafTip(d) {
  const tax = d.tax.map(([k, v]) => `<tr><td>${k}</td><td>${esc(v)}</td></tr>`).join('');
  return `<b>${esc(d.acc)}</b> &nbsp;<span style="color:#666">branch ${d.id}</span><br>${esc(d.name)}<table>
    <tr><td>family</td><td><b>${esc(d.fam)}</b></td></tr><tr><td>lineage</td><td><b>${esc(d.lin)}</b></td></tr>${tax}
    <tr><td>${d.csl || 'cluster size'}</td><td>${d.cs}</td></tr>
    ${d.ti !== undefined ? `<tr><td>tIndex</td><td><b>${d.ti.toFixed(3)}</b></td></tr>` : ''}</table>`;
}
function nodeTip(d) {
  const comp = d.comp.map(([f, c]) => `${esc(f)} ${c}`).join(', ');
  return `<b>${d.id}</b> &nbsp; <b>__SUP__ ${d.sup}</b> &nbsp; ${d.n} leaves<br><span style="color:#666">${comp}</span>`;
}
svg.addEventListener('mouseover', e => {
  const el = e.target.closest('[data-n]'); if (!el) return;
  const i = +el.dataset.n, d = D[i]; clear();
  if (d.leaf) { pathTo(i); tip.innerHTML = leafTip(d); } else { subtree(i); pathTo(D[i].p); mark(i, 'hl'); tip.innerHTML = nodeTip(d); }
  tip.style.display = 'block';
});
svg.addEventListener('mouseout', e => { if (e.target.closest('[data-n]')) { clear(); tip.style.display = 'none'; applySearch(); } });
svg.addEventListener('mousemove', e => {
  const w = tip.offsetWidth, h = tip.offsetHeight;
  tip.style.left = Math.min(e.clientX + 14, innerWidth - w - 8) + 'px';
  tip.style.top = Math.min(e.clientY + 14, innerHeight - h - 8) + 'px';
});
// search
const q = document.getElementById('q'), hits = document.getElementById('hits');
function applySearch() {
  const s = q.value.trim().toLowerCase(); if (!s) { hits.textContent = ''; return; }
  let n = 0;
  D.forEach((d, i) => { if (d.leaf && (d.name + ' ' + d.lin + ' ' + d.tax.map(x => x[1]).join(' ')).toLowerCase().includes(s)) { n++; pathTo(i); } });
  hits.textContent = n + ' leaves';
}
q.addEventListener('input', () => { clear(); applySearch(); });
// zoom and pan
const full = svg.viewBox.baseVal; let vb = [full.x, full.y, full.width, full.height];
const set = () => svg.setAttribute('viewBox', vb.join(' '));
svg.addEventListener('wheel', e => {
  e.preventDefault();
  const r = svg.getBoundingClientRect(), k = e.deltaY > 0 ? 1.15 : 1 / 1.15;
  const s = Math.max(vb[2] / r.width, vb[3] / r.height);
  const px = vb[0] + (e.clientX - r.left - (r.width - vb[2] / s) / 2) * s, py = vb[1] + (e.clientY - r.top - (r.height - vb[3] / s) / 2) * s;
  vb = [px - (px - vb[0]) * k, py - (py - vb[1]) * k, vb[2] * k, vb[3] * k]; set();
}, { passive: false });
let drag = null, moved = 0;
svg.addEventListener('mousedown', e => { drag = [e.clientX, e.clientY]; moved = 0; });
addEventListener('mouseup', () => { drag = null; });
addEventListener('mousemove', e => {
  if (!drag) return;
  const r = svg.getBoundingClientRect(), s = Math.max(vb[2] / r.width, vb[3] / r.height);
  moved += Math.abs(e.clientX - drag[0]) + Math.abs(e.clientY - drag[1]);
  vb[0] -= (e.clientX - drag[0]) * s; vb[1] -= (e.clientY - drag[1]) * s; drag = [e.clientX, e.clientY]; set();
});
// legend filter: clicking family or lineage entries keeps the selected ones and dims the rest;
// a leaf stays visible if it matches the selected families (if any) and the selected lineages (if any)
const sel = { fam: new Set(), lin: new Set() }, clearSel = document.getElementById('clearSel');
function applySelection() {
  const any = sel.fam.size + sel.lin.size > 0, keep = D.map(() => false);
  D.forEach((d, i) => {
    if (!d.leaf || (sel.fam.size && !sel.fam.has(d.fam)) || (sel.lin.size && !sel.lin.has(d.lin))) return;
    for (let k = i; k >= 0; k = D[k].p) { if (keep[k]) break; keep[k] = true; }
  });
  Object.entries(byNode).forEach(([i, els]) => els.forEach(el => el.classList.toggle('dim', !keep[i])));
  svg.querySelectorAll('.lg').forEach(g => {
    const S = sel[g.dataset.key], on = S.size > 0 && S.has(g.dataset.val);
    g.classList.toggle('on', on); g.classList.toggle('off', S.size > 0 && !on);
  });
  clearSel.style.display = any ? '' : 'none';
}
svg.querySelectorAll('.lg').forEach(g => g.addEventListener('click', () => {
  const S = sel[g.dataset.key]; S.has(g.dataset.val) ? S.delete(g.dataset.val) : S.add(g.dataset.val); applySelection();
}));
clearSel.addEventListener('click', () => { sel.fam.clear(); sel.lin.clear(); applySelection(); });
// branch marking: click a branch to mark or unmark it; the panel lists and copies the marked IDs
const KEY = 'marks:' + document.title, marked = new Set();
function summary(i) {
  const d = D[i];
  if (d.leaf) return `${d.id} | leaf ${d.acc} | ${d.fam} | ${d.lin}`;
  return `${d.id} | __SUP__ ${d.sup} | ${d.n} leaves | ` + d.comp.map(([f, c]) => `${f} ${c}`).join(', ');
}
function renderMarks() {
  svg.querySelectorAll('.edges path.mark').forEach(el => el.classList.remove('mark'));
  marked.forEach(i => (byNode[i] || []).forEach(el => { if (el.closest('.edges')) el.classList.add('mark'); }));
  document.getElementById('mcount').textContent = marked.size;
  document.getElementById('mlist').innerHTML = [...marked].map(i => `<li>${esc(summary(i))}</li>`).join('');
  document.getElementById('marks').style.display = marked.size ? 'block' : 'none';
  try { localStorage.setItem(KEY, JSON.stringify([...marked])); } catch (e) {}
}
try { JSON.parse(localStorage.getItem(KEY) || '[]').forEach(i => { if (D[i]) marked.add(i); }); } catch (e) {}
renderMarks();
svg.addEventListener('click', e => {
  if (moved > 4) return;
  const el = e.target.closest('.hit [data-n], .edges [data-n]'); if (!el) return;
  const i = +el.dataset.n; marked.has(i) ? marked.delete(i) : marked.add(i); renderMarks();
});
document.getElementById('mclear').addEventListener('click', () => { marked.clear(); renderMarks(); });
document.getElementById('mcopy').addEventListener('click', () => {
  const text = document.title + '\\n' + [...marked].map(summary).join('\\n');
  const fallback = () => { const ta = document.createElement('textarea'); ta.value = text; document.body.appendChild(ta);
                           ta.select(); document.execCommand('copy'); ta.remove(); };
  navigator.clipboard ? navigator.clipboard.writeText(text).catch(fallback) : fallback();
});
svg.addEventListener('dblclick', () => { vb = [full.x, full.y, full.width, full.height]; set(); });
</script></body></html>
"""

for p in sorted(TREE_DIR.glob("*.families.nwk")):
    tag, TAX = load_tax(p.name)
    # A _bsTBE tree carries bootstrap TBE; a plain consensus carries the percentage of ML
    # trees that hold the split.
    PHYLO = "_MRE" not in p.name
    SUP = "Bipartitions" if p.name.endswith("_MRE.raxml.families.nwk") else "TBE"
    BIPARTITIONS = load_bipartitions(p)
    TIX = load_tindex(p.name)
    EVERYONE = frozenset(l.split("|", 1)[1] for l in re.findall(r"(?<=[(,])([^(),:;]+)", p.read_text()))
    NOTE = ("best ML tree, rooted on the class A + GPR1 + fungal cAMP clade" if PHYLO else
            "MRE consensus drawn as a cladogram, rooted on the class A + GPR1 + fungal cAMP clade")
    svg, data, ids = build(p)
    BUILD.mkdir(parents=True, exist_ok=True)
    out = BUILD / p.name
    with open(out.with_suffix(".branch_ids.tsv"), "w", encoding="utf-8") as fh:
        fh.write("branch_id\tis_leaf\tbranch_length\tTBE\tn_leaves\tfamilies\tleaves\n")
        for bid, node in ids:
            comp = Counter(l.name.split("|")[0] for l in node).most_common()
            fh.write(f"{bid}\t{int(node.is_leaf())}\t{node.dist:.6f}\t{'' if node.is_leaf() else node.support}\t"
                     f"{len(node)}\t{', '.join(f'{f} {c}' for f, c in comp)}\t{','.join(node.get_leaf_names())}\n")
    out.with_suffix(".circular.svg").write_text(svg, encoding="utf-8")
    title = p.name.replace(".raxml.families.nwk", "")
    page = (HTML.replace("__TITLE__", html.escape(title)).replace("__SVG__", svg)
            .replace("__SUP__", SUP)
            .replace("__DATA__", json.dumps(data, separators=(",", ":"))))
    out.with_suffix(".circular.html").write_text(page, encoding="utf-8")
    leaves_n = sum(1 for d in data if d.get("leaf"))
    size = sum(c for _, c in OUTGROUP_COMP)
    print(f"{p.name}: {leaves_n} leaves, lineages from {tag}; outgroup {size} leaves: "
          + ", ".join(f"{f} {c}" for f, c in OUTGROUP_COMP))
