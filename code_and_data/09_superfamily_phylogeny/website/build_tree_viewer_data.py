"""Tree data for the website, one file per final tree in 09_superfamily_phylogeny/trees/ (read by
build_superfamily_trees_data.py).

Nothing about the trees is recomputed here. Every value is lifted from the tree's
<tree>.families.circular.html, which plot_circular_trees.py already built: the node table
(parents, cladogram depths, support, family, lineage), the two ring colours of each leaf,
and the family/lineage legends. Both circular layouts are computed in the viewer page.

Each tree becomes tree_viewer_data/<tree>.js, which sets window.TREE_DATA, and
tree_viewer_data/index.js lists them for the viewer's dropdown. They are scripts rather than
JSON so the viewer can load them with a <script> tag and open straight from disk, where
fetch() of a neighbouring file is blocked.

Usage:
    python build_tree_viewer_data.py
"""
import csv
import json
import re
from pathlib import Path

import sys

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE.parents[1]))
from gpcr_paths import stage  # noqa: E402

HERE = stage("09_superfamily_phylogeny")
# local build folder (not part of the archive): plot_circular_trees.py writes its input there
OUT = CODE / "build" / "tree_viewer_data"
TREE_DIR = CODE / "build" / "circular"

# The 8 final trees (README.md): best ML trees of the three Oct08 rep sets and of the
# literature set, each aligned with ginsi_ep0.123 and with plain ginsi.
FINAL = [f"representatives_{s}_{a}_kpicgappy.raxml"
         for s in ("all", "no_ste2", "no_ste2ste3classCV2", "and_literature")
         for a in ("ginsi_ep0.123", "ginsi")]

RING = re.compile(r'<path data-n="(\d+)"[^>]*fill="(#[0-9a-fA-F]+)"')


def extract(circular):
    s = circular.read_text(encoding="utf-8")
    D = json.loads(re.search(r"^const D = (.*);$", s, re.M).group(1))
    ring = re.search(r'<g class="ring".*?</g>', s, re.S).group(0)
    fills = {}
    for n, c in RING.findall(ring):
        fills.setdefault(int(n), []).append(c)
    for i, d in enumerate(D):
        if d.get("leaf"):
            # family, lineage; with a tIndex ring the plotter writes that wedge first,
            # so it is moved to the end: family, lineage, tIndex
            d["rc"] = fills[i][1:] + fills[i][:1] if len(fills[i]) == 3 else fills[i]
    # Legends and colour bars sit between the family brackets and the filename footer. The
    # lineage legend is the right-hand column (x=2270 on the 2600 canvas), the rest the left.
    lines = s.split("\n")
    start = next(k for k, l in enumerate(lines) if l.startswith('<g class="brackets">'))
    start = next(k for k in range(start, len(lines)) if lines[k] == "</g>") + 1
    left, right = [], []
    for l in lines[start:]:
        if ".nwk</text>" in l:
            break
        (right if l.startswith('<text x="2270"') else left).append(l)
    tis = [d["ti"] for d in D if d.get("leaf") and "ti" in d]
    return {"D": D, "left": "\n".join(left), "right": "\n".join(right),
            "tixMax": max(tis) if tis else None}


def mean_tindex():
    """accession -> mean tIndex over the final trees, and the lowest STE2 mean (the cutoff the
    viewer's "Dim unstable" switch uses: STE2 has no detectable homology to the other families, so a
    sequence placed at least as unstably as the most stable STE2 sequence is not interpreted)."""
    path = HERE / "booster_stats" / "tindex_mean_across_final_trees.tsv"
    if not path.exists():
        return {}, None
    rows = list(csv.DictReader(path.open(), delimiter="\t"))
    ste2 = [float(r["mean_tIndex"]) for r in rows if r["family"] == "STE2"]
    return {r["accession"]: float(r["mean_tIndex"]) for r in rows}, (min(ste2) if ste2 else None)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    MTI, CUT = mean_tindex()
    names = []
    for circular in sorted(TREE_DIR.glob("*.families.circular.html")):
        stem = circular.name[:-len(".families.circular.html")]
        if stem not in FINAL:
            continue
        data = extract(circular)
        data["name"] = stem
        for d in data["D"]:
            if d.get("leaf") and d["acc"] in MTI:
                d["mti"] = MTI[d["acc"]]          # mean tIndex across the final trees
        data["ste2Cut"] = CUT
        # the run is the rep set and its filtering; each run has a best ML tree and may have an MRE
        data["run"] = stem.split("_ginsi")[0]
        # Runs are aligned with ginsi_ep0.123 unless the name says otherwise; another alignment
        # of the same rep set (the literature trees' plain ginsi) is kept apart by its name.
        aln = re.search(r"_(ginsi[^_]*(?:_ep[\d.]+)?)_", stem)
        if aln and aln.group(1) != "ginsi_ep0.123":
            data["run"] += "_" + aln.group(1)
        data["kind"] = "MRE" if "_MRE" in stem else "Best tree"
        data["phylo"] = data["kind"] == "Best tree"      # drawn with branch lengths
        # TBE on a best tree or a _bsTBE consensus, else the share of ML trees holding the split
        data["supLabel"] = "Bipartitions" if stem.endswith("_MRE.raxml") else "TBE"
        (OUT / f"{stem}.js").write_text(
            "window.TREE_DATA = " + json.dumps(data, separators=(",", ":")) + ";\n",
            encoding="utf-8")
        names.append({"run": data["run"], "kind": data["kind"], "file": stem})
        print(f"{stem}: {sum(1 for d in data['D'] if d.get('leaf'))} leaves")
    (OUT / "index.js").write_text("window.TREE_INDEX = " + json.dumps(names, indent=1) + ";\n",
                                  encoding="utf-8")


if __name__ == "__main__":
    main()
