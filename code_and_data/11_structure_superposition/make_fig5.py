"""Figure 5 (structural overlap of the shared motifs) as one SVG, replacing the Illustrator assembly.

Layout of the earlier hand-made figure: the three lineages on top (coloured bar, lineage name, its
families), below them six helix pills with the per-position overlap renders placed on them.
All 18 lineage families are listed; STE2 is not shown.

Panel places, sizes, angles and the GPCRdb number circles come from fig4_layout.json, saved in the
layout tool (python build_fig5_layout_tool.py -> fig4_layout_tool.html); defaults otherwise.

Usage: python make_fig5.py   (writes fig5_structure.svg)
"""
import csv
import json
import statistics
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import stage  # noqa: E402

HERE = stage("11_structure_superposition")
OUT = HERE / "fig5_structure.svg"

W, H = 975, 683          # the earlier figure's size
FONT = "Arial"
TITLE_PT, NAME_PT = 18, 17
# lineage colours of the earlier figure (bar, and later the residues)
COLOURS = {"Rhodopsin": "#F99696", "cAMP": "#BCBCFA", "Glutamate": "#E2BB94"}
# lineage -> (title, family rows in display names, centre x)
HEADER = [
    ("Rhodopsin", "Rhodopsin-related", [["Class A", "Olfactory", "Taste-2"], ["Vomeronasal 1", "Nematode"]], 130),
    ("cAMP", "cAMP-related", [["cAMP", "GPR1", "GP157", "GP143", "STE3", "TM116"],
                              ["Class F", "FSL", "Class B1", "Class B2", "Mth"]], 500),
    ("Glutamate", "Glutamate-related", [["Class C", "Vomeronasal 2"]], 850),
]
NAME_GAP = 12            # space between family names in a row
# six helix pills, as in the earlier figure
PILL_X0, PILL_W, PILL_STEP = 58, 104, 150.4
PILL_Y0, PILL_Y1 = 172, 664
PILL_FILL = "#EDE6E6"
# overlap renders (render_fig4_panels.py), linked from fig4_panels/. Their place, size and in-plane
# angle, and the GPCRdb number circles, come from fig4_layout.json (saved in fig4_layout_tool.html);
# without it, the defaults below (panel centre x, y: the user's layout of Oct 15)
PANEL_DIR = "fig4_panels"
LAYOUT = HERE / "fig4_layout.json"
PANEL_SIZE = 190
PLACES = {"7x49+7x53": (105, 455), "1x57": (139.5, 641), "2x42+2x46+2x50": (250.5, 384.5), "3x25+45x50": (500, 215), "3x39": (362.5, 452), "3x49+3x51+3x54": (400.5, 591), "4x50+4x53": (534, 385.5), "P5x50": (704.5, 404), "5x58": (727, 556), "5x65": (719.5, 640), "6x48": (829.5, 355.5), "6x44": (834, 434.5)}
# GPCRdb number circles per panel: (number, scheme subscript); class A numbers carry none, as in
# the earlier figure; other schemes as the text names them
BADGES = {"7x49+7x53": [("7x49", ""), ("7x53", "")], "1x57": [("1x61", "C")],
          "2x42+2x46+2x50": [("2x38", "C"), ("2x46", ""), ("2x50", "")],
          "3x25+45x50": [("3x25", ""), ("45x50", "")], "3x39": [("3x39", "")],
          "3x49+3x51+3x54": [("3x49", ""), ("3x51", ""), ("3x54", "")], "4x50+4x53": [("4x50", ""), ("4x53", "B"), ("4x43", "C")],
          "P5x50": [("5x50", ""), ("5x42", "B")],
          "5x58": [("5x58", "")], "5x65": [("5x65", "")], "6x48": [("6x48", "")], "6x44": [("6x44", "")]}
BADGE_R, BADGE_PT, BADGE_SUB_PT = 17, 11, 8
# distance labels, one per circle, placed on their own in the layout tool, bold (user, Oct 15): median C-alpha distance of the carriers' residues
# from the centroid of the others (residue_dev of local_backbone_check), as the Results text gives it.
# Circle (number, subscript) -> (position in the table, lineages counted; None: all carriers).
# 4x53B: the cAMP-related prolines; 4x43C: the plant Class C proline's distance from them (text, Oct 15)
DIST_TSV = HERE / "local_backbone_check_gp157.tsv"
DIST_POS = {("1x61", "C"): ("1x57", None), ("2x38", "C"): ("2x42", None), ("4x53", "B"): ("4x53", {"cAMP"}),
            ("4x43", "C"): ("4x53", {"Glutamate"}), ("5x50", ""): ("P5x50", None), ("5x42", "B"): ("P5x461", None)}
DIST_PT, DIST_FILL = 9, "#000"
DIST_BG_OPACITY, DIST_PAD = 0.6, 2# white box behind each distance label (user, Oct 15); padding in pt
# motif names, as in the earlier figure (bold, one shared size): panel they start next to, text, angle
LABELS = [("3x25+45x50", "S-S Bridge", 0), ("2x42+2x46+2x50", "Na+", 0), ("7x49+7x53", "NPxxY", -90),
          ("3x49+3x51+3x54", "DRY", 0), ("4x50+4x53", "Cholesterol", 0), ("6x48", "Toggle Switch", 0),
          ("6x44", "PIF", 0), ("3x39", "Na+", 0, 2)]   # 4th field: layout version the label was added in
LABEL_PT = 17
# double-headed distance arrows (user, Oct 19), horizontal, label centred just above the line; place,
# length and label text set in the layout tool. 3x54-5x65:
# closest heavy-atom distance in the GP157 anchor (I108-I196); median 3.8 A over the 15 structures
# carrying both motifs
ARROWS = [{"text": "4.0 Å", "x1": 477.5, "y1": 627, "x2": 628.5, "y2": 627}]
ARROW_W, ARROW_HEAD = 2.25, 8   # line width, head length
ARROW_GAP = 1.5                 # between the line and the label's box
LAYOUT_VERSION = 2   # saved layouts of an older version get the labels added since


def default_layout():
    return {"panels": {p: {"x": x, "y": y, "size": PANEL_SIZE, "angle": 0, "hidden": False}
                       for p, (x, y) in PLACES.items()},
            "badges": [{"text": t, "sub": sub, "x": PLACES[p][0] + 60, "y": PLACES[p][1] - 30 + 40 * i}
                       for p, lst in BADGES.items() for i, (t, sub) in enumerate(lst)],
            "badge_r": BADGE_R, "badge_pt": BADGE_PT,
            "labels": [{"text": t, "x": PLACES[p][0] - 70 if a else PLACES[p][0],
                        "y": PLACES[p][1] if a else PLACES[p][1] - 85, "angle": a, "since": since[0] if since else 1}
                       for p, t, a, *since in LABELS],
            "label_pt": LABEL_PT, "version": LAYOUT_VERSION}


def layout():
    """The saved layout (fig4_layout.json), with defaults for anything it lacks (new panels, labels)."""
    lay = default_layout()
    if LAYOUT.exists():
        saved = json.loads(LAYOUT.read_text(encoding="utf-8"))
        panels = {**lay["panels"], **saved.get("panels", {})}
        if "labels" in saved:
            saved["labels"] += [t for t in lay["labels"] if t["since"] > saved.get("version", 1)]
        saved["version"] = LAYOUT_VERSION
        lay.update(saved)
        # deleted in the layout tool: kept in the file (so the defaults do not bring them back), not drawn
        lay["panels"] = {k: v for k, v in panels.items() if k in PLACES and not v.get("deleted")}
    if "dists" not in lay:
        lay["dists"] = default_dists(lay["badges"], lay.get("badge_r", BADGE_R))
    lay.setdefault("dist_pt", DIST_PT)
    lay.setdefault("arrows", [dict(a) for a in ARROWS])
    return lay


def width(text, size, bold=False):
    prop = FontProperties(family=FONT, weight="bold" if bold else "normal")
    return TextPath((0, 0), text, size=size, prop=prop).get_extents().width


def text(x, y, s, size, bold=False):
    weight = ' font-weight="bold"' if bold else ""
    return (f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="middle" font-family="{FONT}" '
            f'font-size="{size}"{weight}>{s}</text>')


def header():
    out = []
    for lineage, title, rows, cx in HEADER:
        out.append(text(cx, 18, title, TITLE_PT, bold=True))
        # the bar spans the lineage's widest line (title or family row)
        row_w = [sum(width(n, NAME_PT) for n in r) + NAME_GAP * (len(r) - 1) for r in rows]
        half = max([width(title, TITLE_PT, bold=True)] + row_w) / 2 + 4
        out.append(f'<line x1="{cx - half:.1f}" y1="27" x2="{cx + half:.1f}" y2="27" '
                   f'stroke="{COLOURS[lineage]}" stroke-width="6" stroke-linecap="round"/>')
        for i, (r, rw) in enumerate(zip(rows, row_w)):
            x = cx - rw / 2
            for n in r:
                w = width(n, NAME_PT)
                out.append(text(x + w / 2, 50 + i * 21, n, NAME_PT))
                x += w + NAME_GAP
    return out


def pills():
    return [f'<rect x="{PILL_X0 + i * PILL_STEP:.1f}" y="{PILL_Y0}" width="{PILL_W}" '
            f'height="{PILL_Y1 - PILL_Y0}" rx="{PILL_W / 2}" fill="{PILL_FILL}"/>' for i in range(6)]


def panels(lay):
    out = []
    for p, v in lay["panels"].items():
        if v.get("hidden") or not (HERE / PANEL_DIR / f"{p}.png").exists():
            continue
        x, y, s = v["x"], v["y"], v["size"]
        out.append(f'<image href="{PANEL_DIR}/{p}.png" x="{x - s / 2:.1f}" y="{y - s / 2:.1f}" width="{s:.1f}" '
                   f'height="{s:.1f}" transform="rotate({v.get("angle", 0):.1f} {x:.1f} {y:.1f})"/>')
    return out


def badge(b, r, pt):
    """A GPCRdb number circle: white disc, dark ring, the number with its scheme subscript."""
    sub = (f'<tspan dy="{pt * 0.3:.1f}" font-size="{pt * BADGE_SUB_PT / BADGE_PT:.1f}">{b["sub"]}</tspan>'
           if b.get("sub") else "")
    return (f'<circle cx="{b["x"]:.1f}" cy="{b["y"]:.1f}" r="{r}" fill="white" stroke="#222" stroke-width="1.2"/>'
            f'<text x="{b["x"]:.1f}" y="{b["y"] + pt * 0.35:.1f}" text-anchor="middle" font-family="{FONT}" '
            f'font-size="{pt}">{b["text"]}{sub}</text>')


def distances():
    """Circle "number|subscript" -> its distance label: the median distance (A) of the carriers'
    residues from the centroid of the others (residue_dev), as the Results text reports it."""
    rows = [r for r in csv.DictReader(open(DIST_TSV), delimiter="\t") if r["role"] == "carrier"]
    out = {}
    for text_, sub in {(t, s) for lst in BADGES.values() for t, s in lst}:
        pos, lins = DIST_POS.get((text_, sub), (text_, None))
        d = [float(r["residue_dev"]) for r in rows
             if r["position"] == pos and (lins is None or r["lineage"] in lins)]
        if d:
            out[f"{text_}|{sub}"] = f"{statistics.median(d):.1f} Å"
    return out


def badges(lay):
    return [badge(b, lay.get("badge_r", BADGE_R), lay.get("badge_pt", BADGE_PT)) for b in lay["badges"]]


def default_dists(badges_, r):
    """Distance labels placed under their circles (layouts saved before they existed)."""
    dist = distances()
    return [{"key": k, "x": b["x"], "y": b["y"] + r + DIST_PT * 0.7} for b in badges_
            if (k := f"{b['text']}|{b.get('sub', '')}") in dist]


def dists(lay):
    """The distance labels, centred on (x, y), bold, on a faint white box; the value comes from
    distances() by its circle key."""
    pt, dist = lay.get("dist_pt", DIST_PT), distances()
    return [dist_label(d["x"], d["y"], dist[d["key"]], pt) for d in lay.get("dists", []) if d["key"] in dist]


def dist_label(x, y, s, pt):
    """A distance label centred on (x, y), bold, on a faint white box."""
    w, h = width(s, pt, bold=True) + 2 * DIST_PAD, 0.8 * pt + 2 * DIST_PAD
    return (f'<rect x="{x - w / 2:.1f}" y="{y - h / 2:.1f}" width="{w:.1f}" height="{h:.1f}" '
            f'fill="white" fill-opacity="{DIST_BG_OPACITY}"/>'
            f'<text x="{x:.1f}" y="{y + 0.35 * pt:.1f}" text-anchor="middle" font-family="{FONT}" '
            f'font-size="{pt}" font-weight="bold" fill="{DIST_FILL}">{s}</text>')


def arrow_shape(a):
    """A double-headed arrow from (x1, y1) to (x2, y2): the shaft between the two heads, and the heads."""
    x1, y1, x2, y2 = a["x1"], a["y1"], a["x2"], a["y2"]
    n = max(np.hypot(x2 - x1, y2 - y1), 1e-6)
    ux, uy = (x2 - x1) / n, (y2 - y1) / n
    hw = ARROW_HEAD * 0.45
    heads, bases = [], []
    for x, y, sg in ((x1, y1, 1), (x2, y2, -1)):
        bx, by = x + sg * ux * ARROW_HEAD, y + sg * uy * ARROW_HEAD
        heads.append([(x, y), (bx - uy * hw, by + ux * hw), (bx + uy * hw, by - ux * hw)])
        bases.append((bx, by))
    return bases, heads


def arrow_label_at(a, pt):
    """Centre of an arrow's label: the middle of the arrow, its box just above the line."""
    return (a["x1"] + a["x2"]) / 2, a["y1"] - ARROW_W / 2 - ARROW_GAP - (0.8 * pt + 2 * DIST_PAD) / 2


def arrows(lay):
    out = []
    for a in lay.get("arrows", []):
        (b1, b2), heads = arrow_shape(a)
        out.append(f'<line x1="{b1[0]:.1f}" y1="{b1[1]:.1f}" x2="{b2[0]:.1f}" y2="{b2[1]:.1f}" '
                   f'stroke="{DIST_FILL}" stroke-width="{ARROW_W}"/>')
        out += [f'<polygon points="{" ".join(f"{x:.1f},{y:.1f}" for x, y in h)}" fill="{DIST_FILL}"/>'
                for h in heads]
        if a.get("text"):
            out.append(dist_label(*arrow_label_at(a, lay.get("dist_pt", DIST_PT)), a["text"], lay.get("dist_pt", DIST_PT)))
    return out


def labels(lay):
    pt = lay.get("label_pt", LABEL_PT)
    # centred on (x, y): baseline 0.35 em below the centre (Arial cap height ~0.7 em), turned about (x, y)
    return [f'<text x="{t["x"]:.1f}" y="{t["y"] + 0.35 * pt:.1f}" text-anchor="middle" '
            f'font-family="{FONT}" font-size="{pt}" font-weight="bold" '
            f'transform="rotate({t.get("angle", 0):.1f} {t["x"]:.1f} {t["y"]:.1f})">{t["text"]}</text>'
            for t in lay.get("labels", [])]


def turned(points, x, y, angle):
    t = np.radians(angle)
    c, sn = np.cos(t), np.sin(t)
    return [(x + (px - x) * c - (py - y) * sn, y + (px - x) * sn + (py - y) * c) for px, py in points]


def ink(png):
    """The drawn part of a render (alpha > 0.05) as fractions (left, top, right, bottom) of the image."""
    a = plt.imread(png)[..., 3]
    rows, cols = (a > 0.05).any(1).nonzero()[0], (a > 0.05).any(0).nonzero()[0]
    h, w = a.shape
    return cols[0] / w, rows[0] / h, (cols[-1] + 1) / w, (rows[-1] + 1) / h


def extent(lay):
    """Bounding box of everything drawn: the template, the panels' drawn parts (after size and
    rotation), the number circles and the motif names (after rotation)."""
    pts = [(0, 0), (W, H)]
    for p, v in lay["panels"].items():
        png = HERE / PANEL_DIR / f"{p}.png"
        if v.get("hidden") or not png.exists():
            continue
        x, y, sz = v["x"], v["y"], v["size"]
        l, t, r, b = ink(png)
        x0, y0 = x - sz / 2, y - sz / 2
        pts += turned([(x0 + l * sz, y0 + t * sz), (x0 + r * sz, y0 + t * sz), (x0 + l * sz, y0 + b * sz),
                       (x0 + r * sz, y0 + b * sz)], x, y, v.get("angle", 0))
    r = lay.get("badge_r", BADGE_R)
    for bd in lay["badges"]:
        pts += [(bd["x"] - r - 1, bd["y"] - r - 1), (bd["x"] + r + 1, bd["y"] + r + 1)]
    dpt, dist = lay.get("dist_pt", DIST_PT), distances()
    for d in lay.get("dists", []):
        if d["key"] in dist:
            hw = width(dist[d["key"]], dpt, bold=True) / 2
            pts += [(d["x"] - hw, d["y"] - dpt / 2), (d["x"] + hw, d["y"] + dpt / 2)]
    for a in lay.get("arrows", []):
        pts += [(a["x1"], a["y1"]), (a["x2"], a["y2"])]
        if a.get("text"):
            hw, (lx, ly) = width(a["text"], dpt, bold=True) / 2, arrow_label_at(a, dpt)
            pts += [(lx - hw, ly - dpt / 2), (lx + hw, ly + dpt / 2)]
    pt = lay.get("label_pt", LABEL_PT)
    for lb in lay.get("labels", []):
        hw = width(lb["text"], pt, bold=True) / 2
        pts += turned([(lb["x"] - hw, lb["y"] - pt / 2), (lb["x"] + hw, lb["y"] - pt / 2),
                       (lb["x"] - hw, lb["y"] + pt / 2), (lb["x"] + hw, lb["y"] + pt / 2)],
                      lb["x"], lb["y"], lb.get("angle", 0))
    xs, ys = [p[0] for p in pts], [p[1] for p in pts]
    m = 4   # margin
    return min(xs) - m, min(ys) - m, max(xs) - min(xs) + 2 * m, max(ys) - min(ys) + 2 * m


def main():
    lay = layout()
    body = "".join(pills() + panels(lay) + arrows(lay) + badges(lay) + dists(lay) + labels(lay) + header())
    x, y, w, h = extent(lay)
    OUT.write_text(f'<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" '
                   f'width="{w:.0f}" height="{h:.0f}" viewBox="{x:.1f} {y:.1f} {w:.1f} {h:.1f}">'
                   f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" fill="white"/>{body}</svg>',
                   encoding="utf-8")
    print(f"wrote {OUT}; viewBox {x:.0f} {y:.0f} {w:.0f} {h:.0f} (template 0 0 {W} {H})")


if __name__ == "__main__":
    main()
