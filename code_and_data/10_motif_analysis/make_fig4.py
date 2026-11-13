"""Figure 3 (shared sequence motifs) as one SVG, replacing the Illustrator assembly.

a. Every superfamily column: number of the 18 counted families conserved at 60% (y) against the
   Shannon entropy of the families' most frequent residues (x), the rule of step9
   superfamily_motif_positions.py, on the ML-tree (treein) G-INS-i --ep 0.123 mapping JSONs. The
   selected positions (>= MIN_FAMILIES families) sit in a rounded dashed box that points to b.
b. The structure render of those positions (render_fig3b_structure.py output).
c. The logo exported from the GPCRevolution superfamily-logo page (Download SVG; ML-tree guide tree,
   60% / >= 12 families, cAMP and Class C per kingdom), inserted unchanged. Motif names (MOTIFS) are
   added over its columns, one plain line per motif.

Usage: python make_fig3.py   (writes fig3_motifs.svg; the a+b part is also kept as fig3ab.svg)
"""
import math
import re
import sys
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

HERE = Path(__file__).resolve().parent

import superfamily_motif_positions as smp  # noqa: E402


FIG = smp.stage("10_motif_analysis", "figure4")
STRUCTURE = FIG / "fig3b_motif_positions_2rh1.png"
LOGO = FIG / "fig3c_logo_website.svg"          # the website export
TOP = FIG / "fig3ab.svg"
OUT = FIG / "fig4_motifs.svg"
MIN_FAMILIES = 12

# motif names over their columns: (name, scheme subscript, family whose GPCRdb numbers define it, numbers)
MOTIFS = [
    ("S-S", "ABCF", "classA", ["3x25", "45x50"]),
    ("Na+", "A", "classA", ["2x50", "7x45"]),       # 7x46 dropped (user, Oct 14): no direct ion contact, variable
    ("PIF", "A", "classA", ["3x40", "6x44"]),
    ("DRY", "A", "classA", ["3x49", "3x50", "3x51"]),
    ("Chol.", "ABF", "classA", ["4x50"]),           # no cholesterol site in Class C (user, Oct 24)
    ("CWxP", "A", "classA", ["6x48"]),
    ("NPxxY", "A", "classA", ["7x49", "7x50", "7x53"]),
    ("HETx", "B", "classB1", ["3x50"]),            # H2.50 E3.50 T6.42 (Wootten 2013); Y7.57 is VxxxY, not HETx
    ("NPGQ", "B", "classB1", ["7x49"]),            # N5.50 P6.47 G6.50 Q7.49; G7.50 is not part of it
    ("Ionic lock", "C", "classC", ["3x50"]),
    ("E. switch", "F", "classF", ["3x43", "6x36", "6x40"]),   # extended switch W3.43-Y6.40 (Bous 2024); W7.55 is the R6.32-W7.55 switch
    ("YPERPI", "F", "classA", ["2x40", "2x42"]),   # placed by the user (Oct 14), at these Class A columns
    ("YHW", "F", "classF", ["4x50"]),              # at the cholesterol site (user, Oct 14)
]

W = 941                  # width of the logo export, used for the whole figure
TOP_H = 456              # height of the a + b block
GAP = 6                  # clearance between the a + b block and the tallest motif label
FONT = "Arial"
# one type scale for the whole figure: MAIN for names and titles (the logo's family names are 12),
# SMALL for ticks, numbers and TM labels; motif names are MAIN, SUBSUP for their sub/superscripts
LABEL_PT, SMALL, SUBSUP, LETTER = 12, 10, 9, 20
AXIS_LW = 1.6


def scatter_points(fams):
    pts = []
    for c in range(len(fams["classA"])):
        present = [fams[f][c] for f in smp.COUNTED if fams[f][c]]
        if not present:
            continue
        n = sum(p["cons"] >= smp.THRESHOLD for p in present)
        modes = Counter(p["mode"] for p in present)
        k = sum(modes.values())
        pts.append((c, n, -sum(v / k * math.log2(v / k) for v in modes.values())))
    return pts


def draw_top(pts, ncol, right):
    """a (scatter, box, arrow) and b (structure), W x TOP_H points, as an SVG. Returns b's left edge."""
    plt.rcParams.update({"font.family": "Arial", "svg.fonttype": "none", "font.size": LABEL_PT,
                         "axes.linewidth": AXIS_LW, "xtick.major.width": AXIS_LW, "ytick.major.width": AXIS_LW,
                         "xtick.major.size": 5, "ytick.major.size": 5,
                         "xtick.labelsize": SMALL, "ytick.labelsize": SMALL})
    fig = plt.figure(figsize=(W / 72, TOP_H / 72))
    frac = lambda x, y, w, h: [x / W, y / TOP_H, w / W, h / TOP_H]   # axes from a point rectangle
    axa = fig.add_axes(frac(75, 65 + (TOP_H - 380) / 2, 339, 274))
    sel = {c for c, n, _ in pts if n >= MIN_FAMILIES}
    lo = [(e, n) for c, n, e in pts if c not in sel]
    hi = [(e, n) for c, n, e in pts if c in sel]
    axa.scatter([p[0] for p in lo], [p[1] for p in lo], s=14, color="#c4c8d4", lw=0)
    axa.scatter([p[0] for p in hi], [p[1] for p in hi], s=20, color="#3d4a6b", lw=0)
    axa.set_xlabel("Shannon entropy\n(sequence divergence across families)")
    axa.set_ylabel("Number of families with\n60% conservation")
    axa.set_xlim(0, 4.5)
    axa.set_ylim(0, 17)
    for s in ("top", "right"):
        axa.spines[s].set_visible(False)
    ex, ey = [p[0] for p in hi], [p[1] for p in hi]
    bx0, bx1, by0, by1 = min(ex) - 0.18, max(ex) + 0.18, MIN_FAMILIES - 0.6, max(ey) + 0.7
    axa.add_patch(FancyBboxPatch((bx0, by0), bx1 - bx0, by1 - by0, boxstyle="round,pad=0,rounding_size=0.25",
                                 mutation_aspect=4, fill=False, ls=(0, (4, 3)), lw=AXIS_LW, ec="#222"))
    axa.text((bx0 + bx1) / 2, by1 + 0.4, f"Top {ncol} Positions", ha="center", va="bottom", fontsize=LABEL_PT)
    # b, 20% larger than the first layout (429 x 438 points), right-aligned
    # b: the render cropped to the protein (plus PAD_PX), at the scale that fitted the whole render
    # into 429 x 438 points (20% above the first layout), its right edge on the logo's (right)
    full = plt.imread(STRUCTURE)
    rgb = full[..., :3] if full.ndim == 3 else full
    ink = (rgb.max(axis=-1) if rgb.ndim == 3 else rgb) < 0.97
    ic, ir = ink.any(axis=0).nonzero()[0], ink.any(axis=1).nonzero()[0]
    PAD_PX = 6
    r0, r1 = max(0, ir[0] - PAD_PX), min(full.shape[0], ir[-1] + PAD_PX + 1)
    c0, c1 = max(0, ic[0] - PAD_PX), min(full.shape[1], ic[-1] + PAD_PX + 1)
    img = full[r0:r1, c0:c1]
    scale = min(429 / full.shape[1], 438 / full.shape[0])
    bw, bh = img.shape[1] * scale, img.shape[0] * scale
    axb = fig.add_axes(frac(right - bw, (TOP_H - bh) / 2, bw, bh))   # right edge on the logo's
    axb.imshow(img)
    axb.axis("off")
    # the arrow runs from the box to the left edge of the drawn protein; both ends in figure
    # fractions, which survive the change of resolution between drawing and saving the SVG
    fig.canvas.draw()
    cols, rows = ic - c0, ir - r0             # the protein's extent in the cropped image
    to_fig = fig.transFigure.inverted()
    start = to_fig.transform(axa.transData.transform((bx1, (by0 + by1) / 2)))
    end = to_fig.transform(axb.transData.transform((cols[0] - 25, 0)))
    end[1] = start[1]                       # horizontal, at the height of the box
    fig.add_artist(FancyArrowPatch(tuple(start), tuple(end), transform=fig.transFigure,
                                   arrowstyle="-|>", mutation_scale=30, lw=2 * AXIS_LW, color="#222",
                                   shrinkA=12, shrinkB=12))   # body 2x the axis lines, head 3x size; gaps at both ends
    fig.text((start[0] + end[0]) / 2, start[1] + 5 / TOP_H, "Structural visualization",
             ha="center", va="bottom", fontsize=16)
    b_left = to_fig.transform(axb.get_window_extent().p0)[0] * W
    # lowest drawn point of a and b, in points from the top of the block: the scatter's tight box
    # and the lowest pixel of the protein (the render has a white margin)
    r = fig.canvas.get_renderer()
    low_a = to_fig.transform(axa.get_tightbbox(r).p0)[1]
    low_b = to_fig.transform(axb.transData.transform((0, rows[-1])))[1]
    content_bottom = (1 - min(low_a, low_b)) * TOP_H
    fig.savefig(TOP)
    plt.close(fig)
    return b_left, content_bottom


def label_length(name, sub):
    """Length in points of a motif label as drawn (Arial, LABEL_PT, subscript SUBSUP)."""
    from matplotlib.font_manager import FontProperties
    from matplotlib.textpath import TextPath
    w = lambda t, size: TextPath((0, 0), t, size=size, prop=FontProperties(family="Arial")).get_extents().width
    return w("Na" if name == "Na+" else name, LABEL_PT) + w(sub, SUBSUP) + (w("+", SUBSUP) if name == "Na+" else 0) + 2


def group_reference_rows(logo):
    """Merge the export's GPCRdb numbering rows that are identical at every column into one row
    named after all their families, restack the rows with the export's alternating greys, and
    move the protein-region bar and the panel height up by the rows saved."""
    num = re.compile(r'<text[^>]*transform="translate\(([\d.]+), ([\d.]+)\) rotate\(-90\)"[^>]*>(?:\d+x\d+|\d+)</text>')
    name = re.compile(r'<text text-anchor="end" x="162" y="(\d+)"[^>]*>([^<]*)</text>')
    bg = re.compile(r'<rect x="0" y="(\d+)" width="[\d.]+" height="30" fill="#000000" fill-opacity="(0\.0\d)"/>')
    rows = {}
    for m in num.finditer(logo):
        rows.setdefault(float(m.group(2)), []).append(m.group(0))
    if not rows:
        return logo, 0
    ys = sorted(rows)
    names = {float(m.group(1)) - 4: m.group(2) for m in name.finditer(logo) if float(m.group(1)) - 4 in rows}
    key = lambda y: tuple(sorted((m.group(1), re.sub(r"<[^>]+>", "", t)) for t in rows[y] for m in [num.match(t)]))
    groups = []
    for y in ys:
        for g in groups:
            if key(g[0]) == key(y):
                g.append(y)
                break
        else:
            groups.append([y])
    step, top = ys[1] - ys[0], ys[0]
    fills = [m.group(2) for m in bg.finditer(logo)]
    rect_w = re.search(r'<rect x="0" y="\d+" width="([\d.]+)" height="30" fill="#000000"', logo).group(1)
    # drop the old rows (numbers, names, backgrounds) and write the grouped ones in their place
    for y in ys:
        for t in rows[y]:
            logo = logo.replace(t, "", 1)
    logo = name.sub(lambda m: "" if float(m.group(1)) - 4 in rows else m.group(0), logo)
    logo = bg.sub("", logo)
    new_bg, new_text, new_names = [], [], []
    for k, g in enumerate(groups):
        y = top + k * step
        new_bg.append(f'<rect x="0" y="{y - step / 2:.0f}" width="{rect_w}" height="{step:.0f}" fill="#000000" '
                      f'fill-opacity="{fills[k % 2]}"/>')
        new_text += [t.replace(f", {g[0]:g}) rotate(-90)", f", {y:g}) rotate(-90)") for t in rows[g[0]]]
        new_names.append(f'<text text-anchor="end" x="162" y="{y + 4:g}" style="font-size: 12px; font-family: Arial;">'
                         f'{", ".join(names.get(v, "") for v in g)}</text>')
    saved = (len(ys) - len(groups)) * step
    # the names go in the label panel (the first nested svg), backgrounds and numbers in the plot panel
    label_end = logo.index("</svg>")
    logo = logo[:label_end] + "".join(new_names) + logo[label_end:]
    pos = logo.index('<rect x="0" y="', logo.index("</svg>"))
    logo = logo[:pos] + "".join(new_bg) + "".join(new_text) + logo[pos:]
    # everything below the old rows (the protein-region bar) moves up by the rows saved
    bottom = ys[-1] + step / 2
    def lift(m):
        v = float(m.group(2))
        return f'{m.group(1)}"{v - saved:g}"' if v >= bottom else m.group(0)
    logo = re.sub(r'(\by=)"([\d.]+)"', lift, logo)
    logo = re.sub(r'height="(\d+)"', lambda m: f'height="{int(m.group(1)) - int(saved)}"'
                  if int(m.group(1)) > bottom else m.group(0), logo)
    # the root's viewBox must shrink with its height, or the viewer scales the whole logo down to fit
    logo = re.sub(r'viewBox="0 0 ([\d.]+) ([\d.]+)"',
                  lambda m: f'viewBox="0 0 {m.group(1)} {float(m.group(2)) - saved:g}"', logo, count=1)
    return logo, saved


def label_ecl2(logo):
    """The region bar's unlabeled box between TM4 and TM5 is ECL2: drop the box, write ECL2 at its
    centre in the TM labels' style, re-alternate the remaining boxes' greys, and remove the bar's
    background. All the bar's labels are then centred on their boxes."""
    box = re.compile(r'<rect x="([\d.]+)" y="([\d.]+)" width="([\d.]+)" height="25" fill="(#[0-9A-F]{6})" '
                     r'stroke="#999999"[^>]*/>(<text[^>]*>TM\d</text>)?')
    found = list(box.finditer(logo))
    greys = [found[0].group(4), found[1].group(4)]
    out, k, last = [], 0, 0
    for m in found:
        out.append(logo[last:m.start()])
        last = m.end()
        x, y, w = (float(m.group(i)) for i in (1, 2, 3))
        if m.group(5) is None:     # the unlabeled box
            out.append(f'<text x="{x + w / 2:g}" y="{y + 12.5:g}" text-anchor="middle" dominant-baseline="middle" '
                       f'fill="#000000" font-size="11px" font-family="Helvetica" font-weight="bold">ECL2</text>')
            continue
        g = greys[k % 2]
        rgb = ", ".join(str(int(g[i:i + 2], 16)) for i in (1, 3, 5))
        out.append(re.sub(r'fill: rgb\([^)]*\)', f'fill: rgb({rgb})', m.group(0).replace(m.group(4), g, 1), count=1))
        k += 1
    # drop the bar's shared light background, so the gaps between the boxes are white
    logo = re.sub(r'<rect x="0" y="[\d.]+" width="[\d.]+" height="25" fill="#F7F7F7" fill-opacity="1"/>', "",
                  "".join(out) + logo[last:])
    # and centre the labels on their boxes with an explicit baseline (half the ink height of the
    # bold label at its final size, SMALL, below the box centre), not dominant-baseline="middle",
    # which SVG viewers treat differently
    from matplotlib.font_manager import FontProperties
    from matplotlib.textpath import TextPath
    ink = TextPath((0, 0), "TM1ECL2", size=SMALL, prop=FontProperties(family="Arial", weight="bold")).get_extents()
    drop = (ink.y0 + ink.y1) / 2
    return re.sub(r'<text x="([\d.]+)" y="([\d.]+)" text-anchor="middle" dominant-baseline="middle"([^>]*>(?:TM\d|ECL2)</text>)',
                  lambda m: f'<text x="{m.group(1)}" y="{float(m.group(2)) + drop:.2f}" text-anchor="middle"{m.group(3)}',
                  logo)


def inner(svg_text):
    """An SVG document without its XML declaration, so it can be nested."""
    return re.sub(r"<\?xml[^>]*\?>", "", re.sub(r"<!DOCTYPE[^>]*>", "", svg_text)).strip()


SUB_DROP, SUP_RISE = 3, 4      # points the subscript drops and the Na superscript rises


def motif_label(name, sub, x, y):
    """Vertical label reading upward from (x, y). Its ink (name plus the dropped subscript, and the
    raised + of Na) is centred on x: the baseline is offset by the middle of that ink's extent, which
    the rotation turns into a sideways shift. Explicit offsets, not dominant-baseline or
    baseline-shift="sub", which SVG viewers treat differently."""
    from matplotlib.font_manager import FontProperties
    from matplotlib.textpath import TextPath
    ext = lambda t, size: TextPath((0, 0), t, size=size, prop=FontProperties(family="Arial")).get_extents()
    main_t = "Na" if name == "Na+" else name
    spans = [(ext(main_t, LABEL_PT).y0, ext(main_t, LABEL_PT).y1),
             (ext(sub, SUBSUP).y0 - SUB_DROP, ext(sub, SUBSUP).y1 - SUB_DROP)]
    if name == "Na+":
        spans.append((ext("+", SUBSUP).y0 + SUP_RISE, ext("+", SUBSUP).y1 + SUP_RISE))
    mid = (min(a for a, _ in spans) + max(b for _, b in spans)) / 2        # above the baseline
    head = (f'Na<tspan dy="{-SUP_RISE}" font-size="{SUBSUP}">+</tspan><tspan dy="{SUP_RISE}"></tspan>'
            if name == "Na+" else name)
    return (f'<text transform="translate({x:.2f},{y:.2f}) rotate(-90)" font-family="{FONT}" font-size="{LABEL_PT}" '
            f'text-anchor="start" dy="{mid:.2f}">{head}<tspan dy="{SUB_DROP}" font-size="{SUBSUP}">{sub}</tspan></text>')


def main():
    fams = smp.load_all()
    pts = scatter_points(fams)
    chosen = sorted(c for c, n, _ in pts if n >= MIN_FAMILIES)

    logo, saved = group_reference_rows(LOGO.read_text(encoding="utf-8"))
    logo = label_ecl2(logo)
    lw, lh = (float(v) for v in re.search(r'<svg width="([\d.]+)" height="([\d.]+)"', logo).groups())
    # the export nests a label panel and the plot panel; columns are in the second, offset by its x
    plot_x = float(re.findall(r'<svg width="[\d.]+" height="[\d.]+" x="([\d.]+)" y="0"', logo)[1])
    # the logo's right edge: its widest row background (the reference rows and the region bar)
    right = plot_x + max(float(w) for w in re.findall(r'<rect x="0" y="[\d.]+" width="([\d.]+)"', logo))
    b_left, content_bottom = draw_top(pts, len(chosen), right)
    # column centres: the first reference row's labels, translate(x, y) rotate(-90), one per column
    labels = re.findall(r'transform="translate\(([\d.]+), ([\d.]+)\) rotate\(-90\)"[^>]*>(\d+x\d+|\d+)</text>', logo)
    first_y = labels[0][1]
    xs = [plot_x + float(x) for x, y, _ in labels if y == first_y]
    if len(xs) != len(chosen):
        raise SystemExit(f"the logo has {len(xs)} columns but the selection has {len(chosen)}; re-export the logo")
    logo_top = float(re.search(r'<g transform="translate\([\d.]+, ([\d.]+)\)" class="axis"', logo).group(1))

    # motif lines, one per motif; several at a column stand side by side
    tags = {}
    for name, sub, fam, nums in MOTIFS:
        for j, c in enumerate(chosen):
            if str((fams[fam][c] or {}).get("label")) in nums:
                tags.setdefault(j, []).append((name, sub))
    # the logo sits as high as the tallest motif label allows: that label's top just clears the
    # lowest point of a and b (GAP points)
    tallest = max((label_length(n, s) for lst in tags.values() for n, s in lst), default=0)
    y0 = content_bottom + GAP + tallest + 33 - logo_top   # 33: label start below the logo's top row
    # a single motif gets a straight line; several at one column share a trunk that forks sideways,
    # so their labels stand apart (SPREAD points between neighbouring labels)
    base, fork, tip, SPREAD = y0 + logo_top - 4, y0 + logo_top - 14, y0 + logo_top - 30, 14
    marks = []
    for j, lst in sorted(tags.items()):
        k = len(lst)
        if k > 1:
            marks.append(f'<line x1="{xs[j]:.2f}" y1="{base:.2f}" x2="{xs[j]:.2f}" y2="{fork:.2f}" stroke="#555" stroke-width="1"/>')
        for i, (name, sub) in enumerate(lst):
            x = xs[j] + (i - (k - 1) / 2) * SPREAD
            path = (f"M{xs[j]:.2f},{base:.2f} V{tip:.2f}" if k == 1 else
                    f"M{xs[j]:.2f},{fork:.2f} L{x:.2f},{fork - 8:.2f} V{tip:.2f}")
            marks.append(f'<path d="{path}" fill="none" stroke="#555" stroke-width="1"/>')
            marks.append(motif_label(name, sub, x, tip - 3))

    total_h = y0 + lh
    letters = "".join(f'<text x="{x}" y="{y}" font-family="{FONT}" font-size="{LETTER}">{t}.</text>'
                      for t, x, y in (("a", 4, 22), ("b", round(b_left) - 20, 22), ("c", 4, content_bottom + GAP + 22)))
    top = inner(TOP.read_text(encoding="utf-8"))
    top = re.sub(r'<svg ', f'<svg x="0" y="0" ', top, count=1)
    top = re.sub(r'width="[\d.]+pt" height="[\d.]+pt"', f'width="{W}" height="{TOP_H}"', top, count=1)
    body = inner(logo)
    # Arial everywhere: the website draws in Helvetica and its axes in sans-serif
    body = (body.replace("font-family: Helvetica", "font-family: Arial")
            .replace('font-family="Helvetica"', 'font-family="Arial"')
            .replace('font-family="sans-serif"', 'font-family="Arial"')
            # and onto the figure's type scale: lineage labels (14) to MAIN, TM labels (11) to SMALL
            .replace("font-size: 14px", f"font-size: {LABEL_PT}px")
            .replace("font-size: 11px", f"font-size: {SMALL}px")
            .replace('font-size="11px"', f'font-size="{SMALL}px"'))
    body = re.sub(r'<svg width="[\d.]+" height="[\d.]+"', f'<svg x="0" y="{y0}" width="{lw:.0f}" height="{lh:.0f}"', body, count=1)
    OUT.write_text(f'<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" '
                   f'width="{W}" height="{total_h:.0f}" viewBox="0 0 {W} {total_h:.0f}">'
                   f'<rect width="100%" height="100%" fill="white"/>{top}{body}{"".join(marks)}{letters}</svg>',
                   encoding="utf-8")
    print(f"{len(chosen)} positions, {len(xs)} logo columns; wrote {OUT}")
    for j, lst in sorted(tags.items()):
        print(" ", (fams["classA"][chosen[j]] or {}).get("label"), lst)


if __name__ == "__main__":
    main()
