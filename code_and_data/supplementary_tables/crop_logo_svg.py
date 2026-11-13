"""Crop a superfamily-logo SVG export (GPCRevolution) to a range of Class A GPCRdb positions.

Keeps the family-name panel on the left and the plot columns between the two Class A positions (inclusive,
with half a column of margin); drops logo letters and labels outside the window, so the file stays small.

    python crop_logo_svg.py <export.svg> <out.svg> 5x41 5x70
"""
import re
import sys

src, out, first, last = sys.argv[1:5]
s = open(src, encoding="utf-8").read()

# Class A reference row: the first row of rotated GPCRdb labels under the logo
rows = {}
for x, y, t in re.findall(r'transform="translate\(([\d.]+), ([\d.]+)\) rotate\(-90\)"[^>]*>(\d+x\d+)</text>', s):
    rows.setdefault(float(y), {})[t] = float(x)
class_a = rows[min(rows)]
x0, x1 = class_a[first], class_a[last]
step = min(b - a for a, b in zip(sorted(class_a.values()), sorted(class_a.values())[1:]))
lo, hi = x0 - step / 2 - 2, x1 + step / 2 + 2

plot_open = re.search(r'<svg width="([\d.]+)" height="([\d.]+)" x="([\d.]+)" y="0">', s[200:]).group(0)
start = s.index(plot_open)
label_panel, plot = s[:start], s[start + len(plot_open):]
height = re.search(r'height="([\d.]+)"', plot_open).group(1)
px = float(re.search(r' x="([\d.]+)"', plot_open).group(1))
width = hi - lo


def inside(x, w=0.0):
    return x + w >= lo and x <= hi


# drop letter glyph boxes (<svg x=.. width=..>...</svg>, no nesting) outside the window
plot = re.sub(r'<svg x="([\d.]+)"[^>]*?width="([\d.]+)"[^>]*>.*?</svg>',
              lambda m: m.group(0) if inside(float(m.group(1)), float(m.group(2))) else "", plot, flags=re.S)
# drop rotated position labels outside the window
plot = re.sub(r'<text[^>]*transform="translate\(([\d.]+), [\d.]+\) rotate\(-90\)"[^>]*>[^<]*</text>',
              lambda m: m.group(0) if inside(float(m.group(1))) else "", plot)

total = px + width
head = re.match(r'<\?xml[^>]*>\s*<svg[^>]*>', label_panel).group(0)
head_new = re.sub(r'width="[\d.]+"', f'width="{total:.0f}"', head, count=1)
head_new = re.sub(r'viewBox="[^"]*"', f'viewBox="0 0 {total:.0f} {height}"', head_new, count=1)
label_panel = head_new + label_panel[len(head):]
new_open = f'<svg width="{width:.2f}" height="{height}" x="{px:.0f}" y="0" viewBox="{lo:.2f} 0 {width:.2f} {height}">'
open(out, "w", encoding="utf-8").write(label_panel + new_open + plot)
print(f"window {first}-{last}: x {lo:.1f}-{hi:.1f} ({width:.0f} px of plot); "
      f"{len(s) / 1e6:.1f} MB -> {len(label_panel + plot) / 1e6:.1f} MB")
