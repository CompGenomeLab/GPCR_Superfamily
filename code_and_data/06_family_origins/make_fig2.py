"""Figure 2: lineage distribution and family origins, as a standalone SVG and a PNG for the Word document.

Input: 06_family_origins/lineage_origin_figure.html of the archive (build_lineage_origin_figure.py).
Steps: render the page with headless Chrome and take its SVG, carry over the page styles the SVG relies on,
then rasterize to PNG at 3x.

    python make_fig2.py
"""
import os
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import stage  # noqa: E402

SRC = Path(stage("06_family_origins", "lineage_origin_figure.html"))
OUT = SRC.with_name("Figure_2_lineage_origins.svg")
CHROME = os.environ.get("CHROME", "chrome")  # path of a Chrome/Chromium executable

# the page styles the figure through CSS classes, which the standalone SVG does not carry
STYLE = ("<style>text { fill: #1d1d1f; } .muted { fill: #6e6e73; } "
         ".clade line { stroke: #8e8e93; stroke-width: 1; } .clade .lab { font-size: 10.5px; fill: #6e6e73; } "
         ".clade .mark { fill: #1f5fa8; } .rowbg { fill: none; }</style>")

dom = subprocess.run([CHROME, "--headless=new", "--disable-gpu", "--allow-file-access-from-files", "--dump-dom",
                      SRC.as_uri()], check=True, capture_output=True, text=True, encoding="utf-8").stdout
svg = re.search(r"<svg.*</svg>", dom, re.S).group(0)
w, h = map(float, re.search(r'<svg[^>]*width="([\d.]+)" height="([\d.]+)"', svg).groups())
if "viewBox" not in svg[:500]:
    svg = svg.replace("<svg ", f'<svg viewBox="0 0 {w:g} {h:g}" ', 1)
svg = svg.replace("<svg ", '<svg font-family="Segoe UI, Arial, sans-serif" ', 1)
svg = re.sub(r"(<svg[^>]*>)", r"\1" + STYLE, svg, count=1)
OUT.write_text('<?xml version="1.0" encoding="UTF-8"?>\n' + svg, encoding="utf-8")

png = OUT.with_suffix(".png")
subprocess.run([CHROME, "--headless=new", "--disable-gpu", "--hide-scrollbars", "--default-background-color=ffffffff",
                "--force-device-scale-factor=3", f"--window-size={int(w) + 1},{int(h)}", f"--screenshot={png}",
                OUT.as_uri()], check=True, capture_output=True)
print(f"wrote {OUT.name} ({w:g} x {h:g}) and {png.name}")
