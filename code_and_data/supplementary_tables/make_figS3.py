"""Figure S3: shared-motif matrix (sequence conservation and structural overlap per family), from the Fig. 7 draft.

Input: ../11_structure_superposition/figures/figS3_motif_evidence.html, which draws the matrix from motif_evidence_data.js (Tables S8 and S10b; rerun
build_fig7_data.py there after the tables change).
Steps: render the page with headless Chrome and take its SVG, drop the in-figure title and the hover note (the legend
carries both), then rasterize to PNG at 3x for the Word document.

    python make_figS3.py
"""
import re
import subprocess
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import stage  # noqa: E402
HERE = Path(__file__).resolve().parent
FIG = stage("supplementary_figures")
SRC = HERE.parent / "11_structure_superposition" / "figures" / "figS3_motif_evidence.html"
OUT = FIG / "Figure_S3_shared_motif_evidence.svg"
CHROME = os.environ.get("CHROME", "chrome")  # path of a Chrome/Chromium executable

dom = subprocess.run([CHROME, "--headless=new", "--disable-gpu", "--allow-file-access-from-files", "--dump-dom",
                      SRC.as_uri()], check=True, capture_output=True, text=True, encoding="utf-8").stdout
svg = re.search(r"<svg.*</svg>", dom, re.S).group(0)
svg = re.sub(r'<text[^>]*>Shared motifs that define the GPCR superfamily</text>', "", svg, count=1)
svg = re.sub(r'<text[^>]*>Sequence: a residue of the motif[^<]*</text>', "", svg, count=1)
# trim the space the title left above the motif group bands (y = 36) and the hover note left below the legend
TOP = 28
w, h = map(float, re.search(r'<svg[^>]*width="([\d.]+)" height="([\d.]+)"', svg).groups())
h -= 20 + TOP
svg = re.sub(r'(<svg[^>]*?)height="[\d.]+"', rf'\1height="{h:g}"', svg, count=1)
svg = re.sub(r'viewBox="0 0 ([\d.]+) [\d.]+"', rf'viewBox="0 {TOP} \1 {h:g}"', svg, count=1)
# the page set Arial through CSS, which the standalone SVG does not carry
svg = svg.replace("<svg ", '<svg font-family="Arial, sans-serif" ', 1)
OUT.write_text('<?xml version="1.0" encoding="UTF-8"?>\n' + svg, encoding="utf-8")

png = OUT.with_suffix(".png")
subprocess.run([CHROME, "--headless=new", "--disable-gpu", "--hide-scrollbars", "--default-background-color=ffffffff",
                "--force-device-scale-factor=3", f"--window-size={int(w) + 1},{int(h)}", f"--screenshot={png}",
                OUT.as_uri()], check=True, capture_output=True)
print(f"wrote {OUT.name} ({w:g} x {h:g}) and {png.name}")
