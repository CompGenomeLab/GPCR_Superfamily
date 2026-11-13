"""Figure S2: TM5 region of the superfamily logo built with the default (MAFFT) guide tree.

Input: the superfamily-logo export (guide tree "MAFFT default"), figures/superfamily_logo_default_guide_tree_full.svg.
Steps: crop to Class A 5x41-5x70 (crop_logo_svg.py), then merge the GPCRdb reference rows that are identical
over that window into one row named after all their families (make_fig3.group_reference_rows, as in Fig. 4c),
then rasterize to PNG at 3x with headless Chrome for the Word document.

    python make_figS2.py
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
SRC = FIG / "superfamily_logo_default_guide_tree_full.svg"
CROP = FIG / "Figure_S2_TM5_prolines_default_guide_tree_crop.svg"
OUT = FIG / "Figure_S2_TM5_prolines_default_guide_tree.svg"
CHROME = os.environ.get("CHROME", "chrome")  # path of a Chrome/Chromium executable
sys.path.insert(0, str(HERE.parent / "10_motif_analysis"))
import make_fig4 as make_fig3  # noqa: E402

subprocess.run([sys.executable, str(HERE / "crop_logo_svg.py"), str(SRC), str(CROP), "5x41", "5x70"], check=True)
svg = CROP.read_text(encoding="utf-8")
before = re.search(r'viewBox="([\d.]+) 0 ([\d.]+) ([\d.]+)"', svg[200:])  # the cropped plot panel
svg, saved = make_fig3.group_reference_rows(svg)
# the cropped plot panel has its own viewBox (starting at the window's x), which must shrink like the root's
svg = svg.replace(before.group(0),
                  f'viewBox="{before.group(1)} 0 {before.group(2)} {float(before.group(3)) - saved:g}"', 1)
OUT.write_text(svg, encoding="utf-8")
CROP.unlink()
w, h = map(float, re.search(r'<svg width="([\d.]+)" height="([\d.]+)"', svg).groups())
png = FIG / "Figure_S2_TM5_prolines_default_guide_tree.png"
subprocess.run([CHROME, "--headless=new", "--disable-gpu", "--hide-scrollbars", "--default-background-color=ffffffff",
                "--force-device-scale-factor=3", f"--window-size={int(w) + 1},{int(h)}", f"--screenshot={png}",
                OUT.as_uri()], check=True, capture_output=True)
names = re.findall(r'<text text-anchor="end" x="162" y="[\d.]+" style="font-size: 12px; font-family: Arial;">([^<]*)</text>', svg)
print(f"saved {saved:g} px of reference rows; grouped rows: {names}")
print("wrote", OUT.name, "and", png.name)
