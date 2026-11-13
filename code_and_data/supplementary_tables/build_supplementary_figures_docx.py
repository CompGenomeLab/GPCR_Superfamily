"""Supplementary_Figures.docx: every supplementary figure on its own page with its legend below.

Figures: figures/Figure_S<n>_*.png (rasterized from the exported SVGs). Legends: the bold-titled paragraphs
of figures/supplementary_figure_legends.md ("**Figure S<n>: title.** text"). A figure whose PNG is missing
gets a placeholder, so the document can be rebuilt as figures arrive.

    python build_supplementary_figures_docx.py
"""
import re
import os
import sys
from pathlib import Path

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK
from docx.oxml.ns import qn
from docx.shared import Inches, Pt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import stage  # noqa: E402
HERE = Path(__file__).resolve().parent
FIG = stage("supplementary_figures")
FIGURES = ["S1", "S2", "S3", "S4"]


def legends():
    text = (FIG / "supplementary_figure_legends.md").read_text(encoding="utf-8")
    out = {}
    for m in re.finditer(r"\*\*Figure (S\d+): (.+?)\*\*\s*(.+?)(?=\n\s*\n)", text, re.S):
        out[m.group(1)] = (m.group(2).strip(), " ".join(m.group(3).split()))
    return out


def main():
    leg = legends()
    doc = Document()
    sec = doc.sections[0]
    sec.page_width, sec.page_height = Inches(8.5), Inches(11)
    for side in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
        setattr(sec, side, Inches(1))
    for name in ("Normal", "Heading 1"):
        st = doc.styles[name]
        st.font.name = "Arial"
        rfonts = st.element.get_or_add_rPr().get_or_add_rFonts()
        for attr in ("w:ascii", "w:hAnsi", "w:eastAsia", "w:cs"):
            rfonts.set(qn(attr), "Arial")
        for attr in ("w:asciiTheme", "w:hAnsiTheme", "w:eastAsiaTheme", "w:cstheme"):
            rfonts.attrib.pop(qn(attr), None)
    doc.styles["Normal"].font.size = Pt(11)
    doc.add_heading("Supplementary Figures", level=1)
    for i, fid in enumerate(FIGURES):
        if i:
            doc.add_paragraph().add_run().add_break(WD_BREAK.PAGE)
        png = next(iter(sorted(FIG.glob(f"Figure_{fid}_*.png"))), None)
        if png:
            doc.add_picture(str(png), width=Inches(6.5))
        else:
            doc.add_paragraph(f"[Figure {fid} to be added]")
        title, body = leg.get(fid, ("", "[legend to be added]"))
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        p.add_run(f"Figure {fid}: {title}" + (" " if title else "")).bold = True
        p.add_run(body)
        print(fid, "figure:", png.name if png else "missing", "| legend:", "yes" if fid in leg else "missing")
    out = stage("supplementary_figures", "Supplementary_Figures.docx")
    doc.save(out)
    print("wrote", out)


if __name__ == "__main__":
    main()
