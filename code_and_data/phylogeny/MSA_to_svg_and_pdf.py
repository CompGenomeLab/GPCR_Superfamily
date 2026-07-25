from __future__ import annotations
from typing import Dict, List, Tuple, Optional
import os

# ===== Your color groups (gaps as holes) =====
GROUPS: Dict[str, Dict[str, object]] = {
    "aromatic":    {"residues": ['W', 'Y', 'H', 'F'],          "color": "#FCB315"},
    "polar":       {"residues": ['S', 'T', 'Q', 'N'],          "color": "#7D2985"},
    "small":       {"residues": ['P', 'G', 'A'],               "color": "#231F20"},
    "acidic":      {"residues": ['E', 'D'],                    "color": "#DD6030"},
    "basic":       {"residues": ['R', 'K'],                    "color": "#7CAEC4"},
    "hydrophobic": {"residues": ['V', 'C', 'I', 'M', 'L'],     "color": "#B4B4B4"},
}
GAP_CHAR = "-"
GAP_AS_HOLE = None        # None = don't draw gaps (true holes, fewer SVG elements)
DEFAULT_HEX = "#808080"   # fallback color if an unexpected char appears


# ---------- FASTA ----------
def parse_fasta_aligned_from_path(path: str) -> List[Tuple[str, str]]:
    """Read aligned FASTA from a file path. Returns [(header, sequence), ...]."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"FASTA not found: {path}")
    with open(path, "r", encoding="utf-8") as f:
        lines = [ln.rstrip("\r\n") for ln in f]
    records: List[Tuple[str, str]] = []
    header: Optional[str] = None
    buf: List[str] = []
    for ln in lines:
        if not ln.strip():
            continue
        if ln.startswith(">"):
            if header is not None:
                records.append((header, "".join(buf)))
            header = ln[1:].strip()
            buf = []
        else:
            buf.append(ln.strip())
    if header is not None:
        records.append((header, "".join(buf)))
    if not records:
        raise ValueError("No FASTA records found.")
    lengths = {len(seq) for _, seq in records}
    if len(lengths) != 1:
        raise ValueError(f"Not aligned: sequences have differing lengths {sorted(lengths)}")
    return records


# ---------- Colors ----------
def build_color_map_from_groups(
    groups: Dict[str, Dict[str, object]],
    gap_color: Optional[str] = GAP_AS_HOLE
) -> Dict[str, Optional[str]]:
    """
    Build {char: '#RRGGBB' or None}. If None -> skip drawing (gap).
    Mirrors case automatically unless explicitly overridden.
    """
    cmap: Dict[str, Optional[str]] = {}
    for obj in groups.values():
        residues = obj["residues"]
        color = str(obj["color"])
        for r in residues:
            cmap[r] = color
            if r.isupper() and r.lower() not in cmap:
                cmap[r.lower()] = color
            if r.islower() and r.upper() not in cmap:
                cmap[r.upper()] = color
    cmap.setdefault(GAP_CHAR, gap_color)  # gaps -> None (hole) or a color (e.g., "#FFFFFF")
    return cmap


# ---------- Layout ----------
def compute_cell_dims(
    n_rows: int, n_cols: int,
    width: Optional[int], height: Optional[int],
    x_pad: int, y_pad: int
) -> Tuple[int, int, int, int]:
    """Compute cell_w, cell_h, total_w, total_h. Squares if only one dimension given."""
    def fit_1d(n: int, target: Optional[int], pad: int) -> Optional[int]:
        if target is None: return None
        return max(1, (target - (n - 1) * pad) // n)

    cw = fit_1d(n_cols, width, x_pad)
    ch = fit_1d(n_rows, height, y_pad)

    if width and height:
        cell_w = cw or 1
        cell_h = ch or 1
    elif width:
        cell_w = cw or 1
        cell_h = cell_w
    elif height:
        cell_h = ch or 1
        cell_w = cell_h
    else:
        cell_w = cell_h = 20

    total_w = n_cols * cell_w + (n_cols - 1) * x_pad
    total_h = n_rows * cell_h + (n_rows - 1) * y_pad
    return cell_w, cell_h, total_w, total_h


# ---------- Run-length merge (horizontal + vertical) ----------
def _row_runs(seq: str, cmap: Dict[str, Optional[str]]) -> List[Tuple[int, int, str]]:
    """
    Compress one row into runs of same non-gap color.
    Returns list of (c_start, c_end_inclusive, hex_color).
    """
    runs: List[Tuple[int, int, str]] = []
    c = 0
    n = len(seq)
    while c < n:
        ch = seq[c]
        color = cmap.get(ch, cmap.get(ch.upper(), cmap.get(ch.lower(), DEFAULT_HEX)))
        if color is None:  # gap -> skip
            c += 1
            continue
        c_start = c
        c += 1
        while c < n:
            ch2 = seq[c]
            color2 = cmap.get(ch2, cmap.get(ch2.upper(), cmap.get(ch2.lower(), DEFAULT_HEX)))
            if color2 != color or color2 is None:
                break
            c += 1
        runs.append((c_start, c - 1, color))
    return runs


def _merge_rows_vertically(all_row_runs: List[List[Tuple[int, int, str]]]) -> List[Tuple[int, int, int, int, str]]:
    """
    Given per-row horizontal runs, merge vertically when a run with the same (c_start, c_end, color)
    appears on consecutive rows. Returns rectangles as (r0, r1_inclusive, c0, c1_inclusive, color).
    """
    merged: List[Tuple[int, int, int, int, str]] = []
    prev: dict[Tuple[int, int, str], Tuple[int, int, int, int, str]] = {}  # key->(r0,r1,c0,c1,color)

    for r, runs in enumerate(all_row_runs):
        curr: dict[Tuple[int, int, str], Tuple[int, int, int, int, str]] = {}
        for (c0, c1, color) in runs:
            key = (c0, c1, color)
            if key in prev:
                r0, r1, pc0, pc1, col = prev[key]
                curr[key] = (r0, r, pc0, pc1, col)  # extend downward
            else:
                curr[key] = (r, r, c0, c1, color)   # start new rect
        # finalize anything in prev that didn't continue
        for key, rect in prev.items():
            if key not in curr:
                merged.append(rect)
        prev = curr
    # flush remaining
    merged.extend(prev.values())
    return merged


# ---------- SVG (rectangles only, merged) ----------
def generate_alignment_svg_merged_from_fasta(
    fasta_path: str,
    out_svg: str,
    width: Optional[int] = 800,
    height: Optional[int] = None,
    x_pad: int = 0,
    y_pad: int = 0,
    corner_radius: int = 0,
    color_map: Optional[Dict[str, Optional[str]]] = None,
) -> str:
    """
    Render an aligned FASTA to SVG with **merged rectangles**:
      - gaps not drawn (holes) if color_map maps '-' to None
      - same-color contiguous cells merged horizontally and vertically
    """
    records = parse_fasta_aligned_from_path(fasta_path)
    n_rows = len(records)
    n_cols = len(records[0][1])

    cmap = color_map or build_color_map_from_groups(GROUPS, gap_color=GAP_AS_HOLE)
    cell_w, cell_h, total_w, total_h = compute_cell_dims(n_rows, n_cols, width, height, x_pad, y_pad)

    # Collect horizontal runs for each row
    all_runs: List[List[Tuple[int, int, str]]] = []
    for _, seq in records:
        all_runs.append(_row_runs(seq, cmap))

    # Merge vertically into maximal rectangles
    rects = _merge_rows_vertically(all_runs)

    # If only one dimension provided, set the other to natural size
    if width is None:  width = total_w
    if height is None: height = total_h

    # Transparent background; gaps remain holes
    lines: List[str] = []
    lines.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {total_w} {total_h}">')
    for (r0, r1, c0, c1, color) in rects:
        x = c0 * (cell_w + x_pad)
        y = r0 * (cell_h + y_pad)
        w = (c1 - c0 + 1) * cell_w + (c1 - c0) * x_pad
        h = (r1 - r0 + 1) * cell_h + (r1 - r0) * y_pad
        lines.append(
            f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{corner_radius}" ry="{corner_radius}" fill="{color}"/>'
        )
    lines.append("</svg>")
    svg = "\n".join(lines)
    with open(out_svg, "w", encoding="utf-8") as f:
        f.write(svg)
    return out_svg


# ---------- Optional: SVG Ã¢â€ â€™ PDF ----------
def svg_to_pdf(svg_path: str, pdf_path: str) -> str:
    """
    Convert SVG to PDF (vector). Requires:
      pip install svglib reportlab
    """
    from reportlab.graphics import renderPDF
    from svglib.svglib import svg2rlg
    drawing = svg2rlg(svg_path)
    renderPDF.drawToFile(drawing, pdf_path)
    return pdf_path


# ======================
# === YOUR SETTINGS ====
# ======================

PHYLOGENY_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.environ.get("GPCR_PHYLOGENY_DIR", PHYLOGENY_DIR)
OUTPUT_DIR = os.environ.get("GPCR_MSA_FIGURE_DIR", os.path.join(PHYLOGENY_DIR, "data", "visualization"))
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 1) Path to your aligned FASTA file
fasta_path = os.path.join(
    PROJECT_DIR,
    "data",
    "alignments",
    "All_woC-Ste3_clust0.7_cov1_Oct24_minsize2_linsi_ep0.123_trimends_treein_ginsi_ep0.123_kpicgappy.fasta",
)

# 2) Where you want the SVG saved
out_svg = os.path.join(OUTPUT_DIR, os.path.basename(fasta_path) + ".svg")

# 3) Build the color map
#    IMPORTANT: Using gap_color=None to leave gaps as holes (reduces element count).
#    If you want visible white cells instead, change to gap_color="#FFFFFF".
color_map = build_color_map_from_groups(GROUPS, gap_color=None)

# 4) Generate merged-rectangles SVG
generate_alignment_svg_merged_from_fasta(
    fasta_path=fasta_path,
    out_svg=out_svg,
    width=800,       # set width (px); or set height instead
    height=None,
    x_pad=0,         # keep pads 0 for maximal merging
    y_pad=0,
    corner_radius=0,
    color_map=color_map
)
print("SVG written to:", out_svg)

# 5) Convert SVG Ã¢â€ â€™ PDF
pdf_path = os.path.join(OUTPUT_DIR, os.path.basename(fasta_path) + ".pdf")
try:
    svg_to_pdf(out_svg, pdf_path)
    print("Saved:", out_svg, "and", pdf_path)
except Exception as e:
    print("SVG was saved, but PDF conversion requires 'svglib' and 'reportlab'.")
    print("Install with: pip install svglib reportlab")
    print("Error:", e)
