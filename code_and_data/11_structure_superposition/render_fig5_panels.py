"""Figure 5 panels: one close-up per shared motif position from the saved superposition.

    python render_fig5_panels.py 2x46 4x50 ...   [--tag gp157]   (no positions: all of PANELS)

Each panel shows the carrier residues (local_backbone_check_<tag>.tsv, role carrier; STE2 never)
that lie within SHOW_MAX of the other carriers' centroid (residue_dev) as sticks in the Fig. 4 header's lineage colours, and one TM helix: the carrier's helix stretch
that fits the others best (smallest mean distance to the other carriers' stretches), drawn as the
Fig. 3b cartoon (opaque grey, PyMOL default helix style and lighting), with the stretch of the
neighbouring TM helix the side chains reach towards (neighbour()). The camera keeps the
extracellular side up and turns about the bundle axis to show the two helices side by side with the
residues facing the viewer (camera_axis()), unless an angle was saved for the panel in the 3Dmol
viewer (fig4_views.json, from export_fig4_viewer.py). Transparent background.
Writes fig4_panels/<panel>.png.
"""
import csv
import json
import sys

import numpy as np
from pymol import cmd

import build_superposition as B
from make_fig5 import COLOURS

TAG = sys.argv[sys.argv.index("--tag") + 1] if "--tag" in sys.argv else "gp157"
OUT = B.HERE / "fig4_panels"
VIEWS = B.HERE / "fig4_views.json"   # angles saved in the 3Dmol viewer (fig4_viewer.html)
# panel -> positions of local_backbone_check (class A numbering); the TM5 prolines share a panel
PANELS = {p: [p] for p in ["2x46", "3x25", "45x50", "7x53", "5x65", "6x48",
                           "2x50", "3x39", "3x49", "3x51", "3x54", "4x50", "5x58", "6x44", "7x49",
                           "1x57", "2x42", "4x53"]}   # 12x51 dropped: carriers too divergent to read (user, Oct 15)
PANELS["P5x50"] = ["P5x50", "P5x461"]
# motifs shown together (user, Oct 14): NPxxY, L2x38C with L2x46 and D2x50, the disulfide bridge, DRY
# and the cholesterol W4x50 with the TM4 prolines one turn above it (user, Oct 15)
for together in (["7x49", "7x53"], ["2x42", "2x46", "2x50"], ["3x25", "45x50"], ["3x49", "3x51", "3x54"],
                 ["4x50", "4x53"]):
    for p in together:
        del PANELS[p]
    PANELS["+".join(together)] = together
HELIX = 7          # residues either side of a carrier's residue that make up its helix stretch
SHOW_MAX = 4.0     # only residues within this distance (A) of the other carriers' centroid are shown
# shown whatever their distance (position, structure): plant Class C P4x43C, one residue off P4x53B
ALWAYS_SHOW = {("4x53", "C_plant_A0A8B8NZB9"), ("6x44", "T2_TAS2R14_9IJA")}   # + Taste-2 Y240 at F6x44A
# positions shown with only these residue types (3x54: I/V/L, the text's motif; F carriers hidden)
SHOW_ONLY = {"3x54": "IVL"}
HELIX_TRANSPARENCY = 0.5   # Fig. 3b's cartoon is opaque; lighter here so the residues stand out (user, Oct 15)
SS_MAX = 2.6      # A: SG-SG distance drawn as a disulfide bond
SHOW_NEIGHBOUR = False   # the neighbouring helix only sets the automatic angle; not drawn (user, Oct 14)
SHOWN = 1          # residues either side drawn as helix, for the main and the neighbouring helix
NO_CARTOON = {"45x50"}   # positions without a backbone piece (the ECL2 cysteine: a loop, user Oct 15)
EXTRA_SHOWN = {"3x25": 2}   # position -> extra helix residues either side (TM3 piece of the disulfide panel)
JOIN = 4           # pieces of the same helix at most this many residues apart are drawn as one
DISTANCE = 30.0
FIELD = 18.0       # width of the view (A)
SIZE = 1200        # pixels


def stretch(obj, n, k=None):
    k = HELIX if k is None else k
    return np.array([a.coord for a in cmd.get_model(f"{obj} and resi {n - k}-{n + k} and name CA").atom])


def best_helix(rows):
    """The structure whose helix stretches fit the other carriers' best: per position, the mean
    nearest-CA distance of its stretch to the others'; summed over the panel's positions, among the
    structures that carry all of them (else those of the first position). Returns its rows."""
    def gap(a, b):
        d = np.linalg.norm(a[:, None] - b[None], axis=-1)
        return (d.min(1).mean() + d.min(0).mean()) / 2
    by_pos = {}
    for r in rows:
        by_pos.setdefault(r["position"], []).append(r)
    total = {}
    for prs in by_pos.values():
        pts = [stretch(r["structure"], int(r["residue"][1:])) for r in prs]
        for i, r in enumerate(prs):
            total.setdefault(r["structure"], []).append(
                np.mean([gap(pts[i], b) for j, b in enumerate(pts) if j != i]) if len(pts) > 1 else 0)
    common = [o for o, v in total.items() if len(v) == len(by_pos)]
    pool = common or [r["structure"] for r in next(iter(by_pos.values()))]
    obj = min(pool, key=lambda o: sum(total[o]))
    return [r for r in rows if r["structure"] == obj]


def disulfides(rows):
    """Selections of bonded SG pairs (within SS_MAX) between shown cysteines of the same structure."""
    sg = {}
    for r in rows:
        if r["residue"].startswith("C"):
            sel = f"{r['structure']} and resi {r['residue'][1:]} and name SG"
            if cmd.count_atoms(sel):
                sg.setdefault(r["structure"], []).append(sel)
    pairs = []
    for sels in sg.values():
        for i, a in enumerate(sels):
            for b in sels[i + 1:]:
                if cmd.get_distance(a, b) <= SS_MAX:
                    pairs.append((a, b))
    return pairs


def neighbour(obj, n, tip):
    """The other TM helix of obj that the side chains reach towards: the TM segment (other than the
    one holding residue n) with a C-alpha closest to the side chains' centroid, shown as the
    SHOWN residues either side of that C-alpha. Returns ((first, last), centre of the shown part)."""
    best = None
    for a, b in B.runs(B.TM_RESI[obj]):
        if a - 8 <= n <= b + 8:    # the site's own helix, also when the site lies just past its end
            continue
        for at in cmd.get_model(f"{obj} and resi {a}-{b} and name CA").atom:
            d = np.linalg.norm(np.array(at.coord) - tip)
            if best is None or d < best[0]:
                best = (d, int(at.resi), a, b)
    if best is None:
        return None, tip
    _, m, a, b = best
    span = (max(a, m - SHOWN), min(b, m + SHOWN))
    return span, stretch(obj, m, SHOWN).mean(0)


def camera_axis(up, main, nb, atoms, sc_dir):
    """Viewing direction (towards the viewer), turned only about the bundle axis so the
    extracellular side stays up. Start perpendicular to the line from the main helix to the
    neighbour, so the two stand side by side with the residues between them, on the side the side
    chains point to; then, within +-30 degrees, take the angle with the widest spread of the side
    chain atoms (projected area), keeping the side chains turned towards the viewer."""
    flat = lambda v: (v - (v @ up) * up) / np.linalg.norm(v - (v @ up) * up)
    d = flat(nb - main)
    best = None
    for sign in (1, -1):
        z0 = sign * np.cross(up, d)
        for deg in range(-30, 31, 5):
            t = np.radians(deg)
            z = np.cos(t) * z0 + np.sin(t) * np.cross(up, z0)
            if sc_dir @ z < 0:
                continue
            x = np.cross(up, z)
            p = np.c_[atoms @ x, atoms @ up]
            area = np.ptp(p[:, 0]) * np.ptp(p[:, 1])
            if best is None or area > best[0]:
                best = (area, z)
    return best[1] if best else flat(sc_dir)


def style():
    """The Fig. 3b session's look (render_fig3b_structure.py starts from it)."""
    for k, v in {"cartoon_transparency": HELIX_TRANSPARENCY, "cartoon_oval_width": 0.25, "cartoon_oval_length": 0.7,
                 "cartoon_loop_radius": 0.2, "cartoon_fancy_helices": 0, "cartoon_cylindrical_helices": 0,
                 "ray_trace_mode": 0, "ambient": 0.14, "specular": 1, "shininess": 55, "direct": 0.45,
                 "reflect": 0.45, "light_count": 2, "depth_cue": 1, "antialias": 1, "ray_shadows": 1,
                 "stick_radius": 0.25, "orthoscopic": 1, "ray_opaque_background": 0}.items():
        cmd.set(k, v)
    cmd.bg_color("white")
    cmd.set_color("fig3b_grey", [0.707, 0.707, 0.707])
    for lin, hexcol in COLOURS.items():
        cmd.set_color(f"fig4_{lin}", [int(hexcol[i:i + 2], 16) / 255 for i in (1, 3, 5)])


def view_matrix(x, y, z, origin):
    r = np.array([x, y, z])
    return [*r.T.flatten(), 0.0, 0.0, -DISTANCE, *origin, DISTANCE - 15.0, DISTANCE + 15.0, 1.0]


def scene(panel, positions, all_rows, up, centre):
    """Show the panel's residues and helices in PyMOL and choose its automatic camera.
    Returns a dict (rows, helix, neighbour span, site, camera rows x/y/z) or None."""
    rows = [r for r in all_rows if r["position"] in positions]
    keep = lambda r: ((float(r["residue_dev"]) <= SHOW_MAX or (r["position"], r["structure"]) in ALWAYS_SHOW)
                      and r["residue"][0] in SHOW_ONLY.get(r["position"], r["residue"][0]))
    shown = [r for r in rows if keep(r)]
    # a cysteine past the cutoff stays when it is disulfide-bonded to a shown cysteine (whole bridges)
    def bonded(r):
        if not r["residue"].startswith("C"):
            return False
        sg = f"{r['structure']} and resi {r['residue'][1:]} and name SG"
        return any(o["structure"] == r["structure"] and o["residue"].startswith("C") and cmd.count_atoms(sg)
                   and cmd.get_distance(sg, f"{o['structure']} and resi {o['residue'][1:]} and name SG") <= SS_MAX
                   for o in shown)
    off = [r for r in rows if not keep(r) and not bonded(r)]
    rows = [r for r in rows if keep(r) or bonded(r)]
    for r in off:
        print(f"  not shown: {r['structure']} {r['residue']} ({float(r['residue_dev']):.1f} A)", flush=True)
    if not rows:
        print("no carriers for", panel)
        return None
    cmd.hide("everything")
    cmd.delete("site")
    for obj in cmd.get_names("objects"):
        cmd.unset("cartoon_color", obj)
    own = best_helix(rows)
    obj, n = own[0]["structure"], int(own[0]["residue"][1:])
    centres = [int(r["residue"][1:]) for r in own if r["position"] not in NO_CARTOON]
    extra = {int(r["residue"][1:]): EXTRA_SHOWN.get(r["position"], 0) for r in own}
    # a position the helix structure does not carry (the other lineage's TM5 prolines): extend its
    # helix to its residue closest to that position's carriers, so the helix passes them too
    for pos in {r["position"] for r in rows} - {r["position"] for r in own}:
        target = np.mean([cmd.get_model(f"{r['structure']} and resi {r['residue'][1:]} and name CA").atom[0].coord
                          for r in rows if r["position"] == pos], 0)
        near = min(cmd.get_model(f"{obj} and name CA").atom, key=lambda a: np.linalg.norm(np.array(a.coord) - target))
        centres.append(int(near.resi))
    spans = []
    for c in sorted(centres):    # pieces within JOIN residues joined (both NPxxY, DRY, TM5 prolines)
        k = SHOWN + extra.get(c, 0)
        if spans and c - k <= spans[-1][1] + JOIN:
            spans[-1] = (spans[-1][0], c + k)
        else:
            spans.append((c - k, c + k))
    for a, b in spans:
        cmd.show("cartoon", f"{obj} and resi {a}-{b}")
    cmd.set("cartoon_color", "fig3b_grey", obj)   # stays grey whatever the atoms' colour
    ca, sc_atoms, sc_vec = [], [], []
    for r in rows:
        o, m = r["structure"], int(r["residue"][1:])
        sel = f"{o} and resi {m} and (sidechain or name CA) and not hydro"
        cmd.show("sticks", sel)
        cmd.color(f"fig4_{r['lineage']}", sel)
        ca.append(np.array(cmd.get_model(f"{o} and resi {m} and name CA").atom[0].coord))
        side = [a.coord for a in cmd.get_model(f"{o} and resi {m} and sidechain and not hydro").atom]
        side = np.array(side) if side else ca[-1][None]
        sc_atoms.extend(side)
        sc_vec.append(side.mean(0) - ca[-1])
    for a, b in disulfides(rows):
        cmd.bond(a, b)
        cmd.show("sticks", f"{a} or {b}")
    site = np.mean(ca, 0)
    sc_atoms = np.array(sc_atoms)
    nb = neighbour(obj, n, sc_atoms.mean(0))
    z = camera_axis(up, stretch(obj, n).mean(0), nb[1], sc_atoms, np.mean(sc_vec, 0))
    if nb[0] and SHOW_NEIGHBOUR:
        cmd.show("cartoon", f"{obj} and resi {nb[0][0]}-{nb[0][1]}")
    print(panel, len(rows), "carriers; helix", obj, [r["residue"] for r in own], "; neighbour", nb[0], flush=True)
    return {"rows": rows, "obj": obj, "main": spans, "neighbour": nb[0], "site": site,
            "camera": np.array([np.cross(up, z), up, z])}


def saved_rotation(panel):
    """The rotation saved for this panel in the 3Dmol viewer (fig4_views.json; quaternion x, y, z, w
    of the viewer's model rotation), or None. Applied on top of the automatic camera."""
    if not VIEWS.exists():
        return None
    q = json.loads(VIEWS.read_text()).get(panel)
    if not q:
        return None
    x, y, z, w = q
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                     [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                     [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])


def render(panel, positions, all_rows, up, centre):
    sc = scene(panel, positions, all_rows, up, centre)
    if sc is None:
        return
    cam, site = sc["camera"], sc["site"]
    rot = saved_rotation(panel)
    if rot is not None:          # the viewer shows camera-frame coordinates; its rotation acts on them
        cam = rot @ cam
        print("  saved angle used", flush=True)
    cmd.pseudoatom("site", pos=[float(v) for v in site])   # camera target only, not shown
    cmd.hide("everything", "site")
    cmd.set_view(view_matrix(*cam, site))
    cmd.zoom("site", FIELD / 2, complete=0)
    cmd.set_view(view_matrix(*cam, site)[:9] + list(cmd.get_view())[9:])
    cmd.ray(SIZE, SIZE)
    cmd.png(str(OUT / f"{panel}.png"), dpi=300)
    print("wrote", OUT / f"{panel}.png", flush=True)


def setup():
    """Load the superposition and the carriers; returns (carrier rows, bundle up axis, bundle centre)."""
    cmd.load(str(B.HERE / f"motif_superposition_tmcore_{TAG}.pse"))
    rows = [r for r in csv.DictReader(open(B.HERE / f"local_backbone_check_{TAG}.tsv"), delimiter="	")
            if r["role"] == "carrier"]
    B.assign_tm(B.read_topologies({s[4] for s in B.STRUCTURES}))
    resi, _, core = B.ca_table(B.ANCHORS[TAG], tm_only=True)
    up = -B.bundle_axis(core, resi)      # extracellular up
    style()
    return rows, up, core.mean(0)


def main():
    wanted = [a for a in sys.argv[1:] if not a.startswith("--") and a != TAG] or list(PANELS)
    rows, up, centre = setup()
    OUT.mkdir(exist_ok=True)
    for p in wanted:
        render(p, PANELS[p], rows, up, centre)


if __name__ == "__main__":
    main()
