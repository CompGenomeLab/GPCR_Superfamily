"""PyMOL session in the Figure 4 panel style: every carrier at the figure's positions, on the anchor's TM domain.

    python save_fig4_session.py   [--tag gp157]

Loads motif_superposition_tmcore_<tag>.pse, applies render_fig4_panels.style() (Fig. 3b look, Fig. 4
lineage colours) and shows all carrier residues (local_backbone_check_<tag>.tsv, role carrier) at
the positions of the figure panels (PANELS) as sticks in the lineage colours, with their disulfide
bonds; no distance cutoff, so the divergent carriers the figure leaves out are there too (user, Oct 15).
The anchor's TM domain (as kept in the session, trimmed to the 7TM boundary) is the grey cartoon.
Selections per panel: fig4_<panel> (all carriers) and fig4_<panel>_offfig (carriers the figure
panel does not show, scene()).
Writes motif_superposition_fig4style_<tag>.pse; the source session is not changed.
"""
from pymol import cmd

import render_fig4_panels as R

B = R.B


def residues(rows):
    return " or ".join(f"({r['structure']} and resi {r['residue'][1:]})" for r in rows)


def main():
    rows, up, centre = R.setup()
    panels = {}
    for panel, positions in R.PANELS.items():
        sc = R.scene(panel, positions, rows, up, centre)
        shown = {(r["structure"], r["residue"]) for r in sc["rows"]} if sc else set()
        allr = [r for r in rows if r["position"] in positions]
        panels[panel] = (allr, [r for r in allr if (r["structure"], r["residue"]) not in shown])
    cmd.hide("everything")
    for obj in cmd.get_names("objects"):
        cmd.unset("cartoon_color", obj)
    anchor = B.ANCHORS[R.TAG]
    cmd.show("cartoon", anchor)
    cmd.set("cartoon_color", "fig3b_grey", anchor)
    for panel, (allr, off) in panels.items():
        name = f"fig4_{panel.replace('+', '_')}"
        atoms = "(sidechain or name CA) and not hydro"
        cmd.select(name, f"({residues(allr)}) and {atoms}")
        if off:
            cmd.select(f"{name}_offfig", f"({residues(off)}) and {atoms}")
        for r in allr:
            s = f"{r['structure']} and resi {r['residue'][1:]} and {atoms}"
            cmd.show("sticks", s)
            cmd.color(f"fig4_{r['lineage']}", s)
        for a, b in R.disulfides(allr):
            cmd.bond(a, b)
            cmd.show("sticks", f"{a} or {b}")
        print(f"{panel}: {len(allr)} carriers, {len(off)} not in the figure", flush=True)
    cmd.delete("site")
    cmd.deselect()
    cmd.orient(anchor)
    out = B.HERE / f"motif_superposition_fig4style_{R.TAG}.pse"
    cmd.save(str(out))
    print("wrote", out)


if __name__ == "__main__":
    main()
