"""Figure 4b: the final motif positions mapped on beta2AR (PDB 2RH1).

Starts from the session of the earlier figure (figure4/superfamily_logo_positions_overall_look.pse in the
archive: view, cartoon, colours and sphere style unchanged) and only
replaces the highlighted residues (selection logo_residues, all atoms as spheres) with the final ones.

Positions are the final motif list of stage 11 (motif_position_spread_<tag>.tsv), with the beta2AR
residue at each from the stage 11 mapping of 6PS2 (same UniProt numbering). Left out: P6x50 (STE2 test
only) and E6x35C (Class C only, no class A column). 5x65 is not mapped in 6PS2; it is beta2AR A226,
as in the old session. 45x50 is C191, the cysteine bonded to C106 (3x25).

Writes figure4/fig3b_motif_positions_2rh1.pse and .png in the archive (inputs of make_fig4.py).
Usage: python render_fig4b_structure.py [tag]  (default gp157)
"""
import csv
import sys
from pathlib import Path

from pymol import cmd, stored

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import stage  # noqa: E402

HERE = stage("10_motif_analysis", "figure4")
STEP12 = stage("11_structure_superposition")
OLD = HERE / "superfamily_logo_positions_overall_look.pse"
SKIP = {"P6x50", "E6x35C"}
EXTRA = {"5x65": 226}


def residues(tag):
    pos = [r["position"] for r in csv.DictReader(open(STEP12 / f"motif_position_spread_{tag}.tsv"), delimiter="\t")]
    got = {r["position"]: int(r["residue"][1:]) for r in
           csv.DictReader(open(STEP12 / f"motif_residues_per_structure_{tag}.tsv"), delimiter="\t")
           if r["structure"] == "A_ADRB2_6PS2" and r["residue"] not in ("", "-")}
    got.update(EXTRA)
    return {p: got[p] for p in pos if p not in SKIP and p in got}


def main():
    tag = sys.argv[1] if len(sys.argv) > 1 else "gp157"
    res = residues(tag)
    cmd.load(str(OLD))
    stored.c = []
    cmd.iterate("logo_residues", "stored.c.append(color)")
    colour = max(set(stored.c), key=stored.c.count)  # the old spheres' colour
    old = sorted({int(a.resi) for a in cmd.get_model("logo_residues and name CA").atom})

    cmd.hide("spheres", "2rh1")
    new = sorted(res.values())
    cmd.select("logo_residues", "2rh1 and resi " + "+".join(map(str, new)))
    cmd.show("spheres", "logo_residues")
    cmd.color(colour, "logo_residues")
    for p, n in sorted(res.items(), key=lambda x: x[1]):
        print(p, cmd.get_model(f"2rh1 and resi {n} and name CA").atom[0].resn + str(n))
    print("removed", sorted(set(old) - set(new)), "added", sorted(set(new) - set(old)), flush=True)

    cmd.save(str(HERE / "fig3b_motif_positions_2rh1.pse"))
    cmd.ray(2400, 2400)
    cmd.png(str(HERE / "fig3b_motif_positions_2rh1.png"), dpi=300)
    print("wrote", HERE / "fig3b_motif_positions_2rh1.png")


if __name__ == "__main__":
    main()
