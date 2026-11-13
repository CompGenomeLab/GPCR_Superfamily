"""Per-family calls at the final superfamily motif positions.

Positions are columns of the ML-tree-guided superfamily alignment (motif_columns.SUP_FASTA): the 28
positions conserved in >= 12 of 18 families (60%, similarity groups; Oct 11) plus the positions the
motif Results text adds (12x51, the sodium pocket 3x39, W3x46B = 3x42A, W45x51B, E6x35C, the two
TM5 prolines P5x50A and P5x42B = 5x461A, and the aromatic 5x58).

Where the text names a residue class, a family carries the motif when that class makes up >= 60%
of its sequences at the column (conserved), or when its most common residue belongs to the class
but falls short of 60% (consensus). At the other positions (accepted = None) families conserve
different residues, and each family is judged on its own most common residue's similarity group.
Fractions are over all family sequences, gaps included, as on the website.

Writes final_motif_list_per_family.tsv (families and the cAMP and Class C kingdom rows).
"""
import json

from motif_columns import HERE, MAPS

AR, HY = "FWYH", "VILM"
GROUPS = ["ED", "WYHF", "RK", "QN", "VI", "ML"]  # as in the position selection
FAMILY = "family consensus"
# label, superfamily column (1-based), class name, accepted residues (None: the family's own
# similarity group), role in the text
POSITIONS = [
    ("1x46", 68, "S", "S", "cAMP-specific (S1x50B)"),
    ("1x50", 77, "N", "N", "Rho-specific"),
    ("1x57", 91, "aromatic", AR, "Glu+cAMP (F1x61C)"),
    ("12x51", 181, "R/K", "RK", "Glu+cAMP (R/K12x51C)"),
    ("2x40", 214, FAMILY, None, "top 28"),
    ("2x42", 221, "hydrophobic", HY, "Glu+cAMP (L2x38C)"),
    ("2x46", 245, "L", "L", "all three"),
    ("2x50", 256, "D/E", "DE", "cAMP+Rho sodium pocket"),
    ("3x25", 441, "C", "C", "all three (disulfide)"),
    ("3x39", 468, "D/E/S", "DES", "sodium pocket (plant cAMP S, fungal cAMP D, Olfactory E)"),
    ("3x40", 469, FAMILY, None, "top 28"),
    ("3x42", 471, "aromatic", AR, "cAMP-specific (W3x46B)"),
    ("3x43", 472, FAMILY, None, "top 28"),
    ("3x46", 475, FAMILY, None, "top 28"),
    ("3x49", 495, "D/E", "DE", "cAMP+Rho (DRY)"),
    ("3x50", 496, "R/K", "RK", "Rho-specific (DRY)"),
    ("3x51", 498, "aromatic", AR, "cAMP+Rho (DRY)"),
    ("3x54", 507, "hydrophobic", HY + "F", "all three; Glu-specific F3x58C"),
    ("4x50", 655, "aromatic (W)", AR, "cAMP+Rho (cholesterol)"),
    ("4x53", 665, "P", "P", "cAMP + plant Class C (P4x43C)"),
    ("45x50", 875, "C", "C", "all three (disulfide)"),
    ("45x51", 876, "aromatic", AR, "cAMP-specific (W45x51B)"),
    # TM5 prolines: one helical turn apart in the main alignment, the same column in several
    # alternative alignments; Rhodopsin-related P5x50A and cAMP-related P5x42B (5x461A)
    ("P5x461", 1010, "P", "P", "cAMP TM5 proline (P5x42B)"),
    ("P5x50", 1024, "P", "P", "Rho TM5 proline"),
    ("5x58", 1051, "aromatic", AR, "cAMP+Rho"),
    ("5x61", 1058, FAMILY, None, "top 28"),
    ("5x65", 1063, "hydrophobic", HY + "F", "all three (G protein interface)"),
    ("E6x35C", 1144, "E", "E", "Glu-specific"),
    ("6x40", 1878, FAMILY, None, "top 28"),
    ("6x44", 1883, "aromatic", AR, "cAMP+Rho (PIF)"),
    ("6x48", 1897, "aromatic", AR, "all three (toggle switch)"),
    ("7x45", 2100, FAMILY, None, "top 28"),
    ("7x46", 2102, FAMILY, None, "top 28"),
    ("7x49", 2107, "N/D", "ND", "cAMP+Rho (NPxxY)"),
    ("7x50", 2108, FAMILY, None, "top 28"),
    ("7x53", 2111, "aromatic", AR, "all three (NPxxY)"),
]
LINEAGES = {
    "Glutamate": ["classC", "Vomeronasal2"],
    "cAMP": ["cAMP", "classB1", "classB2", "Mth", "classF", "FSLB", "GP143", "GP157", "GPR1",
             "STE3", "TM116"],
    "Rhodopsin": ["classA", "Olfactory", "Vomeronasal1", "classT", "Nematode"],
    "background": ["STE2"],
}
# kingdom rows the text discusses (plant cAMP sodium pocket, plant Class C P4x43C)
KINGDOM_ROWS = {"cAMP_Plants": "cAMP", "classC_Plants": "classC"}


def family_accept(accepted, position):
    """Accepted residues at a column: the given class, or (None) the similarity group of the
    family's most common residue there."""
    if accepted is not None:
        return accepted
    counts = position.get("residueCounts") or {}
    if not counts:
        return ""
    top = max(counts, key=counts.get)
    return next((g for g in GROUPS if top in g), top)


def call(position, accepted):
    counts = position.get("residueCounts") or {}
    total = position.get("totalSequences") or 0
    if not counts or total <= 0:
        return "gap", "-", 0.0
    accepted = family_accept(accepted, position)
    frac = 100 * sum(counts.get(a, 0) for a in accepted) / total
    top = max(counts, key=counts.get)
    gaps = total - sum(counts.values())
    if frac >= 60:
        return "conserved", top, frac
    if top in accepted and counts[top] > gaps:
        return "consensus", top, frac
    return "absent", top, frac


def main():
    rows = []
    groups = dict(LINEAGES, kingdom=list(KINGDOM_ROWS))
    for lineage, fams in groups.items():
        for fam in fams:
            pos = json.load(open(MAPS / f"{fam}.json"))["positions"]
            for label, col, cls, acc, pattern in POSITIONS:
                status, top, frac = call(pos[col - 1], acc)
                rows.append((lineage, fam, label, col, cls, pattern, status, top, frac))
    with open(HERE / "final_motif_list_per_family.tsv", "w") as f:
        f.write("lineage\tfamily\tposition\tsuperfamily_column\tmotif_class\tlineage_pattern\t"
                "status\tmost_common_residue\tclass_percent\n")
        for r in rows:
            f.write("\t".join(map(str, r[:8])) + f"\t{r[8]:.1f}\n")
    for lineage, fams in groups.items():
        print(f"\n## {lineage}")
        for fam in fams:
            mine = [r for r in rows if r[1] == fam]
            has = [f"{r[2]} {r[7]}{int(r[8])}{'' if r[6] == 'conserved' else '(c)'}"
                   for r in mine if r[6] in ("conserved", "consensus")]
            lacks = [f"{r[2]} {r[7]}" for r in mine if r[6] in ("absent", "gap")]
            print(f"{fam}: HAS {', '.join(has)}\n    LACKS {', '.join(lacks)}")


if __name__ == "__main__":
    main()
