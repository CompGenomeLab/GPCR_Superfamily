"""Superimpose the motif-figure structures on their 7TM core and mark the sequence motifs.

Structures are loaded whole (receptor chain only, nothing trimmed), but only their TM residues
enter the superposition. A structure's TM residues are the M/W positions of its protein's
DeepTMHMM topology (step2 oriented 3line), carried onto the structure by aligning the topology's
sequence to the chain's own sequence. Residue correspondence between structures is never taken
from sequence, and the motif residues play no part:

  1. The reference core is the anchor's TM C-alphas: beta2AR (6PS2) or GP157 (Q5UAW9),
     chosen on the command line (b2ar | gp157).
  2. Each structure's TM C-alphas are placed on that core by CE (sequence-independent).
  3. Register search: starting poses shifted by 0, +-1, +-2 helical turns along the bundle axis
     are each refined, and the pose with the best spatial overlap is kept.
  4. Joint refinement: each structure's TM C-alphas are matched to the mean core by an
     order-preserving spatial alignment, fitted (Kabsch) on pairs within a distance cutoff, the
     mean core is rebuilt from all structures, and this repeats. The anchor stays fixed.

Only then are the sequence-derived motif positions (superfamily alignment columns) mapped onto
each structure, by aligning the family alignment row to the structure's own sequence, and
measured residue by residue. Selections: <position> holds every residue that carries the motif
at the alignment column, whether or not it lines up in 3D; <position>_struct holds, for families
that carry the motif but whose column residue lies off-position, the nearby motif-class residue
that does line up. Residues that do not carry the motif are not selected.

Writes, per anchor tag, motif_superposition_tmcore_<tag>.pse, superposition_core_fit_<tag>.tsv,
motif_position_spread_<tag>.tsv, motif_residue_deviation_<tag>.tsv and
motif_residues_per_structure_<tag>.tsv.

Usage: python build_superposition.py [b2ar|gp157|gp157_hier]
gp157_hier: the GP157-anchored fit, then Olfactory and class T refitted onto class A (beta2AR).
"""
import csv
import json
import sys
from itertools import combinations

import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices
from pymol import cmd

from final_motif_list import FAMILY, LINEAGES, POSITIONS, call, family_accept
from motif_columns import (HERE, open_text, stage, MAPS, SUP_FASTA, accession, anchors, family_motif_columns, read_fasta,
                           trimmed_alignment)

STRUCT = HERE / "structures"
# object name, family, file, chain, family alignment row used for the column mapping, source
STRUCTURES = [
    ("A_ADRB2_6PS2", "classA", "6PS2.cif", "A", "P07550", "PDB"),
    ("Olf_consOR2_8UY0", "Olfactory", "8UY0.cif", "A", "A0A8C9GQD3", "PDB"),
    ("T2_TAS2R14_9IJA", "classT", "9IJA.cif", "R", "Q9NYV8", "PDB"),
    ("T2_TAS2R46_7XP4", "classT", "7XP4.cif", "R", "P59540", "PDB"),
    ("V1R_A0A1U8CRZ8", "Vomeronasal1", "AF-A0A1U8CRZ8-F1-model_v6.cif", "A", "A0A1U8CRZ8", "AF"),
    ("Nem_A8XV99", "Nematode", "AF-A8XV99-F1-model_v6.cif", "A", "A8XV99", "AF"),
    ("GPR1_A0A8H3U2U3", "GPR1", "AF-A0A8H3U2U3-F1-model_v6.cif", "A", "A0A8H3U2U3", "AF"),
    ("cAMP_A0A094FDE0", "cAMP", "AF-A0A094FDE0-F1-model_v6.cif", "A", "A0A094FDE0", "AF"),
    ("cAMP_plant_A0A0A0L5G2", "cAMP", "AF-A0A0A0L5G2-F1-model_v6.cif", "A", "A0A0A0L5G2", "AF"),
    ("F_FZD7_9EPO", "classF", "9EPO.cif", "D", "O75084", "PDB"),
    ("FSLB_A0A428UQK6", "FSLB", "AF-A0A428UQK6-F1-model_v6.cif", "A", "A0A428UQK6", "AF"),
    ("GP143_A0A091JTC0", "GP143", "AF-A0A091JTC0-F1-model_v6.cif", "A", "A0A091JTC0", "AF"),
    ("GP157_Q5UAW9", "GP157", "AF-Q5UAW9-F1-model_v6.cif", "A", "Q5UAW9", "AF"),
    ("TM116_A0A1S3WIB0", "TM116", "AF-A0A1S3WIB0-F1-model_v6.cif", "A", "A0A1S3WIB0", "AF"),
    ("Mth_V3ZTX9", "Mth", "AF-V3ZTX9-F1-model_v6.cif", "A", "V3ZTX9", "AF"),
    ("B2_ADGRL3_8JMT", "classB2", "8JMT.cif", "A", "Q9HAR2", "PDB"),
    ("B1_PTH1R_6FJ3", "classB1", "6FJ3.cif", "A", "Q03431", "PDB"),
    ("STE3_A0A4Y7PPB0", "STE3", "AF-A0A4Y7PPB0-F1-model_v6.cif", "A", "A0A4Y7PPB0", "AF"),
    ("C_CASR_7M3G", "classC", "7M3G.cif", "A", "P41180", "PDB"),
    ("C_plant_A0A8B8NZB9", "classC", "AF-A0A8B8NZB9-F1-model_v6.cif", "A", "A0A8B8NZB9", "AF"),
    ("V2R_A0A1L8FMT1", "Vomeronasal2", "AF-A0A1L8FMT1-F1-model_v6.cif", "A", "A0A1L8FMT1", "AF"),
    ("STE2_Ste2_7AD3", "STE2", "7AD3.cif", "A", "D6VTK4", "PDB"),
    ("STE2_A0A9W9H1E1", "STE2", "AF-A0A9W9H1E1-F1-model_v6.cif", "A", "A0A9W9H1E1", "AF"),
]
ANCHORS = {"b2ar": "A_ADRB2_6PS2", "gp157": "GP157_Q5UAW9", "gp157_hier": "GP157_Q5UAW9"}
# hierarchical variant: after the GP157-anchored fit, these are refitted onto their class A parent
HIERARCHY = {"gp157_hier": ("A_ADRB2_6PS2",
                            ["Olf_consOR2_8UY0", "T2_TAS2R14_9IJA", "T2_TAS2R46_7XP4"])}
ANCHOR = ANCHORS["b2ar"]  # set from the command line in main()
TAG = "b2ar"
TOPOLOGY = stage("03_topology_filter", "deeptmhmm_predictions_oriented.3line")
# exact motif residue numbers from the full family alignments (map_motif_residues.py, run on the
# cluster); without them the older, less reliable reference-based mapping is used
MAPPING = HERE / "motif_residue_numbers.tsv"
MEMBER_SEQS = HERE / "member_sequences.tsv"
CUTOFFS = [6.0, 5.0, 4.0] + [3.5] * 17
LINEAGE_OF = {f: lin for lin, fams in LINEAGES.items() for f in fams}
# one colour per lineage; STE2 (background, no lineage) stays grey
COLOURS = {"Glutamate": "tv_orange", "cAMP": "tv_blue", "Rhodopsin": "forest",
           "background": "grey60"}
# extra test position requested on top of the final list: the class A TM6 proline
EXTRA = [("P6x50", 1901, "P", "P", "test: STE2 TM6 proline")]
EXTRA_FAMILIES = {"P6x50": {"classA", "STE2"}}  # the P6x50 test compares Ste2 with class A only
# structures judged on a kingdom row of their family's logo (plant cAMP sodium pocket, plant
# Class C P4x43C): structure name -> logo key
SITE_KEY = {"cAMP_plant_A0A0A0L5G2": "cAMP_Plants", "C_plant_A0A8B8NZB9": "classC_Plants"}
# family-specific motif residues: class T carries the charge-reversed R at the sodium-pocket D2x50
FAMILY_ACCEPT = {("2x50", "classT"): "RK"}
# residues removed from a position on inspection of the superposition (position, structure):
EXCLUDE = set()  # was V2R F667 at 3x49 in the Sep27 set (model A0A674BXB5, since replaced)
# structural equivalent: in a family that carries the motif, a residue of the motif class within
# STRUCT_WINDOW residues of the column residue whose C-alpha lies within STRUCT_MAX of the
# cAMP/Rhodopsin carriers' C-alpha centroid, used when the column residue lies off-position;
# the family must also carry the motif class (conserved or consensus) at the neighbouring column
STRUCT_WINDOW, STRUCT_MAX = 1, 3.0
# canonical class A motif names, appended to the PyMOL selection names
MOTIF_NAMES = {"2x50": "Na", "3x39": "Na", "3x25": "disulfide", "45x50": "disulfide", "3x49": "DRY",
               "3x50": "DRY", "3x51": "DRY", "4x50": "cholesterol", "P5x50": "PIF",
               "6x44": "PIF", "6x48": "CWxP", "P6x50": "CWxP", "7x49": "NPxxY", "7x50": "NPxxY",
               "7x53": "NPxxY"}
# expected residue prefixed to the selection name; broad motif classes get a short tag
CLASS_TAGS = {FAMILY: "Fc", "hydrophobic": "Hy", "aromatic": "Ar", "hydrophobic/aromatic": "HyAr",
              "D/E or aromatic": "DEAr", "aromatic (W)": "W"}


def selection_name(label, cls):
    sel = label.replace("+", "")
    if not sel[0].isalpha():  # labels such as P5x50 already carry their residue
        sel = CLASS_TAGS.get(cls, cls.replace("/", "")) + sel
    return f"{sel}_{MOTIF_NAMES[label]}" if label in MOTIF_NAMES else sel


def disulfide_partner(obj, cys):
    """Residue number of the cysteine after `cys` (ECL2 side) whose SG lies within 2.5 A of residue
    `cys`'s SG, or None. Earlier cysteines (e.g. an N-terminal partner, GP143 C10) are not ECL2."""
    sg = cmd.get_model(f"{obj} and resi {cys} and resn CYS and name SG and not alt B").atom
    if not sg:
        return None
    ref = np.array(sg[0].coord)
    near = [(float(np.linalg.norm(np.array(a.coord) - ref)), int(a.resi))
            for a in cmd.get_model(f"{obj} and resn CYS and name SG and not alt B and not resi {cys}").atom]
    near = [x for x in near if x[0] < 2.5 and x[1] > cys]
    return min(near)[1] if near else None


def family_carries_at(positions, col, offset, accept):
    """Whether the family carries the motif class at the column `offset` residues from `col`
    along its reference sequence (conserved or consensus in the family logo)."""
    ref = positions[col - 1].get("refResidue")
    if ref is None:
        return False
    nb = next((p for p in positions if p.get("refResidue") == ref + offset), None)
    return nb is not None and call(nb, accept)[0] in ("conserved", "consensus")
BACKBONE = {"N", "CA", "C", "O", "OXT"}


ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
       "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
       "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
TM_RESI = {}  # object -> set of TM residue numbers
KEEP_RESI = {}  # object -> receptor residues within the 7TM boundary, kept in the saved session
MARGIN = 5  # residues kept beyond the first and last TM residue


def ca_table(obj, tm_only=False):
    m = cmd.get_model(f"{obj} and name CA and not alt B")
    atoms = [a for a in m.atom if a.resi.lstrip("-").isdigit()
             and (not tm_only or int(a.resi) in TM_RESI[obj])]
    return ([int(a.resi) for a in atoms], [a.resn for a in atoms],
            np.array([a.coord for a in atoms]))


def read_topologies(accs):
    """Accession -> (sequence, oriented topology) from the step2 3line file."""
    out, lines = {}, open_text(TOPOLOGY).read().splitlines()
    for k in range(0, len(lines) - 2, 3):
        acc = accession(lines[k][1:].split(" |")[0])
        if acc in accs:
            out[acc] = (lines[k + 1].strip(), lines[k + 2].strip())
    return out


def assign_tm(topologies):
    """Mark each structure's TM residues: DeepTMHMM M/W positions carried onto the chain."""
    al = PairwiseAligner(mode="global", open_gap_score=-10, extend_gap_score=-0.5,
                         end_gap_score=0)
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    for name, _, _, _, acc, _ in STRUCTURES:
        seq, topo = topologies[acc]
        resi, aa3, _ = ca_table(name)
        target = "".join(ONE.get(a, "X") for a in aa3)
        aln = al.align(target, seq)[0]
        tm, receptor = set(), set()
        for (t0, t1), (q0, q1) in zip(*aln.aligned):
            for k in range(t1 - t0):
                receptor.add(resi[t0 + k])
                if topo[q0 + k] in "MW":
                    tm.add(resi[t0 + k])
        TM_RESI[name] = tm
        # residues of the receptor itself (fusion partners and tags do not align to it) between
        # the first and last TM residue, with a small margin
        lo, hi = min(tm) - MARGIN, max(tm) + MARGIN
        KEEP_RESI[name] = {r for r in receptor if lo <= r <= hi}


def dp_align(core, xyz, d0=3.0, gap=-0.6):
    """Order-preserving alignment of core points to C-alphas -> list of (core i, ca j)."""
    near = np.where(np.min(np.linalg.norm(xyz[:, None] - core[None], axis=2), axis=1) < 8)[0]
    if not len(near):
        return []
    x = xyz[near]
    s = 1.0 / (1.0 + (np.linalg.norm(core[:, None] - x[None], axis=2) / d0) ** 2)
    m, n = s.shape
    val = np.zeros((m + 1, n + 1))
    path = np.zeros((m + 1, n + 1), dtype=np.int8)
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            d = val[i - 1, j - 1] + s[i - 1, j - 1]
            u = val[i - 1, j] + (gap if path[i - 1, j] == 1 else 0)
            l = val[i, j - 1] + (gap if path[i, j - 1] == 1 else 0)
            val[i, j], path[i, j] = max((d, 1), (u, 2), (l, 3))
    pairs, i, j = [], m, n
    while i > 0 and j > 0:
        if path[i, j] == 1:
            pairs.append((i - 1, int(near[j - 1])))
            i, j = i - 1, j - 1
        elif path[i, j] == 2:
            i -= 1
        else:
            j -= 1
    return pairs[::-1]


def kabsch(p, q):
    pc, qc = p.mean(0), q.mean(0)
    u, _, vt = np.linalg.svd((p - pc).T @ (q - qc))
    r = vt.T @ np.diag([1, 1, np.sign(np.linalg.det(vt.T @ u.T))]) @ u.T
    return r, qc - pc @ r.T


def apply(obj, r, t):
    m = [*r[0], t[0], *r[1], t[1], *r[2], t[2], 0, 0, 0, 1]
    cmd.transform_selection(obj, [float(x) for x in m], homogenous=1)


def load():
    for name, fam, fn, chain, _, _ in STRUCTURES:
        cmd.load(str(STRUCT / fn), name)
        cmd.remove(f"{name} and not chain {chain}")
        cmd.remove(f"{name} and (solvent or hydro)")
        cmd.alter(name, "segi=''")


def overlap_score(core, xyz, d0=3.0):
    """TM-score-like overlap of the core with a structure: sum of 1/(1+(d/d0)^2) over the
    order-preserving alignment, divided by the core length."""
    pairs = dp_align(core, xyz, d0)
    d = np.array([np.linalg.norm(core[i] - xyz[j]) for i, j in pairs]) if pairs else np.array([99])
    return float(np.sum(1 / (1 + (d / d0) ** 2)) / len(core))


def fit_pairwise(core, xyz):
    """Iterative spatial alignment + Kabsch of xyz onto core; returns R, t (for xyz)."""
    R, t = np.eye(3), np.zeros(3)
    cur = xyz.copy()
    for cut in CUTOFFS:
        pairs = [(i, j) for i, j in dp_align(core, cur) if np.linalg.norm(core[i] - cur[j]) < cut]
        if len(pairs) < 30:
            break
        ii, jj = map(list, zip(*pairs))
        r, tt = kabsch(cur[jj], core[ii])
        cur = cur @ r.T + tt
        R, t = r @ R, t @ r.T + tt
    return R, t


def runs(numbers):
    """Sorted residue numbers -> list of contiguous (start, end) runs."""
    out = []
    for r in sorted(numbers):
        if out and r == out[-1][1] + 1:
            out[-1][1] = r
        else:
            out.append([r, r])
    return [tuple(x) for x in out]


def tm_selection(obj):
    return f"{obj} and name CA and ({' or '.join(f'resi {a}-{b}' for a, b in runs(TM_RESI[obj]))})"


def bundle_axis(core_xyz, resi):
    """Mean direction of the TM helices (contiguous TM runs), alternating sign so all point the
    same way across the membrane."""
    v = np.zeros(3)
    for k, (a, b) in enumerate(runs(resi)):
        idx = [n for n, r in enumerate(resi) if a <= r <= b]
        d = core_xyz[idx[-1]] - core_xyz[idx[0]]
        v += d if k % 2 == 0 else -d
    return v / np.linalg.norm(v)


def superpose():
    cmd.select("ref_core", tm_selection(ANCHOR))
    ce = {}
    for name, *_ in STRUCTURES:
        if name != ANCHOR:
            ce[name] = cmd.cealign("ref_core", tm_selection(name))
    resi, _, ref_core = ca_table(ANCHOR, tm_only=True)
    axis = bundle_axis(ref_core, resi)
    tables = {n: ca_table(n, tm_only=True) for n, *_ in STRUCTURES}

    # register search: CE can land a whole helical turn off, which a proximity fit then keeps.
    # Start from the current pose shifted by 0, +-1 and +-2 turns (5.4 A) along the bundle axis,
    # refine each against the target core, keep the pose with the best overlap score. The first
    # pass uses the anchor core; the second uses the mean core of all structures, so a distant
    # family is judged against the whole set rather than against class A alone.
    def register_pass(target):
        result = {}
        for name, *_ in STRUCTURES:
            if name == ANCHOR:
                continue
            xyz0 = tables[name][2]
            trials = []
            for turns in (-2, -1, 0, 1, 2):
                R, t = fit_pairwise(target, xyz0 + turns * 5.4 * axis)
                t = (turns * 5.4 * axis) @ R.T + t
                trials.append((overlap_score(target, xyz0 @ R.T + t), turns, R, t))
            trials.sort(key=lambda x: -x[0])
            apply(name, trials[0][2], trials[0][3])
            tables[name] = ca_table(name, tm_only=True)
            result[name] = (trials[0][1], trials[0][0], trials[1][1], trials[1][0])
        return result

    register_pass(ref_core)
    mean_core, matches = joint_refinement(ref_core, tables)
    register = register_pass(mean_core)
    mean_core, matches = joint_refinement(mean_core, tables)
    # the anchor is never moved; the mean core only defines the targets
    return ce, register, matches, tables, mean_core


def refit_onto(parent, names, tables, axis):
    """Refit each named structure onto one parent's TM core (TM residues only), with the same
    register search as the main fit; the parent does not move."""
    target = tables[parent][2]
    for name in names:
        xyz0 = tables[name][2]
        trials = []
        for turns in (-2, -1, 0, 1, 2):
            R, t = fit_pairwise(target, xyz0 + turns * 5.4 * axis)
            t = (turns * 5.4 * axis) @ R.T + t
            trials.append((overlap_score(target, xyz0 @ R.T + t), R, t))
        best = max(trials, key=lambda x: x[0])
        apply(name, best[1], best[2])
        tables[name] = ca_table(name, tm_only=True)


def joint_refinement(mean_core, tables):
    """Fit every structure to the mean core, rebuild the mean from all of them, repeat."""
    mean_core = mean_core.copy()
    matches = {}
    for cut in CUTOFFS[2:]:
        for name, *_ in STRUCTURES:
            if name == ANCHOR:
                continue
            for _ in range(2):
                xyz = tables[name][2]
                pairs = [(i, j) for i, j in dp_align(mean_core, xyz)
                         if np.linalg.norm(mean_core[i] - xyz[j]) < cut]
                if len(pairs) < 30:
                    break
                ii, jj = map(list, zip(*pairs))
                r, t = kabsch(xyz[jj], mean_core[ii])
                apply(name, r, t)
                tables[name] = ca_table(name, tm_only=True)
        for name, *_ in STRUCTURES:
            xyz = tables[name][2]
            matches[name] = [(i, j) for i, j in dp_align(mean_core, xyz)
                             if np.linalg.norm(mean_core[i] - xyz[j]) < cut]
        acc = [[] for _ in mean_core]
        for name, pairs in matches.items():
            for i, j in pairs:
                acc[i].append(tables[name][2][j])
        mean_core = np.array([np.mean(a, 0) if a else mean_core[k] for k, a in enumerate(acc)])
    return mean_core, matches


def residue_map(obj, row):
    """Alignment column -> structure residue number, via the family row aligned to the chain."""
    resi, aa3, _ = ca_table(obj)
    one = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
           "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
           "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
    target = "".join(one.get(a, "X") for a in aa3)
    cols = [i for i, a in enumerate(row) if a != "-"]
    query = "".join(row[i].upper() for i in cols)
    al = PairwiseAligner(mode="global", open_gap_score=-10, extend_gap_score=-0.5,
                         end_gap_score=0)
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aln = al.align(target, query)[0]
    out = {}
    for (t0, t1), (q0, q1) in zip(*aln.aligned):
        for k in range(t1 - t0):
            out[cols[q0 + k]] = (resi[t0 + k], target[t0 + k])
    return out


def full_alignment_numbers():
    """Exact member residue numbers from the full family alignments (map_motif_residues.py).

    Returns {(structure, position): UniProt number or None for a gap} and
    {accession: (start, sequence)}; both empty when the mapping has not been produced."""
    if not (MAPPING.exists() and MEMBER_SEQS.exists()):
        return {}, {}
    numbers = {}
    for r in csv.DictReader(open(MAPPING), delimiter="\t"):
        res = r["member_residue"]
        numbers[(r["structure"], r["position"])] = int(res[1:]) if res[1:].isdigit() else None
        if r["status"] != "ok":
            print(f"mapping note: {r['structure']} {r['position']}: {r['status']}")
    seqs = {r["accession"]: (int(r["start"]), r["sequence"])
            for r in csv.DictReader(open(MEMBER_SEQS), delimiter="\t")}
    return numbers, seqs


def member_to_structure(obj, start, seq):
    """Member residue number -> (structure residue number, letter), by aligning the member's
    sequence from the full alignment to the chain. For the same protein this only absorbs
    numbering differences (isoforms, fusions, tags); for consOR2 it aligns to its homologue."""
    resi, aa3, _ = ca_table(obj)
    target = "".join(ONE.get(a, "X") for a in aa3)
    al = PairwiseAligner(mode="global", open_gap_score=-10, extend_gap_score=-0.5,
                         end_gap_score=0)
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aln = al.align(target, seq)[0]
    out = {}
    for (t0, t1), (q0, q1) in zip(*aln.aligned):
        for k in range(t1 - t0):
            out[start + q0 + k] = (resi[t0 + k], target[t0 + k])
    return out


def reference_sequence(fam):
    """Full sequence of the family's reference (acc2 where it has one, else acc1), from the
    family's search hits in the searched UniProt release (zenodo_data), not live UniProt."""
    import gzip
    d = json.load(open(MAPS / f"{fam}.json"))
    acc = d["reference"]["accession"]
    path = stage("02_clade_extraction", "clade_sequences", f"{fam}.fasta.gz")
    seq = next(("".join(s) for h, s in _fasta_chunks(gzip.open(path, "rt")) if accession(h) == acc), "")
    if not seq:
        raise SystemExit(f"reference {acc} of {fam} not found in {path}")
    return acc, seq, d["positions"]


def _fasta_chunks(handle):
    name, seq = None, []
    for line in handle:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name is not None:
                yield name, seq
            name, seq = line[1:], []
        else:
            seq.append(line.strip())
    if name is not None:
        yield name, seq


def reference_map(obj, fam):
    """Superfamily column -> structure residue, through the family reference residue.

    The website records each column's reference residue number; references keep every residue
    in the trimmed alignments, so this path has no deleted residues to misplace. The structure's
    sequence is aligned to the reference's full sequence to reach the structure residue."""
    acc, ref, positions = reference_sequence(fam)
    resi, aa3, _ = ca_table(obj)
    one = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
           "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
           "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
    target = "".join(one.get(a, "X") for a in aa3)
    al = PairwiseAligner(mode="global", open_gap_score=-10, extend_gap_score=-0.5,
                         end_gap_score=0)
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aln = al.align(target, ref)[0]
    ref_to_struct = {}
    for (t0, t1), (q0, q1) in zip(*aln.aligned):
        for k in range(t1 - t0):
            ref_to_struct[q0 + k + 1] = (resi[t0 + k], target[t0 + k])
    out = {}
    for col, p in enumerate(positions, start=1):
        r = p.get("refResidue") if p else None
        if r is not None and r in ref_to_struct:
            out[col] = ref_to_struct[r]
    return acc, out


def main():
    global ANCHOR, TAG
    TAG = sys.argv[1] if len(sys.argv) > 1 else "b2ar"
    ANCHOR = ANCHORS[TAG]
    cmd.reinitialize()
    load()
    topologies = read_topologies({s[4] for s in STRUCTURES})
    missing = [s[4] for s in STRUCTURES if s[4] not in topologies]
    if missing:
        raise SystemExit(f"no DeepTMHMM topology for: {', '.join(missing)}")
    assign_tm(topologies)
    ce, register, matches, tables, mean_core = superpose()
    if TAG in HIERARCHY:
        parent, names = HIERARCHY[TAG]
        resi, _, xyz = tables[ANCHOR]
        refit_onto(parent, names, tables, bundle_axis(xyz, resi))
        for name in names:  # core statistics for the refitted structures, against the mean core
            xyz = tables[name][2]
            matches[name] = [(i, j) for i, j in dp_align(mean_core, xyz)
                             if np.linalg.norm(mean_core[i] - xyz[j]) < CUTOFFS[-1]]

    # core fit statistics: per structure RMSD to the mean core on matched pairs
    with open(HERE / f"superposition_core_fit_{TAG}.tsv", "w") as f:
        f.write("structure\tfamily\tsource\tTM_residues\tTM_segments\tCE_RMSD\tCE_aligned\t"
                "turns_from_CE\toverlap_score\trunner_up_turns\trunner_up_score\tcore_pairs\t"
                "core_RMSD\toverlap_final_mean\n")
        for name, fam, _, _, _, src in STRUCTURES:
            p = matches[name]
            d = [np.linalg.norm(tables[name][2][j] - mean_core[i]) for i, j in p]
            c = ce.get(name, {"RMSD": 0.0, "alignment_length": len(p)})
            g = register.get(name, (0, 1.0, 0, 0.0))
            f.write(f"{name}\t{fam}\t{src}\t{len(TM_RESI[name])}\t{len(runs(TM_RESI[name]))}\t"
                    f"{c['RMSD']:.2f}\t{c['alignment_length']}\t{g[0]}\t"
                    f"{g[1]:.3f}\t{g[2]}\t{g[3]:.3f}\t{len(p)}\t"
                    f"{np.sqrt(np.mean(np.square(d))):.2f}\t"
                    f"{overlap_score(mean_core, tables[name][2]):.3f}\n")  # every structure, anchor too

    # map motif columns onto each structure
    positions = POSITIONS + EXTRA
    sup = list(read_fasta(SUP_FASTA))
    acc1s = anchors()
    colnum = {p[0]: p[1] for p in positions}
    site = {}
    placed = {}  # position -> list of (structure, resi, aa, carrier)
    tm3_cys = {}  # structure -> its C3x25 residue number
    carrying = {}  # position -> structures whose family carries the motif
    full, member_seqs = full_alignment_numbers()
    if not full:
        print(f"WARNING: {MAPPING.name} not found; mapping motif columns through the family "
              f"reference (less reliable). Run map_motif_residues.py on the full alignments.")
    disagree = open(HERE / f"mapping_disagreements_{TAG}.tsv", "w")
    disagree.write("position\tstructure\tused\tfull_alignment\tvia_reference\tvia_family_row\n")
    for name, fam, _, _, rowacc, _ in STRUCTURES:
        fcols = family_motif_columns(fam, acc1s[fam], sup, colnum)
        row = next(s for h, s in read_fasta(trimmed_alignment(fam)) if accession(h) == rowacc)
        rmap = residue_map(name, row)
        _, refmap = reference_map(name, fam)
        to_struct = member_to_structure(name, *member_seqs[rowacc]) if rowacc in member_seqs else {}
        skey = SITE_KEY.get(name, fam)
        if skey not in site:
            site[skey] = json.load(open(MAPS / f"{skey}.json"))["positions"]
        for label, col, _, accept, _ in positions:
            c = fcols[label]
            via_row = rmap.get(c) if c is not None else None
            via_ref = refmap.get(col)
            # full family alignment (exact numbers) when available; otherwise the older routes
            via_full = None
            if (name, label) in full:
                number = full[(name, label)]
                via_full = to_struct.get(number) if number else None
            use = via_full if full else (via_ref or via_row)
            fmt = lambda v: f"{v[1]}{v[0]}" if v else "-"
            if len({fmt(v) for v in (via_full, via_ref, via_row) if v}) > 1:
                disagree.write(f"{label}\t{name}\t{fmt(use)}\t{fmt(via_full)}\t{fmt(via_ref)}\t"
                               f"{fmt(via_row)}\n")
            # 45x50: the cysteine that forms the disulfide with the structure's C3x25, wherever the
            # alignment column falls (ECL2 aligns poorly); the column residue where there is none
            bonded = None
            if label == "45x50" and name in tm3_cys:
                bonded = disulfide_partner(name, tm3_cys[name])
                if bonded is not None:
                    if not use or use[0] != bonded:
                        print(f"45x50 {name}: column {fmt(use)}, disulfide partner C{bonded} used")
                    use = (bonded, "C")
            if use is None:
                continue
            resi, aa = use
            if label == "3x25" and aa == "C":
                tm3_cys[name] = resi
            if label == "2x46":
                accept = "L"
            if label in EXTRA_FAMILIES and fam not in EXTRA_FAMILIES[label]:
                continue
            if (label, name) in EXCLUDE:
                continue
            accept = family_accept(FAMILY_ACCEPT.get((label, fam), accept), site[skey][col - 1])
            status = call(site[skey][col - 1], accept)[0]
            family_carries = status in ("conserved", "consensus") or bonded is not None
            carrier = family_carries and aa in accept
            placed.setdefault(label, []).append((name, resi, aa, carrier))
            if family_carries:
                carrying.setdefault(label, set()).add(name)
    disagree.close()

    # selections and per-residue geometry
    def points(name, resi):
        # get_model, not get_coords: this PyMOL build's numpy bridge fails under numpy 2
        atoms = cmd.get_model(f"{name} and resi {resi} and not hydro and not alt B").atom
        ca = next((np.array(a.coord) for a in atoms if a.name == "CA"), None)
        cb = next((np.array(a.coord) for a in atoms if a.name == "CB"), None)
        sc = [a.coord for a in atoms if a.name not in BACKBONE]
        return ca, (np.mean(sc, 0) if sc else ca), cb

    # side-chain stem direction (C-alpha -> C-beta) relative to the receptor core: the facing angle
    # is 0 deg pointing at the bundle axis, 90 deg tangential, 180 deg pointing away from the axis
    # (geometry only: it does not say what the side chain contacts)
    a_resi, _, a_xyz = ca_table(ANCHOR, tm_only=True)
    core_centre, core_axis = a_xyz.mean(0), bundle_axis(a_xyz, a_resi)

    def stem(ca, cb):
        return None if cb is None else (cb - ca) / np.linalg.norm(cb - ca)

    def facing(ca, cb):
        v = stem(ca, cb)
        if v is None:
            return None
        r = ca - core_centre
        r = r - (r @ core_axis) * core_axis
        return float(np.degrees(np.arccos(np.clip(v @ (-r / np.linalg.norm(r)), -1, 1))))

    def angle(v, w):
        return float(np.degrees(np.arccos(np.clip(v @ w, -1, 1))))

    fam_of = {s[0]: s[1] for s in STRUCTURES}
    src_of = {s[0]: s[5] for s in STRUCTURES}
    spread_rows, dev_rows, struct_rows = [], [], []
    rms = lambda a, m: float(np.sqrt(np.mean(np.sum((a - m) ** 2, 1))))
    pair = lambda a: float(np.mean([np.linalg.norm(x - y) for x, y in combinations(a, 2)]))
    for label, col, cls, accept, pattern in positions:
        if label == "2x46":
            accept = "L"
        entries = placed.get(label, [])
        sel = selection_name(label, cls)
        carriers = [e for e in entries if e[3]]
        others = [e for e in entries if not e[3]]
        if carriers:  # every residue carrying the motif at the column, whether or not it aligns
            cmd.select(sel, " or ".join(f"({n} and resi {r})" for n, r, _, _ in carriers))
            cmd.deselect()
        pts = {n: points(n, r) for n, r, _, _ in carriers}
        pts = {n: p for n, p in pts.items() if p[0] is not None}

        # structural equivalents: families that carry the motif but whose column residue lies
        # off-position (or is not the motif residue) get the nearby motif-class residue that
        # sits at the cAMP/Rhodopsin carriers' position, if one lies within STRUCT_MAX
        ref_names = [n for n in pts if LINEAGE_OF[fam_of[n]] in ("cAMP", "Rhodopsin")]
        struct = []
        if len(ref_names) >= 2:
            centre = np.mean([pts[n][0] for n in ref_names], 0)
            for n in sorted(carrying.get(label, ())):
                entry = next((e for e in entries if e[0] == n), None)
                if entry is None:
                    continue
                d_col = float(np.linalg.norm(pts[n][0] - centre)) if n in pts else None
                if d_col is not None and d_col <= STRUCT_MAX:
                    continue
                around = entry[1]
                skey = SITE_KEY.get(n, fam_of[n])
                fam_accept = family_accept(FAMILY_ACCEPT.get((label, fam_of[n]), accept),
                                           site[skey][col - 1])
                cands = [(float(np.linalg.norm(np.array(a.coord) - centre)), a)
                         for a in cmd.get_model(f"{n} and name CA and resi "
                                                f"{around - STRUCT_WINDOW}-{around + STRUCT_WINDOW}").atom
                         if ONE.get(a.resn, "X") in fam_accept
                         and family_carries_at(site[skey], col, int(a.resi) - around, fam_accept)]
                if cands:
                    d, a = min(cands, key=lambda x: x[0])
                    if d <= STRUCT_MAX and int(a.resi) != around:
                        struct.append((n, int(a.resi), ONE[a.resn], entry, d_col, d))
        # <sel>_struct: structural equivalents, shown alongside the column residues in <sel>
        if struct:
            cmd.select(f"{sel}_struct", " or ".join(f"({n} and resi {r})" for n, r, *_ in struct))
            cmd.deselect()
            for n, r, aa, entry, d_col, d in struct:
                struct_rows.append((label, n, fam_of[n], f"{entry[2]}{entry[1]}",
                                    "" if d_col is None else f"{d_col:.2f}", f"{aa}{r}", f"{d:.2f}"))
        with_struct = dict(pts)
        for n, r, *_ in struct:
            with_struct[n] = points(n, r)

        if len(pts) >= 2:
            ca = np.array([p[0] for p in pts.values()])
            sc = np.array([p[1] for p in pts.values()])
            ca_m, sc_m = ca.mean(0), sc.mean(0)
            ca2 = np.array([p[0] for p in with_struct.values()])
            sc2 = np.array([p[1] for p in with_struct.values()])
            by_lin = {}
            for n in pts:
                by_lin.setdefault(LINEAGE_OF[fam_of[n]], []).append(n)
            # stems (with structural equivalents where found)
            stems = {n: stem(p[0], p[2]) for n, p in with_struct.items()}
            stems = {n: v for n, v in stems.items() if v is not None}
            faces = {n: facing(p[0], p[2]) for n, p in with_struct.items() if p[2] is not None}
            if len(stems) >= 2:
                m = np.sum(list(stems.values()), 0)
                m = m / np.linalg.norm(m)
                to_mean = {n: angle(v, m) for n, v in stems.items()}
                disp = float(np.mean(list(to_mean.values())))
            else:
                to_mean, disp = {}, float("nan")
            fvals = list(faces.values())
            spread_rows.append((label, col, cls, pattern, len(pts),
                                sum(src_of[n] == "PDB" for n in pts), rms(ca, ca_m), pair(ca),
                                rms(sc, sc_m), pair(sc),
                                " ".join(f"{k}:{len(v)}" for k, v in by_lin.items()),
                                len(struct), rms(ca2, ca2.mean(0)), rms(sc2, sc2.mean(0)),
                                float(np.mean(fvals)) if fvals else float("nan"),
                                float(np.std(fvals)) if fvals else float("nan"), disp))
            struct_of = {n: (r, aa) for n, r, aa, *_ in struct}
            for n, p in with_struct.items():
                if n in struct_of:
                    resid = f"{struct_of[n][1]}{struct_of[n][0]} (structural)"
                else:
                    res = next(e for e in carriers if e[0] == n)
                    resid = f"{res[2]}{res[1]}"
                f = faces.get(n)
                dev_rows.append((label, n, fam_of[n], LINEAGE_OF[fam_of[n]], src_of[n], resid,
                                 np.linalg.norm(p[0] - ca2.mean(0)),
                                 np.linalg.norm(p[1] - sc2.mean(0)),
                                 "" if f is None else f"{f:.0f}",
                                 "" if n not in to_mean else f"{to_mean[n]:.0f}"))
    with open(HERE / f"motif_position_spread_{TAG}.tsv", "w") as f:
        f.write("position\tsuperfamily_column\tmotif_class\tlineage_pattern\tn_structures\t"
                "n_experimental\tCA_spread_RMSD\tCA_mean_pairwise\tSC_spread_RMSD\t"
                "SC_mean_pairwise\tstructures_per_lineage\tn_structural_equivalents\t"
                "CA_spread_RMSD_with_equivalents\tSC_spread_RMSD_with_equivalents\t"
                "stem_angle_to_axis_mean_deg\tstem_angle_to_axis_sd_deg\tstem_divergence_deg\n")
        for r in sorted(spread_rows, key=lambda r: r[12]):  # by C-alpha RMSD with equivalents
            f.write("\t".join(map(str, r[:6])) + "\t" + "\t".join(f"{x:.2f}" for x in r[6:10])
                    + f"\t{r[10]}\t{r[11]}\t{r[12]:.2f}\t{r[13]:.2f}\t{r[14]:.0f}\t{r[15]:.0f}\t"
                    f"{r[16]:.0f}\n")
    with open(HERE / f"structural_equivalents_{TAG}.tsv", "w") as f:
        f.write("position\tstructure\tfamily\tcolumn_residue\tcolumn_CA_to_centre\t"
                "structural_equivalent\tequivalent_CA_to_centre\n")
        for r in struct_rows:
            f.write("\t".join(r) + "\n")
    with open(HERE / f"motif_residue_deviation_{TAG}.tsv", "w") as f:
        f.write("position\tstructure\tfamily\tlineage\tsource\tresidue\tCA_dev_from_mean\t"
                "SC_dev_from_mean\tstem_angle_to_axis_deg\tstem_angle_to_mean_deg\n")
        for r in dev_rows:
            f.write("\t".join(map(str, r[:6])) + f"\t{r[6]:.2f}\t{r[7]:.2f}\t{r[8]}\t{r[9]}\n")
    with open(HERE / f"motif_residues_per_structure_{TAG}.tsv", "w") as f:
        f.write("position\tstructure\tfamily\tresidue\tcarries_motif\n")
        for label, entries in placed.items():
            for n, r, aa, carrier in entries:
                f.write(f"{label}\t{n}\t{fam_of[n]}\t{aa}{r}\t{carrier}\n")

    # presentation: trim to the 7TM boundary (the fit above used TM residues only, so trimming
    # here changes nothing measured); motif residues are always kept
    for name, *_ in STRUCTURES:
        keep = set(KEEP_RESI[name])
        for entries in placed.values():
            keep |= {r for n, r, _, _ in entries if n == name}
        keep |= {int(r[5][1:]) for r in struct_rows if r[1] == name}
        spans = " or ".join(f"resi {a}-{b}" for a, b in runs(keep))
        cmd.remove(f"{name} and not ({spans})")
    cmd.bg_color("white")
    cmd.hide("everything")
    cmd.show("cartoon", ANCHOR)  # only the template in cartoon; the other structures show motifs only
    cmd.set("cartoon_transparency", 0.6)
    for name, fam, *_ in STRUCTURES:
        cmd.color(COLOURS[LINEAGE_OF[fam]], name)
    cmd.util.cnc("all")  # side chains: lineage-coloured carbons, element-coloured N/O/S
    for label, _, cls, *_ in positions:
        sel = selection_name(label, cls)
        if sel in cmd.get_names("selections"):
            cmd.show("sticks", f"{sel} and not name N+C+O")
        if f"{sel}_struct" in cmd.get_names("selections"):
            cmd.show("sticks", f"{sel}_struct and not name N+C+O")
    cmd.delete("ref_core")
    cmd.orient(ANCHOR)
    cmd.save(str(HERE / f"motif_superposition_tmcore_{TAG}.pse"))


if __name__ == "__main__":
    main()
