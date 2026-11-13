"""Assemble the supplementary tables as Excel workbooks, one per table, with sheets for the lettered parts.

Reads the Table_S*.tsv files the build_S*.py scripts write to supplementary_tables/ in the archive, replaces internal family
names with display names (classA -> Class A, FSLB -> FSL, ...) and the literature "tight / loose" labels with
"matched / loose match", and puts the table title and legend in the first two rows of every sheet (table from row 4). Numbers stay numbers.

    python build_xlsx.py
"""
import csv
import re
from pathlib import Path

from openpyxl import Workbook
from openpyxl.styles import Alignment, Font
from openpyxl.utils import get_column_letter

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import stage  # noqa: E402

HERE = stage("supplementary_tables")
OUT = HERE / "xlsx"

NAMES = {"classA": "Class A", "classB1": "Class B1", "classB2": "Class B2", "classC": "Class C",
         "classF": "Class F", "classT": "Taste-2", "FSLB": "FSL", "Vomeronasal1": "Vomeronasal 1",
         "Vomeronasal2": "Vomeronasal 2"}
NAME_RE = re.compile(r"(?<![A-Za-z0-9_])(" + "|".join(NAMES) + r")(?![A-Za-z0-9_])")
HEADER = [(re.compile(r"\btight\b"), "matched"), (re.compile(r"\bloose\b"), "loose match"),
          (re.compile(r"\bpct\b"), "%")]

TABLES = {
    "S1": ("Search and filtering per family", [
        ("S1", "Table_S1_search_and_filtering.tsv",
         "Per family and control: search method and query (Pfam / KEGG profile HMM or iterative MMseqs2 seed), "
         "families sharing one search, sequences retrieved, the gene tree used for clade extraction (full-length or "
         "HMM-matched region), number of extracted clades, and sequence counts after clade extraction, the single "
         "7TM core filter, the topology filter (score >= 0.8) and the overlap assignment.")]),
    "S2": ("Sequences claimed by two families", [
        ("S2", "Table_S2_contested_sequences.tsv",
         "The 382 sequences that remained in two families after the topology filter, the families claiming each, "
         "the family it was assigned to and how: by its clade in the assignment tree or by individual review.")]),
    "S3": ("Comparison with published GPCR surveys", [
        ("S3a receptors", "Table_S3a_literature_receptors.tsv",
         "Every receptor reported in 15 published surveys and its fate in our pipeline: identifier or sequence "
         "match (same species or genus, >= 90% identity), clade extraction, 7TM and topology filters, final family, "
         "and the nearest family of the best MMseqs2 hit where no match was found (indicator only)."),
        ("S3b by study", "Table_S3b_literature_agreement_by_study.tsv",
         "Family agreement per study for matched receptors and for loose matches (best MMseqs2 hit)."),
        ("S3c by family", "Table_S3c_literature_agreement_by_family.tsv",
         "Family agreement per published family, with the families that disagreeing receptors fall in."),
        ("S3d deeper origins", "Table_S3d_deeper_origin_candidates.tsv",
         "Literature receptors from lineages outside our inferred origin of their family, and whether each was "
         "tested in the trees.")]),
    "S4": ("Taxonomic distribution of GPCR families across eukaryotic lineages", [
        ("S4", "Table_S4_lineage_distribution.tsv",
         "The data behind Figure 2. For each family and control and each of the lineages of the figure: genomes "
         "sampled and genomes carrying the family, genera sampled and genera carrying it, and how the presence is "
         "drawn (present, sporadic, or attributed to contamination or horizontal gene transfer; the number of "
         "genomes of each kind). The last columns give the Dollo origin (k = 2, NCBI taxonomy of all searched "
         "genomes) and our interpretation of the origin.")]),
    "S5": ("Representative sequences", [
        ("S5", "Table_S5_representatives.tsv",
         "The 590 family representatives (family, taxonomy, topology score) and the 18 literature sequences, with "
         "the tree sets each belongs to.")]),
    "S6": ("Support for the superfamily hierarchy", [
        ("S6a tree runs", "Table_S6a_tree_runs.tsv",
         "The eight maximum likelihood inferences: sequence set, alignment, columns after trimming, model, "
         "log-likelihood, ML searches and bootstrap trees."),
        ("S6b family stability", "Table_S6b_family_stability.tsv",
         "Per family: median of the mean transfer index, sequences at or above the STE2 cutoff, clade F1, stem "
         "length and sister group in the trees with all families."),
        ("S6c family stems", "Table_S6c_family_stems.tsv",
         "Stem length of each family across the 100 ML searches of each run, and how often it is among the three "
         "longest."),
        ("S6d groupings", "Table_S6d_groupings.tsv",
         "Each family grouping: TBE in the best tree, recovery in the 100 ML searches, and median site concordance "
         "factor, for each of the eight runs."),
        ("S6e placements", "Table_S6e_placements.tsv",
         "Placement of named sequences and tests (non-animal family members, order within Class B and Class F, "
         "split dating) in the best tree and the ML searches."),
        ("S6f literature placement", "Table_S6f_literature_placement.tsv",
         "Placement of the 18 literature sequences relative to the clade of their claimed family."),
        ("S6g splits sCF main", "Table_S6g_family_splits_sCF_ep0.123.tsv",
         "Every family-level split: recovery in the ML searches and median sCF, main alignment (G-INS-i --ep 0.123)."),
        ("S6h splits sCF plain", "Table_S6h_family_splits_sCF_ginsi.tsv",
         "As S6g, alignment with plain G-INS-i.")]),
    "S7": ("Placement stability of every representative", [
        ("S7", "Table_S7_sequence_tindex.tsv",
         "Transfer index of every representative in each tree it was in, its mean, and whether the mean is at or "
         "above that of the most stable STE2 sequence (0.692), the cutoff for exclusion from interpretation.")]),
    "S8": ("Conservation at the superfamily motif positions", [
        ("S8", "Table_S8_motif_conservation.tsv",
         "Superfamily alignment built with the maximum likelihood guide tree (G-INS-i --ep 0.123). Rows: the 28 "
         "positions conserved in >= 12 of 18 families (Fig. 4) and every other position named in the text. Cells: "
         "most frequent residue and the share of sequences carrying it or a similar residue; * marks >= 60%. cAMP "
         "and Class C are also given per kingdom.")]),
    "S9": ("Conserved residues shared by STE3 and STE2 with other families", [
        ("S9", "Table_S9_STE2_STE3_shared_residues.tsv",
         "Number of alignment columns where both families have a consensus residue reaching the threshold "
         "(>= 40% or >= 60% of sequences) and it is the identical residue, in both MSAs built with the maximum "
         "likelihood guide tree.")]),
    "S10": ("Structural overlap of the shared motifs", [
        ("S10a structures", "Table_S10a_structures.tsv",
         "Structures superimposed: family, experimental structure or AlphaFold model, TM residues used and core "
         "RMSD after refinement against the mean TM core."),
        ("S10b motif residues", "Table_S10b_motif_residue_distances.tsv",
         "Per motif position and structure: residue, distance of its CA from the centroid of the other carrying "
         "residues, variation of the local backbone (+-3 residues), and the call; STE2 residues at each site are "
         "given as the control.")]),
}


def value(v):
    v = NAME_RE.sub(lambda m: NAMES[m.group(1)], v)
    try:
        f = float(v)
        return int(f) if re.fullmatch(r"-?\d+", v) else f
    except ValueError:
        return v


def header(h):
    h = NAME_RE.sub(lambda m: NAMES[m.group(1)], h)
    for pat, rep in HEADER:
        h = pat.sub(rep, h)
    h = h.replace("_", " ")
    return h[:1].upper() + h[1:]


def main():
    OUT.mkdir(exist_ok=True)
    bold = Font(bold=True)
    for tid, (title, sheets) in TABLES.items():
        wb = Workbook()
        wb.remove(wb.active)
        for name, src, desc in sheets:
            ws = wb.create_sheet(name[:31])
            with open(HERE / src, encoding="utf-8") as fh:
                rows = list(csv.reader(fh, delimiter="	"))
            ncol = len(rows[0])
            part = name.split()[0]
            ws.append([f"Supplementary Table {part}. {title}" if part == tid else
                       f"Supplementary Table {part}. {title}: {name.split(' ', 1)[1]}"])
            ws.append([desc])
            ws.append([])
            ws["A1"].font = Font(bold=True, size=12)
            ws.merge_cells(start_row=2, start_column=1, end_row=2, end_column=max(ncol, 8))
            ws["A2"].alignment = Alignment(wrap_text=True, vertical="top")
            ws.row_dimensions[2].height = 15 * (1 + len(desc) // 140)
            ws.append([header(h) for h in rows[0]])
            for r in rows[1:]:
                ws.append([value(v) for v in r])
            for c in ws[4]:
                c.font = bold
                c.alignment = Alignment(wrap_text=True, vertical="top")
            ws.freeze_panes = "B5"
            for i, h in enumerate(rows[0], 1):
                longest = max([len(h) // 2] + [len(r[i - 1]) for r in rows[1:200] if len(r) >= i])
                ws.column_dimensions[get_column_letter(i)].width = min(max(10, longest + 2), 60)
        path = OUT / f"Supplementary_Table_{tid}.xlsx"
        try:
            wb.save(path)
            print("wrote", path.name, [s[0] for s in sheets])
        except PermissionError:
            print("SKIPPED", path.name, "- open in Excel; close it and rerun")


if __name__ == "__main__":
    main()
