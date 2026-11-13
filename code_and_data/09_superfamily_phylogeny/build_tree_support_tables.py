"""Supplementary tables for the hierarchy section, from the Oct08 readouts (no new tree analysis).

  S_tree_runs.tsv             written separately (run summary from the RAxML logs)
  S_family_stability.tsv      per family: clade F1 and stem in the two "all" trees, median mean
                              tIndex, sequences at/above the STE2 cutoff
  S_sequence_tindex.tsv       per sequence: family, taxonomy, tIndex in each tree, mean, flag
  S_family_stems.tsv          per family: mean stem over the searches of the 4 runs with divergent
                              families, top-3 counts, and pooled over all 4 runs
  S_groupings.tsv             per grouping and run: TBE in the best tree, ML searches recovering it,
                              median sCF (ginsi_ep0.123 is the main alignment)
  S_placements.tsv            named non-animal sequences and lineage tests: best tree + searches
  S_literature.tsv            literature sequences: in claimed family clade, distance, tIndex
  S_family_splits_sCF_<aln>.tsv   every family-level split with recovery and median sCF (step 10)

Usage: python build_tree_support_tables.py
"""
import csv
import shutil
from collections import defaultdict
from pathlib import Path

import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import stage  # noqa: E402
S9 = stage("09_superfamily_phylogeny")
OUT = S9 / "tree_support_tables"
IHA = S9 / "readouts"
TI = S9 / "readouts"
S10 = S9 / "sCF"
TAXA = stage("07_representative_selection", "representatives_all_families.tsv")
TINDEX = S9 / "booster_stats" / "tindex_mean_across_final_trees.tsv"
CUTOFF = 0.6919        # most stable STE2 sequence

NICE = {"classA": "Class A", "classB1": "Class B1", "classB2": "Class B2", "classC": "Class C",
        "classF": "Class F", "classT": "Taste-2", "Vomeronasal1": "Vomeronasal 1",
        "Vomeronasal2": "Vomeronasal 2", "Mth": "Methuselah", "GP143": "GPR143", "GP157": "GPR157",
        "TM116": "TMEM116", "FSLB": "FSL", "LIT_GPR1like": "GPR1-like (lit.)", "LIT_GP157": "GPR157 (lit.)",
        "LIT_classT": "Taste-2 (lit.)", "LIT_Vomeronasal1": "Vomeronasal 1 (lit.)"}
nice = lambda f: NICE.get(f, f)
RUN_ORDER = [f"{s}, {a}" for s in ("all", "no STE2", "core", "core+lit") for a in ("ep0.123", "ginsi")]
DATASET = {"all": "all", "no STE2": "no_ste2", "core": "no_ste2ste3classCV2", "core+lit": "and_literature"}


def read(path):
    return list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))


def write(name, header, rows):
    with open(OUT / name, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)
    print(f"{name}: {len(rows)} rows")


tax = {r["accession"]: r for r in read(TAXA)}
ti_rows = read(TINDEX)

# --- per sequence tIndex
trees = [c for c in ti_rows[0] if c not in ("accession", "family", "n_trees", "mean_tIndex", "max_tIndex", "leaf")]
rows = []
for r in ti_rows:
    t = tax.get(r["accession"], {})
    rows.append([r["accession"], nice(r["family"]), t.get("kingdom", ""), t.get("phylum", ""),
                 t.get("species", r["leaf"].split("|")[3].replace("_", " ") if r["leaf"].startswith("LIT") else ""),
                 r["n_trees"], r["mean_tIndex"], "yes" if float(r["mean_tIndex"]) >= CUTOFF else "no"]
                + [r[c] for c in trees])
write("S_sequence_tindex.tsv", ["Accession", "Family", "Kingdom", "Phylum", "Species", "Trees",
                               "Mean tIndex", f"At/above STE2 cutoff ({CUTOFF:.3f})"] + trees, rows)

# --- per family stability (units: nested lineages merged, as in tree_metrics.py)
metrics = {(r["run"], r["unit"]): r for r in read(IHA / "tree_unit_metrics.tsv")}
UNIT = {"classA": "Class A-like", "Olfactory": "Class A-like", "Vomeronasal1": "Class A-like",
        "classT": "Class A-like", "Nematode": "Class A-like", "GPR1": "GPR1", "classB1": "Class B",
        "classB2": "Class B", "Mth": "Class B", "classF": "Class F/FSL", "FSLB": "Class F/FSL",
        "GP143": "GP143", "TM116": "TM116", "GP157": "GP157", "cAMP": "cAMP", "classC": "Class C/V2",
        "Vomeronasal2": "Class C/V2", "STE2": "STE2", "STE3": "STE3"}
by_unit = defaultdict(list)
for r in ti_rows:
    if r["family"] in UNIT:
        by_unit[UNIT[r["family"]]].append(float(r["mean_tIndex"]))
rows = []
for u, v in sorted(by_unit.items(), key=lambda x: -sorted(x[1])[len(x[1]) // 2]):
    v = sorted(v)
    med = v[len(v) // 2] if len(v) % 2 else (v[len(v) // 2 - 1] + v[len(v) // 2]) / 2
    row = [u, len(v), f"{med:.3f}", sum(x >= CUTOFF for x in v)]
    for run in ("all, ep0.123", "all, ginsi"):
        m = metrics.get((run, u), {})
        row += [m.get("F1", ""), m.get("stem", ""), m.get("sister_in_best_tree", "")]
    rows.append(row)
write("S_family_stability.tsv",
      ["Family (nested lineages merged)", "Sequences", "Median of mean tIndex", "At/above STE2 cutoff",
       "Clade F1, all ep0.123", "Stem length, all ep0.123", "Sister, all ep0.123",
       "Clade F1, all ginsi", "Stem length, all ginsi", "Sister, all ginsi"], rows)

# --- stems: per run over the 100 searches, and pooled (GP143 + TM116 merged; stem_ranks_Oct08.py)
rows = [[r["run"], r["unit"], r["mean_stem"], r["median_stem"], r["in_top3"], r["longest"], r["of"]]
        for r in read(IHA / "stem_ranks.tsv")]
write("S_family_stems.tsv", ["Run", "Family (nested lineages merged; GP143 + TM116 as one)",
                             "Mean stem", "Median stem", "Searches with a top-3 stem",
                             "Searches with the longest stem", "Searches"], rows)

# --- groupings: TBE, recovery, sCF
SPLIT = {"Taste-2 + V1": "classT+Vomeronasal1", "GP143 + TM116": "GP143+TM116",
         "Class B1 + B2": "classB1+classB2", "Class B (B1+B2+Mth)": "classB1+classB2+Mth",
         "Class F + FSL": "classF+FSLB", "Class A-like": "classA+classT+Nematode+Olfactory+Vomeronasal1",
         "Class A-like + GPR1": "classA+classT+GPR1+Nematode+Olfactory+Vomeronasal1",
         "Class B + Class F/FSL": "classB1+classB2+classF+FSLB+Mth",
         "Class B + F/FSL + GP143 + TM116": "classB1+classB2+classF+FSLB+GP143+Mth+TM116",
         "cAMP (non-fungal) + GP157": "cAMP+GP157"}
scf = {}
for aln in ("ginsi_ep0.123", "ginsi"):
    short = "ep0.123" if aln == "ginsi_ep0.123" else "ginsi"
    for r in read(S10 / aln / "grouping_discovery_divergent_excluded_long.tsv"):
        lab = {"all": "all", "no STE2": "no STE2", "core (no STE2/STE3/C/V2)": "core",
               "core + literature": "core+lit"}[r["label"]]
        for side in (r["side_A"], r["side_B"]):
            scf[(f"{lab}, {short}", side)] = r["median_sCF"]
tbe = {(r["run"], r["grouping"]): r for r in read(IHA / "core_groupings_tbe_recovery.tsv")}
rows = []
for g in dict.fromkeys(r["grouping"] for r in read(IHA / "core_groupings_tbe_recovery.tsv")):
    row = [g]
    for run in RUN_ORDER:
        t = tbe[(run, g)]
        row += [t["TBE_best_tree"], t["ml_trees_recovered_of_100"],
                scf.get((run, SPLIT[g]), "") if g in SPLIT else "n/a"]
    rows.append(row)
hdr = ["Grouping"] + [f"{run} {k}" for run in RUN_ORDER for k in ("TBE", "recovered/100", "median sCF")]
write("S_groupings.tsv", hdr, rows)

# --- placements of named sequences and lineage tests
rows = []
for r in read(TI / "named_sequences.tsv"):
    rows.append([r["run"], f"{nice(r['family'])}: {r['species']} ({r['accession']})",
                 "inside the family clade" if r["family"] != "Mth" else "inside the animal Mth clade",
                 r["in_family_clade_best_tree"], r["ml_trees"], r["of"], r["mean_tIndex"], r["unstable"]])
for r in read(TI / "nonanimal_classA.tsv"):
    rows.append([r["run"], f"Class A, {r['lineage']}: {r['species']} ({r['accession']})",
                 "inside animal Class A-like clade / inside GPR1 clade / neither",
                 r["best_tree"], f"{r['ml_in_ClassA_like']} / {r['ml_in_GPR1']} / {r['ml_neither']}",
                 r["of"], r["mean_tIndex"], r["unstable"]])
for r in read(TI / "internal_order.tsv"):
    rows.append([r["run"], "stable leaves: " + r["test"], "test holds", r["best_tree"], r["ml_trees"],
                 r["of"], "", ""])
for r in read(IHA / "split_dating.tsv"):
    if r["trees"] != "ML searches":
        continue
    rows.append([r["run"], f"{r['split']}: {r['test_group']}",
                 "mean leaves on first side / second side / outside per search", "",
                 f"{r['mean_on_first_side']} / {r['mean_on_second_side']} / {r['mean_outside']}",
                 r["n_trees"], "", ""])
write("S_placements.tsv", ["Run", "Sequence or test", "Measure", "Best tree", "ML searches", "Of",
                           "Mean tIndex", "At/above STE2 cutoff"], rows)

# --- literature
rows = []
for r in read(S9 / "readouts" / "literature_placement.tsv"):
    rows.append([r["alignment"], r["accession"], nice(r["candidate"]), r["species"].replace("_", " "),
                 r["family_clade"], r["in_clade_best_tree"], r["in_clade_ml_trees"], r["of"],
                 r["foreign_best_tree"], r["foreign_ml_median"], r["nearest_best_tree"],
                 r["mean_tIndex"], "yes" if r["mean_tIndex"] and float(r["mean_tIndex"]) >= CUTOFF else "no"])
write("S_literature.tsv", ["Alignment", "ID", "Claimed family", "Species", "Family clade tested",
                           "In clade, best tree", "In clade, ML searches", "Of",
                           "Other-family leaves to join the clade, best tree",
                           "Other-family leaves to join the clade, median over searches",
                           "Nearest families, best tree", "Mean tIndex", "At/above STE2 cutoff"], rows)

# --- full split tables from step 10
for aln in ("ginsi_ep0.123", "ginsi"):
    shutil.copy(S10 / aln / "grouping_discovery_divergent_excluded_wide.tsv",
                OUT / f"S_family_splits_sCF_{aln}.tsv")
    print(f"S_family_splits_sCF_{aln}.tsv: copied")
