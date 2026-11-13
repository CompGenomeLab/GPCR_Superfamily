"""Supplementary tables of agreement between published family labels and ours, per study and per family.

Tight: receptors matched by identifier or by >= 0.9 identity in the same species or genus that are
final in one of our families; agreement = final in the published family.
Loose: receptors without such a match; agreement = closest sequence in our data is in the published family.
Receptors of Nematode and Olfactory that end in Class A count as agreement (both families are nested in Class A);
how many of the agreeing receptors were placed in Class A is reported in its own column.
Only study families with a counterpart among ours are included.
"""
from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import stage  # noqa: E402

OUT = stage("05_literature_comparison")
TIGHT_ROUTES = ["ID match", "species match", "genus match"]
AGREE = {"final in claimed family", "final in its family"}
NESTED = {"Nematode", "Olfactory"}  # nested in Class A
DISPLAY = {"classA": "Class A", "classB1": "Class B1", "classB2": "Class B2", "classC": "Class C",
           "classF": "Class F", "classT": "Taste-2", "FSLB": "FSL", "Vomeronasal1": "Vomeronasal 1",
           "Vomeronasal2": "Vomeronasal 2"}


def show(fams):
    return ", ".join(DISPLAY.get(f, f) for f in fams.split(",") if f)


outcomes = pd.read_csv(OUT / "literature_all_outcomes.tsv", sep="\t", dtype=str).fillna("")
supp = pd.read_csv(OUT / "Supplementary_Table_literature_receptors.tsv", sep="\t", dtype=str).fillna("")
assert len(outcomes) == len(supp)
outcomes["final"] = supp["step8_final_family"].values

# family of each hit protein, for the closest-sequence (loose) placements
hit_family = {}
for name in ["literature_mmseqs_besthit_all.tsv", "unmatched_literature_mmseqs_besthit.tsv",
             "remaining_literature_mmseqs_besthit.tsv"]:
    b = pd.read_csv(OUT / name, sep="\t", dtype=str).fillna("")
    for hit, fam in list(zip(b.best_hit, b.best_family)) + list(zip(b.best_gpcr_hit, b.best_gpcr_family)):
        if hit and fam:
            hit_family.setdefault(hit, fam)

# study families linked to ours by where their receptors end (agree by construction)
unmapped = pd.read_csv(OUT / "literature_unmapped_families.tsv", sep="\t", dtype=str).fillna("")
linked = unmapped[unmapped.status.str.startswith("mapped to")]
by_outcome = set(zip(linked.study, linked.their_family))

d = outcomes[outcomes.our_families != ""].copy()
d["linked_by_outcome"] = [(s, f) in by_outcome for s, f in zip(d.study, d.their_family)]

tight = d[d.route.isin(TIGHT_ROUTES)].copy()
tight = tight[tight.outcome.isin(AGREE) | tight.outcome.str.startswith("assigned to another")
              | tight.outcome.isin(["final in another family", "final in our classC"])]
tight["placed"] = tight.final
loose = d[d.route == "indicator"].copy()
loose["placed"] = loose.our_accession.map(hit_family).fillna("")
loose.loc[loose.outcome == "nearest: claimed family", "placed"] = "(claimed)"
loose.loc[loose.outcome == "nearest: another family", "placed"] = loose.placed.replace("", "family not resolved")
loose.loc[~loose.outcome.isin(["nearest: claimed family", "nearest: another family"]), "placed"] = "none of our families"

for x, ok in [(tight, tight.outcome.isin(AGREE)), (loose, loose.outcome == "nearest: claimed family")]:
    x["agree"] = ok
    x["nested"] = ~ok & x.our_families.isin(NESTED) & (x.placed == "classA")
    x["agree"] = ok | x.nested
    x["mismatch"] = ~x.agree


def pct(a, n):
    return round(100 * a / n, 1) if n else None


def mismatch_targets(x):
    m = x[x.mismatch & (x.placed != "none of our families")].placed
    return "; ".join(f"{DISPLAY.get(f, f)} {n}" for f, n in m.value_counts().items())


def summarise(key):
    rows = []
    for k in sorted(set(tight[key]) | set(loose[key])):
        t, l = tight[tight[key] == k], loose[loose[key] == k]
        r = {key: k}
        for tag, x in [("tight", t), ("loose", l)]:
            r[f"{tag}_compared"] = len(x)
            r[f"{tag}_agree"] = int(x.agree.sum())
            r[f"{tag}_agree_placed_in_Class_A"] = int(x.nested.sum())
            r[f"{tag}_agree_pct"] = pct(x.agree.sum(), len(x))
        if key == "our_families":
            r["tight_mismatches_go_to"] = mismatch_targets(t)
            r["loose_mismatches_go_to"] = mismatch_targets(l)
        r["includes_family_linked_by_outcome"] = bool(t.linked_by_outcome.any() or l.linked_by_outcome.any())
        rows.append(r)
    return pd.DataFrame(rows)


study = summarise("study").sort_values("tight_compared", ascending=False)
family = summarise("our_families")
names = d.groupby("our_families").their_family.apply(lambda s: "; ".join(sorted(set(s))))
family.insert(0, "published_family", family.our_families.map(names))
family["our_families"] = family.our_families.map(show)
family = family.sort_values(["tight_agree_pct", "loose_agree_pct"], na_position="last")

notes = pd.DataFrame({"note": [
    "Only study families that correspond to one of our families are included.",
    "Tight: matched by identifier, or by MMseqs2 at >= 0.9 identity over the longer sequence to a protein "
    "from the same species or genus, and final in one of our families. Agreement: final in the published family.",
    "Loose (MMseqs2 best hit): no such match (mostly genomes outside our set). Agreement: the closest sequence in our data "
    "(E <= 1e-5, query coverage >= 0.5) is in the published family.",
    "Agreement includes Nematode or Olfactory receptors placed in Class A, within which both families are nested; "
    "agree_placed_in_Class_A gives how many of the agreeing receptors these are (Krishnan 2014 srw, Niimura 2009 group eta and amphioxus ORs).",
    "Mismatch columns list families of ours receiving the mismatched receptors; loose mismatches with no family "
    "of ours are not listed, and 'family not resolved' marks closest hits whose family is not in the best-hit tables.",
    "includes_family_linked_by_outcome: the row contains a study family linked to ours by where its receptors "
    "end in our data (Liu et al. C7, C15, C16, C17; Nordstrom 2011 7tm_2_other); these agree by construction.",
    "Source: literature_all_outcomes.tsv, built by oct05_scripts/build_agreement_tables.py.",
]})

def readable(df):
    return df.rename(columns=lambda c: c.replace("_", " "))


study, family = readable(study), readable(family)
study.to_csv(OUT / "Supplementary_Table_literature_agreement_by_study.tsv", sep="\t", index=False)
family.to_csv(OUT / "Supplementary_Table_literature_agreement_by_family.tsv", sep="\t", index=False)
with pd.ExcelWriter(OUT / "Supplementary_Table_literature_agreement.xlsx") as xw:
    study.to_excel(xw, sheet_name="by study", index=False)
    family.to_excel(xw, sheet_name="by family", index=False)
    notes.to_excel(xw, sheet_name="notes", index=False)
print(study.to_string(index=False))
print(family.to_string(index=False))
