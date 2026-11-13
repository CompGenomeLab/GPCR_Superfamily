"""Quantify cross-family accession overlap before and after the 7TM filter, then
score the remaining overlaps by median-normalized search score.

The baseline is the Step-1 extraction set, one accession list per family. The 7TM
filter reduces it. Whatever accessions are still claimed by more than one family are
scored against each claiming family's own post-filter member set, so that families of
very different size stay comparable.

E-values are deliberately not used. A HMMER E-value depends on the profile's length
and information content and on the target database size, and the seed searches report
MMseqs2 E-values from a different statistical framework, so the values are not
comparable across families. The score used here never leaves the search that produced
it.

Contested accessions are scored by `median_ratios`: the search bit score divided by
the median bit score of that family's own post-filter members, so the typical member
of every family is 1.0 and higher is a better match. The value is internal to a single
search, so it never compares raw scores across profiles, and unlike an ordinal rank it
distinguishes a decisive win from a narrow one.

No deduplicated table is written. The diagnostics name the family each contested
accession scores best in, but the input membership is left intact and the decision is
taken downstream.
"""

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
from statistics import median

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import exists, open_text, stage  # noqa: E402

# Families whose accessions are ranked against a broader shared search rather than a
# search of their own.
SHARED_SEARCH = {
    "classB1": "classB",
    "classB2": "classB",
    "Vomeronasal1": "classT",
}

# Families whose search is an MMseqs2 seed search, given as the directory holding it.
M8_SEARCH = {
    "GPR1": "GPR1_YEAST",
    "FSLB": "FSLB_DICDI",
}

# Taxonomic fields from broadest to narrowest, with Domain deliberately absent: it is
# "Eukaryota" for every taxon in this study and carries no information. Kingdom is the
# broadest informative rank, and the Clade fields cover the protist supergroups that
# NCBI leaves without a kingdom.
BROADEST_TAXON_FIELDS = (
    ["Kingdom", "Clade"]
    + [f"Clade {index}" for index in range(1, 22)]
    + ["Subkingdom", "Phylum", "Subphylum", "Superclass", "Class", "Superorder",
       "Order", "Family", "Genus", "Species"]
)

SERIES_BEFORE = "#2a78d6"   # categorical slot 1
SERIES_AFTER = "#eb6834"    # categorical slot 2
SURFACE = "#fcfcfb"
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
GRID = "#d9d8d4"


def accession_of(token):
    """tr|A0A8H5TWS9|A0A8H5TWS9_FUSHE -> A0A8H5TWS9; anything else passes through."""
    fields = token.split("|")
    return fields[1].strip() if len(fields) >= 3 else token.strip()


def read_wide_table(path):
    """Read a family-per-column accession table into {family: [accessions]}."""
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    header = rows[0]
    table = {name: [] for name in header}
    for row in rows[1:]:
        for name, value in zip(header, row):
            value = value.strip()
            if value:
                table[name].append(value)
    return table


def read_baseline(before_tsv, families):
    """The clade membership (stage 02) of each family in the post-filter table."""
    table = read_wide_table(before_tsv)
    return {family: set(table.get(family, ())) for family in families}


def hmm_search_hits(path):
    """(accession, score) in the order hmmsearch reported them, best-scoring first.

    The score is the full-sequence bit score, the second column of the per-sequence
    table, not the best-single-domain score and not the E-value.
    """
    hits = []
    seen = set()
    in_scores = header_seen = reading = False
    with open_text(path) as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("Domain annotation for each sequence:"):
                break
            if not in_scores:
                if line.startswith("Scores for complete sequences"):
                    in_scores = True
                continue
            if not header_seen:
                if line.lstrip().startswith("E-value"):
                    header_seen = True
                continue
            if not reading:
                if set(line.strip()) <= {"-", " "} and "-" in line:
                    reading = True
                continue
            parts = line.split()
            if len(parts) < 9:
                continue
            accession = accession_of(parts[8])
            if accession in seen:
                continue
            try:
                score = float(parts[1])
            except ValueError:
                continue
            seen.add(accession)
            hits.append((accession, score))
    return hits


def m8_search_hits(path):
    """(accession, score) in the order the MMseqs2 result file lists them.

    The layout is query, target, theader, pident, fident, nident, alnlen, qlen, tlen,
    mismatch, raw, bits, qcov, tcov, evalue, so the bit score is column 11. Using bits
    rather than the raw score keeps the value on the same footing as HMMER's.
    """
    hits = []
    seen = set()
    with open_text(path) as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            accession = accession_of(parts[1])
            if not accession or accession in seen:
                continue
            try:
                score = float(parts[11])
            except ValueError:
                continue
            seen.add(accession)
            hits.append((accession, score))
    return hits


def search_hits_for(family, search_dir):
    """(accession, score) for the search this family's members were found by."""
    search_dir = Path(search_dir)
    if family in M8_SEARCH:
        folder = search_dir / M8_SEARCH[family]
        paths = sorted(list(folder.glob("*.m8")) + list(folder.glob("*.m8.gz")))
        if not paths:
            sys.exit(f"no .m8 file for {family} in {folder}")
        return m8_search_hits(paths[0])
    search = SHARED_SEARCH.get(family, family)
    path = search_dir / search / f"{search}_hmmsearch_E1000.txt"
    if not exists(path):
        sys.exit(f"missing search output for {family}: {path}")
    return hmm_search_hits(path)


def median_ratios(members, scores):
    """Divide each member's search score by the median score of that family's members.

    Dividing by the family's own median is what makes families of very different size
    comparable, because the typical member of every family is 1.0 by construction, so
    a small family is not penalised against a large one. Unlike an ordinal rank it
    also preserves how far above or below typical a sequence actually scores, which
    is what separates a decisive win from a hair's-breadth one.

    Higher is better.
    """
    present = {a: scores[a] for a in members if a in scores}
    missing = [a for a in members if a not in scores]
    if not present:
        return {}, missing
    middle = median(present.values())
    if middle <= 0:
        return {a: 0.0 for a in present}, missing
    return {a: value / middle for a, value in present.items()}, missing


def pairwise_overlaps(table):
    """Count accessions shared by each family pair."""
    counts = Counter()
    for family_a, family_b in combinations(sorted(table), 2):
        shared = len(set(table[family_a]) & set(table[family_b]))
        if shared:
            counts[(family_a, family_b)] = shared
    return counts


def load_taxonomy(path):
    """taxID -> (broadest informative taxon, species).

    Domain is skipped because it is "Eukaryota" throughout, and Kingdom is empty for
    the protist supergroups, so the walk falls through to the Clade fields rather than
    leaving the column blank.
    """
    with open(path, encoding="utf-8") as handle:
        records = json.load(handle)
    lookup = {}
    for record in records:
        broad = "unassigned"
        for field in BROADEST_TAXON_FIELDS:
            value = str(record.get(field, "NA")).strip()
            if value and value != "NA":
                broad = value
                break
        lookup[str(record.get("taxID"))] = (broad, str(record.get("Species", "NA")))
    return lookup


def read_acc_taxid(path):
    """accession -> taxID, from a long-format mapping table."""
    mapping = {}
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            mapping.setdefault(row["accession"], row["taxID"])
    return mapping


def write_overlap_table(path, before, after, baseline=None, reinstated=None):
    """Write one row per family pair.

    `from_lifted_curation` counts how many of the pair's before-overlap accessions
    exist only because a Step-1 exclusion was lifted, so a pair that the curation had
    zeroed is visible as such rather than being confused with a pair that merely
    involves a family which gained sequences elsewhere.
    """
    baseline = baseline or {}
    reinstated = reinstated or {}
    pairs = sorted(set(before) | set(after),
                   key=lambda pair: (-before.get(pair, 0), -after.get(pair, 0), pair))
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["family_a", "family_b", "overlap_before", "overlap_after",
                         "removed_by_7TM", "from_lifted_curation"])
        for family_a, family_b in pairs:
            n_before = before.get((family_a, family_b), 0)
            n_after = after.get((family_a, family_b), 0)
            lifted = len((reinstated.get(family_a, set()) & baseline.get(family_b, set()))
                         | (reinstated.get(family_b, set()) & baseline.get(family_a, set())))
            writer.writerow([family_a, family_b, n_before, n_after,
                             n_before - n_after, lifted])
    return pairs


def plot_overlaps(path, pairs, before, after, top, lifted_counts=None):
    """Interactive grouped bar chart, one pair of bars per family pair."""
    import plotly.graph_objects as go

    shown = pairs[:top]
    labels = [f"{a} / {b}" for a, b in shown]
    before_values = [before.get(pair, 0) for pair in shown]
    after_values = [after.get(pair, 0) for pair in shown]
    removed = [b - a for b, a in zip(before_values, after_values)]
    lifted = [lifted_counts.get(pair, 0) for pair in shown] if lifted_counts else [0] * len(shown)
    kept = [f"{100 * a / b:.1f}%" if b else "n/a"
            for a, b in zip(after_values, before_values)]

    hover = ("<b>%{y}</b><br>%{fullData.name}<br>"
             "shared accessions: %{x:,}<br>"
             "removed by 7TM: %{customdata[0]:,}<br>"
             "surviving: %{customdata[1]}<br>"
             "from lifted curation: %{customdata[2]:,}<extra></extra>")
    custom = list(zip(removed, kept, lifted))

    figure = go.Figure([
        go.Bar(y=labels, x=before_values, name="Before 7TM filter", orientation="h",
               marker_color=SERIES_BEFORE, customdata=custom, hovertemplate=hover),
        go.Bar(y=labels, x=after_values, name="After 7TM filter", orientation="h",
               marker_color=SERIES_AFTER, customdata=custom, hovertemplate=hover),
    ])
    figure.update_layout(
        barmode="group",
        bargap=0.28,
        bargroupgap=0.08,
        title=dict(text="Cross-family accession overlap before and after the 7TM filter",
                   x=0, xanchor="left", font=dict(size=16, color=INK_PRIMARY)),
        xaxis=dict(title="Accessions shared by the pair", gridcolor=GRID,
                   zerolinecolor=GRID, tickfont=dict(color=INK_SECONDARY)),
        yaxis=dict(autorange="reversed", tickfont=dict(color=INK_PRIMARY)),
        legend=dict(orientation="h", yanchor="bottom", y=1.01, xanchor="left", x=0,
                    font=dict(color=INK_PRIMARY)),
        plot_bgcolor=SURFACE,
        paper_bgcolor=SURFACE,
        font=dict(family="system-ui, -apple-system, Segoe UI, sans-serif", size=12),
        height=max(420, 26 * len(shown) + 180),
        margin=dict(l=10, r=40, t=90, b=60),
        hoverlabel=dict(font_size=12),
    )
    figure.write_html(path, include_plotlyjs="cdn")


def main():
    parser = argparse.ArgumentParser(
        description="Overlap before/after the 7TM filter, and median-normalized scoring "
                    "of the accessions more than one family claims")
    parser.add_argument("--after-tsv", default=str(stage("03_topology_filter", "topology_membership.tsv")),
                        help="wide family-by-accession table after the topology filter")
    parser.add_argument("--before-tsv", default=str(stage("02_clade_extraction", "clade_membership.tsv")),
                        help="wide family-by-accession table of the clade membership")
    parser.add_argument("--search-dir", default=str(stage("01_homology_search")),
                        help="the search outputs, one folder per search")
    parser.add_argument("--acc-taxid", default=None,
                        help="optional accession/taxID table, to label each contested sequence")
    parser.add_argument("--taxonomy", default=None,
                        help="optional taxonomy JSON (06_family_origins/taxonomy/"
                             "taxonomy_ncbi_searched.json), to label each contested sequence")
    parser.add_argument("--output-dir", default=str(stage("04_overlap_assignment")))
    parser.add_argument("--top-pairs", type=int, default=40,
                        help="family pairs drawn in the bar plot")
    parser.add_argument("--baseline-json", default=None,
                        help="JSON of {family: [accessions]} built with the Step-1 "
                             "EXCLUDE_ITEMS lifted. When given, the accessions the "
                             "exclusion had removed are reinstated on both sides so "
                             "those pairs are measured rather than zeroed by fiat.")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    after = read_wide_table(args.after_tsv)
    families = sorted(after)
    before = read_baseline(args.before_tsv, families)

    reinstated = {family: set() for family in families}
    if args.baseline_json:
        with open(args.baseline_json, encoding="utf-8") as handle:
            unexcluded = json.load(handle)
        # A sequence's 7TM verdict does not depend on which family claimed it, so an
        # accession the exclusion removed from one family can be reinstated there if
        # any family's DeepTMHMM run passed it.
        passed_7tm = set().union(*(set(v) for v in after.values()))
        for family in families:
            added = set(unexcluded.get(family, ())) - before[family]
            if added:
                reinstated[family] = added
                before[family] = before[family] | added
                after[family] = sorted(set(after[family]) | (added & passed_7tm))
                print(f"  reinstated into {family}: {len(added):,} baseline, "
                      f"{len(added & passed_7tm):,} of them 7TM")

    print(f"families: {len(families)}")
    print(f"baseline accessions: {sum(len(v) for v in before.values()):,}")
    print(f"after 7TM:           {sum(len(v) for v in after.values()):,}")

    overlaps_before = pairwise_overlaps(before)
    overlaps_after = pairwise_overlaps(after)
    pairs = write_overlap_table(output_dir / "overlap_counts_before_after.tsv",
                                overlaps_before, overlaps_after, before, reinstated)
    lifted_counts = {
        (a, b): len((reinstated.get(a, set()) & before.get(b, set()))
                    | (reinstated.get(b, set()) & before.get(a, set())))
        for a, b in pairs}
    print(f"\noverlapping pairs: {len(overlaps_before)} before, "
          f"{len(overlaps_after)} after")
    print(f"overlapping assignments: {sum(overlaps_before.values()):,} before, "
          f"{sum(overlaps_after.values()):,} after")

    # Median-normalized score of every accession inside its own family's post-filter
    # member set.
    ratios = {}
    for family in families:
        scores = dict(search_hits_for(family, args.search_dir))
        ratios[family], missing = median_ratios(after[family], scores)
        note = f", {len(missing)} not found in search" if missing else ""
        print(f"  scored {family}: {len(ratios[family]):,} of {len(after[family]):,}{note}")

    claimed = defaultdict(list)
    for family in families:
        for accession in after[family]:
            claimed[accession].append(family)
    contested = {a: f for a, f in claimed.items() if len(f) > 1}
    print(f"\ncontested accessions: {len(contested):,}")

    acc_taxid = read_acc_taxid(args.acc_taxid) if args.acc_taxid else {}
    taxonomy = load_taxonomy(args.taxonomy) if args.taxonomy else {}

    rows = []
    for accession, family_list in sorted(contested.items()):
        candidates = [f for f in family_list if accession in ratios[f]]
        if not candidates:
            continue
        ordered = sorted(candidates, key=lambda f: -ratios[f][accession])
        taxid = acc_taxid.get(accession, "NA")
        broadest, species = taxonomy.get(taxid, ("unassigned", "NA"))
        # Families and ratios run in the same order, highest ratio first, so the two
        # columns read as parallel lists and the assignment is the leading entry.
        rows.append({
            "accession": accession,
            "n_families": len(family_list),
            "families": ",".join(ordered),
            "median_ratios": ",".join(f"{ratios[f][accession]:.4f}" for f in ordered),
            "assigned_family": ordered[0],
            "taxID": taxid,
            "broadest_taxon": broadest,
            "species": species,
        })

    diagnostics = output_dir / "overlap_diagnostics.tsv"
    with diagnostics.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    plot_overlaps(output_dir / "overlap_before_after.html", pairs,
                  overlaps_before, overlaps_after, args.top_pairs, lifted_counts)

    # No deduplicated table is written. The contested accessions are reported in the
    # diagnostics with the family their median-normalized score favours, and the
    # decision about which family actually keeps them is made downstream rather than
    # applied here, so no membership is silently dropped from the input table.
    print(f"contested accessions scored: {len(rows):,}")
    print(f"\nwrote {output_dir}")


if __name__ == "__main__":
    main()
