# -*- coding: utf-8 -*-
"""Select representatives per family by taxonomic diversity, preferring the sequences
that align best against their family's topology consensus.

This replaces the mmseqs2 clustering that `select_representatives.py` uses, and changes
nothing else about how the taxonomic spread is built.

Why the clustering went. Clustering supplied two things: a candidate pool of cluster
representatives, and `cluster_size` as the measure of which candidate to prefer. A
`min_cluster_size` of 2 then discarded every singleton cluster, which is 63% of all
cluster representatives across the superfamily. Singletons are by construction the most
isolated sequences at 70% identity, so the filter removed lineages in rough proportion
to how divergent they were. In class A it left 39 of 241 fungal cluster representatives,
erasing eight fungal classes and seventeen fungal orders from the candidate pool
outright, among them Spizellomycetales with thirteen sequences in thirteen distinct
clusters. Those lineages are what a representative set exists to carry.

What replaces it. The candidate pool is the family's complete final membership, and the
preference among candidates is the topology helix distance score from
`03_topology_filter/score_all_family_topologies.py`, the same score the
topology viewer ranks by. No minimum is applied. The sequences that align badly enough
to be worth removing were already removed upstream by the topology filter, whose
exclusions are what produced this membership in the first place, so a second cutoff here
would only shave the surviving tail, and the tail it would shave is the divergent one.
A lineage with exactly one candidate therefore contributes that candidate whatever it
scores.

What is unchanged. `RANKS_PRIORITY`, the choice of the primary rank as the first rank
holding more than one known value, the round robin across primary bins taking the least
filled first, the greedy walk that opens a new deeper rank prefix inside a bin, the
`max_per_species` cap in both the diversification pass and the top-up, the per-family
targets in `TARGETS`, and the three pooled output groups. Wherever the old code sorted
by `cluster_size` descending, this one sorts by `score` descending.

The inner loop is written against numpy arrays with a monotone pointer per bin and
depth, because the old row-by-row `DataFrame.apply` inside the round robin is not
workable once class A's pool is 209,604 rows rather than 14,147. The selection order is
the same; only its bookkeeping differs. A row is skipped forever once it is picked, once
its species reaches the cap, or once its prefix at that depth has been opened, and all
three are one-way, so a pointer that only moves forward cannot miss a candidate.

Usage:
    python select_representatives_by_topology.py

    python select_representatives_by_topology.py \
        --accession-tsv ../step5_family_assignment/family_tables_final/acc_taxid_map_final.tsv \
        --score-tsv ../sequence_structure_conversion/topology_helix_scores.tsv \
        --taxonomy-json ../../zenodo_data/families/reference_data/taxonomy_eukaryotes_filtered.json \
        --output-dir . --tag Sep21
"""
from __future__ import annotations
import os
import sys
import json
import time
import argparse
import logging
from collections import OrderedDict, defaultdict
from typing import Dict, Iterable, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gpcr_paths import stage  # noqa: E402

logging.basicConfig(
    level=os.environ.get("GPCR_LOG_LEVEL", "INFO").upper(),
    format="%(asctime)s | %(levelname)-7s | %(message)s",
)

UNK = "__unk__"


# ---------- 0) Canonicalization ----------
def _canonicalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Lowercase column names and coalesce any that collide once lowercased, taking the
    first non-empty value left to right. The taxonomy JSON carries both 'Species' and
    'species', which is exactly this case."""
    g = df.copy()
    g.columns = [str(c).strip().lower() for c in g.columns]
    lower = list(g.columns)
    out = OrderedDict()
    for name in list(OrderedDict.fromkeys(lower)):
        idxs = [i for i, c in enumerate(lower) if c == name]
        block = g.iloc[:, idxs]
        if block.shape[1] == 1:
            out[name] = block.iloc[:, 0]
        else:
            out[name] = block.replace("", np.nan).bfill(axis=1).iloc[:, 0]
    return pd.DataFrame(out, index=g.index)


def _normalize_str(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip()


def _norm_rank(s: pd.Series) -> pd.Series:
    return (s.astype(str).str.strip()
             .replace({"": UNK, "na": UNK, "nan": UNK, "none": UNK, "NA": UNK, "None": UNK})
             .fillna(UNK))


def _require_file(path: str, kind: str) -> None:
    if not os.path.exists(path):
        raise FileNotFoundError(f"{kind}: file not found -> {path}")


# ---------- 1) Inputs ----------
def load_membership(tsv_path: str) -> pd.DataFrame:
    """The family membership, which is also the candidate pool: every accession that
    survived the 7TM filter, the topology filter and overlap resolution."""
    _require_file(tsv_path, "accession_tsv")
    df = _canonicalize_columns(pd.read_csv(tsv_path, sep="\t", dtype=str, keep_default_na=False))
    for c in ["family", "accession", "taxid"]:
        if c not in df.columns:
            raise ValueError(f"accession_tsv: missing column '{c}'; got {list(df.columns)}")
        df[c] = _normalize_str(df[c])
    df = df.loc[:, ["family", "accession", "taxid"]].drop_duplicates()
    df = df[(df["accession"] != "") & (df["taxid"] != "")]
    logging.info("membership: %d accessions across %d families",
                 len(df), df["family"].nunique())
    return df.rename(columns={"family": "gpcrfamily"})


def load_scores(tsv_path: str, rank_by: str = "mean_jaccard") -> pd.DataFrame:
    """The topology score, one row per family and accession.

    `rank_by` names the score_tsv column the selection ranks on, written to `score`.
    The default is `mean_jaccard`, the topology score used everywhere else in the
    manuscript; `score`, the Euclidean norm of the per-helix Jaccards, is the legacy
    ranking that produced the Sep18 to Sep22 sets. `mean_jaccard` is also carried
    through on its own because `--min-jaccard` filters on it."""
    _require_file(tsv_path, "score_tsv")
    df = _canonicalize_columns(pd.read_csv(tsv_path, sep="\t", dtype=str, keep_default_na=False))
    for c in ["family", "accession", rank_by]:
        if c not in df.columns:
            raise ValueError(f"score_tsv: missing column '{c}'; got {list(df.columns)}")
    out = pd.DataFrame({
        "gpcrfamily": _normalize_str(df["family"]),
        "accession": _normalize_str(df["accession"]),
        "score": pd.to_numeric(df[rank_by], errors="coerce"),
    })
    if "mean_jaccard" in df.columns:
        out["mean_jaccard"] = pd.to_numeric(df["mean_jaccard"], errors="coerce")
    if "coverage" in df.columns:
        out["coverage"] = pd.to_numeric(df["coverage"], errors="coerce")
    out = out.dropna(subset=["score"]).drop_duplicates(subset=["gpcrfamily", "accession"])
    logging.info("scores: %d rows across %d families", len(out), out["gpcrfamily"].nunique())
    return out


def load_taxonomy(json_path: str, ranks: Iterable[str]) -> pd.DataFrame:
    _require_file(json_path, "taxonomy_json")
    with open(json_path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    tax = _canonicalize_columns(pd.DataFrame(data))
    if "family" in tax.columns:
        tax = tax.rename(columns={"family": "family_tax"})
    if "taxid" not in tax.columns:
        raise ValueError("taxonomy_json: missing 'taxid' column")
    tax["taxid"] = _normalize_str(tax["taxid"])
    keep = ["taxid"] + [r for r in ranks if r in tax.columns]
    missing = [r for r in ranks if r not in tax.columns]
    if missing:
        logging.warning("taxonomy JSON has no column for %s; those ranks are skipped", missing)
    tax = tax.loc[:, keep].drop_duplicates(subset=["taxid"])
    for r in keep[1:]:
        tax[r] = _norm_rank(tax[r])
    logging.info("taxonomy: %d taxids, ranks %s", len(tax), keep[1:])
    return tax


# ---------- 2) Selection ----------
def _first_diverse_rank(df: pd.DataFrame, ranks: List[str]) -> Optional[str]:
    for r in ranks:
        if r in df.columns:
            known = df[r][df[r] != UNK]
            if known.nunique() > 1:
                return r
    return None


def hierarchical_diverse_topk_by_score(
    family_df: pd.DataFrame,
    k: int,
    ranks_priority: Iterable[str],
    max_per_species: Optional[int] = 2,
    species_rank: str = "species",
) -> pd.DataFrame:
    """Up to k rows maximizing taxonomic spread, preferring high `score` within each bin
    and allowing no more than `max_per_species` picks to share one known species.

    The primary rank is the first in `ranks_priority` with more than one known value.
    Bins of that rank are visited least-filled first; inside a bin the best-scoring
    candidate that opens a rank prefix not yet used in that bin is taken, trying
    shallower prefixes before deeper ones, and falling back to the best-scoring
    candidate left in the bin when every prefix is already open.
    """
    fam = str(family_df["gpcrfamily"].iloc[0]) if len(family_df) else "?"
    if k <= 0 or family_df.empty:
        logging.warning("[%s] k<=0 or empty candidate pool; nothing selected.", fam)
        return family_df.head(0)

    ranks = [r.lower() for r in ranks_priority]
    species_rank = species_rank.lower()
    df = family_df.copy()
    for r in ranks:
        if r in df.columns:
            df[r] = _norm_rank(df[r])
    df["score"] = pd.to_numeric(df["score"], errors="coerce").fillna(-np.inf)

    # The one sort the whole selection reads. Ties break on accession so a rerun on the
    # same inputs gives the same set.
    df = df.sort_values(["score", "accession"], ascending=[False, True]).reset_index(drop=True)
    n = len(df)

    primary = _first_diverse_rank(df, ranks)
    species = (df[species_rank].to_numpy() if species_rank in df.columns
               else np.full(n, UNK, dtype=object))

    picked = np.zeros(n, dtype=bool)
    species_count: Dict[str, int] = defaultdict(int)
    cap = max_per_species if max_per_species else None

    def eligible(i: int) -> bool:
        if picked[i]:
            return False
        if cap is None:
            return True
        sp = species[i]
        # An unknown species is never capped; there is nothing to diversify against.
        return sp == UNK or species_count[sp] < cap

    def take(i: int) -> None:
        picked[i] = True
        if cap is not None and species[i] != UNK:
            species_count[species[i]] += 1

    if primary is None:
        logging.info("[%s] no rank has more than one known value; taking the best scoring "
                     "under the species cap.", fam)
        chosen = []
        for i in range(n):
            if len(chosen) >= k:
                break
            if eligible(i):
                take(i)
                chosen.append(i)
        return df.iloc[chosen].reset_index(drop=True)

    ranks_ordered = [primary] + [r for r in ranks if r != primary and r in df.columns]
    max_depth = len(ranks_ordered)
    rank_cols = [df[r].to_numpy() for r in ranks_ordered]

    def prefix(i: int, depth: int) -> Tuple:
        return tuple(rank_cols[d][i] for d in range(depth))

    primary_values = df[primary].to_numpy()
    bins: Dict[str, np.ndarray] = {}
    for pv in sorted(set(primary_values.tolist())):
        bins[pv] = np.flatnonzero(primary_values == pv)

    counts_by_primary: Dict[str, int] = {pv: 0 for pv in bins}
    used_prefixes: Dict[str, List[Set[Tuple]]] = {
        pv: [set() for _ in range(max_depth + 1)] for pv in bins
    }
    # Monotone scan pointers. Everything a pointer steps over is disqualified for good,
    # so neither pointer ever needs to move back.
    ptr_depth: Dict[str, List[int]] = {pv: [0] * (max_depth + 1) for pv in bins}
    ptr_best: Dict[str, int] = {pv: 0 for pv in bins}

    def choose(pv: str) -> Optional[int]:
        rows = bins[pv]
        for depth in range(2, max_depth + 1):
            p = ptr_depth[pv][depth]
            seen = used_prefixes[pv][depth]
            while p < rows.size:
                i = int(rows[p])
                if not eligible(i) or prefix(i, depth) in seen:
                    p += 1
                    continue
                break
            ptr_depth[pv][depth] = p
            if p < rows.size:
                return int(rows[p])
        p = ptr_best[pv]
        while p < rows.size and not eligible(int(rows[p])):
            p += 1
        ptr_best[pv] = p
        return int(rows[p]) if p < rows.size else None

    chosen: List[int] = []
    progressing = True
    while len(chosen) < k and progressing:
        progressing = False
        for pv in sorted(bins, key=lambda b: (counts_by_primary[b], b)):
            if len(chosen) >= k:
                break
            i = choose(pv)
            if i is None:
                continue
            take(i)
            chosen.append(i)
            counts_by_primary[pv] += 1
            for depth in range(1, max_depth + 1):
                used_prefixes[pv][depth].add(prefix(i, depth))
            progressing = True

    if len(chosen) < k:
        for i in range(n):
            if len(chosen) >= k:
                break
            if eligible(i):
                take(i)
                chosen.append(i)

    out = df.iloc[chosen].reset_index(drop=True)
    logging.info("[%s] selected %d/%d from %d candidates (balanced by %s; max %s per species; "
                 "score %.3f-%.3f)", fam, len(out), k, n, primary, max_per_species,
                 out["score"].min(), out["score"].max())
    return out


# ---------- 3) Reporting ----------
def print_kingdom_summary(df: pd.DataFrame, label: str) -> None:
    if "kingdom" not in df.columns:
        return
    print(f"=== {label}: kingdoms per family ===")
    counts = df.groupby(["gpcrfamily", "kingdom"]).size().rename("count").reset_index()
    for fam, grp in counts.groupby("gpcrfamily", sort=True):
        grp = grp.sort_values("kingdom")
        distinct = grp[grp["kingdom"] != UNK]["kingdom"].nunique()
        parts = [f"{k}:{int(c)}" for k, c in zip(grp["kingdom"], grp["count"])]
        print(f"{fam}\t{distinct}\t[{', '.join(parts)}]")


def print_species_cap_summary(df: pd.DataFrame, species_rank: str = "species") -> None:
    if species_rank not in df.columns:
        return
    print(f"=== species with more than one representative (rank={species_rank}) ===")
    for fam, grp in df.groupby("gpcrfamily", sort=True):
        counts = grp[grp[species_rank] != UNK][species_rank].value_counts()
        repeated = counts[counts > 1]
        if not repeated.empty:
            print(f"{fam}\t{', '.join(f'{s}:{int(c)}' for s, c in repeated.items())}")


EXCLUSION_DETAIL_RANKS = ["kingdom", "phylum", "class", "order", "family_tax"]


def lineage_exclusion_report(pool_before: pd.DataFrame, min_jaccard: float,
                             ranks: Iterable[str] = EXCLUSION_DETAIL_RANKS) -> pd.DataFrame:
    """Every lineage, at each rank in `ranks`, that the `--min-jaccard` floor removes
    entirely from a family's candidate pool -- every member of that taxon scored below
    the floor (or had no score at all), so no sequence from it can become a
    representative, however the taxonomic-diversity pass would otherwise have weighted
    it. A lineage that keeps even one eligible member is not reported here; it is still a
    candidate, just a smaller one.
    """
    below = pd.to_numeric(pool_before["mean_jaccard"], errors="coerce") < min_jaccard
    below = below | pool_before["mean_jaccard"].isna()
    rows = []
    for fam, g in pool_before.groupby("gpcrfamily", sort=True):
        for r in [x for x in ranks if x in g.columns]:
            known = g[g[r] != UNK]
            if known.empty:
                continue
            all_taxa = set(known[r])
            eligible_taxa = set(known.loc[~below.loc[known.index], r])
            lost = all_taxa - eligible_taxa
            if not lost:
                continue
            for taxon in sorted(lost):
                sub = known[known[r] == taxon]
                rows.append({
                    "gpcrfamily": fam,
                    "rank": r,
                    "lineage": taxon,
                    "n_sequences_excluded": len(sub),
                    "best_mean_jaccard": round(float(pd.to_numeric(
                        sub["mean_jaccard"], errors="coerce").max()), 5)
                        if sub["mean_jaccard"].notna().any() else float("nan"),
                })
    if not rows:
        return pd.DataFrame(columns=["gpcrfamily", "rank", "lineage",
                                     "n_sequences_excluded", "best_mean_jaccard"])
    out = pd.DataFrame(rows).sort_values(
        ["gpcrfamily", "rank", "n_sequences_excluded"], ascending=[True, True, False])
    return out.reset_index(drop=True)


def print_score_summary(chosen: pd.DataFrame, pool: pd.DataFrame) -> None:
    """The selected scores against the pool they came from, so it is visible whether
    preferring the score actually moved the set."""
    print("=== topology score, selection against candidate pool ===")
    print(f"{'family':<14}{'pool':>9}{'poolMed':>9}{'selMed':>9}{'selMin':>9}{'selMax':>9}")
    pool_med = pool.groupby("gpcrfamily")["score"].median()
    pool_n = pool.groupby("gpcrfamily")["score"].size()
    for fam, grp in chosen.groupby("gpcrfamily", sort=True):
        print(f"{fam:<14}{int(pool_n.get(fam, 0)):>9}{pool_med.get(fam, float('nan')):>9.3f}"
              f"{grp['score'].median():>9.3f}{grp['score'].min():>9.3f}{grp['score'].max():>9.3f}")


# ---------- 4) Orchestration ----------
def select_all_families(
    accession_tsv: str,
    score_tsv: str,
    taxonomy_json: str,
    targets_per_family: Dict[str, int],
    ranks_priority: Iterable[str],
    max_per_species: Optional[int] = 2,
    min_jaccard: Optional[float] = None,
    rank_by: str = "mean_jaccard",
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Returns (chosen, pool, exclusions). `pool` is the candidate pool the selection
    actually ran on -- after the `min_jaccard` floor, when one is given. `exclusions` is
    the lineage-exclusion report for that floor, empty when `min_jaccard` is None."""
    ranks = [r.lower() for r in ranks_priority]

    members = load_membership(accession_tsv)
    scores = load_scores(score_tsv, rank_by)
    tax = load_taxonomy(taxonomy_json, ranks)

    pool = members.merge(scores, on=["gpcrfamily", "accession"], how="left")
    unscored = int(pool["score"].isna().sum())
    if unscored:
        # Worst-ranked rather than disqualified: the pool has already passed the only
        # quality cutoff this pipeline applies, so a missing score must not remove a
        # sequence, only push it behind everything that has one.
        logging.warning("%d members have no topology score; ranking them last", unscored)
        pool["score"] = pool["score"].fillna(-np.inf)

    before = len(pool)
    pool = pool.merge(tax, on="taxid", how="inner")
    if len(pool) < before:
        logging.warning("taxonomy join dropped %d of %d members", before - len(pool), before)

    exclusions = pd.DataFrame()
    if min_jaccard is not None:
        if "mean_jaccard" not in pool.columns:
            raise ValueError("--min-jaccard was given but score_tsv has no 'mean_jaccard' column")
        exclusions = lineage_exclusion_report(pool, min_jaccard)
        eligible = pd.to_numeric(pool["mean_jaccard"], errors="coerce") >= min_jaccard
        n_dropped = int((~eligible).sum())
        logging.info("min_jaccard=%.2f: %d/%d candidates below the floor (or unscored) "
                     "across %d families", min_jaccard, n_dropped, len(pool),
                     exclusions["gpcrfamily"].nunique() if not exclusions.empty else 0)
        pool = pool[eligible].copy()

    known = set(pool["gpcrfamily"].unique())
    for fam in targets_per_family:
        if fam not in known and fam.lower() not in {f.lower() for f in known}:
            logging.warning("TARGETS names '%s', which the membership does not contain", fam)

    selections = []
    for fam, fam_df in pool.groupby("gpcrfamily", sort=True):
        k = int(targets_per_family.get(fam, targets_per_family.get(str(fam).lower(), 0)) or 0)
        if k <= 0:
            logging.info("[%s] no target; skipped.", fam)
            continue
        sel = hierarchical_diverse_topk_by_score(
            family_df=fam_df, k=k, ranks_priority=ranks, max_per_species=max_per_species)
        if sel.empty:
            logging.warning("[%s] selection returned nothing (k=%d).", fam, k)
            continue
        selections.append(sel)

    if not selections:
        return pd.DataFrame(), pool, exclusions

    keep = ["gpcrfamily", "accession", "score", "taxid"]
    if "coverage" in pool.columns:
        keep.append("coverage")
    keep += [r for r in ranks if r in pool.columns]
    if "kingdom" in pool.columns and "kingdom" not in keep:
        keep.append("kingdom")
    result = pd.concat(selections, ignore_index=True)
    result = result.loc[:, [c for c in keep if c in result.columns]]
    # Unknown ranks go out as "NA", which is what the step 8 scripts and the earlier
    # selection table both use. The internal sentinel never reaches the file.
    rank_cols = [c for c in result.columns if c not in ("gpcrfamily", "accession", "score",
                                                        "taxid", "coverage")]
    result[rank_cols] = result[rank_cols].replace(UNK, "NA")
    return result, pool, exclusions


def write_accession_list(df: pd.DataFrame, out_path: str) -> None:
    (df["accession"].dropna().astype(str).str.strip()
       .replace("", pd.NA).dropna().drop_duplicates()
       .to_csv(out_path, index=False, header=False))


# Each group is a family exclusion over the one set of per-family selections, so a
# family's representatives are identical in every group that contains it. The names are
# the ones the step 8 alignments and trees already use.
GROUPS: "OrderedDict[str, Set[str]]" = OrderedDict([
    ("all", set()),
    ("no_ste2", {"ste2"}),
    ("no_ste2ste3classCV2", {"classc", "vomeronasal2", "ste3", "ste2"}),
])


def write_groups(chosen: pd.DataFrame, output_dir: str, tag: str) -> Dict[str, str]:
    paths: Dict[str, str] = {}
    for name, excluded in GROUPS.items():
        subset = chosen[~chosen["gpcrfamily"].str.lower().isin(excluded)]
        if subset.empty:
            logging.warning("[%s] excludes every selected family; nothing written.", name)
            continue
        out_path = os.path.join(output_dir, f"representatives_{name}{tag}.txt")
        write_accession_list(subset, out_path)
        paths[name] = out_path
        logging.info("[%s] %d accessions across %d families -> %s", name,
                     subset["accession"].nunique(), subset["gpcrfamily"].nunique(), out_path)
    return paths


# Unchanged from select_representatives.py.
TARGETS = {
    "camp": 100, "classa": 100, "classb1": 20, "classb2": 20, "classf": 20,
    "classt": 20, "gp143": 25, "mth": 20, "nematode": 20, "olfactory": 20,
    "vomeronasal1": 20, "gpr1": 25, "fslb": 20, "gp157": 25, "vomeronasal2": 20,
    "classc": 50, "ste3": 25, "ste2": 20, "tm116": 20,
}

RANKS_PRIORITY = [
    "kingdom", "phylum", "subphylum", "superclass", "class", "superorder", "order",
    "suborder", "infraorder", "parvorder", "superfamily", "family_tax", "subfamily",
    "genus", "species",
]


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--accession-tsv",
                   default=str(stage("04_overlap_assignment", "acc_taxid_map.tsv")),
                   help="the final membership, which is also the candidate pool")
    p.add_argument("--score-tsv",
                   default=str(stage("03_topology_filter", "topology_helix_scores.tsv")),
                   help="output of score_all_family_topologies.py")
    p.add_argument("--taxonomy-json",
                   default=str(stage("07_representative_selection", "taxonomy_eukaryotes_filtered.json")))
    p.add_argument("--output-dir", default=str(stage("07_representative_selection")))
    p.add_argument("--tag", default="", help="optional filename suffix")
    p.add_argument("--max-per-species", type=int, default=2,
                   help="cap on representatives sharing one species inside one family's "
                        "own selection. 0 disables it.")
    p.add_argument("--min-jaccard", type=float, default=None,
                   help="minimum mean_jaccard a candidate must clear to be eligible as a "
                        "representative (e.g. 0.85 for 85%%). Unset by default, matching "
                        "the no-floor selection.")
    p.add_argument("--rank-by", default="mean_jaccard", choices=["mean_jaccard", "score"],
                   help="score_tsv column candidates are ranked by (default mean_jaccard; "
                        "'score' is the legacy Euclidean ranking of the Sep18-Sep22 sets)")
    return p


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    tag = ("_" + args.tag) if args.tag else ""
    os.makedirs(args.output_dir, exist_ok=True)

    chosen, pool, exclusions = select_all_families(
        accession_tsv=args.accession_tsv,
        score_tsv=args.score_tsv,
        taxonomy_json=args.taxonomy_json,
        targets_per_family=TARGETS,
        ranks_priority=RANKS_PRIORITY,
        max_per_species=(args.max_per_species or None),
        min_jaccard=args.min_jaccard,
        rank_by=args.rank_by,
    )

    if chosen.empty:
        raise SystemExit("No selections produced; see the log for the stage that emptied.")

    print_kingdom_summary(chosen, "post-selection")
    print_species_cap_summary(chosen)
    print_score_summary(chosen, pool)

    if args.min_jaccard is not None and not exclusions.empty:
        print(f"\n=== lineages excluded entirely by --min-jaccard {args.min_jaccard:.2f} ===")
        for r in EXCLUSION_DETAIL_RANKS:
            sub = exclusions[exclusions["rank"] == r]
            if sub.empty:
                continue
            print(f"\n-- rank: {r} --")
            for _, row in sub.iterrows():
                print(f"{row['gpcrfamily']:<14}{row['lineage']:<28}"
                      f"{row['n_sequences_excluded']:>6} seqs excluded"
                      f"  best mean_jaccard {row['best_mean_jaccard']:.3f}")
        excl_out = os.path.join(args.output_dir, f"lineage_exclusions_minjaccard{tag}.tsv")
        exclusions.to_csv(excl_out, sep="\t", index=False)
        logging.info("wrote the lineage-exclusion report -> %s", excl_out)

    full_out = os.path.join(args.output_dir, f"representatives_all_families{tag}.tsv")
    chosen.to_csv(full_out, sep="\t", index=False)
    logging.info("wrote the full per-family selection -> %s", full_out)

    for name, path in write_groups(chosen, args.output_dir, tag).items():
        logging.info("group '%s' -> %s", name, path)
