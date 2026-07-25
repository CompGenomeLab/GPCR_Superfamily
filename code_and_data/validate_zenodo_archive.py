#!/usr/bin/env python3
"""Validate the filenames and pipeline handoffs in an unpacked Zenodo archive."""

from __future__ import annotations

import gzip
import csv
import filecmp
import os
import re
import sys
from pathlib import Path


HMM_SEARCHES = (
    "STE3", "STE2", "Vomeronasal1", "Vomeronasal2", "Nematode", "Mth",
    "classB", "classC", "classF", "classT", "classA", "cAMP", "GP143",
    "Galpha", "Ggamma", "Gbeta", "Olfactory", "GP157", "TBP", "RBP1",
)
SEED_SEARCHES = ("GPR1_YEAST", "FSLB_DICDI", "LRLA_DICDI")
FAMILIES = (
    "Olfactory", "Vomeronasal1", "Vomeronasal2", "Nematode", "STE2",
    "STE3", "Mth", "GP143", "classF", "classC", "classT", "classB1",
    "classB2", "GPR1", "FSLB", "classA", "cAMP", "GP157",
)
TRIMMED_FAMILIES = tuple(family for family in FAMILIES if family != "STE2")
SUPERFAMILY_RUNS = (
    "All_woC-Ste3_clust0.7_cov1_Oct24_minsize2",
    "All_clust0.7_cov1_Oct25_minsize2",
)
HMM_REGION_SELECTIONS = {
    "STE3", "STE2", "Vomeronasal2", "classB", "classC", "classA",
    "cAMP", "Ggamma", "GP157", "TBP",
}


def fasta_accessions(path: Path) -> set[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    accessions = set()
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            identifier = line[1:].split(None, 1)[0]
            parts = identifier.split("|")
            accessions.add(parts[1] if len(parts) >= 3 else identifier)
    return accessions


def text_ids(path: Path) -> set[str]:
    with path.open(encoding="utf-8", errors="replace") as handle:
        return {line.strip().split()[0] for line in handle if line.strip()}


def fasta_identifiers(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        return [line[1:].split(None, 1)[0] for line in handle if line.startswith(">")]


def newick_leaf_identifiers(path: Path) -> set[str]:
    text = path.read_text(encoding="utf-8", errors="replace")
    return {
        match.group(1).strip().strip("'\"")
        for match in re.finditer(r"(?<=[(,])\s*([^():,;]+?)\s*(?=:)", text)
    }


def file_or_gzip(path: Path) -> bool:
    """Return True when a file is retained either directly or gzip-compressed."""
    return resolve_file_or_gzip(path) is not None


def resolve_file_or_gzip(path: Path) -> Path | None:
    """Return the retained plain or gzip-compressed path, when present."""
    if path.is_file():
        return path
    compressed = Path(f"{path}.gz")
    return compressed if compressed.is_file() else None


def read_count_table(path: Path) -> tuple[list[str], dict[str, tuple[str, ...]]]:
    """Read a small tax-count table into taxID-keyed immutable rows."""
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = reader.fieldnames or []
        return columns, {
            row["taxID"]: tuple(row[column] for column in columns)
            for row in reader
        }


def main() -> int:
    repository = Path(__file__).resolve().parents[1]
    archive = Path(os.environ.get("GPCR_ZENODO_DIR", repository.parent / "zenodo_data"))
    missing: list[str] = []
    mismatched: list[str] = []
    warnings: list[str] = []

    expected_top = {"search_results", "extracted_clades", "families"}
    actual_top = {path.name for path in archive.iterdir() if path.is_dir()} if archive.is_dir() else set()
    if actual_top != expected_top:
        mismatched.append(
            f"top-level directories: expected {sorted(expected_top)}, found {sorted(actual_top)}"
        )

    search_root = archive / "search_results"
    expected_search_dirs = set(HMM_SEARCHES + SEED_SEARCHES)
    actual_search_dirs = {path.name for path in search_root.iterdir() if path.is_dir()} if search_root.is_dir() else set()
    for name in sorted(expected_search_dirs - actual_search_dirs):
        missing.append(f"search_results/{name}/")
    for name in sorted(actual_search_dirs - expected_search_dirs):
        mismatched.append(f"unexpected search directory: search_results/{name}/")

    for family in HMM_SEARCHES:
        folder = search_root / family
        stem = f"{family}_hmm_E1000_db"
        required = (
            f"{family}_hmmsearch_E1000.txt",
            f"{family}_hmmsearch_E1000_headers.txt",
            f"{stem}.fasta",
            f"{stem}_FAMSA_gt0.1.fasta",
            f"{stem}_FAMSA_gt0.1.tree",
            f"{stem}_FAMSA_gt0.1_hmmordered.tree",
            f"{stem}_FAMSA_gt0.1_treeordered.fasta",
            f"{stem}_FAMSA_gt0.1_hmmordered_score_lineplot.html",
            f"{stem}_hmmregions.fasta",
            f"{stem}_hmmregions_FAMSA_gt0.1.fasta",
            f"{stem}_hmmregions_FAMSA_gt0.1.tree",
            f"{stem}_hmmregions_FAMSA_gt0.1_hmmordered.tree",
            f"{stem}_hmmregions_FAMSA_gt0.1_treeordered.fasta",
            f"{stem}_hmmregions_FAMSA_gt0.1_hmmordered_score_lineplot.html",
        )
        for name in required:
            path = folder / name
            present = file_or_gzip(path) if name.endswith((".fasta", ".txt")) else path.is_file()
            if not present:
                missing.append(f"search_results/{family}/{name}")
        for initial_alignment in (
            f"{stem}_FAMSA.fasta",
            f"{stem}_hmmregions_FAMSA.fasta",
        ):
            if not file_or_gzip(folder / initial_alignment):
                missing.append(f"search_results/{family}/{initial_alignment}[.gz]")

    for family in SEED_SEARCHES:
        folder = search_root / family
        result = f"{family}_UniProt_s7.5_filter0_eprofile0.1_db"
        stem = f"{family}_db_FAMSA_gt0.1"
        for name in (
            f"{result}.m8", f"{result}_headers10000.txt", f"{family}_db.fasta",
            f"{stem}.fasta", f"{stem}.tree", f"{stem}_m8ordered.tree",
            f"{stem}_treeordered.fasta", f"{stem}_m8ordered_score_lineplot.html",
        ):
            path = folder / name
            present = file_or_gzip(path) if name.endswith((".fasta", ".txt")) else path.is_file()
            if not present:
                missing.append(f"search_results/{family}/{name}")

    family_root = archive / "families"
    expected_family_dirs = set(FAMILIES) | {
        "accessions", "clusters", "tables",
        "reference_data", "superfamily", "trimmed_msas",
    }
    actual_family_dirs = {path.name for path in family_root.iterdir() if path.is_dir()} if family_root.is_dir() else set()
    for name in sorted(actual_family_dirs - expected_family_dirs):
        mismatched.append(f"unexpected family directory: families/{name}/")

    expected_clade_alignments = {
        f"{family}_hmm_E1000_db_FAMSA_gt0.1_treeordered_extracted_clade.fasta"
        for family in HMM_SEARCHES
    }
    expected_clade_alignments.remove("classB_hmm_E1000_db_FAMSA_gt0.1_treeordered_extracted_clade.fasta")
    expected_clade_alignments.add("classB1_hmm_E1000_db_FAMSA_gt0.1_treeordered_extracted_clade.fasta")
    expected_clade_alignments.add("classB2_hmm_E1000_db_FAMSA_gt0.1_treeordered_extracted_clade.fasta")
    expected_clade_alignments.update(
        f"{family}_db_FAMSA_gt0.1_treeordered_extracted_clade.fasta"
        for family in SEED_SEARCHES
    )
    extracted_root = archive / "extracted_clades"
    alignment_root = extracted_root / "alignments"
    actual_clade_alignments = {
        path.name[:-3] if path.name.endswith(".gz") else path.name
        for path in alignment_root.glob("*_extracted_clade.fasta*")
    }
    for name in sorted(expected_clade_alignments - actual_clade_alignments):
        missing.append(f"extracted_clades/alignments/{name}")
    for name in sorted(actual_clade_alignments - expected_clade_alignments):
        mismatched.append(f"unexpected extracted-clade alignment: extracted_clades/alignments/{name}")
    expected_table_files = {
        "accessions_by_family_prefilter_May22.tsv",
        "acc_taxid_map_unfiltered.tsv",
        "tax_counts_per_family_unfiltered.tsv",
        "acc_taxid_map_e0.1_BestEval.tsv",
        "tax_counts_per_family_e0.1_BestEval.tsv",
        "tax_counts_per_family_e0.1_BestEval_allControls.tsv",
        "overlap_plot_unfiltered.html",
        "overlap_plot_e0.1_BestEval.html",
        "overlap_plot_e0.1_BestEval_allControls.html",
        "missing_acc_taxid_e0.1_BestEval.tsv",
        "acc_taxid_map_BestEval_May22_e0.1.tsv",
        "acc_taxid_map_prefilter_Oct24_e1000_all.tsv",
        "missing_acc_taxid_BestEval_May22_e0.1.tsv",
        "tax_counts_per_family_BestEval_May22_e0.1.tsv",
        "overlap_plot_prefilter_May22.html",
        "tax_counts_per_family_filtered_Oct24_e0.1_BestEval.tsv",
        "tax_counts_per_family_prefilter_Oct24_e0.1_all.tsv",
        "tax_counts_per_family_filtered_Oct24_e0.1_BestEval_withControls_Filtered.tsv",
    }
    table_root = family_root / "tables"
    actual_table_files = {path.name for path in table_root.iterdir() if path.is_file()}
    for name in sorted(actual_table_files - expected_table_files):
        mismatched.append(f"unexpected family-table artifact: families/tables/{name}")

    expected_clade_trees = set()
    for family in HMM_SEARCHES:
        region = "hmmregions_" if family in HMM_REGION_SELECTIONS else ""
        expected_clade_trees.add(
            f"{family}_hmm_E1000_db_{region}FAMSA_gt0.1_hmmordered_extracted_clade.tree"
        )
    expected_clade_trees.remove(
        "classB_hmm_E1000_db_hmmregions_FAMSA_gt0.1_hmmordered_extracted_clade.tree"
    )
    expected_clade_trees.update({
        "classB1_hmm_E1000_db_hmmregions_FAMSA_gt0.1_hmmordered_extracted_clade.tree",
        "classB2_hmm_E1000_db_hmmregions_FAMSA_gt0.1_hmmordered_extracted_clade.tree",
    })
    expected_clade_trees.update(
        f"{family}_db_FAMSA_gt0.1_m8ordered_extracted_clade.tree"
        for family in SEED_SEARCHES
    )
    tree_root = extracted_root / "trees"
    actual_clade_trees = {path.name for path in tree_root.glob("*_extracted_clade.tree")}
    for name in sorted(expected_clade_trees - actual_clade_trees):
        missing.append(f"extracted_clades/trees/{name}")
    for name in sorted(actual_clade_trees - expected_clade_trees):
        mismatched.append(f"unexpected extracted-clade tree: extracted_clades/trees/{name}")

    clade_outputs = []
    for family in HMM_SEARCHES:
        output_family = "classB1" if family == "classB" else family
        region = "hmmregions_" if family in HMM_REGION_SELECTIONS else ""
        clade_outputs.append((
            output_family,
            alignment_root / f"{output_family}_hmm_E1000_db_FAMSA_gt0.1_treeordered_extracted_clade.fasta",
            tree_root / f"{output_family}_hmm_E1000_db_{region}FAMSA_gt0.1_hmmordered_extracted_clade.tree",
        ))
    clade_outputs.append((
        "classB2",
        alignment_root / "classB2_hmm_E1000_db_FAMSA_gt0.1_treeordered_extracted_clade.fasta",
        tree_root / "classB2_hmm_E1000_db_hmmregions_FAMSA_gt0.1_hmmordered_extracted_clade.tree",
    ))
    clade_outputs.extend(
        (
            family,
            alignment_root / f"{family}_db_FAMSA_gt0.1_treeordered_extracted_clade.fasta",
            tree_root / f"{family}_db_FAMSA_gt0.1_m8ordered_extracted_clade.tree",
        )
        for family in SEED_SEARCHES
    )
    for family, fasta_path, tree_path in clade_outputs:
        fasta_path = resolve_file_or_gzip(fasta_path)
        if fasta_path is None or not tree_path.is_file():
            continue
        fasta_ids = fasta_identifiers(fasta_path)
        tree_ids = newick_leaf_identifiers(tree_path)
        if len(fasta_ids) != len(set(fasta_ids)):
            mismatched.append(f"{family}: duplicate identifiers in extracted-clade alignment")
        if set(fasta_ids) != tree_ids:
            mismatched.append(
                f"{family}: extracted FASTA has {len(set(fasta_ids))} IDs and pruned tree "
                f"has {len(tree_ids)} IDs ({len(set(fasta_ids) - tree_ids)} missing from tree, "
                f"{len(tree_ids - set(fasta_ids))} extra in tree)"
            )

    for family in FAMILIES:
        folder = family_root / family
        if not folder.is_dir():
            missing.append(f"families/{family}/")
            continue
        for name in (
            f"{family}_genes_db.fasta",
            f"{family}_genes_db_FAMSA_gt0.1.fasta",
            f"{family}_genes_db_FAMSA_gt0.1_filtered.fasta",
            f"{family}_genes_db_FAMSA_gt0.1_filtered_ids.txt",
            f"{family}_genes_filtered_db.fasta",
            f"{family}_genes_filtered_db_clust0.7_clust0_cov1.tsv",
            f"{family}_genes_filtered_db_clust0.7_clust0_cov1_reps.tsv",
            f"{family}_genes_filtered_reps_db.fasta",
        ):
            path = folder / name
            present = file_or_gzip(path) if name.endswith(".fasta") else path.is_file()
            if not present:
                missing.append(f"families/{family}/{name}")
        aligned = folder / f"{family}_genes_filtered_db_FAMSA.fasta"
        if not file_or_gzip(aligned):
            missing.append(f"families/{family}/{aligned.name}[.gz]")
        initial_family_alignment = folder / f"{family}_genes_db_FAMSA.fasta"
        if not file_or_gzip(initial_family_alignment):
            missing.append(
                f"families/{family}/{initial_family_alignment.name}[.gz]"
            )

        cluster_name = f"{family}_genes_filtered_db_clust0.7_clust0_cov1.tsv"
        family_cluster = folder / cluster_name
        staged_cluster = family_root / "clusters" / cluster_name
        if not staged_cluster.is_file():
            missing.append(f"families/clusters/{cluster_name}")
        elif family_cluster.is_file() and not filecmp.cmp(
            family_cluster, staged_cluster, shallow=False
        ):
            mismatched.append(
                f"{family}: families/clusters/{cluster_name} is not byte-identical "
                f"to families/{family}/{cluster_name}"
            )

        accession_file = family_root / "accessions" / f"{family}_accessions.txt"
        database_fasta = resolve_file_or_gzip(folder / f"{family}_genes_db.fasta")
        if not accession_file.is_file():
            missing.append(f"families/accessions/{family}_accessions.txt")
        elif database_fasta is not None:
            expected = text_ids(accession_file)
            observed = fasta_accessions(database_fasta)
            if expected != observed:
                mismatched.append(
                    f"{family}: accession list has {len(expected)} IDs, database FASTA has "
                    f"{len(observed)} IDs ({len(expected - observed)} missing, "
                    f"{len(observed - expected)} extra)"
                )

    for family in TRIMMED_FAMILIES:
        name = f"{family}_genes_filtered_db_FAMSA.ref_trimmed.fasta"
        if not file_or_gzip(family_root / "trimmed_msas" / name):
            missing.append(f"families/trimmed_msas/{name}")

    for relative in (
        "families/tables/accessions_by_family_prefilter_May22.tsv",
        "families/tables/acc_taxid_map_unfiltered.tsv",
        "families/tables/tax_counts_per_family_unfiltered.tsv",
        "families/tables/acc_taxid_map_e0.1_BestEval.tsv",
        "families/tables/tax_counts_per_family_e0.1_BestEval.tsv",
        "families/tables/tax_counts_per_family_e0.1_BestEval_allControls.tsv",
        "families/tables/overlap_plot_unfiltered.html",
        "families/tables/overlap_plot_e0.1_BestEval.html",
        "families/tables/overlap_plot_e0.1_BestEval_allControls.html",
        "families/tables/acc_taxid_map_BestEval_May22_e0.1.tsv",
        "families/tables/acc_taxid_map_prefilter_Oct24_e1000_all.tsv",
        "families/reference_data/taxonomy_eukaryotes_filtered.json",
    ):
        if not (archive / relative).is_file():
            missing.append(relative)

    canonical_counts = [
        table_root / "tax_counts_per_family_unfiltered.tsv",
        table_root / "tax_counts_per_family_e0.1_BestEval.tsv",
        table_root / "tax_counts_per_family_e0.1_BestEval_allControls.tsv",
    ]
    if all(path.is_file() for path in canonical_counts):
        unfiltered_columns, unfiltered_rows = read_count_table(canonical_counts[0])
        best_columns, best_rows = read_count_table(canonical_counts[1])
        complete_columns, complete_rows = read_count_table(canonical_counts[2])
        if unfiltered_columns != best_columns or best_columns != complete_columns:
            mismatched.append("canonical tax-count tables do not share identical columns")
        if not set(best_rows).issubset(unfiltered_rows):
            mismatched.append("BestEval taxIDs are not a subset of unfiltered taxIDs")
        control_indexes = [
            best_columns.index(name)
            for name in ("Galpha_count", "Gbeta_count", "Ggamma_count", "RBP1_count", "TBP_count")
        ]
        expected_complete = {
            taxid: row
            for taxid, row in best_rows.items()
            if all(int(row[index]) >= 1 for index in control_indexes)
        }
        if complete_rows != expected_complete:
            mismatched.append(
                "all-controls tax-count table is not the exact five-controls-present subset of BestEval"
            )

    superfamily_root = family_root / "superfamily"
    actual_runs = {path.name for path in superfamily_root.iterdir() if path.is_dir()} if superfamily_root.is_dir() else set()
    if actual_runs != set(SUPERFAMILY_RUNS):
        mismatched.append(
            f"superfamily runs: expected {sorted(SUPERFAMILY_RUNS)}, found {sorted(actual_runs)}"
    )
    for run in SUPERFAMILY_RUNS:
        if not list((superfamily_root / run).glob("*trimends*.fasta*")):
            missing.append(f"families/superfamily/{run}/*trimends*.fasta")

    obsolete = [
        str(path.relative_to(archive))
        for path in archive.rglob("*")
        if "first_round" in path.name.lower() or "second_round" in path.name.lower()
    ]
    if obsolete:
        mismatched.extend(f"obsolete round label: {path}" for path in obsolete)

    print(f"Archive root: {archive}")
    print(f"Search folders checked: {len(expected_search_dirs)}")
    print(f"Family datasets checked: {len(FAMILIES)}")
    print(f"Warnings: {len(warnings)}")
    for item in warnings:
        print(f"WARNING: {item}")
    if missing or mismatched:
        print(f"Missing required items: {len(missing)}")
        for item in missing:
            print(f"MISSING: {item}")
        print(f"Naming/content mismatches: {len(mismatched)}")
        for item in mismatched:
            print(f"MISMATCH: {item}")
        return 1

    print("All required pipeline handoffs and accession contents match.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
