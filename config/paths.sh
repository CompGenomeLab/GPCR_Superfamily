#!/usr/bin/env bash
# Shared data-path resolver. Source this file before running a pipeline stage.

GPCR_REPOSITORY_DIR="${GPCR_REPOSITORY_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
GPCR_ZENODO_DIR="${GPCR_ZENODO_DIR:-${GPCR_REPOSITORY_DIR}/../zenodo_data}"

export GPCR_REPOSITORY_DIR GPCR_ZENODO_DIR
export GPCR_SEARCH_RESULTS_DIR="${GPCR_SEARCH_RESULTS_DIR:-${GPCR_ZENODO_DIR}/search_results}"
export GPCR_EXTRACTED_CLADE_DIR="${GPCR_EXTRACTED_CLADE_DIR:-${GPCR_ZENODO_DIR}/extracted_clades}"
export GPCR_FAMILY_DIR="${GPCR_FAMILY_DIR:-${GPCR_ZENODO_DIR}/families}"
export GPCR_ANALYSIS_DIR="${GPCR_ANALYSIS_DIR:-${GPCR_FAMILY_DIR}}"
export GPCR_FAMILY_TABLE_DIR="${GPCR_FAMILY_TABLE_DIR:-${GPCR_FAMILY_DIR}/tables}"
export GPCR_ACCESSION_TXT_DIR="${GPCR_ACCESSION_TXT_DIR:-${GPCR_FAMILY_DIR}/accessions}"
export GPCR_CLUSTER_DIR="${GPCR_CLUSTER_DIR:-${GPCR_FAMILY_DIR}/clusters}"
export GPCR_TRIMMED_MSA_DIR="${GPCR_TRIMMED_MSA_DIR:-${GPCR_FAMILY_DIR}/trimmed_msas}"
export GPCR_TAXONOMY_JSON="${GPCR_TAXONOMY_JSON:-${GPCR_FAMILY_DIR}/reference_data/taxonomy_eukaryotes_filtered.json}"
export GPCR_ACCESSION_TSV="${GPCR_ACCESSION_TSV:-${GPCR_FAMILY_TABLE_DIR}/acc_taxid_map_BestEval_May22_e0.1.tsv}"
export GPCR_FAMILY_DB_ACCESSION_TSV="${GPCR_FAMILY_DB_ACCESSION_TSV:-${GPCR_FAMILY_TABLE_DIR}/acc_taxid_map_prefilter_Oct24_e1000_all.tsv}"
export GPCR_CLASSA_CLUSTER_TSV="${GPCR_CLASSA_CLUSTER_TSV:-${GPCR_SEARCH_RESULTS_DIR}/classA/classA_hmm_E1000_db_clust0.8.tsv}"
