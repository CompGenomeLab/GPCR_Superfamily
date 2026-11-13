#!/usr/bin/env bash
# Shared data-path resolver for the shell (SLURM) scripts. Source it before running one:
#     source config/paths.sh
# The Zenodo archive is expected beside the repository; set GPCR_ZENODO_DIR first if it is
# elsewhere. Python scripts resolve the same locations through code_and_data/gpcr_paths.py.

GPCR_REPOSITORY_DIR="${GPCR_REPOSITORY_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
GPCR_ZENODO_DIR="${GPCR_ZENODO_DIR:-${GPCR_REPOSITORY_DIR}/../zenodo_data}"
GPCR_CODE_DIR="${GPCR_REPOSITORY_DIR}/code_and_data"
export GPCR_REPOSITORY_DIR GPCR_ZENODO_DIR GPCR_CODE_DIR

# External inputs that are not part of the archive (only needed to rerun the searches).
export UNIPROT_FASTA="${UNIPROT_FASTA:-/path/to/Uniprot_Ref_Eukaryota_Proteome.fasta}"
export UNIPROT_MMSEQS_DB="${UNIPROT_MMSEQS_DB:-/path/to/UniprotRefEuk}"
