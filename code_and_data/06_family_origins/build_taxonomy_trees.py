"""Build the NCBI reference taxonomy of every searched proteome.

Stage 06. The taxonomy is built from every taxID in the searched database, UniProt 2025_02
Uniprot_Ref_Eukaryota_Proteome (taxonomy/searched_taxids_2025_02.tsv: one row per OX= taxID with
its protein count, taken from the FASTA headers), so a searched genome with no hit in any family
or control is a leaf scored absent everywhere rather than missing from the tree.

  ncbi_taxonomy_create.py            taxonomy JSON from nodes.dmp and names.dmp, the eukaryotic
                                     lineage of each taxID, with Unikonta, Archaeplastida, Holozoa
                                     and Holomycota added (CUSTOM_GROUPS)
  generate-ncbi-taxonomy-newick.mjs  the same hierarchy as a newick tree, every tip at the same depth
  ete3 prune()                       single-child ranks merged into their child, lengths summed

Both builders are the ones the GPCRevolution phyletic-distribution page uses. The taxdump used
is the copy of 7 September 2025 (06_family_origins/taxonomy/ncbi_taxdump/ in the archive).

Usage (needs node):
    gunzip -k $GPCR_ZENODO_DIR/06_family_origins/taxonomy/ncbi_taxdump/*.dmp.gz
    python build_taxonomy_trees.py --nodes .../ncbi_taxdump/nodes.dmp --names .../ncbi_taxdump/names.dmp
"""

import argparse
import csv
import json
import os
import subprocess
import sys

from ete3 import Tree

from key_nodes import TAXONOMY_DIR

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SEARCHED_TAXIDS = os.path.join(TAXONOMY_DIR, "searched_taxids_2025_02.tsv")
NCBI_JSON = os.path.join(TAXONOMY_DIR, "taxonomy_ncbi_searched.json")
NCBI_STRETCHED = os.path.join(TAXONOMY_DIR, "taxonomy_ncbi_searched.stretched.nwk")
NCBI_ORDER = os.path.join(TAXONOMY_DIR, "taxonomy_ncbi_searched.tree-order.json")
NCBI_TREE = os.path.join(TAXONOMY_DIR, "taxonomy_ncbi_searched.rankdepth.nwk")

# Searched taxIDs that NCBI has since merged into another, so the taxdump no longer has
# them. The NCBI tree places each at its current taxID and the leaf keeps the searched
# one, which is what the count tables use.
MERGED_TAXIDS = {"1450536": "856835"}  # Guyanagaster necrorhizus


def read_searched_taxids():
    with open(SEARCHED_TAXIDS, newline="", encoding="utf-8") as handle:
        return [row[0] for row in csv.reader(handle, delimiter="	") if row and row[0].isdigit()]


def run(cmd, cwd=None):
    print("$ " + " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=cwd)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--nodes", required=True, help="NCBI nodes.dmp (decompressed)")
    parser.add_argument("--names", required=True, help="NCBI names.dmp (decompressed)")
    args = parser.parse_args()

    taxids = read_searched_taxids()
    print("%d searched taxIDs from %s\n" % (len(taxids), os.path.basename(SEARCHED_TAXIDS)))

    ncbi_taxids = os.path.join(TAXONOMY_DIR, "searched_taxids_2025_02.ncbi_current.tsv")
    with open(ncbi_taxids, "w", newline="", encoding="utf-8") as handle:
        handle.write("\n".join(MERGED_TAXIDS.get(t, t) for t in taxids) + "\n")
    run([sys.executable, os.path.join(SCRIPT_DIR, "ncbi_taxonomy_create.py"),
         "--nodes", args.nodes, "--names", args.names,
         "--taxids", ncbi_taxids, "--output", NCBI_JSON])
    os.remove(ncbi_taxids)
    current_to_searched = {int(v): int(k) for k, v in MERGED_TAXIDS.items()}
    with open(NCBI_JSON, encoding="utf-8") as handle:
        records = json.load(handle)
    for record in records:
        record["taxID"] = current_to_searched.get(record["taxID"], record["taxID"])
    with open(NCBI_JSON, "w", encoding="utf-8") as handle:
        json.dump(records, handle, ensure_ascii=False)
    run(["node", os.path.join(SCRIPT_DIR, "generate-ncbi-taxonomy-newick.mjs"),
         "--input=" + NCBI_JSON, "--output=" + NCBI_STRETCHED, "--order-output=" + NCBI_ORDER])
    tree = Tree(NCBI_STRETCHED, format=1)
    tree.prune(tree.get_leaf_names(), preserve_branch_length=True)
    tree.write(format=1, outfile=NCBI_TREE)
    missing = sorted(set(taxids) - set(tree.get_leaf_names()), key=int)
    print("NCBI: %d leaves -> %s; searched but not placed: %s\n"
          % (len(tree), os.path.basename(NCBI_TREE), ", ".join(missing) or "none"))


if __name__ == "__main__":
    main()
