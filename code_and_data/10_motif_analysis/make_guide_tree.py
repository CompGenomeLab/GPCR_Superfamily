"""Guide tree for a tree-guided MAFFT G-INS-i superfamily alignment, from an ML best tree:
(STE2, (Glutamate-related, rest)).

STE2 is taken out of the tree and joined back at the root, so it is merged last, as one block,
against the whole superfamily profile instead of wherever the ML tree places it. The remaining
tree is rooted between the Glutamate-related families (Class C, Vomeronasal 2) and the rest.
Family membership is read from 07_representative_selection/representatives_all_families.tsv; leaves not
in it (literature sequences) count as "rest".

A non-STE2 leaf inside the STE2 clade stays in the rest and is reported. A tree without STE2
gives (Glutamate, rest); a tree without Glutamate-related families stops, since the root for
those sets is not decided. Support values and internal labels are dropped (MAFFT merge steps
cannot carry them). The STE2 subtree keeps the stem length it had in the ML tree; the root
splits the Glutamate stem in half.

python make_guide_tree.py <best tree> <out.nwk>
Then prepare_mafft_treein.py (GPCR_Superfamily/code_and_data/superfamily_alignment) with the
unaligned FASTA of the same sequences, --tree <out.nwk>.
"""
import copy
import csv
import sys
from pathlib import Path

from Bio import Phylo

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import stage  # noqa: E402

MEMBERSHIP = stage("07_representative_selection", "representatives_all_families.tsv")
GLUTAMATE = {"classC", "Vomeronasal2"}


def family_of(leaf, fam):
    parts = leaf.name.split("|")
    return fam.get(parts[1] if len(parts) > 2 else leaf.name.split("/")[0], "rest")


def clade_of(tree, members, leaves):
    """Clade spanned by members when the tree is rooted on a leaf outside them."""
    tree.root_with_outgroup(next(x for x in leaves if x not in members))
    return tree.common_ancestor(members)


def main(tree_path, out_path):
    fam = {r["accession"]: r["gpcrfamily"] for r in csv.DictReader(open(MEMBERSHIP), delimiter="\t")}
    tree = Phylo.read(tree_path, "newick")
    leaves = tree.get_terminals()
    ste2 = [x for x in leaves if family_of(x, fam) == "STE2"]
    glu = [x for x in leaves if family_of(x, fam) in GLUTAMATE]
    print(f"{len(leaves)} leaves: STE2 {len(ste2)}, Glutamate-related {len(glu)}, rest "
          f"{len(leaves) - len(ste2) - len(glu)}")
    if not glu:
        sys.exit("no Glutamate-related leaves: root for this set not decided")

    ste2_tree = None
    if ste2:
        clade = clade_of(tree, ste2, leaves)
        extra = [x.name for x in clade.get_terminals() if x not in ste2]
        if extra:
            print(f"  non-STE2 leaves inside the STE2 clade, kept in the rest: {extra}")
        stem = clade.branch_length or 0.0
        ste2_tree = copy.deepcopy(clade)
        for name in extra:
            ste2_tree.prune(name)
        ste2_tree.branch_length = stem
        for leaf in ste2:
            tree.prune(leaf)

    leaves = tree.get_terminals()
    glu_clade = clade_of(tree, glu, leaves)
    inside = [x.name for x in glu_clade.get_terminals() if x not in glu]
    if inside:
        sys.exit(f"Glutamate-related families are not one clade; also inside: {inside}")
    stem = glu_clade.branch_length or 0.0
    tree.root_with_outgroup(glu_clade, outgroup_branch_length=stem / 2)
    if ste2_tree is not None:
        rest = tree.root
        rest.branch_length = 0.0
        tree = Phylo.BaseTree.Tree(root=Phylo.BaseTree.Clade(clades=[ste2_tree, rest]), rooted=True)
    for c in tree.find_clades():
        c.confidence = None
        if not c.is_terminal():
            c.name = None
    Phylo.write(tree, out_path, "newick", format_branch_length="%.6f")
    top = [len(c.get_terminals()) for c in tree.root.clades]
    print(f"  wrote {out_path}: root split {top}, {len(tree.get_terminals())} leaves")


if __name__ == "__main__":
    main(*sys.argv[1:3])
