"""Write the residue-number request for the full (untrimmed) family alignments.

For each structure in the motif superposition and each motif position, gives the full family
alignment column (1-based) that the consensus column mapping (motif_columns.COLUMN_MAP) ties to the
superfamily column. map_motif_residues.py reads the member's residue and exact number there.

Writes residue_number_request.tsv.
"""
from build_superposition import EXTRA, STRUCTURES
from final_motif_list import POSITIONS
from motif_columns import HERE, full_family_columns


def main():
    colnum = {p[0]: p[1] for p in POSITIONS + EXTRA}
    with open(HERE / "residue_number_request.tsv", "w") as f:
        f.write("structure\tfamily\tmember_accession\tposition\tsuperfamily_column\t"
                "family_column\tmember_residue\n")
        cols = {}
        for name, fam, _, _, rowacc, _ in STRUCTURES:
            if fam not in cols:
                cols[fam] = full_family_columns(fam, colnum)
            # consOR2 (8UY0) has no UniProt entry; its closest olfactory member stands in, and
            # build_superposition.py carries that member's numbers onto consOR2 by alignment
            for pos, col in colnum.items():
                fc = cols[fam][pos]
                f.write(f"{name}\t{fam}\t{rowacc}\t{pos}\t{col}\t{fc if fc else 'unmapped'}\t\n")
    print("wrote residue_number_request.tsv")


if __name__ == "__main__":
    main()
