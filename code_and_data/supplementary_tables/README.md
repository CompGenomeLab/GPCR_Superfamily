# Supplementary tables and figures

Each `build_S*.py` writes one supplementary table as a TSV in `supplementary_tables/` of the archive,
from the outputs of the stage that produces it; `build_xlsx.py` turns them into the Excel workbooks
that are distributed (title in row 1, legend in row 2, table from row 4, a sheet per lettered part,
display names for the families).

| Table | Contents | Built by | Stage |
|---|---|---|---|
| S1 | per family and control: search, sequences retrieved, gene tree used, extracted clades, counts after each filter | `build_S1_search_and_filtering.py` | 01-04 |
| S2 | the 382 sequences claimed by two families and how each was assigned | `build_S2_contested_sequences.py` | 04 |
| S3a-d | literature receptors, agreement by study and by family, deeper-origin candidates | `build_S3_literature.py` | 05 |
| S4 | the distribution behind Figure 2: genomes and genera sampled and carrying each family per lineage, presence call, Dollo origin and our interpretation | `build_S4_lineage_distribution.py` | 06 |
| S5 | the 590 representatives and 18 literature sequences | `build_S5_representatives.py` | 07 |
| S6a-h, S7 | tree runs, family stability and stems, groupings, placements, literature placement, family-level splits with sCF; transfer index of every sequence | `build_S6_S7_tree_support.py` | 09 |
| S8 | conservation at the 28 motif positions and every other position named in the text | `build_S8_motif_conservation.py` | 10 |
| S9 | conserved residues STE3 and STE2 share with each family | `build_S9_ste_similarity.py` | 10 |
| S10a-b | structures superimposed; distance of every motif residue from the centroid and its local backbone | `build_S10_structures.py` | 11 |

```bash
for b in build_S*.py; do python "$b"; done
python build_xlsx.py
```

Supplementary figures (`supplementary_figures/` in the archive, legends in
`supplementary_figure_legends.md`): Figs. S1 and S4 are exports of gpcrevolution.org/taxonomic-distribution,
before and after the reassignment of Dikarya Class A sequences to GPR1;
Fig. S2 is the superfamily-logo export with the default guide tree, cropped to the TM5 prolines by
`make_figS2.py` (`crop_logo_svg.py`); Fig. S3 is `../11_structure_superposition/figures/figS3_motif_evidence.html`
exported by `make_figS3.py`. `build_supplementary_figures_docx.py` assembles `Supplementary_Figures.docx`.
The two export scripts rasterize SVGs with headless Chrome (set `CHROME` to its path).
