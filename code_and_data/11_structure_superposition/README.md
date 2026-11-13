# 11 · Structural superimposition of the shared motifs

Methods: *Superimposition of Identified Motifs*; Figures 5 and 6; Table S10.

One structure per family is used: the experimental structure where one exists, otherwise the AlphaFold
model of the family member carrying the most motif residues (`model_candidates.py`, `af_confidence.py`),
plus AlphaFold models of a plant cAMP receptor, a plant Class C receptor and a second STE2.
`final_motif_list.py` lists the motif positions and the residue each family carries there (from stage
10), `residue_number_request.py` and `map_motif_residues.py` (cluster) give the residue number of each
motif residue in each structure, and `build_superposition.py` superimposes the structures on their
DeepTMHMM transmembrane residues only (first on the GP157 model, then refined against the mean
transmembrane core; register shifts of up to two helical turns tested), so no sequence correspondence
enters the fit. For each motif, the distance of each carrying residue's Cα from the centroid of the other
carrying residues and the same measure for its local backbone (±3 residues) are computed
(`local_backbone_check.py`, `per_family_overlap.py`). `render_fig5_panels.py` renders one PyMOL close-up
per motif (`fig4_panels/`), `build_fig5_layout_tool.py` writes a page for placing the panels
(`fig4_layout.json`), `make_fig5.py` assembles Figure 5 from them, and `save_fig4_session.py` saves the
PyMOL session of the figure.

`figures/build_motif_evidence.py` applies the fixed rules of Figure 6 (sequence: conserved by
the 60% rule, I+V+L together at 3x54 and 5x65, any kingdom for cAMP and Class C; structure: within 2 Å
of the centroid of the other carriers or at most 1 Å beyond the local backbone displacement, any carrier
at C45x50; STE2 a control only) to Tables S8 and S10b, and writes `motif_evidence_data.js`, which
`fig6_hierarchy.html` (Figure 6) draws. Open it in a browser; it has a "Download SVG" button.

```bash
python final_motif_list.py
python residue_number_request.py          # then, where the family alignments of stage 10 are:
#   python map_motif_residues.py --request residue_number_request.tsv \n#       --families-dir $GPCR_ZENODO_DIR/10_motif_analysis/family_alignments/families \n#       --out motif_residue_numbers.tsv --sequences-out member_sequences.tsv
python build_superposition.py gp157       # needs PyMOL
python local_backbone_check.py gp157
python per_family_overlap.py gp157
python render_fig5_panels.py             # needs PyMOL; writes fig4_panels/
python make_fig5.py                       # fig4_layout.json + fig4_panels/ -> fig5_structure.svg
python figures/build_motif_evidence.py    # after the supplementary tables
```

Verified: rerun from the archive, the motif list, the core fits and every distance table are identical
to those of the paper, and `make_fig5.py` rebuilds Figure 5 byte-identical.

**Archive** (`11_structure_superposition/`): `structures/` (PDB entries) and `af_models/` (AlphaFold
models), the motif residue numbers, the per-family candidate tables, the superposition and distance
tables, the PyMOL sessions (`motif_superposition_tmcore_gp157.pse`), the motif panels (`fig4_panels/`), the
panel angles (`fig4_views.json`) and layout (`fig4_layout.json`) and Figure 5 (`fig5_structure.svg`).
