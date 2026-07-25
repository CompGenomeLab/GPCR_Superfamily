# Analysis code and data

This directory follows the order of the GPCR superfamily analysis. The top-level [repository README](../README.md) gives the general overview and Zenodo data policy; the stage READMEs document the files and settings needed to run each part.

| Step | Directory | Main result |
|---|---|---|
| 1 | `homolog_search_pipeline/` | Broad HMMER/MMseqs2 search results, exploratory MSAs, and gene trees |
| 2 | `post_processing/` | Score-ordered trees/alignments and manually extracted homolog sets |
| 3 | `family_processing/` | Accession preparation, E-value filtering, taxonomic mapping, clustering, and family-wide MSAs |
| 4 | `phylogeny/` | Selected representatives, primary alignments and trees, and the final MSA visualization |
| 5 | `superfamily_alignment/` | Extended 550-sequence alignment and motif/logo mappings |
| 6 | `structure_visualization/` | Structural validation resources |

The MSA visualization generator is stored directly under `phylogeny/` because it renders the final alignment used for tree inference.

The workflow includes GPCR families, STE2, the heterotrimeric G-protein subunits, and the TBP/RBP1 controls. RAMP, RGS, Arrestin, and Dvl are excluded from the final configured analysis.
