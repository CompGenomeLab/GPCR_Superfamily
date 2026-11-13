# Literature sources

The receptors and claims of the published GPCR surveys, one table per study in `parsed/` (`parsed/<study>.tsv`:
one row per claim or receptor, as the study's reader in `../literature_readers.py` takes it from the study's
supplementary files; columns kind, study, their_family, genome, status, count, source, evidence, seq_id,
sequence, ids). The comparison (`../compare_literature.py`) reads these tables, so it reruns without the
supplementary files, which are not distributed here (copyright).

To read the supplementary files themselves (to rebuild `parsed/` with `../write_parsed_literature.py`, or to
check the readers against the papers with `../verify_literature_comparison.py`), download the files below from
each study's supplementary material into a folder laid out like this one (one subfolder per study, the names in
the table) and set `GPCR_LITERATURE_ORIGINALS` to it. Files no reader uses are listed in `FILES_NOT_READ` of
`../literature_readers.py`. Our own files are here: `curated_claims.tsv` holds claims stated only in a paper's
text, `genome_name_decisions.tsv` the genome names resolved by hand, `uniprot_accession_taxa.json` the species of
UniProt accessions the studies list without one, and `id_cache.json` the identifier cross-reference table.

| Folder | Study (as in Table S3) | Files |
|---|---|---|
| `behrens_2023/` | Behrens et al. 2023 | `behrens-et-al-2023-a-singular-shark-bitter-taste-receptor-provides-insights-into-the-evolution-of-bitter-taste.pdf`, `pnas.2310347120.sapp.pdf` |
| `cardoso_2014/` | Cardoso et al. 2014 | `pone.0092220.s006.pdf`, `pone.0092220.s007.pdf`, `pone.0092220.s014.fasta` |
| `de_mendoza_2014/` | de Mendoza et al. 2014 | `supp_evu038_SupplementaryTable1.xls` |
| `de_mendoza_2016/` | de Mendoza et al. 2016 | `srep21801-s1.pdf` |
| `fernandez_figueroa_2026/` | Fernandez-Figueroa et al. 2026 | `media-1_SuppTables1-2.xlsx` |
| `garcia_de_las_bayonas_king_2025/` | Garcia De Las Bayonas & King 2025 | `Adhesion_Receptors.txt`, `Frizzled_Receptors.txt`, `Glutamate_Receptors.txt`, `Rhodopsin_Receptors.txt`, `SF1_proteome_accessions.xlsx`, `SF25_GPR157_EukProt_hits.xlsx`, `SF25_GPR157_homologs.txt`, `SF4_choanoflagellate_GPCRs.txt` |
| `gotkhindikar_2025/` | Gotkhindikar et al. 2025 | `kiaf057_supplementary_data.pdf` |
| `grus_zhang_2009/` | Grus & Zhang 2009 | `msn262_Supplementary_Data.zip`, `msn262_paper.pdf` |
| `itoigawa_2024/` | Itoigawa et al. 2024 | `Tables_S1-S2.docx`, `mmc1_supplement.pdf`, `zenodo_10042321` |
| `kowatschew_korsching_2022/` | Kowatschew & Korsching 2022 | `bjac007_suppl_supplementary_material.pdf` |
| `krishnan_2012/` | Krishnan et al. 2012 | `pone.0029817.s001.txt`, `pone.0029817.s002.xls`, `pone.0029817.s003.doc`, `pone.0029817.s004.xls` |
| `krishnan_2014_nematode/` | Krishnan et al. 2014 (nematode) | `pone.0093048.s001.pdf`, `pone.0093048.s002.pdf`, `pone.0093048.s003.pdf`, `pone.0093048.s004.pdf`, `pone.0093048.s005.docx`, `pone.0093048.s006.docx`, `pone.0093048.s007.pdf` |
| `krishnan_2014_sponge/` | Krishnan et al. 2014 (sponge) | `12862_2014_270_MOESM1_ESM.pdf`, `12862_2014_270_MOESM2_ESM.xlsx`, `12862_2014_270_MOESM3_ESM.xlsx`, `12862_2014_270_MOESM4_ESM.pdf`, `12862_2014_270_MOESM5_ESM.pdf`, `12862_2014_270_MOESM6_ESM.pdf`, `12862_2014_270_MOESM7_ESM.pdf`, `12862_2014_270_MOESM8_ESM.pdf` |
| `liu_2026/` | Liu et al. 2026 (Journal of Fungi) | `jof-12-00030-s001.zip` |
| `niimura_2009/` | Niimura 2009 | `evp003_Supplementary_Data.zip` |
| `nordstrom_2009/` | Nordstrom et al. 2009 | `msn228_Supplementary_Data.zip` |
| `nordstrom_2011/` | Nordstrom et al. 2011 | `msr061_Supplementary_Data.zip` |
| `prabhu_eichinger_2006/` | Prabhu & Eichinger 2006 | `1-s2.0-S0171933506000720-main.pdf` |

The query sets of the MMseqs2 search that `../build_remaining_fasta.py` and `../fetch_idonly_round2.py` write
also read some of these files; their outputs are in the Zenodo archive (`05_literature_comparison/`).

## References

As in the manuscript's reference list.

- `behrens_2023/`: Behrens, M., Lang, T. & Korsching, S. I. A singular shark bitter taste receptor provides insights into the evolution of bitter taste perception. Proceedings of the National Academy of Sciences 120 (2023). https://doi.org/10.1073/pnas.2310347120
- `cardoso_2014/`: Cardoso, J. C. R., Félix, R. C. & Power, D. M. Nematode and Arthropod Genomes Provide New Insights into the Evolution of Class 2 B1 GPCRs. PLoS ONE 9, e92220 (2014). https://doi.org/10.1371/journal.pone.0092220
- `de_mendoza_2014/`: de Mendoza, A., Sebé-Pedrós, A. & Ruiz-Trillo, I. The evolution of the GPCR signaling system in eukaryotes: modularity, conservation, and the transition to metazoan multicellularity. Genome Biol Evol 6, 606-619 (2014). https://doi.org/10.1093/gbe/evu038
- `de_mendoza_2016/`: De Mendoza, A., Jones, J. W. & Friedrich, M. Methuselah/Methuselah-like G protein-coupled receptors constitute an ancient metazoan gene family. Scientific Reports 6, 21801 (2016). https://doi.org/10.1038/srep21801
- `fernandez_figueroa_2026/`: Fernandez-Figueroa, V., Quercia-Raty, C. A., Gallastegui-Ulloa, J., Robeson, L. & Brauchi, S. E. Canonical G-protein coupled receptors of vascular plants. Plant Physiology and Biochemistry 236, 111444 (2026). https://doi.org/10.1016/j.plaphy.2026.111444
- `garcia_de_las_bayonas_king_2025/`: Garcia De Las Bayonas, A. & King, N. G protein-coupled receptor diversity and evolution in the closest living relatives of Metazoa (eLife Sciences Publications, Ltd, 2025).
- `gotkhindikar_2025/`: Gotkhindikar, A., Chakravorty, D., Sengupta, D., Joshi, M. & Assmann, S. M. Is GCR1 the GPR157 of plants? Plant Physiology 197 (2025). https://doi.org/10.1093/plphys/kiaf057
- `grus_zhang_2009/`: Grus, W. E. & Zhang, J. Origin of the Genetic Components of the Vomeronasal System in the Common Ancestor of all Extant Vertebrates. Molecular Biology and Evolution 26, 407-419 (2009). https://doi.org/10.1093/molbev/msn262
- `itoigawa_2024/`: Itoigawa, A., Toda, Y., Kuraku, S. & Ishimaru, Y. Evolutionary origins of bitter taste receptors in jawed vertebrates. Current Biology 34, R271-R272 (2024). https://doi.org/10.1016/j.cub.2024.02.024
- `kowatschew_korsching_2022/`: Kowatschew, D. & Korsching, S. I. Lamprey possess both V1R and V2R olfactory receptors, but only V1Rs are expressed in olfactory sensory neurons. Chemical Senses 47 (2022). https://doi.org/10.1093/chemse/bjac007
- `krishnan_2012/`: Krishnan, A., Almén, M. S., Fredriksson, R. & Schiöth, H. B. The Origin of GPCRs: Identification of Mammalian like Rhodopsin, Adhesion, Glutamate and Frizzled GPCRs in Fungi. PLoS ONE 7, e29817 (2012). https://doi.org/10.1371/journal.pone.0029817
- `krishnan_2014_nematode/`: Krishnan, A., Almén, M. S., Fredriksson, R. & Schiöth, H. B. Insights into the Origin of Nematode Chemosensory GPCRs: Putative Orthologs of the Srw Family Are Found across Several Phyla of Protostomes. PLoS ONE 9, e93048 (2014). https://doi.org/10.1371/journal.pone.0093048
- `krishnan_2014_sponge/`: Krishnan, A. et al. The GPCR repertoire in the demosponge Amphimedon queenslandica: insights into the GPCR system at the early divergence of animals. BMC Evolutionary Biology 14 (2014). https://doi.org/10.1186/s12862-014-0270-4
- `liu_2026/`: Liu, Z., Salamov, A. & Grigoriev, I. V. Identification and Classification of Fungal GPCR Gene Families. Journal of Fungi 12, 30 (2026).
- `niimura_2009/`: Niimura, Y. On the Origin and Evolution of Vertebrate Olfactory Receptor Genes: Comparative Genome Analysis Among 23 Chordate Species. Genome Biology and Evolution 1, 34-44 (2009). https://doi.org/10.1093/gbe/evp003
- `nordstrom_2009/`: Nordstrom, K. J. V., Lagerstrom, M. C., Waller, L. M. J., Fredriksson, R. & Schioth, H. B. The Secretin GPCRs Descended from the Family of Adhesion GPCRs. Molecular Biology and Evolution 26, 71-84 (2008). https://doi.org/10.1093/molbev/msn228
- `nordstrom_2011/`: Nordstrom, K. J. V., Sallman Almen, M., Edstam, M. M., Fredriksson, R. & Schioth, H. B. Independent HHsearch, Needleman-Wunsch-Based, and Motif Analyses Reveal the Overall Hierarchy for Most of the G Protein-Coupled Receptor Families. Molecular Biology and Evolution 28, 2471-2480 (2011). https://doi.org/10.1093/molbev/msr061
- `prabhu_eichinger_2006/`: Prabhu, Y. & Eichinger, L. The Dictyostelium repertoire of seven transmembrane domain receptors. European Journal of Cell Biology 85, 937-946 (2006). https://doi.org/10.1016/j.ejcb.2006.04.003
