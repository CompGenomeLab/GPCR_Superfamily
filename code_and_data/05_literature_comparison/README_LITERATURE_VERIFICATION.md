# Literature comparison: how it is built and how to verify it

This folder compares every published report of a GPCR family in a genome with our family calls,
and follows every receptor the studies list, of any family, through our pipeline. It is written
for an independent reviewer who will rerun and check it. Results are in
`05_literature_comparison/Supplementary_Table_literature_receptors.tsv` and the two
`Supplementary_Table_literature_agreement_by_*.tsv` tables (built by `build_supplementary_table.py` and
`build_agreement_tables.py`). All MMseqs2 queries (one per unique sequence, including those with no hit) are in
`05_literature_comparison/literature_mmseqs_queries_all.fasta`, their hits in `literature_mmseqs_all.m8`, and
`literature_mmseqs_query_map.tsv` links every receptor query to the query its sequence was searched under.

**The reference is our own data.** A receptor is matched only to the database we searched
(UniProt reference proteomes, release 2025_02) as our pipeline holds it: the full-length sequences
of our clade members and search-tree inputs, and the headers of every search hit. The current
UniProt is never queried; it has changed since the release we searched (entries created, merged
or retired), so it cannot say whether a protein was in our data.

## What to verify, in short

1. Each reader recovers what the supplement actually says (badly formatted PDFs, Word and Excel
   files were parsed; some studies' own files disagree with each other).
2. Every supplementary file is either read or listed as not read, with the reason.
3. Genomes are matched correctly between a study's names and our taxIDs.
4. Receptors are matched to the right protein of ours, and only in the species the study names.
5. Our side is read from the latest stage tables.
6. The comparison categories, the explanation column and the family correspondence follow from
   the data.

## Files

| File | Role |
|---|---|
| `literature_readers.py` | one reader per study; returns *claims* (family present/absent in a genome, with count) and *sequences* (listed receptors with sequence and/or identifiers); `FILES_NOT_READ` lists the supplementary files no reader uses, with the reason |
| `literature/curated_claims.tsv` | claims stated only in a paper's text or a small table, each row with its source and, where no receptor is listed, its `evidence` |
| `literature/genome_name_decisions.tsv` | genome names replaced by hand (a name NCBI resolves to another species than the genome the study used, or not at all), each with its evidence; applied by `literature_readers.read_all` |
| `literature/<study>/` | the supplementary files as downloaded (one folder per study) |
| `compare_literature.py` | builds the tables from the readers, our sequences, our search outputs and our stage tables |
| `verify_literature_comparison.py` | reader totals against the supplements, file coverage, and internal consistency of the outputs (local only) |
| `literature_deeper_origins.py` | literature presences outside our origin node where our data lack the family; every receptor with a sequence (and Itoigawa's TAS2R loci, translated) scored against all 19 of our families; writes `DEEPER_ORIGINS.md` |
| `literature_origin_conflicts.py` | literature claims against our origin calls by the 20 lineages of the origin figure (with ctenophores and non-chromadorean nematodes placed separately) |
| `05_literature_comparison/literature_comparison.tsv` | one row per study x family x genome (families with a counterpart) |
| `05_literature_comparison/literature_sequences.tsv` | one row per listed receptor, of every family: where its sequence comes from (`sequence_source`), how its genome was found (`genome_how`), the protein of ours it matches and how (`match`), `stopped_at` and `explanation` |
| `05_literature_comparison/literature_genome_proxies.tsv` | study genomes not in our set that are compared in a sister genome of ours (same genus), with the receptor evidence for each |
| `05_literature_comparison/literature_summary.tsv` | counts per study, family and outcome; the genome claims of families with no counterpart |
| `05_literature_comparison/literature_unresolved_genomes.tsv` | genome names still not matched to ours, with the closest of our species by name and our species of the same genus (candidates for review, not used as matches) |
| `05_literature_comparison/literature_unmapped_families.tsv` | for each study family with no counterpart (and those mapped by where their receptors end): where its receptors end in our data |

## Inputs (all read by `compare_literature.py`, from the archive)

| What | File |
|---|---|
| Clade membership | `02_clade_extraction/clade_membership.tsv` |
| Single 7TM core | `03_topology_filter/core_membership.tsv` |
| Topology filter (0.8) | `03_topology_filter/topology_membership.tsv` |
| Final families | `04_overlap_assignment/final_membership.tsv` |
| Per-genome counts | `04_overlap_assignment/tax_counts_per_family.tsv` |
| Searched genomes | `06_family_origins/taxonomy/taxonomy_ncbi_searched.json` (2,760 taxIDs of the 2025_02 reference proteomes; one, Guyanagaster necrorhizus 1450536, is missing from the taxdump and is placed by its species name) |
| NCBI taxonomy | `06_family_origins/taxonomy/ncbi_taxdump/` (`names.dmp`, `nodes.dmp`) |
| Our sequences | `02_clade_extraction/clade_sequences/<family>.fasta.gz` (full-length clade members) and `01_homology_search/<search>/*_db.fasta.gz` (full-length hits each search tree was built from) |
| Our searches | `01_homology_search/<search>/<search>_hmmsearch_E1000.txt` (HMMER) and the two MMseqs2 `.m8` tables; their headers carry accession, entry name, `OX=` and `GN=` as in the searched release |
| Search trees | `01_homology_search/<search>/*_ordered.tree` (the trees the clades were read from) |

## How the tables are built

1. **Read** every study (`literature_readers.read_all`). Only what a study reports enters; a genome
   a study does not mention gets no row. Nucleotide sequences (Nordstrom 2009 File 2) are
   translated. Sequences from alignments of genome contigs (Grus & Zhang) keep their gaps, and
   Niimura's pseudogenes keep their frameshifts ("/") as breaks. Page numbers inside sequences
   (sponge PDF) are dropped, and a Word record stored in two pieces (Nordstrom 2011 nV_35931) is
   joined, so no sequence is shortened by the reading: every sequence under 150 aa was checked
   against its record in the study's file. The genome names of `literature/genome_name_decisions.tsv`
   are replaced here.
2. **Match genomes** at NCBI species level (scientific names, synonyms, strain -> species;
   abbreviations like "C. elegans" expanded against our genomes and the studies' full names). A
   name NCBI does not know is matched to one of our species if it is a spelling variant (same
   genus; epithet equal once Latin endings, hyphens and one typo are set aside, or cut at a hyphen;
   a single candidate; "sp." never), e.g. "Aspergillus nomius" = *A. nomiae*; a misspelt genus
   counts only with the epithet exact and a single candidate ("Phyomyces blakesleeanus"). A name NCBI files
   under a species we do not have is matched to our genome when its strain designation (>= 4
   characters with a digit; a strain that is only a number counts only within the same genus,
   since collection numbers repeat across collections: "Trametes maxima CIRM-BRFM 1813" is not
   *Coemansia* sp. RSA 1813) appears word for word in exactly one of our genomes' names, e.g.
   "Penicillium chrysogenum Wisconsin 54-1255" = *P. rubens* Wisconsin 54-1255. Claims for several
   strains of one species are combined: present if any strain has the family, the largest count.
   Names that no rule resolves correctly are decided by hand in
   `literature/genome_name_decisions.tsv` (14 rows, each with its evidence: Krishnan 2012 names
   whose receptors identify the genome they used, e.g. "Magnaporthe grisea" = *Pyricularia oryzae*
   70-15 by its MGG_ locus tags; renamed species; Liu et al. "Peniophora sp.", whose receptors are
   identical to our *Peniophora* sp. CONT). These rows read "manual decision" in `match` and
   `genome_how`.
3. **Match receptors to our data**. First by accession: a UniProt accession the study gives
   (alone, or inside an entry name such as Q8NIR1_NEUCR) that is in our search headers is that
   protein, and the protein's genome is the receptor's genome whatever species name the study uses
   (`match` = "accession ...", with "(sequence differs)" when the study's sequence is neither equal
   to nor contained in ours). A genome name that does not resolve to one of ours takes the genome
   of the study's own matched receptors under that name, when they all fall in one species of ours
   ("Endogone sp FLAS 59071" = Endogone sp. FLAS-F59071); its claims are then compared in that
   genome. A receptor listed with no species and not matched takes the genome of its UniProt
   accession (or UniParc ID) as UniProt records it, from the current entry or, for a retired entry,
   its last version in UniSave (`literature/uniprot_accession_taxa.json`; species codes in entry
   names change when an organism is renamed, accessions do not); failing that, the species of its
   entry-name code as our search headers carry it.
   Species named by the readers where a study's own file gives none: Krishnan 2012 Table S3 and
   Dataset S1 locus tags by their Broad prefix (`KRISHNAN12_LOCUS_SPECIES`; taxID from our search
   headers or NCBI Protein, every species checked against their Table S2; SNOT_/PGUT_ transcript
   forms read as the SNOG_/PGUG_ genes); Liu et al. File S2 portals that Table S3 shortens or spells
   differently (`LIU_PORTAL_ALIASES`), and two it lacks, named by receptors identical to our proteins (`LIU_PORTAL_SPECIES`: Amyenc1, Aspnov1); Garcia headers with no
   species (`GARCIA_HEADER_SPECIES`: TRINITY = *Choanoeca flexa*, all 19 found in its published
   proteome; genemark-NKKH001C = *Barroeca monosierra*; Nk52 = an ichthyosporean outside our set),
   and EukProt lineages without a species name ("MAST-03A sp."). The rest, in the species of the receptor:
   - by sequence (our full-length sequences of that species): exact; else one sequence inside the
     other (>= 150 aa); else every gap-free piece of >= 30 aa inside the same protein (contig
     translations whose gaps can stand for missing exons). When a study names no species, only an
     exact match anywhere counts.
   - by identifier, for receptors not matched by sequence: accession, entry name or `GN=` in the
     headers of our search hits, with simple variants (the accession of an entry name; a locus tag
     without its `T0`, `-PA` or version suffix). When the study names no species, an identifier
     counts only if every protein carrying it is in one species.
   - by a near-identical sequence, for the rest (another gene model of the same protein: edited,
     older or newer): in the named species, the 20 candidates sharing most 5-mers (at least 30% of them), the best accepted
     if the identical residues of a global alignment (BLOSUM62, gap open 11, extend 1, end gaps
     free) cover at least 90% of the **longer** sequence (`match` = "similar sequence (score)").
     A short study sequence (partial gene, pseudogene, alignment row) is therefore not matched to a
     full-length protein of ours, nor a fragment of ours to a full-length study sequence; these are
     left to the MMseqs2 search.
   - by a cross-reference table, last: identifiers our headers do not carry (dictyBase, Ensembl,
     NCBI gene, AGI, ORF names) are translated to accessions with `id_cache.json`, the identifier
     table of the earlier UniProt-based run. The accession must be in our own search headers, in
     the named species. These rows read "identifier ... (cross-reference table)"; set
     `USE_XREF_TABLE = False` in `compare_literature.py` to use no UniProt-derived information.
   - A receptor that matches several proteins (identical paralogs or alleles) is followed through
     the one that goes furthest in the family the study named.
4. **Receptors not placed in their own genome** (genome not in our set, or not found in it): an
   identical protein of ours in any other species, or else a near-identical one (>= 90%) in another
   species of the same genus, is reported in `other_species_match`, `other_species_how` and
   `other_species_shared_taxon` (the nearest taxon the two species share, with its rank). A distant
   one may point to contamination or horizontal transfer and is for review.
   **Sister genomes**: a study genome not in our set is compared in a genome of ours of the same
   genus when at least half of its receptors with a sequence (and at least 2) are identical or
   near-identical to proteins of that one genome (`MIN_PROXY_SHARE`, `MIN_PROXY_RECEPTORS`). Its
   claims are then compared in that genome and its receptors followed through those proteins
   (`genome_how` = "sister genome"; `match` = "sister genome <name>: ..."). These compare a claim
   about one species with a closely related one; every pairing, with its counts of identical and
   near-identical receptors, is in `literature_genome_proxies.tsv`.
5. **Our side**: for the matched protein, the family it belongs to at each stage, the searches that
   hit it, and whether it is a leaf of its family's search tree.
6. **Family correspondence** (`FAMILY_MAP` in `compare_literature.py`). Deliberate choices:
   the studies' broad Rhodopsin = our class A + its nested lineages (Olfactory, Nematode, class T,
   Vomeronasal 1, GPR1); Glutamate = class C + Vomeronasal 2; Adhesion = class B2 + Mth; Frizzled =
   class F + FSLB; Pfam domains (de Mendoza 2014) = every family of ours built on that domain;
   srw = Nematode only; GCR1 = cAMP or GP157.
   A study family with no named counterpart is mapped by where its receptors end (`INFERRED_MAP`):
   to every family of ours that holds at least one of its final receptors, with no threshold. This
   maps Liu et al. class 7 to FSLB, class 16 to class C, class 15 to GPR1 + cAMP + class A, class 17 to
   FSLB + cAMP, and Nordstrom 2011 "7tm_2_other" to class B2 + Mth. Study families with no receptor
   final in any family of ours are not compared, but their receptors are traced and `literature_unmapped_families.tsv` shows where they end.
7. **Comparison**: agree present / agree absent / they report it and we do not (with the reason:
   another family of ours, removed by our filters, in no clade, their receptors not found in our
   data, or no receptor to trace) / they report it absent and we find it / genome not in our set.
   `evidence` says whether receptor sequences or IDs back the claim, or what else does (genomic
   loci, named genes, a count).

`stopped_at` for a receptor: genome not in our set -> not found in our data (no sequence or
identifier match) -> no search hit -> not hit by its family's search -> outside its family's tree
(sequence limit; class A and Olfactory take members from complete lists) -> not in its family's
clade -> 7TM filter -> topology filter (0.8) -> overlap resolution -> final in its family, or
assigned to another family of ours. For a study family with no counterpart: "final in our X" or
the step where it stops.

## Presences that would deepen an origin (`literature_deeper_origins.py`)

Our origin of each family is the interpretation column of the lineage figure. A literature presence
is a candidate when the study reports the family in a species outside that origin node (NCBI
taxonomy; Unikonta = Opisthokonta + Amoebozoa and Holozoa = Opisthokonta outside Fungi and
Rotosphaerida, as in our genome table) and our data do not show the family in that genome. A name
the NCBI taxonomy does not resolve is placed by the species its abbreviation stands for or by the
taxon its first word names; a claim is dropped only when every taxon the name may belong to lies
inside our origin, and is otherwise tested. Its receptors are taken from `literature_sequences.tsv`,
and every one with a sequence is tested (the study's sequence; for a receptor listed by identifier
only and matched to a protein of ours, that protein's, with its fate in our pipeline beside it).
The TAS2R loci of Itoigawa et al. 2024 are read at their coordinates (NCBI contigs; the GWH assembly
of *Potamotrygon leopoldi*, GWHAOTN00000000, downloaded once and set with `GPCR_GWH_GENOME`),
cached in `itoigawa_tas2r_loci.fasta`, and translated (longest reading frame from a methionine if
>= 200 aa, otherwise the frame with fewest stops, read through, so pseudogenes are tested too).
Behrens et al. 2023 give no sequences; a Behrens claim is tested with Itoigawa's locus in the same
species. Loci on unpublished contigs (*Hemitrygon akajei*, ctg_hemak) and Behrens species Itoigawa
did not search (*Ginglymostoma cirratum*, *Hemiscyllium ocellatum*, *Scyliorhinus canicula*) have no
sequence and are listed untested. Outputs: `deeper_origin_candidates.tsv`, `DEEPER_ORIGINS.md`.
The tested sequences (`write_deeper_origin_fasta.py`) are searched with MMseqs2 against the searched
proteomes, and `deeper_origin_best_hits.py` gives each the family of its best hit (E <= 1e-5, query
coverage >= 0.5) among the final members of our families, or no hit (`literature_closest_family.tsv`).

## Running

```bash
cd code_and_data/05_literature_comparison
python compare_literature.py            # a few minutes: reads our sequences and every search output
python verify_literature_comparison.py  # local only; prints PASS/FAIL
python literature_origin_conflicts.py   # prints lists A, B, C by lineage
python literature_deeper_origins.py     # fetches the TAS2R loci from NCBI once
```

Only `literature_deeper_origins.py` uses the network, and only for TAS2R loci not yet in its cache.
Of the `*_cache.json` files left by the earlier, UniProt-based version, only `id_cache.json` is
read (as the cross-reference table above).

## Points to check independently

- **Supplements that disagree with themselves**: Krishnan 2012 Dataset S1 has one Glutamate
  sequence (AMAG_16966) not in Table S1 and repeats two sequences; Table S1 is used. Liu et al. 2026
  Table S3 counts disagree with the File S2 receptor list in 4-8% of genome x class cells; File S2 is
  used where it has the genome. Garcia De Las Bayonas & King Supplementary file 4 holds 905
  sequences against 918 in the text. Krishnan 2012 Table S3 holds 618 accessions against 619 in the
  text, and Table S2 lists 76 fungal species against 79 in the text (58 Ascomycota against 61). de Mendoza 2014 Table S1 has 78 genomes against "75 genomes and five transcriptomes".
  Cardoso 2014 s014 spells some genes differently from s007 ("HMELO14214" for HMEL014214); 108 of
  122 genes are joined to their sequence.
- **Parsing from PDF and Word**: de Mendoza 2016 (sequences in the supplement PDF), Kowatschew
  (SI file 3 in the PDF), Prabhu (Tables 1-2 in the paper PDF), Cardoso (Table S2 PDF), Krishnan 2014
  nematode (s007 PDF, s005 docx), Krishnan 2014 sponge (Additional file 1 PDF), Grus & Zhang (Files 2
  and 9, PDF alignments), Niimura (Dataset S1 .doc), Nordstrom 2011 (Additional text 2, a UTF-16
  .doc), Itoigawa (Tables docx).
- **Curated claims**: compare each row of `literature/curated_claims.tsv` with its stated source.
  Grus & Zhang's sentence on sea lamprey V2R is confirmed in the paper (p. 411; `grus_zhang_2009/msn262_paper.pdf`).
- **Similarity and cross-reference matches**: rows whose `match` is "similar sequence" or ends in
  "(cross-reference table)".
- **Strain-matched genomes**: rows of `literature_comparison.tsv` whose `match` begins with "strain".
- **Hand decisions and sister genomes**: each row of `literature/genome_name_decisions.tsv` against
  its evidence, and each pairing in `literature_genome_proxies.tsv`.
- **Sequence matches by containment and by pieces**: a short literature sequence inside a longer
  protein of the same species is accepted; check a sample of rows whose `match` is "sequence
  contained" or "sequence pieces contained", and of `other_matching_accessions`.
- **Identifier matches**: check a sample of rows whose `match` begins with "identifier".
- **Spelling variants of species names**: rows of `literature_comparison.tsv` whose `match` is
  "name (spelling variant)".
- **Mapped-by-receptor families**: `literature_unmapped_families.tsv`.

## Known limitations

- "Not found in our data" covers receptors whose protein is not in the release we searched (a
  gene missing from the proteome, a different strain's proteome, a newer or retired gene model)
  and receptors that are in it but neither among our clades and search-tree inputs (by sequence)
  nor carrying an identifier our search headers show. Telling these apart needs the searched
  database itself (`Uniprot_Ref_Eukaryota_Proteome.fasta`, on the cluster).
- Liu et al. 2026 File S2 lists 749 receptors under two classes; each listing is traced.
- Near-identical matching can pick a close paralog (e.g. among olfactory receptors); the family is
  then the same, the protein may not be.
- "Not identified" in our data does not mean absent: proteome annotation gaps, tree size limits and
  our filters all remove real receptors.
- Claims without receptors cannot inform an origin: de Mendoza 2014 (Pfam counts), the abstract
  counts of the sponge study (its receptors are traced separately). Itoigawa (genomic loci) and
  Behrens (named genes) list no proteins; their loci are translated (`literature_deeper_origins.py`)
  and placed by their MMseqs2 best hit (`deeper_origin_best_hits.py`).
