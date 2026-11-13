"""Readers for the published GPCR surveys in literature/, one per study.

Step 6, literature validation. Each reader returns what the study itself reports, as two kinds
of record (plain dicts):

    claim     {study, their_family, genome, status, count, source}
              a statement about a genome: the family is present (with a count when the study
              gives one) or absent (a 0 or "not identified" in the study's own table or text)
    sequence  {study, their_family, genome, seq_id, sequence, ids}
              one receptor the study lists, with its sequence and/or identifiers, so it can be
              followed through our pipeline; genome is the species as the study names it, or ""
              when only the identifier tells it. A sequence taken from an alignment of genome
              contigs keeps "-" where the alignment breaks it into pieces.

Nothing a study does not report is produced here. compare_literature.py turns these records into
the comparison table.
"""

import csv
import io
import os
import re
import tempfile
import zipfile
from collections import defaultdict

import pandas as pd

LIT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "literature")
# The studies' supplementary files are not distributed (copyright); what each reader takes from them is,
# as one table per study in literature/parsed/ (written by write_parsed_literature.py). To read the
# originals instead, download them (literature/SOURCES.md) into a folder laid out like literature/
# and point GPCR_LITERATURE_ORIGINALS at it.
ORIGINALS_DIR = os.environ.get("GPCR_LITERATURE_ORIGINALS", LIT_DIR)
PARSED_DIR = os.path.join(LIT_DIR, "parsed")
READ_FILES = set()  # every supplementary file a reader opened (checked by verify_literature_comparison.py)

# Supplementary files no reader uses, and why
FILES_NOT_READ = {
    "grus_zhang_2009/msn262_paper.pdf": "the paper; its sea lamprey V2R statement is in curated_claims.tsv",
    "uniprot_accession_taxa.json":
        "read by compare_literature.py: taxID of each UniProt accession listed with no species (UniProtKB or UniSave)",
    "behrens_2023/behrens-et-al-2023-a-singular-shark-bitter-taste-receptor-provides-insights-into-the-evolution-of-bitter-taste.pdf":
        "the paper; its claims are in curated_claims.tsv",
    "behrens_2023/pnas.2310347120.sapp.pdf":
        "T2R sequences appear only as an unlabelled alignment; the genes are in curated_claims.tsv",
    "cardoso_2014/pone.0092220.s006.pdf": "Table S1, the s007 genes split by subfamily; no family-level information beyond s007",
    "garcia_de_las_bayonas_king_2025/SF1_proteome_accessions.xlsx": "the proteomes they searched, not receptors",
    "garcia_de_las_bayonas_king_2025/SF25_GPR157_EukProt_hits.xlsx": "search hits; the receptors are in SF25_GPR157_homologs.txt",
    "gotkhindikar_2025/kiaf057_supplementary_data.pdf": "structural analysis of AtGCR1; the one receptor is read directly",
    "itoigawa_2024/mmc1_supplement.pdf": "figures; the loci are in Tables_S1-S2.docx",
    "krishnan_2014_nematode/pone.0093048.s001.pdf": "tree figure",
    "krishnan_2014_nematode/pone.0093048.s002.pdf": "alignment figure",
    "krishnan_2014_nematode/pone.0093048.s003.pdf": "Pfam domain table",
    "krishnan_2014_nematode/pone.0093048.s004.pdf": "HHpred table",
    "krishnan_2014_nematode/pone.0093048.s006.docx": "Table S2 (vertebrate srsx genes, not a family claim)",
    "krishnan_2014_sponge/12862_2014_270_MOESM2_ESM.xlsx": "BLAST results for the receptors in MOESM1",
    "krishnan_2014_sponge/12862_2014_270_MOESM3_ESM.xlsx": "BLAST results for the receptors in MOESM1",
    "krishnan_2014_sponge/12862_2014_270_MOESM4_ESM.pdf": "tree figure",
    "krishnan_2014_sponge/12862_2014_270_MOESM5_ESM.pdf": "tree figure",
    "krishnan_2014_sponge/12862_2014_270_MOESM6_ESM.pdf": "tree figure",
    "krishnan_2014_sponge/12862_2014_270_MOESM7_ESM.pdf": "tree figure",
    "krishnan_2014_sponge/12862_2014_270_MOESM8_ESM.pdf": "tree figure",
}


def _path(*parts):
    """Our own files (curated claims, genome decisions) are in literature/, a study's file in ORIGINALS_DIR."""
    READ_FILES.add("/".join(parts))
    ours = os.path.join(LIT_DIR, *parts)
    return ours if os.path.exists(ours) else os.path.join(ORIGINALS_DIR, *parts)


def _require(path):
    """Publisher PDFs are not distributed with this repository (copyright); stop with a pointer."""
    if not os.path.exists(path):
        raise SystemExit("%s is not distributed with the repository; download it from the study's "
                         "supplementary material into that folder (literature/SOURCES.md)." % path)
    return path


def _pdf_text(path):
    _require(path)
    import pypdf
    import logging
    logging.getLogger("pypdf").setLevel(logging.ERROR)
    text = "\n".join(page.extract_text() or "" for page in pypdf.PdfReader(path).pages)
    return text.replace("\xa0", " ").replace("\xad", "-")


def _docx_text(data):
    xml = zipfile.ZipFile(io.BytesIO(data)).read("word/document.xml").decode("utf8")
    paragraphs = re.findall(r"<w:p[ >].*?</w:p>", xml, re.S)
    return "\n".join(re.sub(r"<[^>]+>", "", p) for p in paragraphs)


def _docx_tables(data):
    xml = zipfile.ZipFile(io.BytesIO(data)).read("word/document.xml").decode("utf8")
    tables = []
    for table in re.findall(r"<w:tbl>.*?</w:tbl>", xml, re.S):
        rows = []
        for row in re.findall(r"<w:tr[ >].*?</w:tr>", table, re.S):
            rows.append([re.sub(r"<[^>]+>", "", c).strip() for c in re.findall(r"<w:tc>.*?</w:tc>", row, re.S)])
        tables.append(rows)
    return tables


def parse_fasta_text(text, keep_gaps=False):
    """(header, sequence) pairs from FASTA-like text; sequence lines may be wrapped, may end in
    '*' or '.', and may carry alignment gaps (kept as '-' when keep_gaps, else removed)."""
    records, header, seq = [], None, []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq)))
            header, seq = line[1:].strip(), []
        elif header is not None and re.fullmatch(r"[A-Za-z*.\-]+", line):
            seq.append((re.sub(r"[.*]", "", line) if keep_gaps else re.sub(r"[-.*]", "", line)).upper())
        elif header is not None and seq:
            records.append((header, "".join(seq)))
            header, seq = None, []
    if header is not None:
        records.append((header, "".join(seq)))
    return [(h, s) for h, s in records if s.strip("-")]


CODONS = {a + b + c: aa for (a, b, c), aa in zip(
    [(a, b, c) for a in "TCAG" for b in "TCAG" for c in "TCAG"],
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}


def translate(dna):
    """Protein of a coding sequence read from its first base; internal stops (edited or
    pseudogene sequences) become '-' so each stop-free piece is matched on its own."""
    dna = re.sub(r"[^ACGT]", "N", dna.upper())
    protein = "".join(CODONS.get(dna[i:i + 3], "X") for i in range(0, len(dna) - 2, 3))
    return protein.rstrip("*").replace("*", "-")


def _claim(study, family, genome, status, count, source, evidence=""):
    """evidence: what backs a claim that lists no receptor sequence or ID (e.g. genomic loci)."""
    return {"kind": "claim", "study": study, "their_family": family, "genome": genome.strip(),
            "status": status, "count": count, "source": source, "evidence": evidence}


def _seq(study, family, genome, seq_id, sequence="", ids=()):
    return {"kind": "sequence", "study": study, "their_family": family, "genome": (genome or "").strip(),
            "seq_id": seq_id, "sequence": sequence, "ids": [i for i in ids if i]}


# ---------------------------------------------------------------- species from identifiers

PREFIX_SPECIES = {
    "PTSG": "Salpingoeca rosetta", "Monbr1": "Monosiga brevicollis", "Cowc": "Capsaspora owczarzaki",
    "CAOG": "Capsaspora owczarzaki", "Pgem": "Pirum gemmata", "Awhi": "Abeoforma whisleri",
    "Clim": "Corallochytrium limacisporum", "Sarc": "Sphaeroforma arctica", "Mvib": "Ministeria vibrans",
    "Pvie": "Pigoraptor vietnamica", "Pchi": "Pigoraptor chileana", "CFRG": "Creolimax fragrantissima",
    "Colp12": "Syssomonas multiformis", "Partr": "Parvularia atlantis", "Brafl1": "Branchiostoma floridae",
    "AMAG": "Allomyces macrogynus", "SPPG": "Spizellomyces punctatus", "BDEG": "Batrachochytrium dendrobatidis",
    "BDET": "Batrachochytrium dendrobatidis", "RO3G": "Rhizopus delemar",
}


def species_from_header(header):
    """The species a FASTA header names, if any."""
    m = re.search(r"OS=([A-Z][a-z]+[ _][a-z]+)", header)
    if m:
        return m.group(1).replace("_", " ")
    m = re.search(r"\[([A-Z][a-z]+)[ _]([a-z]+)", header)
    if m:
        return "%s %s" % m.groups()
    m = re.match(r"EP\d+_([A-Z][a-z]+)_([a-z]+)_", header)
    if m:
        return "%s %s" % m.groups()
    m = re.match(r"([A-Z][a-z]+)[_ ]([a-z]+)_", header)
    if m:
        return "%s %s" % m.groups()
    m = re.search(r"\|\s*([A-Z][a-z]+ [a-z]+)", header)
    if m:
        return m.group(1)
    for prefix, species in PREFIX_SPECIES.items():
        if re.search(r"(^|[|_ ])%s([|_ ]|$)" % re.escape(prefix), header):
            return species
    for pattern, species in PATTERN_SPECIES:
        if re.match(pattern, header):
            return species
    return ""


# Gene-model names that identify the species by their form
PATTERN_SPECIES = [
    (r"Sarc\d+_g", "Sphaeroforma arctica"), (r"CFRG\d", "Creolimax fragrantissima"),
    (r"ML\d{4,}[a-z]", "Mnemiopsis leidyi"), (r"T_wilhelma_", "Tethya wilhelma"),
    (r"Em\d{4}g\d+", "Ephydatia muelleri"), (r"FBpp\d", "Drosophila melanogaster"),
]

ACCESSION = re.compile(r"\b([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})\b")
ENTRY_NAME = re.compile(r"\b([A-Z0-9]{1,10}_[A-Z][A-Z0-9]{2,4})\b")
LOCUS = re.compile(r"\b([A-Z][A-Z0-9]{2,10}_\d{4,6})(?:T\d)?\b")


def ids_from_header(header):
    """Identifiers in a header that UniProt can resolve: accessions, entry names, locus tags,
    RefSeq proteins."""
    ids = []
    m = re.match(r"(?:tr|sp)\|([A-Z0-9]+)\|", header)
    if m:
        ids.append(m.group(1))
    ids += re.findall(r"\b([NX]P_\d+\.?\d*)", header)
    ids += [t for t in LOCUS.findall(header) if not t.startswith(("XP_", "NP_"))]
    m = re.match(r"jgi\|Monbr1\|(\d+)\|", header)
    if m:
        ids.append("MONBRDRAFT_%s" % m.group(1))
    return list(dict.fromkeys(ids))


# ---------------------------------------------------------------- one reader per study

GARCIA = "Garcia De Las Bayonas & King 2025"
GARCIA_HEADINGS = {"Glutamate Receptors": "Glutamate", "Rhodopsins": "Rhodopsin",
                   "Adhesion Receptors": "Adhesion", "cAMP Receptors": "cAMP"}


# Garcia headers that do not name the species. TRINITY: the Choanoeca flexa proteome (Brunet et al.
# 2019, figshare 8216291; all 19 sequences are in it). genemark-NKKH001C: the Barroeca monosierra
# genome-derived proteome (Hake et al. 2024, figshare 14474214, file NKKH001C_Barroeca_monosierra).
# Nk52: under their Ichthyosporea heading; of their six ichthyosporeans (Table 2) the other four
# with gene models have their own prefixes (Sarc4, Awhi, CFRG, Pgem), so it is Chromosphaera perkinsii
# or Ichthyophonus hoferi, neither in our genome set.
# Garcia entries named only "MTH-1" and "MTH-2" (Adhesion file, Metazoa section): their sequences are
# identical to Caenorhabditis elegans proteins in UniParc (UPI0000225AC9 = Q19931; UPI0000225ACA =
# Q19930; looked up by sequence checksum, October 3, 2026)
GARCIA_UNIPARC = {"MTH-1": ("Caenorhabditis elegans", ["Q19931"]), "MTH-2": ("Caenorhabditis elegans", ["Q19930"])}
GARCIA_HEADER_SPECIES = [
    (r"TRINITY_DN", "Choanoeca flexa"),
    (r"genemark-NKKH001C", "Barroeca monosierra"),
    (r"Nk52_", "Ichthyosporea sp. Nk52 (Chromosphaera perkinsii or Ichthyophonus hoferi)"),
]


def _garcia_file(name, family):
    """Supplementary files 9 and 25 have lineage headings ("Amoebozoa:"); file 4 has family
    headings ("cAMP Receptors:"), read when family is None."""
    records, section, fam = [], "", family
    text = open(_path("garcia_de_las_bayonas_king_2025", name), encoding="utf-8", errors="ignore").read()
    chunks = re.split(r"^([A-Za-z0-9 _/-]+):\s*$", text, flags=re.M)
    for i in range(0, len(chunks)):
        if i % 2 == 1:
            heading = chunks[i].strip()
            if family is None:
                fam = GARCIA_HEADINGS.get(heading, heading)
            else:
                section = heading
            continue
        for header, seq in parse_fasta_text(chunks[i]):
            species = species_from_header(header) or next(
                (sp for pattern, sp in GARCIA_HEADER_SPECIES if re.match(pattern, header)), "")
            m = re.match(r"EP\d{5}_([A-Za-z0-9-]+)_sp_", header)  # EukProt lineages with no species name
            species = species or (m.group(1) + " sp." if m else "")
            ids = ids_from_header(header)
            if header.split()[0] in GARCIA_UNIPARC:
                species, ids = GARCIA_UNIPARC[header.split()[0]]
            records.append(_seq(GARCIA, fam, species, header[:80], seq, ids))
    return records


def read_garcia():
    records = []
    for name, family in [("Rhodopsin_Receptors.txt", "Rhodopsin"), ("Glutamate_Receptors.txt", "Glutamate"),
                         ("Adhesion_Receptors.txt", "Adhesion"), ("Frizzled_Receptors.txt", "Frizzled"),
                         ("SF4_choanoflagellate_GPCRs.txt", None), ("SF25_GPR157_homologs.txt", "GPR157")]:
        records += _garcia_file(name, family)
    # The choanoflagellate receptors appear in both file 4 and file 9; count each once.
    seen, unique = set(), []
    for r in records:
        key = (r["their_family"], r["sequence"])
        if key not in seen:
            seen.add(key)
            unique.append(r)
    return unique


KRISHNAN12 = "Krishnan et al. 2012"


def read_krishnan_2012():
    records = []
    text = open(_path("krishnan_2012", "pone.0029817.s001.txt"), encoding="utf-8", errors="ignore").read()
    # Dataset S1: headings "Rhodopsin (7tm_1)", "Adhesion (7tm_2)", "Glutamate (7tm_3)", "cAMP";
    # wrapped header lines are joined to the header.
    lines, fam, joined = text.splitlines(), None, []
    for line in lines:
        s = line.strip()
        m = re.fullmatch(r"(Rhodopsin|Adhesion|Glutamate|Frizzled|cAMP)( \(7tm_\d\))?", s)
        if m:
            fam = m.group(1)
            joined.append("@@%s" % fam)
        elif s.startswith(">") or re.fullmatch(r"[A-Z*]+", s) or not s:
            joined.append(s)
        elif joined and joined[-1].startswith(">"):
            joined[-1] += " " + s
    # Table S1 (family, species, accession) is the authoritative list: it gives the genome
    # claims, and a Dataset S1 sequence is kept only if Table S1 lists it (Dataset S1 has one
    # Glutamate sequence, AMAG_16966, that Table S1 does not). Table S1 also lists the six
    # Frizzled sequences that Dataset S1 does not contain.
    t1 = pd.read_excel(_path("krishnan_2012", "pone.0029817.s002.xls"), header=None).iloc[4:, :3]
    t1.columns = ["family", "species", "acc"]
    t1[["family", "species"]] = t1[["family", "species"]].ffill()
    t1 = t1.dropna(subset=["acc"])
    t1["family"] = t1["family"].astype(str).str.split("(").str[0].str.strip()
    t1["acc"] = t1["acc"].astype(str).str.strip()
    listed = defaultdict(set)
    for _, row in t1.iterrows():
        listed[row["family"]].add(row["acc"])
    for (family, species), n in t1.groupby(["family", "species"]).size().items():
        records.append(_claim(KRISHNAN12, family, str(species), "present", int(n), "Table S1 (accessions per species)"))
    fam, block = None, []
    for s in joined + ["@@END"]:
        if s.startswith("@@"):
            for header, seq in parse_fasta_text("\n".join(block)):
                if any(acc in header for acc in listed.get(fam, ())):
                    prefix = header.split("_")[0].lstrip(">")
                    records.append(_seq(KRISHNAN12, fam, species_from_header(header) or KRISHNAN12_LOCUS_SPECIES.get(prefix, ""),
                                        header[:80], seq, ids_from_header(header)))
            fam, block = s[2:], []
        else:
            block.append(s)
    # Table S1 accessions with no Dataset S1 sequence (all six Frizzled, and a few whose FASTA
    # header carries another identifier) are traced by the accession alone.
    in_fasta = {(r["their_family"], acc) for r in records if r["kind"] == "sequence"
                for acc in listed.get(r["their_family"], ()) if acc in r["seq_id"]}
    for _, row in t1.iterrows():
        if (row["family"], row["acc"]) not in in_fasta:
            records.append(_seq(KRISHNAN12, row["family"], row["species"], row["acc"], "", [row["acc"]]))
    # Dataset S1 repeats two sequences (AMAG_13197, B6HLF9) and Table S1 one accession; each
    # receptor is traced once, while the claims keep Table S1's counts as published.
    seen, unique = set(), []
    for r in records:
        key = (r["kind"], r["their_family"], r.get("seq_id"), r.get("genome") if r["kind"] == "claim" else "")
        if key not in seen:
            seen.add(key)
            unique.append(r)
    records = unique
    # Table S2: the species they searched. Table S3 lists all fungus-specific GPCRs they found, so a
    # fungal species of Table S2 is searched for each Table S3 family; the claim is "searched" and
    # becomes present or absent once their receptors are placed in genomes (compare_literature.py).
    # (Table S1 lists only their novel sequences, so it gives no absences.)
    for species in krishnan12_fungi_searched():
        for family in ("STE2", "STE3", "Git3", "Nutrient receptor"):
            records.append(_claim(KRISHNAN12, family, species, "searched", "", "Table S2 (species searched)"))
    # Table S3: fungus-specific GPCRs as UniProt entry names or locus tags, by family column.
    t3 = pd.read_excel(_path("krishnan_2012", "pone.0029817.s004.xls"), header=None)
    families = t3.iloc[1].ffill()
    names = {"Ste2": "STE2", "Ste3": "STE3", "GIT3/GIT3C": "Git3", "Nutrient receptors": "Nutrient receptor"}
    for col in range(t3.shape[1]):
        family = names.get(str(families[col]).strip())
        for acc in t3.iloc[3:, col].dropna():
            acc = str(acc).strip()
            if family and acc:
                # Broad locus tags: some cells give the transcript form (SNOT_10426 for gene
                # SNOG_10426); the prefix names the genome (Table S2 lists the species searched)
                gene = re.sub(r"^(SNO|PGU)T_", r"\1G_", acc)
                prefix = gene.split("_")[0] if re.fullmatch(r"[A-Z0-9]+_\d+(T\d)?", gene) else ""
                ids = list(dict.fromkeys([acc, gene]))
                records.append(_seq(KRISHNAN12, family, KRISHNAN12_LOCUS_SPECIES.get(prefix, ""), acc, "", ids))
    return records


def krishnan12_fungi_searched():
    """The fungal species of Krishnan 2012 Table S2 (a Word .doc; text runs read as Latin-1): the
    lines between the "Ascomycota" heading and the first non-fungal group."""
    raw = open(_path("krishnan_2012", "pone.0029817.s003.doc"), "rb").read().decode("latin-1", errors="ignore")
    lines = [r.strip() for r in re.findall(r"[A-Za-z0-9 .,()/_\-]{6,}", raw) if re.search(r"[a-z]{3}", r)]
    groups = {"Ascomycota", "Basidiomycota", "Microspordia", "Microsporidia", "Zygomycota", "Blastocladiomycota",
              "Chytridiomycota"}
    species, inside = [], False
    for line in lines:
        if line in groups:
            inside = True
        elif line in ("Filasterea", "Choanoflagellata", "Alveolata"):
            break
        elif inside and re.fullmatch(r"[A-Z][a-z]+ [a-z][a-z.-]+( .*)?", line):
            species.append(line)
    return species


# Locus-tag prefixes of Krishnan 2012 Table S3 (Broad Institute gene models). Each taxID is the one
# the prefix carries in our search headers (GN=), or, for FGSG, SNOG, SOCG, TEQG and TESG, which our
# headers do not carry, the organism NCBI Protein gives for the locus tag; names are the NCBI
# scientific names. All are species of Krishnan's Table S2 (pone.0029817.s003.doc).
KRISHNAN12_LOCUS_SPECIES = {
    "AMAG": "Allomyces macrogynus ATCC 38327",  # taxID 578462
    "BDEG": "Batrachochytrium dendrobatidis JEL423",  # taxID 403673
    "CIMG": "Coccidioides immitis RS",  # taxID 246410
    "FGSG": "Fusarium graminearum PH-1",  # taxID 229533
    "FOXG": "Fusarium oxysporum f. sp. lycopersici 4287",  # taxID 426428
    "FVEG": "Fusarium verticillioides 7600",  # taxID 334819
    "MGG": "Pyricularia oryzae 70-15",  # taxID 242507
    "MGYG": "Nannizzia gypsea CBS 118893",  # taxID 535722
    "PGUG": "Meyerozyma guilliermondii ATCC 6260",  # taxID 294746
    "RO3G": "Rhizopus delemar RA 99-880",  # taxID 246409
    "SNOG": "Parastagonospora nodorum SN15",  # taxID 321614
    "SOCG": "Schizosaccharomyces octosporus yFS286",  # taxID 483514
    "SPPG": "Spizellomyces punctatus DAOM BR117",  # taxID 645134
    "TEQG": "Trichophyton equinum CBS 127.97",  # taxID 559882
    "TERG": "Trichophyton rubrum CBS 118892",  # taxID 559305
    "TESG": "Trichophyton tonsurans CBS 112818",  # taxID 647933
    "VDAG": "Verticillium dahliae VdLs.17",  # taxID 498257
}


def read_de_mendoza_2014():
    study = "de Mendoza et al. 2014"
    table = pd.read_excel(_path("de_mendoza_2014", "supp_evu038_SupplementaryTable1.xls"), header=None)
    domains = list(table.iloc[0])
    records = []
    for _, row in table.iloc[2:].iterrows():
        species = re.sub(r"\s*\(.*\)", "", str(row[0])).strip()
        if not species or species == "nan":
            continue
        for col, domain in enumerate(domains):
            if col == 0 or str(domain) == "nan":
                continue
            count = int(row[col]) if str(row[col]) not in ("nan", "") else 0
            records.append(_claim(study, str(domain).strip(), species, "present" if count else "absent", count,
                                  "Supplementary Table 1 (PfamScan domain count)"))
    return records


DE_MENDOZA_2016_CODES = {
    "Adig": "Acropora digitifera", "Avag": "Adineta vaga", "Aaeg": "Aedes aegypti", "Aque": "Amphimedon queenslandica",
    "Amel": "Apis mellifera", "Bmor": "Bombyx mori", "Bbel": "Branchiostoma belcheri", "Cele": "Caenorhabditis elegans",
    "Cmil": "Callorhinchus milii", "Ctel": "Capitella teleta", "Cowc": "Capsaspora owczarzaki",
    "Cint": "Ciona intestinalis", "Cgig": "Crassostrea gigas", "Cfra": "Creolimax fragrantissima",
    "Drer": "Danio rerio", "Dpul": "Daphnia pulex", "Dmel": "Drosophila melanogaster", "Hrob": "Helobdella robusta",
    "Hsap": "Homo sapiens", "Hmag": "Hydra magnipapillata", "Lgig": "Lottia gigantea", "Mvib": "Ministeria vibrans",
    "Mlei": "Mnemiopsis leidyi", "Mbre": "Monosiga brevicollis", "Nvec": "Nematostella vectensis",
    "Ocar": "Oscarella carmela", "Lamp": "Petromyzon marinus", "Skow": "Saccoglossus kowalevskii",
    "Sros": "Salpingoeca rosetta", "Sman": "Schistosoma mansoni", "Sarc": "Sphaeroforma arctica",
    "Smar": "Strigamia maritima", "Spur": "Strongylocentrotus purpuratus", "Scil": "Sycon ciliatum",
    "Turt": "Tetranychus urticae", "Tcas": "Tribolium castaneum", "Tspi": "Trichinella spiralis",
    "Tadh": "Trichoplax adhaerens",
    # sequences from species outside Table S1, named in the sequence file by these codes
    "Dvir": "Drosophila virilis", "Agam": "Anopheles gambiae", "Bflo": "Branchiostoma floridae",
    "Mdom": "Musca domestica", "Ccap": "Ceratitis capitata", "Bdor": "Bactrocera dorsalis",
}


def read_de_mendoza_2016():
    study = "de Mendoza et al. 2016"
    text = _pdf_text(_path("de_mendoza_2016", "srep21801-s1.pdf"))
    end = text.find("Table S1:")
    records = []
    for header, seq in parse_fasta_text(text[:end]):
        code = re.split(r"[_ ]", header)[0]
        genome = DE_MENDOZA_2016_CODES.get(code) or species_from_header(header)
        records.append(_seq(study, "Mth", genome, header[:80], seq, ids_from_header(header)))
    return records


def read_nordstrom_2009():
    study = "Nordstrom et al. 2009"
    z = zipfile.ZipFile(_path("nordstrom_2009", "msn228_Supplementary_Data.zip"))
    table = pd.read_excel(io.BytesIO(z.read("mbe-08-0416-File003_msn228.xls")), header=0)
    species_of = dict(zip(table["Transcript"].astype(str), table["Specie"].astype(str)))
    original = dict(zip(table["Transcript"].astype(str), table["Original sequence"].astype(str)))
    records = []
    # File 2 holds coding sequences (nucleotides); they are translated for matching.
    for header, seq in parse_fasta_text(z.read("mbe-08-0416-File002_msn228.fasta").decode("latin-1")):
        seq = translate(seq)
        name = header.split()[0]
        key = re.sub(r"_edited$", "", name)
        ids = ids_from_header(original.get(key, "")) + [t for t in re.findall(r"ENS\w*P\d+", original.get(key, ""))]
        # Monosiga (Mb_) and Drosophila (CG...) sequences are named by gene, not by transcript
        species = species_of.get(key) or ("M. brevicollis" if key.startswith("Mb_") else
                                          "D. melanogaster" if key.startswith("CG") else "")
        if key.startswith("Mb_"):
            ids.append("MONBRDRAFT_%s" % key.split("_")[1])
        if key.startswith("CG"):
            ids.append(key.split("-")[0])
        records.append(_seq(study, "Adhesion/Secretin", species, name, seq, ids))
    return records


def read_nordstrom_2011():
    study = "Nordstrom et al. 2011"
    z = zipfile.ZipFile(_path("nordstrom_2011", "msr061_Supplementary_Data.zip"))
    table = pd.read_excel(io.BytesIO(z.read("Additional_table_1.xls")), header=0)
    counts = pd.crosstab(table["Species"], table["Family"])
    records = []
    for species, row in counts.iterrows():
        for family, count in row.items():
            records.append(_claim(study, family, species, "present" if count else "absent", int(count),
                                  "Additional table 1 (sequences per species and family)"))
    # Additional text 2 (a UTF-16 Word file) holds the sequence of every receptor in the table.
    text = z.read("Additional_Text_2.doc").decode("utf-16le", errors="ignore")
    # Word stores the text in pieces with binary data between them (nV_35931 is split 72 + 191 aa);
    # removing the non-text characters joins each record back together; the binary tail after the
    # last record is cut first
    last = text.rfind(">")
    end = re.search(r"[^\x20-\x7e\r\n]", text[last:])
    text = re.sub(r"[^\x20-\x7e\r\n]+", "", text[:last + end.start()] if end else text)
    sequences = {name: re.sub(r"[^A-Z]", "", seq) for name, seq in
                 re.findall(r">([A-Za-z]{2}_[\w.]+)[^\r\n]*[\r\n]+([A-Z\r\n*-]{20,})", text)}
    for _, row in table.iterrows():
        records.append(_seq(study, row["Family"], row["Species"], row["Sequence"], sequences.get(row["Sequence"], "")))
    return records


def read_krishnan_2014_nematode():
    study = "Krishnan et al. 2014 (nematode)"
    z_text = _docx_text(open(_path("krishnan_2014_nematode", "pone.0093048.s005.docx"), "rb").read())
    prefix_species = {}
    for species, prefix in re.findall(r"([A-Z][a-z]+ [a-z]+)\s+([A-Z][a-z])\s?_", z_text):
        prefix_species.setdefault(prefix, species)
    # Aplysia genes are named by number with their contig in parentheses, e.g. "5a (AASC01117380)".
    srw_ids = set(re.findall(r"([A-Z][a-z]_\w+|\d+[abc] \(AASC\d+\))", z_text))
    records = []
    for header, seq in parse_fasta_text(_pdf_text(_path("krishnan_2014_nematode", "pone.0093048.s007.pdf"))):
        name = re.match(r"(\d+[abc] \(AASC\d+\)|\S+)", header).group(1)
        # s007 also holds the vertebrate srsx-domain sequences of their Table S2; only the srw list
        # of Table S1 is a family claim.
        if name not in srw_ids:
            continue
        species = "Aplysia californica" if "AASC" in name else prefix_species.get(name.split("_")[0], "")
        inner = re.search(r"\(([^)]+)\)", header)
        ids = [re.sub(r"-P[A-Z]$", "", inner.group(1))] if inner and "AASC" not in name else []
        records.append(_seq(study, "srw", species, name, seq, ids))
    return records


SPONGE = "Krishnan et al. 2014 (sponge)"
SPONGE_HEADINGS = {"Glutamate": "Glutamate", "Rhodopsin": "Rhodopsin", "Adhesion": "Adhesion", "Frizzled": "Frizzled",
                   "cAMP": "cAMP", "Ocular albinism (GPR143)": "OA1",
                   "Intimal thickness-related receptor": "Intimal thickness-related receptor",
                   "Lung 7TM receptor": "Lung 7TM receptor (GPR108)"}


def read_krishnan_2014_sponge():
    """Additional file 1: the Amphimedon queenslandica receptors (Aqu1 gene models) as FASTA, under
    family headings ("Glutamate like GPCRs", ...). The genome counts of the abstract stay in
    curated_claims.tsv."""
    text = _pdf_text(_path("krishnan_2014_sponge", "12862_2014_270_MOESM1_ESM.pdf"))
    # page numbers fall inside sequences at page breaks; a line holding only a number is dropped
    text = re.sub(r"\n\s*\d+\s*(?=\n)", "", text)
    pattern = r"\n\s*(%s) like GPCRs[^\n]*\n" % "|".join(re.escape(h) for h in SPONGE_HEADINGS)
    parts = re.split(pattern, text)
    records = []
    for heading, block in zip(parts[1::2], parts[2::2]):
        for header, seq in parse_fasta_text(block):
            gene = re.search(r"gene:(\S+)", header)
            records.append(_seq(SPONGE, SPONGE_HEADINGS[heading], "Amphimedon queenslandica", header.split()[0], seq,
                                [gene.group(1)] if gene else []))
    return records


GRUS = "Grus & Zhang 2009"
# The V1R sequences of Supplementary File 9 named in curated_claims.tsv (lamprey contigs and
# elephant shark WGS contigs); the file's other alignments are V2R and T1R genes of other species.
GRUS_V1R = {"Contig15650.3": "Petromyzon marinus", "Contig26060.4": "Petromyzon marinus",
            "Contig59695.1": "Petromyzon marinus", "AAVX01055762.1": "Callorhinchus milii",
            "AAVX01177037.1": "Callorhinchus milii"}


def _grus_pdf(name):
    z = zipfile.ZipFile(_require(_path("grus_zhang_2009", "msn262_Supplementary_Data.zip")))
    path = os.path.join(tempfile.gettempdir(), name)
    with open(path, "wb") as handle:
        handle.write(z.read(name))
    return _pdf_text(path)


def read_grus_zhang_2009():
    """Translated genome contigs: File 2 (File008 pdf), 32 partial elephant shark V2Rs as an
    alignment; File 9, the V1R alignment (interleaved blocks). Alignment gaps are kept, since a gap
    in a contig translation can stand for a missing exon."""
    records = []
    for header, seq in parse_fasta_text(_grus_pdf("mbe-08-0363-File008_msn262.pdf"), keep_gaps=True):
        records.append(_seq(GRUS, "V2R", "Callorhinchus milii", header.split()[0], seq))
    blocks = defaultdict(str)
    for name, seq in re.findall(r"^(\S+)\s+([A-Z\-]{10,})\s*$", _grus_pdf("mbe-08-0363-File009_msn262.pdf"), re.M):
        if name in GRUS_V1R:
            blocks[name] += seq
    for name, seq in blocks.items():
        records.append(_seq(GRUS, "V1R", GRUS_V1R[name], name, seq))
    return records


def read_cardoso_2014():
    study = "Cardoso et al. 2014"
    text = _pdf_text(_path("cardoso_2014", "pone.0092220.s007.pdf"))
    records, species = [], None
    for line in text.splitlines():
        m = re.match(r"\s*([A-Z]\.\s+[a-zó]+)\s+(\S+)\s+\S+", line)
        if m:
            species = re.sub(r"\s+", " ", m.group(1))
            gene = m.group(2)
        else:
            m = re.match(r"\s+(\S+)\s+\S+\s*$", line)
            if not (m and species):
                continue
            gene = m.group(1)
        if gene in ("Genes",) or not re.search(r"\d", gene):
            continue
        records.append(_seq(study, "Secretin", species, gene, "", [gene]))
    # s014: their alignment of these receptors (transmembrane region), named species code + gene
    # or acronym, with spelling variants ("HMELO14214" for HMEL014214); joined by a normalised key.
    key = lambda name: re.sub(r"[^A-Z0-9]", "", name.upper()).replace("O", "0")
    acronyms = dict((key(a), g) for g, a in re.findall(r"^\s*(?:[A-Z]\.\s+[a-zó]+\s+)?(\S+)\s+(\S+)\s*$", text, re.M))
    by_key = {key(r["seq_id"]): r for r in records}
    for header, seq in parse_fasta_text(open(_path("cardoso_2014", "pone.0092220.s014.fasta"), encoding="utf-8",
                                             errors="ignore").read(), keep_gaps=True):
        name = re.sub(r"^[A-Z][a-z]{2}_?", "", header.split()[0])
        r = by_key.get(key(name)) or by_key.get(key(acronyms.get(key(name), "")))
        if r:
            r["sequence"] = seq
    return records


def read_prabhu_2006():
    study = "Prabhu & Eichinger 2006"
    text = _pdf_text(_path("prabhu_eichinger_2006", "1-s2.0-S0171933506000720-main.pdf"))
    t2 = text.find("Table 2. The Dictyostelium family of Frizzled-like")
    records = []
    for i, family in [(text.find("Table 1. The Dictyostelium family of GABA"), "GABAB"), (t2, "Frizzled")]:
        block = text[i:i + 2500]
        # UniProt indexes dictyBase gene IDs (DDB_G...) but not these DDB0 protein IDs, so the
        # gene name in the next column (grlA, fslA, ...) is what is looked up.
        for ddb, gene in dict.fromkeys(re.findall(r"(DDB\d{7})\*?\s+([a-z]{3}[A-Z])\b", block)):
            records.append(_seq(study, family, "Dictyostelium discoideum", "%s %s" % (ddb, gene), "", [gene]))
    return records


def read_itoigawa_2024():
    study = "Itoigawa et al. 2024"
    tables = _docx_tables(open(_path("itoigawa_2024", "Tables_S1-S2.docx"), "rb").read())
    searched = [r[0] for r in tables[0][1:]]
    loci = {}
    for row in tables[1][1:]:
        loci.setdefault(row[0], []).append(row[2])
    records = []
    for species in searched:
        genes = loci.get(species, [])
        intact = [g for g in genes if not g.endswith("p")]
        source = "Table S1 (searched) and Table S2 (%s)" % (", ".join(genes) if genes else "no TAS2R locus")
        if genes and not intact:
            source += "; pseudogenes only"
        records.append(_claim(study, "T2R", species, "present" if genes else "absent", len(genes), source,
                              "genomic loci (contig and coordinates), no protein" if genes else ""))
    return records


KOWATSCHEW_CODES = {
    "Amra": "Amblyraja radiata", "Cami": "Callorhinchus milii", "Cm": "Callorhinchus milii", "Dare": "Danio rerio",
    "Dr": "Danio rerio", "Gaga": "Gallus gallus", "Hs": "Homo sapiens", "Icpu": "Ictalurus punctatus",
    "Lach": "Latimeria chalumnae", "Leoc": "Lepisosteus oculatus", "Lo": "Lepisosteus oculatus",
    "Leth": "Lethenteron camtschaticum", "Mm": "Mus musculus", "Onki": "Oncorhynchus kisutch",
    "Pema": "Petromyzon marinus", "Pm": "Petromyzon marinus", "Rat": "Rattus norvegicus",
    "Rhin": "Rhincodon typus", "Rh": "Rhincodon typus", "Scfo": "Scleropages formosus",
    "Scan": "Scyliorhinus canicula", "Sqac": "Squalus acanthias", "Brbe": "Branchiostoma belcheri",
    "Bb": "Branchiostoma belcheri", "Brfl": "Branchiostoma floridae", "Bf": "Branchiostoma floridae",
    "Ciin": "Ciona intestinalis", "Ci": "Ciona intestinalis", "Ptfl": "Ptychodera flava", "Pt": "Ptychodera flava",
    "Sako": "Saccoglossus kowalevskii", "Sa": "Saccoglossus kowalevskii", "Stpu": "Strongylocentrotus purpuratus",
    "St": "Strongylocentrotus purpuratus", "Cael": "Caenorhabditis elegans", "Cae": "Caenorhabditis elegans",
    "Drme": "Drosophila melanogaster", "Dm": "Drosophila melanogaster",
}


def read_kowatschew_2022():
    study = "Kowatschew & Korsching 2022"
    text = _pdf_text(_path("kowatschew_korsching_2022", "bjac007_suppl_supplementary_material.pdf"))
    start = text.find("Predicted amino acid sequences for lamprey")
    end = text.find("SI file 4", start)
    records = []
    for header, seq in parse_fasta_text(text[start:end if end > 0 else None]):
        name = header.split()[0]
        code = re.split(r"[-_.]", name)[0]
        upper = name.upper()
        if "V2R" in upper or "OLFC" in upper:
            family = "V2R"
        elif "V1R" in upper or "ORA" in upper:
            family = "V1R"
        elif "T2R" in upper:
            family = "T2R"
        else:
            continue  # outgroup genes (other class A and class C receptors), not a family claim
        codes = {k.lower(): v for k, v in KOWATSCHEW_CODES.items()}
        genome = codes.get(code.lower(), "")
        if genome not in ("Lethenteron camtschaticum", "Petromyzon marinus"):
            continue  # SI file 2 lists these as "reference v1r/ora genes, and outgroup genes", not lamprey receptors
        records.append(_seq(study, family, genome, name, seq, []))
    return records


# A File S2 portal spelled differently in Table S3 (Trichoderma virens Gv29-8, written "G9-8" in
# Table S3; 7 of its receptors are identical to our T. virens proteins)
LIU_PORTAL_ALIASES = {"TriviGv29_8_2": "TriviG9_8_2"}
# File S2 portals named by their receptors, each with receptors identical to proteins of ours in this
# species: Amyenc1 (1), absent from Table S3; Aspnov1 (11), which Table S3 lists shortened as "Aspno"
# (A. novofumigatus IBT 16806, reached first through table_portal; listed here so the species is
# stated whatever Table S3 says - not Aspnovo1, A. novoparasiticus)
LIU_PORTAL_SPECIES = {"Amyenc1": "Amylocarpus encephaloides", "Aspnov1": "Aspergillus novofumigatus"}


def read_liu_2026():
    study = "Liu et al. 2026"
    z = zipfile.ZipFile(_path("liu_2026", "jof-12-00030-s001.zip"))
    table = pd.read_excel(io.BytesIO(z.read("Supplementary file S1.xlsx")), sheet_name="Table S3", header=1)
    # File S2: every identified receptor with its class and sequence ([portal ID]|[protein ID]).
    # Table S3's counts disagree with the receptors of File S2 in 4-8% of genome x class cells,
    # in both directions, so a genome's claims come from its File S2 receptors; Table S3 is used
    # only for the 35 genomes that File S2 does not contain.
    species_of_portal = {str(r["MycoCosm portal ID"]): str(r["Species name"]) for _, r in table.iterrows()}
    seqs = pd.read_excel(io.BytesIO(z.read("Supplementary file S2.xlsx")), header=None).iloc[2:]
    records, in_s2 = [], defaultdict(lambda: defaultdict(int))

    def table_portal(portal):
        """File S2 portal -> Table S3 portal. Table S3 shortens some (Annbo for Annbov1, Tubme1 for
        Tubme1v2): the longest Table S3 portal of >= 4 characters that begins File S2's, if only one
        does; two are spelled differently and are listed in LIU_PORTAL_ALIASES."""
        if portal in species_of_portal:
            return portal
        if portal in LIU_PORTAL_ALIASES:
            return LIU_PORTAL_ALIASES[portal]
        heads = [p for p in species_of_portal if len(p) >= 4 and portal.startswith(p)]
        longest = [p for p in heads if len(p) == max(map(len, heads))] if heads else []
        return longest[0] if len(longest) == 1 else portal

    for seq_id, cls, seq in seqs.itertuples(index=False):
        if not isinstance(seq, str) or not re.fullmatch(r"C?\d+(\.0)?", str(cls).strip()):
            continue  # title and header rows
        portal = table_portal(str(seq_id).split("|")[0])
        cls = str(cls).strip()
        family = cls if cls.startswith("C") else "C%s" % int(float(cls))
        in_s2[portal][family] += 1
        species = species_of_portal.get(portal) or LIU_PORTAL_SPECIES.get(portal, "")
        records.append(_seq(study, family, species, str(seq_id), re.sub(r"[^A-Za-z]", "", seq).upper()))
    for _, row in table.iterrows():
        portal, species = str(row["MycoCosm portal ID"]), str(row["Species name"])
        for c in range(1, 18):
            if portal in in_s2:
                count, source = in_s2[portal]["C%d" % c], "File S2 (receptors listed)"
            else:
                count, source = int(row["C%d" % c]), "Table S3 count; genome not in File S2"
            records.append(_claim(study, "C%d" % c, species, "present" if count else "absent", count, source))
    return records


def read_fernandez_figueroa_2026():
    study = "Fernandez-Figueroa et al. 2026"
    path = _path("fernandez_figueroa_2026", "media-1_SuppTables1-2.xlsx")
    t1 = pd.read_excel(path, sheet_name="Supp Table 1", header=0)
    records = []
    for _, row in t1.dropna(subset=["Uniprot code"]).iterrows():
        kind = str(row["GPCR type"]).strip()
        label = str(row["Tip label"]).strip()
        if kind in ("GPCR plant candidate", "non"):
            family = "GCR1 (canonical GPCR)" if label == "GCR1" else "plant candidate, not a canonical GPCR"
            genome = "Arabidopsis thaliana"
        else:
            family = kind
            genome = ""
        records.append(_seq(study, family, genome, label, "", [str(row["Uniprot code"]).strip()]))
    t2 = pd.read_excel(path, sheet_name="Supp Table 2", header=0)
    for _, row in t2.dropna(subset=["GEN ID"]).iterrows():
        gene = str(row["GEN ID"]).split("/")[0].strip()
        records.append(_seq(study, "GCR1 ortholog", str(row["Species"]), gene, "", [gene]))
    return records


NIIMURA_CODES = {
    "Brfl": "Branchiostoma floridae", "Pema": "Petromyzon marinus", "Cami": "Callorhinchus milii",
    "Dare": "Danio rerio", "Orla": "Oryzias latipes", "Gaac": "Gasterosteus aculeatus", "Taru": "Takifugu rubripes",
    "Teni": "Tetraodon nigroviridis", "Xetr": "Xenopus tropicalis", "Anca": "Anolis carolinensis",
    "Gaga": "Gallus gallus",
}


def read_niimura_2009():
    """Dataset S1: the OR genes behind Table 1 (names like BrflORs100.1P; P = pseudogene,
    T = truncated). Dataset S2 (59 genes of other groups) is not used. The genome counts stay
    those of Table 1 (curated_claims.tsv); these sequences are for tracing. Pseudogene sequences
    carry "/" at frameshifts; it is kept as "-" so each piece is matched on its own."""
    z = zipfile.ZipFile(_path("niimura_2009", "evp003_Supplementary_Data.zip"))
    text = z.read("DatasetS1.doc").decode("latin-1")
    records = []
    for name, seq in re.findall(r">([^\r\n>]{1,80})[\r\n]+([A-Z*/\r\n]{30,})", text):
        code = name[:4]
        records.append(_seq("Niimura 2009", "Olfactory", NIIMURA_CODES.get(code, ""), name.strip(),
                            re.sub(r"[^A-Z-]", "", seq.replace("/", "-"))))
    return records


def read_gotkhindikar_2025():
    return [_seq("Gotkhindikar et al. 2025", "GPR157", "Arabidopsis thaliana", "AtGCR1", "", ["O04714"])]


def read_curated():
    records = []
    with open(_path("curated_claims.tsv"), encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            count = int(row["count"]) if row["count"] else ""
            records.append(_claim(row["study"], row["their_family"], row["genome"], row["status"], count, row["source"],
                                  row.get("evidence") or ""))
    return records


READERS = [read_garcia, read_krishnan_2012, read_de_mendoza_2014, read_de_mendoza_2016, read_nordstrom_2009,
           read_nordstrom_2011, read_krishnan_2014_nematode, read_krishnan_2014_sponge, read_grus_zhang_2009,
           read_cardoso_2014, read_prabhu_2006,
           read_itoigawa_2024, read_kowatschew_2022, read_liu_2026, read_fernandez_figueroa_2026,
           read_gotkhindikar_2025, read_niimura_2009, read_curated]


# Where each study's receptor sequences come from (column sequence_source of literature_sequences.tsv)
SEQUENCE_SOURCE = {
    GARCIA: "Supplementary files 4, 9 and 25 (protein FASTA)",
    KRISHNAN12: "Dataset S1 (protein FASTA)",
    "de Mendoza et al. 2016": "Text file S1 (protein sequences read from the supplement PDF)",
    "Nordstrom et al. 2009": "File 2 (coding sequences, translated)",
    "Nordstrom et al. 2011": "Additional text 2 (protein sequences, Word .doc)",
    "Krishnan et al. 2014 (nematode)": "File s007 (protein sequences read from PDF)",
    SPONGE: "Additional file 1 (protein sequences read from PDF)",
    GRUS: "Files 2 and 9 (translated genome contigs, alignment rows read from PDF; gaps kept)",
    "Cardoso et al. 2014": "File s014 (transmembrane-region alignment rows; gaps kept)",
    "Kowatschew & Korsching 2022": "SI file 3 (protein sequences read from the supplement PDF)",
    "Liu et al. 2026": "File S2 (protein sequences, Excel)",
    "Niimura 2009": "Dataset S1 (protein sequences, Word .doc; frameshifts kept as breaks)",
}


def genome_decisions():
    """Genome names replaced by hand (literature/genome_name_decisions.tsv): names that NCBI resolves
    to another species than the genome the study used, or not at all, with the evidence for each.
    A row with id_contains applies only to receptors whose ID contains it (checked first)."""
    rows = list(csv.DictReader(open(_path("genome_name_decisions.tsv"), encoding="utf-8"), delimiter="\t"))
    return sorted(rows, key=lambda d: not d["id_contains"])


# One parsed table per reader: the records it returns, in its order (read_curated reads our own file)
PARSED_FIELDS = ["kind", "study", "their_family", "genome", "status", "count", "source", "evidence",
                 "seq_id", "sequence", "ids"]
CLAIM_KEYS = ["kind", "study", "their_family", "genome", "status", "count", "source", "evidence"]
SEQ_KEYS = ["kind", "study", "their_family", "genome", "seq_id", "sequence", "ids"]


def parsed_path(reader):
    return os.path.join(PARSED_DIR, reader.__name__[len("read_"):] + ".tsv")


def write_parsed(reader, records):
    os.makedirs(PARSED_DIR, exist_ok=True)
    with open(parsed_path(reader), "w", newline="", encoding="utf-8") as handle:
        w = csv.DictWriter(handle, PARSED_FIELDS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in records:
            row = dict(r)
            if r["kind"] == "sequence":
                row["ids"] = ";".join(r["ids"])
            w.writerow(row)


def read_parsed(reader):
    records = []
    with open(parsed_path(reader), encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["kind"] == "claim":
                r = {k: row[k] for k in CLAIM_KEYS}
                r["count"] = int(r["count"]) if r["count"] else ""
            else:
                r = {k: row[k] for k in SEQ_KEYS}
                r["ids"] = row["ids"].split(";") if row["ids"] else []
            records.append(r)
    return records


def read_all(from_originals=False):
    """Every study's records, from the parsed tables, or with from_originals from the supplementary files."""
    records = []
    for reader in READERS:
        records += reader() if from_originals or reader is read_curated else read_parsed(reader)
    decisions = genome_decisions()
    for r in records:
        for d in decisions:
            if (r["study"], r["genome"]) == (d["study"], d["genome_as_reported"]) and \
                    d["id_contains"] in (r.get("seq_id") or ""):
                r["genome"], r["genome_decision"] = d["our_genome"], d["evidence"]
                break
        if r["kind"] == "sequence":
            r["sequence_source"] = (SEQUENCE_SOURCE[r["study"]] if r["sequence"]
                                    else "no sequence in the study (identifiers only)")
    # The same sequence listed more than once by one study, in one family and genome (under
    # different IDs), is one receptor: kept once. A protein a study lists under two families keeps
    # both listings (two claims), and the same protein listed by two studies stays in both.
    seen, unique, dropped = set(), [], 0
    for r in records:
        if r["kind"] == "sequence" and r["sequence"]:
            key = (r["study"], r["their_family"], r["genome"], re.sub(r"[^A-Z]", "", r["sequence"].upper()))
            if key in seen:
                dropped += 1
                continue
            seen.add(key)
        unique.append(r)
    if dropped:
        print("identical sequences listed again by the same study (same family and genome), kept once: %d" % dropped)
    return unique
