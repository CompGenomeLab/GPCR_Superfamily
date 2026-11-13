"""One FASTA of every literature sequence not yet searched with MMseqs2 (query ids LITX#####), to be
merged with unmatched_literature_sequences.m8. Sources:
  A. identifier-only receptors, sequences fetched and organism-verified (idonly_literature_sequences_r2.fasta)
  B. Behrens et al. 2023 SI Table S1 (cartilaginous fish T2R1, reference T2R, outgroup ORA/V1R)
  C. Itoigawa et al. 2024 Zenodo 10042321 TAS2R nucleotide sequences, translated
  D. reference/outgroup sequences the readers skip: Krishnan 2014 nematode s007 non-srw, Kowatschew &
     Korsching 2022 outgroups, Niimura 2009 Dataset S2, Krishnan 2012 Dataset S1 AMAG_16966
Writes remaining_literature_sequences.fasta and .tsv in 05_literature_comparison/ (in the archive). Reads some
studies' supplementary files directly, which are not distributed: download them (literature/SOURCES.md) and set
GPCR_LITERATURE_ORIGINALS to their folder."""
import csv, os, re, sys, zipfile
from Bio.Seq import Seq

STEP6 = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(STEP6))
from gpcr_paths import stage  # noqa: E402
OUT = str(stage("05_literature_comparison"))
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, STEP6)
import literature_readers as L


def fasta(path):
    out = []
    for l in open(path, encoding="utf-8"):
        l = l.strip()
        if l.startswith(">"):
            out.append([l[1:], ""])
        elif out:
            out[-1][1] += l
    return out


def fewest_stops(nt):
    """Pseudogenes: the reading frame (either strand) with the fewest stops, stops read through as X."""
    best = None
    for s in (nt, str(Seq(nt).reverse_complement())):
        for f in range(3):
            p = str(Seq(s[f:len(s) - (len(s) - f) % 3]).translate()).strip("*")
            if best is None or p.count("*") < best.count("*"):
                best = p
    return best.replace("*", "X")


clean = lambda s: re.sub(r"[^A-Za-z0-9.=_-]+", "_", s).strip("_")
rows = []
# sequences the readers already return (matched directly or in the first MMseqs2 query set) are
# not searched again
ALREADY = {re.sub(r"[^A-Z]", "", r["sequence"].upper()) for r in L.read_all() if r["kind"] == "sequence" and r["sequence"]}
skipped = []


def add(study, family, species, rid, kind, source, seq):
    seq = re.sub(r"[^A-Z]", "", seq.upper().replace("*", "X"))
    if kind.startswith("reference"):
        return  # reference and outgroup sequences are not the studies' claims
    if seq in ALREADY:
        skipped.append((study, rid))
        return
    if len(seq) >= 30:
        rows.append([study, family, species, rid, kind, source, seq])

# A. identifier-only receptors
p = os.path.join(OUT, "idonly_literature_sequences_r2.fasta")
meta = {r["query_id"]: r for r in csv.DictReader(open(os.path.join(OUT, "idonly_literature_sequences_r2.tsv"), encoding="utf-8"), delimiter="\t")}
for h, s in fasta(p):
    m = meta[h.split("|")[0]]
    add(m["study"], m["their_family"], m["genome"], m["id"], "study receptor (identifier only; sequence fetched)", m["route"], s)


# A2. identifier-only receptors recovered after round 2 (study gives no species for CC1G_; cDNA; Ensembl)
EXTRA = {"CC1G_07132": ("cc1g.fasta", "A8NR64", "Coprinopsis cinerea", "UniProt free text (A8NR64; same protein as Liu et al. Table S2)"),
         "CC1G_04180": ("cc1g.fasta", "A8NF49", "Coprinopsis cinerea", "UniProt free text (A8NF49; same protein as Liu et al. Table S2)"),
         "BT110964.1": ("idonly_round3.fasta", "BT110964.1", "Picea glauca", "GenBank cDNA, longest ORF"),
         "DRNTG_19123.1": ("idonly_round3.fasta", "DRNTG_19123", "Dioscorea rotundata", "Ensembl Plants canonical translation")}
r2 = {r["id"]: r for r in csv.DictReader(open(os.path.join(OUT, "idonly_literature_sequences_r2.tsv"), encoding="utf-8"), delimiter="	")}
for rid, (f, key, sp, how) in EXTRA.items():
    s = next(s for h, s in fasta(os.path.join(HERE, "query_inputs", f)) if key in h.split("|")[0] or key in h)
    m = r2[rid]
    add(m["study"], m["their_family"], m["genome"] or sp, rid, "study receptor (identifier only; sequence fetched)", how, s)

# B. Behrens SI Table S1
for h, s in fasta(os.path.join(HERE, "query_inputs", "behrens_si_clean.fasta")):
    name, species, section, src = (h.split("|") + ["", "", "", ""])[:4]
    kind = "study receptor" if section.startswith("cartilaginous") else "reference/outgroup sequence"
    add("Behrens et al. 2023", "T2R" if "T2R" in section else "ORA/V1R (outgroup)", species.replace("_", " "), name, kind,
        "SI Table S1 (%s)" % src.replace("_", " "), s)

# C. Itoigawa Zenodo nucleotide sequences
CODE = {"Leer": "Leucoraja erinacea", "Amra": "Amblyraja radiata", "Prpe": "Pristis pectinata", "Heak": "Hemitrygon akajei",
        "Pole": "Potamotrygon leopoldi", "Scto": "Scyliorhinus torazame", "Spmo": "Sphyrna mokarran", "Stfa": "Stegostoma fasciatum",
        "Rhty": "Rhincodon typus", "Chpu": "Chiloscyllium punctatum", "Chpl": "Chiloscyllium plagiosum"}
for h, nt in fasta(L._path("itoigawa_2024", "zenodo_10042321", "Seq_elasmobranchT2R_nt.fasta")):
    nt = nt.upper()
    pseudo = h.endswith("p")
    prot = fewest_stops(nt) if pseudo else str(Seq(nt[:len(nt) - len(nt) % 3]).translate()).rstrip("*")
    add("Itoigawa et al. 2024", "T2R", CODE[h[:4]], h, "study receptor" + (" (pseudogene)" if pseudo else ""),
        "Zenodo 10042321 Seq_elasmobranchT2R.docx, " + ("fewest-stop frame, stops as X" if pseudo else "frame 1"), prot)

# D1. Krishnan 2014 nematode s007: sequences outside their srw list (vertebrate srsx and other reference genes)
z = L._docx_text(open(L._path("krishnan_2014_nematode", "pone.0093048.s005.docx"), "rb").read())
srw = set(re.findall(r"([A-Z][a-z]_\w+|\d+[abc] \(AASC\d+\))", z))
for header, seq in L.parse_fasta_text(L._pdf_text(L._path("krishnan_2014_nematode", "pone.0093048.s007.pdf"))):
    name = re.match(r"(\d+[abc] \(AASC\d+\)|\S+)", header).group(1)
    if name not in srw:
        add("Krishnan et al. 2014 (nematode)", "reference (not srw)", "", header.strip(), "reference/outgroup sequence", "s007.pdf", seq)

# D2. Kowatschew & Korsching 2022 outgroups
t = L._pdf_text(L._path("kowatschew_korsching_2022", "bjac007_suppl_supplementary_material.pdf"))
a = t.find("Predicted amino acid sequences for lamprey"); b = t.find("SI file 4", a)
codes = {k.lower(): v for k, v in L.KOWATSCHEW_CODES.items()}
for header, seq in L.parse_fasta_text(t[a:b if b > 0 else None]):
    name = header.split()[0]
    if not re.search(r"V2R|OLFC|V1R|ORA|T2R", name.upper()):
        add("Kowatschew & Korsching 2022", "outgroup", codes.get(re.split(r"[-_.]", name)[0].lower(), ""), name,
            "reference/outgroup sequence", "SI file 3", seq)

# D3. Niimura 2009 Dataset S2
zz = zipfile.ZipFile(L._path("niimura_2009", "evp003_Supplementary_Data.zip"))
tx = zz.read("DatasetS2.doc").decode("latin-1")
for name, seq in re.findall(r">([^\r\n>]{1,80})[\r\n]+([A-Z*/\r\n]{30,})", tx):
    add("Niimura 2009", "other OR-like group (Dataset S2)", L.NIIMURA_CODES.get(name[:4], ""), name.strip(),
        "reference/outgroup sequence", "Dataset S2", seq.replace("/", ""))

# D4. Krishnan 2012 Dataset S1 record absent from Table S1
for header, seq in L.parse_fasta_text(open(L._path("krishnan_2012", "pone.0029817.s001.txt"), encoding="utf-8", errors="replace").read()):
    if "AMAG_16966" in header:
        add("Krishnan et al. 2012", "Glutamate (Dataset S1 only)", "Allomyces macrogynus", "AMAG_16966",
            "study sequence not in Table S1", "Dataset S1", seq)

with open(os.path.join(OUT, "remaining_literature_sequences.fasta"), "w", encoding="utf-8") as fa, \
        open(os.path.join(OUT, "remaining_literature_sequences.tsv"), "w", encoding="utf-8", newline="") as tv:
    w = csv.writer(tv, delimiter="\t")
    w.writerow(["query_id", "study", "their_family", "species", "id", "kind", "source", "length"])
    for k, (study, fam, sp, rid, kind, src, seq) in enumerate(rows, 1):
        q = "LITX%05d" % k
        fa.write(">%s|%s|claimed=%s|%s|%s\n%s\n" % (q, clean(study), clean(fam), clean(sp or "no_species"), clean(kind), seq))
        w.writerow([q, study, fam, sp, rid, kind, src, len(seq)])
from collections import Counter
print(len(rows), "sequences;", len(skipped), "skipped as identical to a sequence the readers already return")
for k, v in sorted(Counter((r[0], r[4]) for r in rows).items()):
    print("  %4d  %s | %s" % (v, k[0], k[1]))
