"""Second round for the identifier-only receptors:
1. Fernandez-Figueroa Table 1 receptors by their UniProt code (not the tip label), Table 2 by the
   database the study names (Ensembl REST for Ensembl Plants/Metazoa/Protists/Fungi).
2. Every other receptor not found, or found only through a name search or a UniParc cross-reference,
   retried and checked: UniProt free text, NCBI protein/gene without a species filter, and the organism
   of the record checked against the study's genus.
Writes idonly_literature_sequences_r2.tsv/.fasta (all 434 rows, with a verified flag)."""
import csv, json, os, re, sys, time, urllib.parse, urllib.request
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import openpyxl
from fetch_idonly import get, fasta_records, uniprot_acc, ensembl_canonical, ncbi_gene_proteins, UP, NCBI

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gpcr_paths import stage  # noqa: E402
OUT = str(stage("05_literature_comparison"))
import literature_readers as L  # noqa: E402
# reads the study's supplementary file itself (not distributed; GPCR_LITERATURE_ORIGINALS, literature/SOURCES.md)
FF = L._path("fernandez_figueroa_2026", "media-1_SuppTables1-2.xlsx")
wb = openpyxl.load_workbook(FF, read_only=True)
T1 = {str(r[1]).strip(): str(r[2]).strip() for r in list(wb.worksheets[0].iter_rows(values_only=True))[1:] if r and r[1] and r[2]}
T2 = {str(r[4]).strip(): (str(r[5] or "").strip(), str(r[2] or "").strip()) for r in list(wb.worksheets[1].iter_rows(values_only=True))[1:] if r and r[4]}

SYN = {"Nannizzia": "Arthroderma", "Trichophyton": "Arthroderma", "Microsporum": "Arthroderma", "Fusarium": "Gibberella",
       "Aspergillus": "Emericella", "Histoplasma": "Ajellomyces", "Blastomyces": "Ajellomyces", "Nakaseomyces": "Candida",
       "Plenodomus": "Leptosphaeria", "Leptographium": "Grosmannia", "Pyricularia": "Magnaporthe", "Komagataella": "Pichia",
       "Rhodonia": "Postia", "Ogataea": "Pichia", "Parastagonospora": "Phaeosphaeria"}


def genus_ok(study_genome, organism):
    g = (study_genome.split() or [""])[0]
    o = organism.lower()
    return bool(g) and (g.lower() in o or (g in SYN and SYN[g].lower() in o))


def uniparc_org(upi):
    t = get(UP + "/uniparc/%s?format=json" % upi)
    try:
        d = json.loads(t)
        orgs = {x.get("organism", {}).get("scientificName", "") for x in d.get("uniParcCrossReferences", [])}
        return "; ".join(sorted(o for o in orgs if o))
    except Exception:
        return ""


def unisave_org(acc, header):
    m = re.search(r"OS=(.*?) OX=", header)
    if m:
        return m.group(1)
    t = get(UP + "/unisave/%s?format=txt&versions=1" % acc)
    m = re.search(r"^OS\s+(.*)$", t, re.M)
    return m.group(1).strip(" .") if m else ""


def free_text(i, genome):
    g = genome.split()[0] if genome else ""
    q = '"%s"' % i + (' AND organism_name:"%s"' % g if g else "")
    recs = fasta_records(get(UP + "/uniprotkb/search?format=fasta&size=10&query=" + urllib.parse.quote(q)))
    if recs:
        return recs, "UniProt free text"
    t = get(NCBI + "/esearch.fcgi?db=protein&retmode=json&term=" + urllib.parse.quote(i))
    ids = json.loads(t or "{}").get("esearchresult", {}).get("idlist", [])
    if ids:
        recs = fasta_records(get(NCBI + "/efetch.fcgi?db=protein&rettype=fasta&retmode=text&id=" + ",".join(ids[:10])))
        recs = [r for r in recs if genus_ok(genome, r[0])]
        if recs:
            return recs, "NCBI Protein, identifier without species filter"
    t = get(NCBI + "/esearch.fcgi?db=gene&retmode=json&term=" + urllib.parse.quote(i))
    ids = json.loads(t or "{}").get("esearchresult", {}).get("idlist", [])
    for gid in ids[:3]:
        recs = [r for r in ncbi_gene_proteins("LOC" + gid) if genus_ok(genome, r[0])]
        if recs:
            return [max(recs, key=lambda r: len(r[1]))], "NCBI Gene (identifier as alias) -> RefSeq protein"
    return [], ""

r1 = list(csv.DictReader(open(os.path.join(OUT, "idonly_literature_sequences.tsv"), encoding="utf-8"), delimiter="\t"))
seq1 = {h.split("|")[0]: s for h, s in fasta_records(open(os.path.join(OUT, "idonly_literature_sequences.fasta"), encoding="utf-8").read())}
clean = lambda s: re.sub(r"[^A-Za-z0-9.=_-]+", "_", s).strip("_")
rows_out, n = [], 0
for x in r1:
    i, genome, study = x["id"].strip(), x["genome"], x["study"]
    seq, route, rec, org = seq1.get(x["query_id"], ""), x["route"], x["fetched_record"], x["fetched_organism"]
    redo = False
    if study.startswith("Fernandez") and i in T1:                       # Table 1: UniProt code
        recs, how = uniprot_acc(T1[i])
        seq, route, rec = (recs[0][1], how + " (Table 1 code %s)" % T1[i], recs[0][0][:150]) if recs else ("", "not found", "")
        org = unisave_org(T1[i], rec) if recs else ""
    elif study.startswith("Fernandez") and i in T2 and T2[i][0].startswith("Ensembl") and not route.startswith("UniProt"):
        recs = ensembl_canonical(i)
        if recs:
            seq, route, rec, org = recs[0][1], "Ensembl canonical translation (%s)" % T2[i][0], recs[0][0][:150], genome + " (Ensembl gene ID)"
        else:
            redo = True
    elif route == "not found" or route.startswith(("UniProt gene/ORF name", "NCBI Protein search", "NCBI Protein (GenBank)")) and not genus_ok(genome, rec):
        redo = True
    if redo:
        recs, how = free_text(i, genome)
        if recs:
            b = sorted(recs, key=lambda r: (not r[0].startswith("sp|"), -len(r[1])))[0]
            seq, route, rec, org = b[1], how, b[0][:150], b[0]
        else:
            seq, route, rec, org = "", "not found", "", ""
    # organism check
    if route.startswith("UniParc"):
        org = uniparc_org(rec.split()[0])
    elif route.startswith("UniSave"):
        org = unisave_org(re.sub(r"_.*", "", i), rec)
    elif route.startswith("Ensembl"):
        org = org or genome + " (Ensembl gene ID)"
    elif not org:
        m = re.search(r"OS=(.*?) OX=|\[([^\]]+)\]", rec)
        org = (m.group(1) or m.group(2)) if m else ""
    verified = "yes" if seq and (genus_ok(genome, org) or route.startswith(("Ensembl", "UniSave")) and org) else ("check" if seq else "")
    if seq and verified == "check" and "Table 1 code" in route:
        verified = "yes"  # the study's own UniProt code
    qid = ""
    if seq and verified == "yes":
        n += 1
        qid = "LITI%04d" % n
    rows_out.append([qid, study, x["their_family"], i, genome, x["stopped_at"], route, rec, org, verified, len(seq), seq])
    time.sleep(0.3)

with open(os.path.join(OUT, "idonly_literature_sequences_r2.tsv"), "w", encoding="utf-8", newline="") as fh, \
        open(os.path.join(OUT, "idonly_literature_sequences_r2.fasta"), "w", encoding="utf-8") as fa:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["query_id", "study", "their_family", "id", "genome", "stopped_at", "route", "fetched_record", "fetched_organism", "organism_verified", "length"])
    for r in rows_out:
        w.writerow(r[:-1])
        if r[0]:
            fa.write(">%s|%s|claimed=%s|%s|stopped=%s\n%s\n" % (r[0], clean(r[1]), clean(r[2]), clean(r[4] or "no_species"), clean(r[5]), r[-1]))
from collections import Counter
print("verified sequences:", n, "of", len(rows_out))
print(Counter(r[9] or "none" for r in rows_out))
for r in rows_out:
    if r[9] != "yes":
        print("  ", r[9] or "NOT FOUND", "|", r[1][:18], "|", r[3], "|", r[4], "|", r[6][:40], "|", r[8][:60])
