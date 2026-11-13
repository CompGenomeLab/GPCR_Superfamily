"""Fetch protein sequences for the literature receptors listed by identifier only (idonly_434.tsv), so
they can be searched with MMseqs2 like the others. Sequences come from public databases (UniProt,
UniSave, UniParc, NCBI, Ensembl); family membership is still judged against our own data only.
Every receptor gets a row in the log with the route used, the fetched record, its organism and how
many candidate records the lookup returned."""
import csv, io, json, os, re, sys, time, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from gpcr_paths import stage  # noqa: E402
OUT = str(stage("05_literature_comparison"))
UP = "https://rest.uniprot.org"
NCBI = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def get(url, tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "gpcr-literature-check"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read().decode("utf-8", "replace")
        except urllib.error.HTTPError as e:
            if e.code in (400, 404, 410):
                return ""
            time.sleep(2 * (i + 1))
        except Exception:
            time.sleep(2 * (i + 1))
    return ""


def fasta_records(text):
    recs, h, s = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if h:
                recs.append((h, "".join(s)))
            h, s = line[1:], []
        elif h:
            s.append(line.strip())
    if h:
        recs.append((h, "".join(s)))
    return recs


def organism(header):
    m = re.search(r"OS=(.*?)(?: OX=|$)", header) or re.search(r"\[([^\]]+)\]", header)
    return m.group(1).strip() if m else ""


def best(recs):
    """Reviewed first, then longest."""
    if not recs:
        return None
    return sorted(recs, key=lambda r: (not r[0].startswith("sp|"), -len(r[1])))[0]


def uniprot_search(query):
    return fasta_records(get(UP + "/uniprotkb/search?format=fasta&size=25&query=" + urllib.parse.quote(query)))


def uniprot_acc(acc):
    recs = fasta_records(get(UP + "/uniprotkb/%s.fasta" % acc))
    if recs:
        return recs, "UniProt accession (current)"
    # retired: last UniSave version that carries a sequence
    recs = fasta_records(get(UP + "/unisave/%s?format=fasta" % acc))
    if recs:
        return recs[:1], "UniSave (retired entry, last version)"
    return [], ""


def uniparc_dbid(i):
    return fasta_records(get(UP + "/uniparc/search?format=fasta&size=10&query=" + urllib.parse.quote("dbid:" + i)))


def ncbi_protein(i):
    t = get(NCBI + "/efetch.fcgi?db=protein&rettype=fasta&retmode=text&id=" + urllib.parse.quote(i))
    return fasta_records(t) if t.startswith(">") else []


def ncbi_gene_proteins(loc):
    """NCBI gene LOC id -> all protein isoforms; the longest is kept."""
    gid = loc.replace("LOC", "")
    t = get(NCBI + "/elink.fcgi?dbfrom=gene&db=protein&linkname=gene_protein_refseq&retmode=json&id=" + gid)
    try:
        ids = json.loads(t)["linksets"][0]["linksetdbs"][0]["links"]
    except Exception:
        return []
    return fasta_records(get(NCBI + "/efetch.fcgi?db=protein&rettype=fasta&retmode=text&id=" + ",".join(ids[:50])))


def ensembl_canonical(gene):
    t = get("https://rest.ensembl.org/lookup/id/%s?expand=1;content-type=application/json" % gene)
    try:
        d = json.loads(t)
        tr = [x for x in d["Transcript"] if x.get("is_canonical")] or d["Transcript"]
        pid = tr[0]["Translation"]["id"]
    except Exception:
        return []
    s = get("https://rest.ensembl.org/sequence/id/%s?type=protein;content-type=text/x-fasta" % pid)
    return fasta_records(s)


ACC = re.compile(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$")
GENBANK = re.compile(r"^[A-Z]{3}\d{5}(\.\d+)?$")


def lookup(row):
    i, genome = row["id"].strip(), row["genome"].strip()
    sp = " ".join(genome.split()[:2])
    org = ('organism_id:%s' % row["our_taxid"]) if row["our_taxid"] else ('organism_name:"%s"' % sp)
    tried = []
    # UniProt entry name or accession
    m = re.match(r"^([A-Z0-9]+)_([A-Z0-9]{2,5})$", i)
    if m or ACC.match(i):
        if m:
            recs = uniprot_search("id:" + i)
            if recs:
                return recs, "UniProt entry name (current)", tried
            tried.append("entry name")
        acc = m.group(1) if m else i
        if ACC.match(acc):
            recs, how = uniprot_acc(acc)
            if recs:
                return recs, how, tried
            tried.append("accession/UniSave")
    if GENBANK.match(i) or re.match(r"^[A-Z]{3}\d{5}$", i):
        recs = ncbi_protein(i)
        if recs:
            return recs, "NCBI Protein (GenBank)", tried
        tried.append("NCBI protein")
    if i.startswith("ENS"):
        recs = ensembl_canonical(i)
        if recs:
            return recs, "Ensembl canonical translation", tried
        tried.append("Ensembl")
    if i.startswith("LOC") and i[3:].isdigit():
        recs = ncbi_gene_proteins(i)
        if recs:
            return [max(recs, key=lambda r: len(r[1]))], "NCBI Gene -> RefSeq protein (longest isoform of %d)" % len(recs), tried
        tried.append("NCBI gene")
    # gene names / locus tags (Broad, VectorBase, WormBase, plant names); dictyBase "DDB0... fslJ" -> fslJ
    names = []
    if row["study"].startswith("Prabhu"):
        names = [i.split()[-1]]
    else:
        base = re.sub(r"T0$|-P[A-Z]$|\.\d+$", "", i)
        names = [base, base.replace(" ", ""), base.replace("SNOT_", "SNOG_").replace("PGUT_", "PGUG_")]
    for n in dict.fromkeys(names):
        recs = uniprot_search('(gene_exact:"%s") AND %s' % (n, org)) or uniprot_search('(gene:"%s") AND %s' % (n, org))
        if recs:
            return recs, "UniProt gene/ORF name in species", tried
        tried.append("UniProt gene %s" % n)
    for n in dict.fromkeys(names):
        recs = uniparc_dbid(n)
        if recs:
            return recs, "UniParc cross-reference", tried
        tried.append("UniParc %s" % n)
    for n in dict.fromkeys(names):
        t = get(NCBI + "/esearch.fcgi?db=protein&retmode=json&term=" + urllib.parse.quote('%s AND "%s"[Organism]' % (n, sp)))
        try:
            ids = json.loads(t)["esearchresult"]["idlist"]
        except Exception:
            ids = []
        if ids:
            recs = fasta_records(get(NCBI + "/efetch.fcgi?db=protein&rettype=fasta&retmode=text&id=" + ",".join(ids[:10])))
            if recs:
                return recs, "NCBI Protein search in species", tried
        tried.append("NCBI search %s" % n)
    return [], "", tried


if __name__ == "__main__":
    rows = list(csv.DictReader(open(os.path.join(HERE, "query_inputs", "idonly_434.tsv"), encoding="utf-8"), delimiter="\t"))
    clean = lambda s: re.sub(r"[^A-Za-z0-9.=_-]+", "_", s).strip("_")
    n = 0
    with open(os.path.join(OUT, "idonly_literature_sequences.fasta"), "w", encoding="utf-8") as fa, \
            open(os.path.join(OUT, "idonly_literature_sequences.tsv"), "w", encoding="utf-8", newline="") as lg:
        w = csv.writer(lg, delimiter="\t")
        w.writerow(["query_id", "study", "their_family", "id", "genome", "stopped_at", "route", "n_candidates",
                    "fetched_record", "fetched_organism", "same_species", "length", "routes_failed"])
        for k, row in enumerate(rows, 1):
            recs, how, tried = lookup(row)
            b = best(recs)
            if b:
                n += 1
                qid = "LITI%04d" % n
                org = organism(b[0])
                same = " ".join(org.split()[:2]).lower() == " ".join(row["genome"].split()[:2]).lower()
                fa.write(">%s|%s|claimed=%s|%s|stopped=%s\n%s\n" % (qid, clean(row["study"]), clean(row["their_family"]),
                                                                    clean(row["genome"] or "no_species"), clean(row["stopped_at"]), b[1]))
                w.writerow([qid, row["study"], row["their_family"], row["id"], row["genome"], row["stopped_at"], how, len(recs),
                            b[0][:150], org, "yes" if same else "no", len(b[1]), "; ".join(tried)])
            else:
                w.writerow(["", row["study"], row["their_family"], row["id"], row["genome"], row["stopped_at"], "not found", 0,
                            "", "", "", "", "; ".join(tried)])
            if k % 25 == 0:
                print(k, "done,", n, "fetched", flush=True)
            time.sleep(0.35)
    print("fetched %d of %d" % (n, len(rows)))
