"""Clean Behrens SI set: parsed rows, minus the Lethenteron V1Rs (Kowatschew & Korsching's sequences, already
compared from their own SI) and the Lepisosteus rows the layout parse cannot split; those are taken from
NCBI by the Gene ID or genome coordinates the table gives."""
import json, os, re, time, urllib.request
from Bio.Seq import Seq

NCBI = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
get = lambda u: urllib.request.urlopen(u, timeout=60).read().decode()


def fasta(text):
    out, h = [], None
    for l in text.splitlines():
        if l.startswith(">"):
            out.append([l[1:], ""])
        elif out:
            out[-1][1] += l.strip()
    return out


def gene_protein(gid):
    t = json.loads(get(NCBI + "/elink.fcgi?dbfrom=gene&db=protein&linkname=gene_protein_refseq&retmode=json&id=%s" % gid))
    ids = t["linksets"][0]["linksetdbs"][0]["links"]
    recs = fasta(get(NCBI + "/efetch.fcgi?db=protein&rettype=fasta&retmode=text&id=" + ",".join(ids[:20])))
    h, s = max(recs, key=lambda r: len(r[1]))
    return s, "NCBI Gene %s -> %s (longest of %d)" % (gid, h.split()[0], len(recs))


def locus(acc, a, b):
    nt = "".join(get(NCBI + "/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&strand=2&seq_start=%d&seq_stop=%d&id=%s"
                     % (a, b, acc)).splitlines()[1:])
    best = ""
    for f in range(3):
        prot = str(Seq(nt[f:len(nt) - (len(nt) - f) % 3]).translate())
        for piece in prot.split("*"):
            if "M" in piece and len(piece) - piece.index("M") > len(best):
                best = piece[piece.index("M"):]
    return best, "%s:%d-%d(-) translated, longest ORF" % (acc, a, b)

FETCH = {
    "Leoc-T2R2": lambda: locus("NC_023182.1", 36789856, 36790839),
    "Leoc-T2R3": lambda: locus("NC_023198.1", 9783077, 9784207),
    "Lo-ORA1": lambda: gene_protein(102693077), "Lo-ORA2": lambda: gene_protein(102693803),
    "Lo-ORA3": lambda: gene_protein(102684493), "Lo-ORA4": lambda: gene_protein(102684096),
    "Lo-ORA5": lambda: gene_protein(102682960),
    "Lo-ORA6": lambda: locus("NW_006270091.1", 262491, 263456),
}

recs = fasta(open("behrens_si.fasta", encoding="utf-8").read())
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "query_inputs", "behrens_si_clean.fasta"), "w", encoding="utf-8") as fh:
    kept = 0
    for h, s in recs:
        f = h.split("|")
        name, species = f[0], f[1]
        if species == "Lethenteron_camtschaticum":
            continue
        gene = next((k for k in FETCH if k.replace("-", "") in h.replace("-", "").replace(".", "")), None)
        src = "Behrens SI Table S1"
        if species == "Lepisosteus_oculatus" and gene:
            s, src = FETCH[gene]()
            time.sleep(0.4)
        fh.write(">%s|%s|%s|%s\n%s\n" % (name, species, f[2], src.replace(" ", "_"), s))
        kept += 1
        print(f[2][:14], species, name, gene or "", len(s), src)
print(kept, "kept")
