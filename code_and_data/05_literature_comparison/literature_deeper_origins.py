"""Literature presences that would place a family's origin deeper than ours, and the sequences to test.

Step 6, literature validation. Our origin of each family is the interpretation column of the
lineage figure (06_family_origins/lineage_origin_figure.html in the archive; Dollo origin given beside it). A
literature presence is a candidate when the study reports the family in a species outside that
origin node (NCBI taxonomy) and our data do not show the family in that genome (the genome is not in
our set, or we do not find the family in it). For each candidate the study's receptors are taken
from literature_sequences.tsv, and every one with a sequence is tested: the study's sequence, or,
for a receptor listed by identifier only and matched to a protein of ours, that protein's.
write_deeper_origin_fasta.py writes them as one FASTA; their closest family of ours is their MMseqs2
best hit among our family members (deeper_origin_best_hits.py). The elasmobranch TAS2R loci of Itoigawa et al. 2024 are
read at their coordinates (NCBI contigs, and the GWH assembly of Potamotrygon leopoldi), translated
and tested, pseudogenes included (read through their stops); Behrens et al. 2023 claims are tested
with Itoigawa's locus in the same species. Count-only claims, and named genes with no locus, have
nothing to test.
A species name the NCBI taxonomy does not resolve is placed by the species its abbreviation stands
for, or by the taxon its first word names (genus, or a lineage in EukProt names); a claim is dropped
only when every taxon the name may belong to lies inside our origin, and is otherwise tested.

Writes to 05_literature_comparison/ in the archive:
    deeper_origin_candidates.tsv   one row per candidate claim (study x family x genome)
    DEEPER_ORIGINS.md              the candidates by family (no interpretation)
    itoigawa_tas2r_loci.fasta      cache of the TAS2R loci (nucleotide; NCBI is queried only for loci
                                   not in it)

Usage (run compare_literature.py first):
    python literature_deeper_origins.py
"""

import csv
import glob
import json
import os
import re
import time
import urllib.request
from collections import Counter, defaultdict

import compare_literature as C
import literature_readers as L

OUT = C.OUT_DIR
FIGURE = os.path.join(C.ORIGINS_DIR, "lineage_origin_figure.html")
PRESENCE = os.path.join(C.ORIGINS_DIR, "family_presence_absence_ncbi_all.csv")
TAS2R_CACHE = os.path.join(OUT, "itoigawa_tas2r_loci.fasta")  # nucleotide sequence of each locus
# The GWH assembly of Potamotrygon leopoldi (GWHAOTN00000000), needed only when a locus is
# missing from TAS2R_CACHE; download it once and point GPCR_GWH_GENOME at it.
GWH_GENOME = os.environ.get("GPCR_GWH_GENOME", "")
FAMILIES = ["classA", "Olfactory", "Nematode", "classT", "Vomeronasal1", "GPR1", "classB1", "classB2", "Mth",
            "classC", "Vomeronasal2", "classF", "FSLB", "cAMP", "GP143", "GP157", "TM116", "STE2", "STE3"]
MIN_ORF = 200  # aa; a shorter reading frame means the coordinates do not hold the gene


# ---------------------------------------------------------------- our origins

def our_origins(matcher):
    """family -> (interpretation, Dollo, test: is a species (its ancestor set) inside the interpretation node)."""
    page = open(FIGURE, encoding="utf-8").read()
    data, _ = json.JSONDecoder().raw_decode(page[page.index('{"lineages"'):])
    carriers = defaultdict(list)
    with open(PRESENCE) as handle:
        for row in csv.DictReader(handle):
            for f in FAMILIES:
                if row.get(f) == "1":
                    carriers[f].append(row["taxID"])
    names = ({r["interp"] for r in data["rows"]} - {"LECA"}) | {"Opisthokonta", "Amoebozoa", "Fungi", "Rotosphaerida"}
    by_name = defaultdict(set)
    with C.open_text(os.path.join(C.NCBI_DIR, "names.dmp")) as handle:
        for line in handle:
            f = line.split("\t|\t")
            if f[1] in names and f[3].startswith("scientific name"):
                by_name[f[1]].add(f[0])
    one = lambda name: next(iter(by_name[name]))
    # Unikonta and Holozoa are clades of the newer NCBI taxonomy our genome table was built from
    # (taxonomy_ncbi_searched.json), absent from the taxonomy files on disk; as in that table,
    # Unikonta = Opisthokonta + Amoebozoa, Holozoa = Opisthokonta outside Fungi and Rotosphaerida.
    opi, amo, fun, rot = one("Opisthokonta"), one("Amoebozoa"), one("Fungi"), one("Rotosphaerida")
    special = {"Unikonta": lambda ups: opi in ups or amo in ups,
               "Holozoa": lambda ups: opi in ups and fun not in ups and rot not in ups}
    origins = {}
    for r in data["rows"]:
        if r["family"] not in FAMILIES:
            continue
        if r["interp"] == "LECA":
            inside = lambda ups: True
        elif r["interp"] in special:
            inside = special[r["interp"]]
        else:
            # a node name can belong to more than one taxon; take the one above most of our carriers
            # (set-aside carriers lie outside it by definition)
            ups_all = [set(ancestors(matcher, t)) for t in carriers[r["family"]]]
            node = max(by_name[r["interp"]], key=lambda t: sum(t in u for u in ups_all))
            inside = (lambda node: lambda ups: node in ups)(node)
        origins[r["family"]] = (r["interp"], r["dollo"], inside)
    return origins


def ancestors(matcher, taxid):
    taxid, path = matcher.taxid_alias.get(taxid, taxid), []
    while taxid in matcher.parent and len(path) < 80:
        path.append(taxid)
        if matcher.parent[taxid] == taxid:
            break
        taxid = matcher.parent[taxid]
    return path


# ---------------------------------------------------------------- candidates

def place_unresolved(matcher, labels):
    """label -> (taxIDs it may be, how) for genome names the taxonomy does not resolve to a species:
    the species an abbreviation stands for ("H. contortus": every species with that initial and
    epithet), else the taxon named by the first word (genus, or a higher lineage in EukProt names
    such as "Apusomonadida sp"). A name can belong to several taxa (homonyms); all are returned."""
    wanted, abbrev = defaultdict(set), {}
    for label in labels:
        words = re.sub(r"\b(cf|aff|afn)\.\s*", "", label).split()
        m = re.fullmatch(r"([A-Z])\.\s*(\S+)", " ".join(words[:2])) if len(words) >= 2 else None
        if m:
            abbrev[label] = (m.group(1), C._plain(m.group(2)))
        elif words:
            wanted[words[0]].add(label)
    epithets = {e for _, e in abbrev.values()}
    by_first, by_abbrev = defaultdict(set), defaultdict(set)
    with C.open_text(os.path.join(C.NCBI_DIR, "names.dmp")) as handle:
        for line in handle:
            f = line.split("\t|\t")
            if not f[3].startswith("scientific name"):
                continue
            if f[1] in wanted:
                by_first[f[1]].add(f[0])
            parts = f[1].split()
            if len(parts) == 2 and C._plain(parts[1]) in epithets and matcher.rank.get(f[0]) == "species":
                by_abbrev[(parts[0][0], C._plain(parts[1]))].add(f[0])
    placed = {}
    for label in labels:
        if label in abbrev and by_abbrev.get(abbrev[label]):
            placed[label] = (by_abbrev[abbrev[label]], "abbreviation")
        else:
            first = re.sub(r"\b(cf|aff|afn)\.\s*", "", label).split()[:1]
            if first and by_first.get(first[0]):
                placed[label] = (by_first[first[0]], "by its first word (%s)" % first[0])
    return placed


def candidates(matcher, origins):
    rows = list(csv.DictReader(open(os.path.join(OUT, "literature_comparison.tsv"), encoding="utf-8"),
                               delimiter="\t"))
    seqs = list(csv.DictReader(open(os.path.join(OUT, "literature_sequences.tsv"), encoding="utf-8"),
                               delimiter="\t"))
    receptors = defaultdict(list)
    for s in seqs:
        receptors[(s["study"], s["their_family"], s["genome"])].append(s)
    unresolved = {r["genome_as_reported"] for r in rows if not r["our_taxids"]
                  and not matcher.match(r["genome_as_reported"])[0]}
    placed = place_unresolved(matcher, unresolved)
    out = []
    for r in rows:
        if r["their_status"] != "present" or r["comparison"] == "agree: present" or r["our_families"] == "any":
            continue
        fams = [f for f in r["our_families"].split(",") if f in origins]
        if not fams or all(origins[f][0] == "LECA" for f in fams):
            continue  # no species can lie outside a LECA origin
        species = (matcher.species_of(r["our_taxids"].split(",")[0]) if r["our_taxids"]
                   else matcher.match(r["genome_as_reported"])[0])
        if species:
            options, placement = [species], ""
        elif r["genome_as_reported"] in placed:
            options, how = placed[r["genome_as_reported"]]
            placement = "species name not resolved; placed %s" % how
        else:
            options, placement = [], "species name not resolved; lineage unknown, tested"
        # outside the origin of every family of ours that answers for theirs, under every reading of the name
        verdicts = {all(not origins[f][2](set(ancestors(matcher, t))) for f in fams) for t in options}
        if verdicts == {False}:
            continue
        if len(verdicts) > 1:
            placement += "; the name fits taxa inside and outside our origin, tested"
        # the study's receptors under any of the names combined in this row ("A; B" for several strains)
        names = ({r["genome_matched"]} | set(r["genome_as_reported"].split("; "))) - {""}
        recs = [s for n in sorted(names) for s in receptors.get((r["study"], r["their_family"], n), [])]
        out.append(dict(r, species=species or "", placement=placement, receptors=recs,
                        deeper_than=", ".join("%s (%s)" % (f, origins[f][0]) for f in fams)))
    return out


def our_protein_sequences(cands):
    """Sequences of our proteins matched to candidate receptors that the study lists without one."""
    wanted = {s["our_accession"] for c in cands for s in c["receptors"]
              if s["our_accession"] and s["id"] not in c.get("sequences", {})}
    found = {}
    for path in sorted(p for pattern in C.SEQUENCE_FILES for p in glob.glob(pattern)):
        if len(found) == len(wanted):
            break
        for acc, _, seq in C._fasta_gz(path):
            if acc in wanted and acc not in found:
                found[acc] = seq.replace("*", "")
    return found


def sort_receptors(cands):
    """Per candidate: the receptors in our data (with where they end), and every receptor with a
    sequence to test: the study's own, or, for a receptor listed by identifier only and matched to a
    protein of ours, that protein's."""
    ours = our_protein_sequences(cands)
    queries = []
    for c in cands:
        found = Counter()
        for s in c["receptors"]:
            if s["our_accession"]:
                found[re.sub(r" \(.*", "", s["stopped_at"]) + ("" if not s["our_final"] or s["our_final"] == "-"
                                                               else " (final %s)" % s["our_final"])] += 1
        c["in_our_data"] = "; ".join("%s %d" % kv for kv in found.most_common())
        c["other_species"] = "; ".join(sorted({"%s [%s]" % (s["other_species_match"], s["other_species_shared_taxon"])
                                               for s in c["receptors"] if s.get("other_species_match")}))
        seq_by_id = c.get("sequences", {})
        test = [(s, seq_by_id[s["id"]], "study") if s["id"] in seq_by_id else (s, ours[s["our_accession"]], "ours")
                for s in c["receptors"] if s["id"] in seq_by_id or s["our_accession"] in ours]
        c["to_test"] = len(test)
        c["untested"] = len(c["receptors"]) - len(test)
        for s, seq, source in test:
            queries.append({"candidate": c, "study": c["study"], "their_family": c["their_family"],
                            "genome": c["genome_as_reported"], "id": s["id"], "sequence": seq,
                            "source": "study's sequence" if source == "study" else "our protein %s" % s["our_accession"],
                            "our_pipeline": s["stopped_at"] if s["our_accession"] else ""})
    return queries


def attach_sequences(cands):
    """The study's sequences (literature_sequences.tsv keeps identifiers only)."""
    wanted = {(c["study"], c["their_family"]) for c in cands if c["receptors"]}
    by_key = defaultdict(dict)
    for r in L.read_all():
        if r["kind"] == "sequence" and (r["study"], r["their_family"]) in wanted and r["sequence"]:
            by_key[(r["study"], r["their_family"])][r["seq_id"]] = re.sub(r"[^A-Z]", "", r["sequence"].upper())
    for c in cands:
        c["sequences"] = by_key.get((c["study"], c["their_family"]), {})


def translate_locus(nt):
    """(protein, how) for a TAS2R locus: the longest reading frame from a methionine if it is at least
    MIN_ORF long; otherwise (pseudogenes, frameshifted or partial loci) the frame with the fewest stops,
    read through with each stop as X, so the locus can still be placed."""
    from Bio.Seq import Seq
    frames = [str(Seq(nt[f:len(nt) - (len(nt) - f) % 3]).translate()) for f in range(3)]
    best = ""
    for prot in frames:
        for piece in prot.split("*"):
            if "M" in piece and len(piece) - piece.index("M") > len(best):
                best = piece[piece.index("M"):]
    if len(best) >= MIN_ORF:
        return best, "open reading frame"
    prot = min(frames, key=lambda p: p.count("*")).strip("*")
    return prot.replace("*", "X"), "read through %d stops" % prot.count("*")


def tas2r_queries(cands):
    """Itoigawa et al. 2024 TAS2R loci, intact genes and pseudogenes, translated, for candidate genomes.
    The nucleotide sequence of each locus is cached. Loci on NCBI contigs are fetched from NCBI; loci on
    the GWH assembly of Potamotrygon leopoldi are read from its genome FASTA (GPCR_GWH_GENOME, downloaded
    from the GWH, GWHAOTN00000000) at the same coordinates. Loci on ctg_... contigs (Hemitrygon akajei,
    an unpublished draft) cannot be read and are listed as such. Behrens et al. 2023 name one T2R1 per
    species and give no sequence; Itoigawa's locus in the same species is that gene, so a Behrens claim
    is tested with Itoigawa's sequence where Itoigawa has the species."""
    tables = L._docx_tables(open(L._path("itoigawa_2024", "Tables_S1-S2.docx"), "rb").read())
    cached = {}
    if os.path.exists(TAS2R_CACHE):
        cached = {h: s for h, s in L.parse_fasta_text(open(TAS2R_CACHE).read())}
    loci, gwh, unreadable = defaultdict(list), {}, []
    for species, _, gene, acc, start, end, strand in tables[1][1:]:
        locus = "%s:%s-%s(%s)" % (acc, start, end, strand)
        if locus not in cached and re.fullmatch(r"[A-Z]{1,6}_?\d+\.\d", acc):
            url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=%s&rettype=fasta"
                   "&seq_start=%s&seq_stop=%s&strand=%d" % (acc, start, end, 1 if strand == "+" else 2))
            time.sleep(1)  # NCBI allows three requests per second without an API key
            cached[locus] = "".join(urllib.request.urlopen(url, timeout=60).read().decode().splitlines()[1:]).upper()
        elif locus not in cached and acc.startswith("GWH"):
            gwh[(acc, int(start), int(end), strand)] = locus
        elif locus not in cached:
            unreadable.append((species, gene, locus))
            continue
        loci[species].append((gene, locus))
    if gwh and os.path.exists(GWH_GENOME):
        import gzip
        contigs, name = defaultdict(list), None
        with gzip.open(GWH_GENOME, "rt") as handle:  # keep only the contigs holding a locus
            for line in handle:
                if line.startswith(">"):
                    name = line[1:].split()[0]
                elif name in {k[0] for k in gwh}:
                    contigs[name].append(line.strip().upper())
        for (acc, start, end, strand), locus in gwh.items():
            nt = "".join(contigs.get(acc, []))[start - 1:end]
            if nt:
                cached[locus] = nt if strand == "+" else nt[::-1].translate(str.maketrans("ACGTN", "TGCAN"))
    with open(TAS2R_CACHE, "w") as handle:
        handle.write("".join(">%s\n%s\n" % kv for kv in sorted(cached.items())))
    queries = []
    for c in cands:
        if c["study"] not in ("Itoigawa et al. 2024", "Behrens et al. 2023"):
            continue
        species = c["genome_as_reported"]
        for gene, locus in loci.get(species, []):
            if locus not in cached:
                unreadable.append((species, gene, locus))
                continue
            seq, how = translate_locus(cached[locus])
            queries.append({"candidate": c, "study": c["study"], "their_family": c["their_family"],
                            "genome": species, "id": "%s %s" % (gene, locus), "sequence": seq,
                            "source": "Itoigawa locus translated (%s)" % how, "our_pipeline": ""})
            c["to_test"] = c.get("to_test", 0) + 1
        notes = ["locus %s %s cannot be read (contig not public)" % (g, l) for sp, g, l in unreadable if sp == species]
        if not c.get("to_test"):
            notes.append("no sequence to test")
        if notes:
            c["note"] = "; ".join(dict.fromkeys(notes))
    return queries


# ---------------------------------------------------------------- report

def write(cands, origins, matcher):
    lineage = {}
    for c in cands:
        ups = ancestors(matcher, c["species"]) if c["species"] else []
        named = [t for t in ups if matcher.rank.get(t) in ("phylum", "class", "kingdom")]
        lineage[id(c)] = named
    cols = ["study", "their_family", "our_families", "deeper_than", "genome_as_reported", "genome_matched",
            "species", "evidence", "comparison", "in_our_data", "to_test", "untested", "other_species", "placement", "note"]
    with open(os.path.join(OUT, "deeper_origin_candidates.tsv"), "w", newline="", encoding="utf-8") as handle:
        w = csv.DictWriter(handle, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for c in cands:
            w.writerow({k: c.get(k, "") for k in cols})
    by_fam = defaultdict(list)
    for c in cands:
        by_fam[c["our_families"]].append(c)
    md = ["# Literature presences beyond our origins (October 3, 2026)", "",
          "Built by `../literature_deeper_origins.py` from `literature_comparison.tsv` and",
          "`literature_sequences.tsv`. A case is a study reporting a family in a species outside our origin",
          "node (lineage figure, interpretation column) where our data do not show the family in that genome.",
          "Data readout only. Method in the script's docstring and in `../README_LITERATURE_VERIFICATION.md`.",
          "", "| Our families | Our origin (Dollo) | Cases | Genome in our set | Receptors in our data | Sequences to test |",
          "|---|---|---|---|---|---|"]
    for fams, cs in sorted(by_fam.items()):
        f0 = [f for f in fams.split(",") if f in origins]
        orig = "; ".join("%s %s (%s)" % (f, origins[f][0], origins[f][1]) for f in f0)
        md.append("| %s | %s | %d | %d | %d | %d |" % (
            fams, orig, len(cs), sum(1 for c in cs if c.get("genome_matched")),
            sum(bool(c.get("in_our_data")) for c in cs), sum(c.get("to_test", 0) for c in cs)))
    for fams, cs in sorted(by_fam.items()):
        md += ["", "## %s" % fams, ""]
        md.append("| Study | Their family | Genome | In our set | Evidence | Their receptors in our data | Other-species match | Tested | Not tested (no sequence) |")
        md.append("|---|---|---|---|---|---|---|---|---|")
        for c in sorted(cs, key=lambda c: (c["study"], c["genome_as_reported"])):
            md.append("| %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
                c["study"], c["their_family"], c["genome_as_reported"], "yes" if c.get("genome_matched") else "no",
                c["evidence"], "; ".join(x for x in (c.get("in_our_data", ""), c.get("placement", "")) if x),
                c.get("other_species", ""), c.get("to_test", 0), c.get("untested", 0) or c.get("note", "")))
    open(os.path.join(OUT, "DEEPER_ORIGINS.md"), "w", encoding="utf-8").write("\n".join(md) + "\n")


def main():
    records = L.read_all()
    # the names of the output tables as well (a receptor's genome can come from its UniProt record)
    table_names = {row["genome_as_reported"] for row in csv.DictReader(
        open(os.path.join(OUT, "literature_comparison.tsv"), encoding="utf-8"), delimiter="\t")}
    matcher = C.GenomeMatcher({r["genome"] for r in records if r["genome"]} | table_names, set())
    origins = our_origins(matcher)
    cands = candidates(matcher, origins)
    attach_sequences(cands)
    queries = sort_receptors(cands) + tas2r_queries(cands)
    print("%d candidate claims; %d sequences to test" % (len(cands), len(queries)))
    print(Counter(c["our_families"] for c in cands).most_common())
    print(Counter(q["candidate"]["our_families"] for q in queries).most_common())
    write(cands, origins, matcher)


if __name__ == "__main__":
    main()
