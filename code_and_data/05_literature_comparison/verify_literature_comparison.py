"""Checks on the literature comparison, for an independent reviewer.

Step 6. Two groups of checks, each printed as PASS or FAIL with the numbers behind it.

1. Readers against the papers. Every reader in literature_readers.py must recover the totals the
   supplementary files themselves give (counted directly from the files here, not through the
   reader), so a parsing error in a badly formatted supplement shows up as a mismatch.
   Every supplementary file in literature/ must be read by a reader or listed in
   literature_readers.FILES_NOT_READ with the reason.
2. The output tables against themselves (05_literature_comparison/). Stage counts must
   never increase along clade -> 7TM -> topology -> final, stopped_at must agree with the stage
   columns, rows outside our genome set must carry no call of ours, no genome row may say
   "absent" while listing the study's receptors, and every matched protein must be in one of our
   searched genomes and in the species of its row.

It reads only local files (no network). Run compare_literature.py first.

Usage (from this folder):
    python verify_literature_comparison.py
"""

import collections
import csv
import io
import os
import re
import zipfile

import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gpcr_paths import stage  # noqa: E402

import literature_readers as L

# checks the readers against the papers, so it reads the supplementary files themselves
if not os.path.isdir(os.path.join(L.ORIGINALS_DIR, "garcia_de_las_bayonas_king_2025")):
    raise SystemExit("Needs the studies' supplementary files: download them (literature/SOURCES.md) and set "
                     "GPCR_LITERATURE_ORIGINALS to their folder.")

OUT = str(stage("05_literature_comparison"))
failures = []


def check(name, ok, detail=""):
    print("%s  %s%s" % ("PASS" if ok else "FAIL", name, ("  [%s]" % detail) if detail else ""))
    if not ok:
        failures.append(name)


def count(records, kind=None, family=None, genome=None):
    return sum(1 for r in records if (kind is None or r["kind"] == kind)
               and (family is None or r["their_family"] == family)
               and (genome is None or r["genome"] == genome))


# ---------------------------------------------------------------- 1. readers against the papers
print("== readers against the supplementary files ==")

k12 = L.read_krishnan_2012()
t1 = pd.read_excel(L._path("krishnan_2012", "pone.0029817.s002.xls"), header=None).iloc[4:, :3]
t1.columns = ["f", "s", "a"]
t1[["f", "s"]] = t1[["f", "s"]].ffill()
t1 = t1.dropna(subset=["a"])
t1["f"] = t1["f"].astype(str).str.split("(").str[0].str.strip()
for fam, n in t1.groupby("f")["a"].nunique().items():
    check("Krishnan 2012 %s receptors = unique Table S1 accessions" % fam,
          count(k12, "sequence", fam) == n, "%d vs %d" % (count(k12, "sequence", fam), n))
t3 = pd.read_excel(L._path("krishnan_2012", "pone.0029817.s004.xls"), header=None)
fungal = sum(count(k12, "sequence", f) for f in ("STE2", "STE3", "Git3", "Nutrient receptor"))
check("Krishnan 2012 Table S3 accessions", fungal == int(t3.iloc[3:].notna().sum().sum()),
      "%d vs %d cells (the paper text says 619)" % (fungal, int(t3.iloc[3:].notna().sum().sum())))

dm14 = L.read_de_mendoza_2014()
check("de Mendoza 2014 genomes = Table S1 rows", len({r["genome"] for r in dm14}) == 78, "78 (paper text: 75 + 5)")

dm16 = L.read_de_mendoza_2016()
check("de Mendoza 2016 sequences (Text file S1)", count(dm16, "sequence") == 76, str(count(dm16, "sequence")))
cur = L.read_curated()
check("de Mendoza 2016 Table S1 yes/no rows (curated)", count([r for r in cur if r["study"].startswith("de Mendoza")]) == 39)

n09 = L.read_nordstrom_2009()
check("Nordstrom 2009 sequences (File 2) = table rows (File 3)", len(n09) == 268, str(len(n09)))
check("Nordstrom 2009 coding sequences translated (no nucleotide-only sequence left)",
      not any(set(r["sequence"]) <= set("ACGTN-") for r in n09))
n11 = L.read_nordstrom_2011()
n11c = [r for r in n11 if r["kind"] == "claim"]
n11s = [r for r in n11 if r["kind"] == "sequence"]
check("Nordstrom 2011 claims cover 4 genomes x 9 families", len(n11c) == 36)
check("Nordstrom 2011 sequence total (2,288 rows incl. header)", sum(r["count"] for r in n11c) == 2287,
      str(sum(r["count"] for r in n11c)))
check("Nordstrom 2011 every table receptor has its sequence (Additional text 2)",
      len(n11s) == 2287 and all(r["sequence"] for r in n11s),
      "%d of %d" % (sum(bool(r["sequence"]) for r in n11s), len(n11s)))

sp = L.read_krishnan_2014_sponge()
for fam, n in [("Glutamate", 33), ("Rhodopsin", 126), ("Adhesion", 40), ("Frizzled", 3)]:
    check("Krishnan 2014 (sponge) %s sequences (abstract: %d)" % (fam, n), count(sp, family=fam) == n,
          str(count(sp, family=fam)))
for sheet, fam in [("Adhesion", "Adhesion"), ("Glutamate", "Glutamate"), ("Rhodopsins", "Rhodopsin")]:
    queries = pd.read_excel(L._path("krishnan_2014_sponge", "12862_2014_270_MOESM2_ESM.xlsx"), sheet, header=None)[0]
    n = int(queries.astype(str).str.startswith("PAC:").sum())
    check("Krishnan 2014 (sponge) %s: Additional file 1 sequences = Additional file 2 queries" % fam,
          count(sp, family=fam) == n, "%d vs %d" % (count(sp, family=fam), n))

gz = L.read_grus_zhang_2009()
check("Grus & Zhang 2009 elephant shark partial V2Rs (File 2: 32)", count(gz, family="V2R") == 32,
      str(count(gz, family="V2R")))
check("Grus & Zhang 2009 V1R rows of File 9 (3 lamprey, 2 elephant shark contigs), each >= 250 aa",
      count(gz, family="V1R") == 5 and all(len(r["sequence"].replace("-", "")) >= 250
                                           for r in gz if r["their_family"] == "V1R"))

nem = L.read_krishnan_2014_nematode()
check("Krishnan 2014 (nematode) srw list (Table S1)", len(nem) == 66, str(len(nem)))

car = L.read_cardoso_2014()
s1_counts = {"I. scapularis": 11, "T. castaneum": 9, "N. vitripennis": 7, "C. elegans": 3, "D. melanogaster": 5}
for sp, n in s1_counts.items():
    check("Cardoso 2014 %s (Table S1)" % sp, count(car, genome=sp) == n, str(count(car, genome=sp)))

pr = L.read_prabhu_2006()
check("Prabhu 2006 GABAB-like (Table 1)", count(pr, family="GABAB") == 17)
check("Prabhu 2006 Frizzled-like (Table 2)", count(pr, family="Frizzled") == 25)

it = L.read_itoigawa_2024()
check("Itoigawa 2024 genomes searched (Table S1)", len(it) == 21)
check("Itoigawa 2024 genomes with a TAS2R locus (Table S2: 14 loci in 11 genomes)",
      sum(r["status"] == "present" for r in it) == 11 and sum(r["count"] for r in it) == 14)

kow = L.read_kowatschew_2022()
for sp, fam, n in [("Lethenteron camtschaticum", "V1R", 6), ("Petromyzon marinus", "V1R", 6)]:
    got = sum(1 for r in kow if r["genome"] == sp and r["their_family"] == fam)
    check("Kowatschew 2022 %s %s (abstract: 6 each)" % (sp, fam), got == n, str(got))
v2 = [sum(1 for r in kow if r["genome"] == sp and r["their_family"] == "V2R") for sp in
      ("Lethenteron camtschaticum", "Petromyzon marinus")]
check("Kowatschew 2022 lamprey V2R (abstract: 1-2)", all(1 <= n <= 2 for n in v2), str(v2))

jf = L.read_liu_2026()
check("Liu et al. 2026 genomes (Table S3 rows, 17 classes each)", count(jf, "claim") == 1357 * 17)
check("Liu et al. 2026 receptors (File S2)", count(jf, "sequence") == 28294, str(count(jf, "sequence")))
z = zipfile.ZipFile(L._path("liu_2026", "jof-12-00030-s001.zip"))
s3 = pd.read_excel(io.BytesIO(z.read("Supplementary file S1.xlsx")), sheet_name="Table S3", header=1)
check("Liu et al. 2026 Table S3 rows sum to Total",
      bool((s3[["C%d" % i for i in range(1, 18)]].sum(axis=1) == s3["Total"]).all()))

ni = L.read_niimura_2009()
ni_table1 = {r["genome"]: r["count"] for r in L.read_curated() if r["study"] == "Niimura 2009" and r["status"] == "present"}
check("Niimura 2009 OR sequences per genome (Dataset S1) = Table 1",
      collections.Counter(r["genome"] for r in ni) == collections.Counter(ni_table1),
      "%d read, Table 1 %d" % (len(ni), sum(ni_table1.values())))
check("Niimura 2009 B. floridae ORs (Table 1: 43)", count(ni, genome="Branchiostoma floridae") == 43)

fe = L.read_fernandez_figueroa_2026()
check("Fernandez-Figueroa 2026 Table S1 rows", sum(1 for r in fe if r["their_family"] != "GCR1 ortholog") == 50)
check("Fernandez-Figueroa 2026 Table S2 orthologs", count(fe, family="GCR1 ortholog") == 157)

for name, fam, n in [("Rhodopsin_Receptors.txt", "Rhodopsin", 6137), ("Glutamate_Receptors.txt", "Glutamate", 1762),
                     ("Adhesion_Receptors.txt", "Adhesion", 1965), ("Frizzled_Receptors.txt", "Frizzled", 235),
                     ("SF25_GPR157_homologs.txt", "GPR157", 250)]:
    got = len(L._garcia_file(name, fam))
    raw = open(L._path("garcia_de_las_bayonas_king_2025", name), encoding="utf-8", errors="ignore").read().count(">")
    check("Garcia 2025 %s records = FASTA headers" % name, got == raw == n, "%d / %d" % (got, raw))
got = len(L._garcia_file("SF4_choanoflagellate_GPCRs.txt", None))
check("Garcia 2025 Supplementary file 4 records", got == 905, "%d (the paper says 918)" % got)

cs = [r for r in L.read_cardoso_2014() if r["sequence"]]
print("INFO  Cardoso 2014: %d of the s007 genes carry their s014 sequence" % len(cs))

L.READ_FILES.clear()  # count only what the readers open, not the checks above
everything = L.read_all(from_originals=True)
check("the parsed tables equal what the readers take from the files", everything == L.read_all())
on_disk = {os.path.relpath(os.path.join(d, f), L.ORIGINALS_DIR).replace(os.sep, "/")
           for d, _, fs in os.walk(L.ORIGINALS_DIR) for f in fs}
on_disk |= {os.path.relpath(os.path.join(d, f), L.LIT_DIR).replace(os.sep, "/")
            for d, _, fs in os.walk(L.LIT_DIR) for f in fs if not d.startswith(L.PARSED_DIR)}
unused = on_disk - L.READ_FILES - set(L.FILES_NOT_READ)
check("every supplementary file is read or listed as not read (with the reason)", not unused, ", ".join(sorted(unused)))
both = L.READ_FILES & set(L.FILES_NOT_READ)
check("no file is both read and listed as not read", not both, ", ".join(sorted(both)))
nucleotide = [r for r in everything if r["kind"] == "sequence" and r["sequence"] and set(r["sequence"]) <= set("ACGTN-")]
check("no reader returns a nucleotide sequence", not nucleotide, str(len(nucleotide)))

# ---------------------------------------------------------------- 2. output tables
print("\n== output tables ==")
comp = list(csv.DictReader(open(os.path.join(OUT, "literature_comparison.tsv"), encoding="utf-8"), delimiter="\t"))
seqs = list(csv.DictReader(open(os.path.join(OUT, "literature_sequences.tsv"), encoding="utf-8"), delimiter="\t"))

steps = ("in_our_clade", "in_our_7TM", "in_our_topology", "in_our_final")
bad = [r for r in comp if r["receptors_listed"] and
       not all(int(r[a]) >= int(r[b]) for a, b in zip(steps, steps[1:]))]
check("comparison: stage counts never increase", not bad, "%d rows break it" % len(bad))
bad = [r for r in comp if r["comparison"] == "genome not in our set" and (r["our_taxids"] or r["our_status"])]
check("comparison: rows outside our set carry no call of ours", not bad, str(len(bad)))
bad = [r for r in comp if r["evidence"] == "count only" and r["receptors_listed"]]
check("comparison: 'count only' rows list no receptors", not bad, str(len(bad)))

bad = [r for r in comp if r["their_status"] == "absent" and r["receptors_listed"] not in ("", "0")]
check("comparison: no row says absent while listing the study's receptors", not bad, str(len(bad)))
bad = [r for r in comp if r["comparison"] == "agree: present" and r["our_status"] != "present"]
check("comparison: 'agree: present' rows are present on our side", not bad, str(len(bad)))

bad = [r for r in seqs if r["stopped_at"] in ("not found in our data", "genome not in our set") and (
    r["our_accession"] or any(r["our_" + s] for s in ("clade", "7TM", "topology", "final")))]
check("receptors: not found, or outside our genomes -> no protein and no stage columns", not bad, str(len(bad)))
bad = [r for r in seqs if r["stopped_at"] == "final in its family" and not (
    set(r["our_final"].split(",")) & set(r["our_families"].split(",")))]
check("receptors: 'final in its family' is final in a family mapped to theirs", not bad, str(len(bad)))
families = lambda value: set(value.split(",")) - {"", "-", "SLT_TRNS"}  # SLT_TRNS is a background
bad = [r for r in seqs if r["stopped_at"].startswith("removed by the 7TM") and (not families(r["our_clade"]) or families(r["our_7TM"]))]
check("receptors: '7TM filter' stops are in a clade and not past 7TM", not bad, str(len(bad)))
bad = [r for r in seqs if r["stopped_at"].startswith("removed by the topology") and (not families(r["our_7TM"]) or families(r["our_topology"]))]
check("receptors: 'topology filter' stops passed 7TM and not topology", not bad, str(len(bad)))

# every matched protein is from our searched genomes, in the species its row names
from compare_literature import GenomeMatcher
matched = [r for r in seqs if r["our_accession"]]
M = GenomeMatcher(set(), {r["our_taxid"] for r in matched})
ours = {t for ids in M.ours.values() for t in ids}
no_hit = sum(r["stopped_at"] == "no search hit" for r in matched)
check("receptors: matched proteins are hit by our searches (our sequences come from them)",
      no_hit <= 0.01 * len(matched), "%d of %d matched with no search hit" % (no_hit, len(matched)))
bad = [r for r in matched if r["our_taxid"] not in ours]
check("receptors: every matched protein is from one of our searched genomes", not bad, str(len(bad)))
bad = [r for r in matched if M.sci.get(M.species_of(r["our_taxid"])) != r["genome"]]
check("receptors: every matched protein is in the species of its row", not bad,
      "%d, e.g. %s" % (len(bad), "; ".join("%s %s" % (r["id"][:30], r["genome"]) for r in bad[:3])))
print("INFO  matched by:", dict(collections.Counter(re.sub(r"^(identifier|accession|similar sequence|sister genome) .*", r"\1", r["match"]) for r in matched)))

print("\nstopped_at:", dict(collections.Counter(re.sub(r" \(.*", "", r["stopped_at"]) for r in seqs).most_common()))
print("\n%d checks failed" % len(failures) if failures else "\nall checks passed")
