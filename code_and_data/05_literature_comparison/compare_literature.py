"""Published GPCR family distributions against ours, in one table.

Stage 05, comparison with published classifications. literature_readers.py reads what each study reports: a family
present in a genome (with a count when given) or absent from it (a 0 or "not identified" in the
study's own table or text), and the receptors it lists. This script sets every such report
beside our call for the same family in the same genome, and follows every listed receptor, of
any family, through our pipeline. Nothing a study does not report enters the tables.

Genomes are matched by name through the NCBI taxonomy (names.dmp: scientific names and
synonyms; nodes.dmp: every strain or subspecies is taken up to its species), so a study's
"Allomyces macrogynus" or "A. macrogynus" meets our strain-level taxID. Abbreviated names
("C. elegans") are expanded against our genomes and the full names in the studies. Claims a
study makes for several strains of one species are combined (present if any strain has it).

Receptors are matched to our own data only, never to the current UniProt: the database we
searched (UniProt reference proteomes, release 2025_02) is the reference.
    by accession   first: a UniProt accession the study gives that is in our search headers is
                   that protein, and its genome is the receptor's genome whatever the species name
    by sequence    against the full-length sequences of our clade members and of the sequences
                   each search tree was built from (the archive), in the species the study names:
                   exact; else one sequence inside the other (>= 150 aa); else, for alignment or
                   contig sequences, every gap-free piece of >= 30 aa inside the same protein
    by identifier  against the headers of every hit our searches report (accession, entry name,
                   GN=), in the species the study names
and are then followed through our stages:

    clade      in a family's search-tree clade          (02_clade_extraction/clade_membership.tsv)
    7TM        a single GPCR-like 7TM core              (03_topology_filter/core_membership.tsv)
    topology   kept by the topology filter at 0.8       (03_topology_filter/topology_membership.tsv)
    final      final family membership                  (04_overlap_assignment/final_membership.tsv)

Writes to 05_literature_comparison/ in the archive:
    literature_comparison.tsv        one row per study, family and genome
    literature_sequences.tsv         one row per listed receptor and where it ends in our pipeline
    literature_summary.tsv           rows per study, family and comparison outcome
    literature_unmapped_families.tsv the study families with no counterpart: where their receptors end
    literature_unresolved_genomes.tsv genome names not matched to ours, with our closest species (review only)
    literature_genome_proxies.tsv    study genomes not in our set compared in a sister genome of ours (same
                                     genus), with the receptor evidence for each

Usage:
    python compare_literature.py
"""

import csv
import glob
import gzip
import hashlib
import json
import os
import sys
import re
import unicodedata
from collections import Counter, defaultdict
from difflib import SequenceMatcher, get_close_matches

from literature_readers import read_all

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gpcr_paths import open_text, stage  # noqa: E402

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = str(stage("05_literature_comparison"))
# Stage 06 holds the taxonomy (the searched genomes and the NCBI taxdump), the presence table and
# the lineage figure this stage reads.
ORIGINS_DIR = str(stage("06_family_origins"))
NCBI_DIR = str(stage("06_family_origins", "taxonomy", "ncbi_taxdump"))
SEARCH_DIR = str(stage("01_homology_search"))
# Full-length sequences of our clade members (02_clade_extraction/clade_sequences/<family>.fasta.gz)
# and of the hits each search tree was built from (01_homology_search/<search>/<search>_*db.fasta.gz).
SEQUENCE_FILES = (str(stage("02_clade_extraction", "clade_sequences", "*.fasta.gz")),
                  str(stage("01_homology_search", "*", "*_db.fasta.gz")))
CONTROLS = ("Galpha", "Gbeta", "Ggamma", "RBP1", "TBP")

STAGES = [
    ("clade", str(stage("02_clade_extraction", "clade_membership.tsv"))),
    ("7TM", str(stage("03_topology_filter", "core_membership.tsv"))),
    ("topology", str(stage("03_topology_filter", "topology_membership.tsv"))),
    ("final", str(stage("04_overlap_assignment", "final_membership.tsv"))),
]
# The searches themselves (stage 01): HMMER output at E <= 1000, and the MMseqs2 tables of the
# families searched with a seed sequence. classB1/classB2 share one 7tm_2 search and
# Vomeronasal1 uses the class T (TAS2R) search, so each search is listed once.
SEARCHES = [
    ("classA (7tm_1)", "classA/classA_hmmsearch_E1000.txt"), ("Olfactory", "Olfactory/Olfactory_hmmsearch_E1000.txt"),
    ("Nematode", "Nematode/Nematode_hmmsearch_E1000.txt"), ("classB1/B2 (7tm_2)", "classB/classB_hmmsearch_E1000.txt"),
    ("Mth", "Mth/Mth_hmmsearch_E1000.txt"), ("classC", "classC/classC_hmmsearch_E1000.txt"),
    ("Vomeronasal2", "Vomeronasal2/Vomeronasal2_hmmsearch_E1000.txt"), ("classF", "classF/classF_hmmsearch_E1000.txt"),
    ("classT/Vomeronasal1 (TAS2R)", "classT/classT_hmmsearch_E1000.txt"), ("cAMP", "cAMP/cAMP_hmmsearch_E1000.txt"),
    ("GP143", "GP143/GP143_hmmsearch_E1000.txt"), ("GP157", "GP157/GP157_hmmsearch_E1000.txt"),
    ("TM116", "TM116/TM116_hmmsearch_E1000.txt"), ("STE2", "STE2/STE2_hmmsearch_E1000.txt"),
    ("STE3", "STE3/STE3_hmmsearch_E1000.txt"),
    ("FSLB (MMseqs2)", "FSLB_DICDI/FSLB_DICDI_UniProt_s7.5_filter0_eprofile0.1_db.m8"),
    ("GPR1 (MMseqs2)", "GPR1_YEAST/GPR1_YEAST_UniProt_s7.5_filter0_eprofile0.1_db.m8"),
]
# The families each search feeds; a receptor hit by a search but in none of its families' clades
# is "hit but not assigned" for that search.
SEARCH_FAMILIES = {
    "classA (7tm_1)": ["classA"], "Olfactory": ["Olfactory"], "Nematode": ["Nematode"],
    "classB1/B2 (7tm_2)": ["classB1", "classB2"], "Mth": ["Mth"], "classC": ["classC"],
    "Vomeronasal2": ["Vomeronasal2"], "classF": ["classF"], "classT/Vomeronasal1 (TAS2R)": ["classT", "Vomeronasal1"],
    "cAMP": ["cAMP"], "GP143": ["GP143"], "GP157": ["GP157"], "TM116": ["TM116"], "STE2": ["STE2"],
    "STE3": ["STE3"], "FSLB (MMseqs2)": ["FSLB"], "GPR1 (MMseqs2)": ["GPR1"],
}
# The search tree each family's clade is read from (step 1, as resolved by
# clade_figures/build_figure_data.py). A hit becomes a family member only if it is a leaf of
# this tree, so a hit missing from it was cut by the tree's sequence limit. Class A and
# Olfactory are sampled into their trees and take their members from complete lists instead.
TREE_FILES = {
    "STE3": "STE3_hmm_E1000_db_hmmregions_FAMSA_gt0.1_hmmordered.tree",
    "STE2": "STE2_hmm_E1000_db_hmmregions_FAMSA_gt0.1_hmmordered.tree",
    "Vomeronasal1": "classT_hmm_E1000_db_FAMSA_gt0.1_hmmordered.tree",
    "Vomeronasal2": "Vomeronasal2_hmm_E1000_db_hmmregions_FAMSA_gt0.1_hmmordered.tree",
    "Nematode": "Nematode_hmm_E1000_db_FAMSA_gt0.1_hmmordered.tree",
    "Mth": "Mth_hmm_E1000_db_FAMSA_gt0.1_hmmordered.tree",
    "GPR1": "GPR1_YEAST_db_FAMSA_gt0.1_m8ordered.tree",
    "classB1": "classB_hmm_E1000_db_hmmregions_FAMSA_gt0.1_hmmordered.tree",
    "classB2": "classB_hmm_E1000_db_hmmregions_FAMSA_gt0.1_hmmordered.tree",
    "classC": "classC_hmm_E1000_db_hmmregions_FAMSA_gt0.1_hmmordered.tree",
    "classF": "classF_hmm_E1000_db_FAMSA_gt0.1_hmmordered.tree",
    "classT": "classT_hmm_E1000_db_FAMSA_gt0.1_hmmordered.tree",
    "cAMP": "cAMP_hmm_E1000_db_hmmregions_FAMSA_gt0.1_hmmordered.tree",
    "GP143": "GP143_hmm_E1000_db_FAMSA_gt0.1_hmmordered.tree",
    "TM116": "TM116_hmm_E1000_db_hmmregions_FAMSA_gt0.1_hmmordered.tree",
    "GP157": "GP157_hmm_E1000_db_hmmregions_FAMSA_gt0.1_hmmordered.tree",
    "FSLB": "FSLB_DICDI_db_FAMSA_gt0.1_m8ordered.tree",
}
SAMPLED_TREES = {"classA", "Olfactory"}
COUNTS = str(stage("04_overlap_assignment", "tax_counts_per_family.tsv"))
TAXONOMY = str(stage("06_family_origins", "taxonomy", "taxonomy_ncbi_searched.json"))

# Nested families are removed from the broader family as in step 1 (write_baselines.py); the
# 7TM alignments are built before this removal.
NESTED_MODEL_EXCLUSIONS = {"classA": ["GPR1", "Vomeronasal1"], "classC": ["Vomeronasal2"]}

# A study's family name -> our families that answer for it. The merged families count as one
# (GPR1 with class A, Mth with class B2, FSLB with class F). Pfam domains (de Mendoza 2014) take
# every family of ours built on that domain. A study family with no counterpart among ours is not
# compared and is listed in the summary.
FAMILY_MAP = {
    # The studies' Rhodopsin and Glutamate are the broad GRAFS families; our lineages nested in
    # class A (Olfactory, Nematode, class T, Vomeronasal 1) and class C (Vomeronasal 2) answer for them.
    "Rhodopsin": ["classA", "GPR1", "Olfactory", "Nematode", "classT", "Vomeronasal1"],
    "Class A - Rhodopsin": ["classA", "GPR1", "Olfactory", "Nematode", "classT", "Vomeronasal1"],
    "7tm1": ["classA", "GPR1", "Olfactory"],
    "Adhesion": ["classB2", "Mth"], "Secretin": ["classB1"],
    "Adhesion/Secretin": ["classB1", "classB2", "Mth"], "Secretin-like": ["classB1", "classB2", "Mth"],
    "Class B - Secretin like/Adhesion": ["classB1", "classB2", "Mth"], "7tm2": ["classB1", "classB2", "Mth"],
    "Glutamate": ["classC", "Vomeronasal2"], "GABAB": ["classC"],
    "Class C - Metabotropic Glutamate": ["classC", "Vomeronasal2"],
    "7tm3": ["classC", "Vomeronasal2"],
    "Frizzled": ["classF", "FSLB"], "Class F - Frizzled": ["classF", "FSLB"],
    "Class F - Smoothened": ["classF", "FSLB"],
    "cAMP": ["cAMP"], "Class E": ["cAMP"], "C5": ["cAMP"],
    "GPR157": ["GP157"], "OA1": ["GP143"], "Occular albinism": ["GP143"],
    "Git3": ["GPR1"], "C3": ["GPR1"],
    "STE2": ["STE2"], "C1": ["STE2"], "Class D": ["STE2"],
    "STE3": ["STE3"], "C2": ["STE3"],
    "Mth": ["Mth"], "srw": ["Nematode"], "Olfactory": ["Olfactory"],
    "V1R": ["Vomeronasal1"], "V2R": ["Vomeronasal2"], "T2R": ["classT"],
    "GCR1 (canonical GPCR)": ["cAMP", "GP157"], "GCR1 ortholog": ["cAMP", "GP157"],
    # their claim is that these are not GPCRs: any family of ours holding one disagrees
    "plant candidate, not a canonical GPCR": ["*"],
}
NOT_A_GPCR = "plant candidate, not a canonical GPCR"
# Study families with no named counterpart, mapped by where their receptors end: every family of
# ours that holds at least one of their final receptors (no threshold). Study families with no
# receptor final in any family of ours stay uncompared; literature_unmapped_families.tsv shows where
# the receptors of each end.
INFERRED_MAP = {
    # (final receptors per family of ours in the October 3, 2026 build; literature_unmapped_families.tsv)
    "C7": ["FSLB"],                      # Liu et al. class 7 (GprM group, their Table S2): FSL 664
    "C16": ["classC"],                   # Liu et al. class 16 (no reference protein): class C 337
    "C15": ["GPR1", "cAMP", "classA"],   # Liu et al. class 15 (no reference protein): GPR1 220, cAMP 148, class A 34
    "C17": ["FSLB", "cAMP"],             # Liu et al. class 17 (no reference protein): FSL 69, cAMP 18
    "7tm_2_other": ["classB2", "Mth"],   # Nordstrom 2011: Mth 26, class B2 7
}
FAMILY_MAP.update(INFERRED_MAP)


# ---------------------------------------------------------------- genomes

def _plain(text):
    return unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode().lower().strip()


class GenomeMatcher:
    """Matches a study's genome name or a taxID to our searched genomes at species level."""

    def __init__(self, labels, taxids):
        self.parent, self.rank = {}, {}
        with open_text(os.path.join(NCBI_DIR, "nodes.dmp")) as handle:
            for line in handle:
                f = line.split("\t|\t", 3)
                self.parent[f[0]] = f[1]
                self.rank[f[0]] = f[2]
        records = json.load(open(TAXONOMY))
        self.ours = defaultdict(list)
        # a taxID of ours that the taxonomy files no longer hold (Guyanagaster necrorhizus 1450536)
        # is placed by the species name our taxonomy table gives it
        self.taxid_alias, lost = {}, {}
        for r in records:
            species = self.species_of(str(r["taxID"]))
            if species:
                self.ours[species].append(str(r["taxID"]))
            elif r.get("Species") not in (None, "", "NA"):
                lost[str(r["taxID"])] = _plain(r["Species"])
        our_ids = {t for ids in self.ours.values() for t in ids}
        wanted_names = set(lost.values())
        for label in labels:
            for form in self._forms(label):
                wanted_names.add(form)
        wanted_ids = set(taxids) | set(self.ours) | our_ids
        self.by_name, self.sci = defaultdict(set), {}
        with open_text(os.path.join(NCBI_DIR, "names.dmp")) as handle:
            for line in handle:
                f = line.split("\t|\t")
                name = _plain(f[1])
                if name in wanted_names:
                    self.by_name[name].add(f[0])
                if (f[0] in wanted_ids or name in lost.values()) and f[3].startswith("scientific name"):
                    self.sci[f[0]] = f[1]
        for taxid, name in lost.items():
            species = {self.species_of(t) for t in self.by_name.get(name, ())} - {None}
            if len(species) == 1:
                s = species.pop()
                self.taxid_alias[taxid] = s
                self.ours[s].append(taxid)
                self.sci.setdefault(taxid, self.sci.get(s, name))
                print("our taxID %s is not in the taxonomy files; placed by name in %s (%s)" % (taxid, s, self.sci.get(s)))
        # full binomials to expand abbreviations ("C. elegans") against
        self.binomials = {self.sci[s] for s in self.ours if s in self.sci}
        self.binomials |= {l for l in labels if re.fullmatch(r"[A-Z][a-z]+ [a-z]+", l)}
        # an abbreviation can expand to one of our species that no study spells out in full
        self.our_species_by_name = {_plain(self.sci[s]): s for s in self.ours if s in self.sci}

    def species_of(self, taxid):
        if taxid in getattr(self, "taxid_alias", {}):
            return self.taxid_alias[taxid]
        seen = 0
        while taxid in self.rank and seen < 60:
            if self.rank[taxid] == "species":
                return taxid
            if taxid == self.parent.get(taxid):
                return None
            taxid, seen = self.parent[taxid], seen + 1
        return None

    @staticmethod
    def _forms(label):
        label = re.sub(r"\s+", " ", label.replace("_", " ")).strip()
        inline = re.sub(r"\s*\([^)]*\)", "", label)  # "Blastobotrys (Arxula) adeninivorans"
        forms = [_plain(label), _plain(inline), _plain(re.sub(r"\s*\(.*", "", label))]
        words = inline.split()
        if len(words) >= 2:
            forms.append(_plain(" ".join(words[:2])))
        return [f for f in dict.fromkeys(forms) if f]

    def _expand(self, label):
        m = re.fullmatch(r"([A-Z])\.\s*([a-zA-ZÀ-ſ]+)", label.strip())
        if not m:
            return label
        initial, epithet = m.group(1), _plain(m.group(2))
        options = {b for b in self.binomials if b.startswith(initial) and _plain(b.split()[1]) == epithet}
        if not options:
            near = {b for b in self.binomials if b.startswith(initial)}
            close = get_close_matches(epithet, [_plain(b.split()[1]) for b in near], n=1, cutoff=0.85)
            options = {b for b in near if close and _plain(b.split()[1]) == close[0]}
        return options.pop() if len(options) == 1 else label

    def match(self, label="", taxid=""):
        """(species taxID or None, our taxIDs, how it was matched, display name)."""
        if taxid:
            species = self.species_of(taxid)
            if species:
                return species, self.ours.get(species, []), "taxID", self.sci.get(species, taxid)
        if not label:
            return None, [], "unresolved", ""
        expanded = self._expand(label)
        for form in self._forms(expanded):
            species = {self.species_of(t) for t in self.by_name.get(form, ())} - {None}
            if not species and form in self.our_species_by_name:
                species = {self.our_species_by_name[form]}
            if len(species) == 1:
                s = species.pop()
                if not self.ours.get(s):
                    strain = self._strain_match(label)
                    if strain:
                        return strain, self.ours[strain], "strain (NCBI names the species differently)", self.sci[strain]
                how = "abbreviation" if expanded != label else "name"
                return s, self.ours.get(s, []), how, self.sci.get(s, expanded)
        s = self._strain_match(label)
        if s:
            return s, self.ours[s], "strain (NCBI names the species differently)", self.sci[s]
        s = self._spelling_variant(expanded)
        if s:
            return s, self.ours.get(s, []), "name (spelling variant)", self.sci[s]
        return None, [], "unresolved", label

    def _strain_match(self, label):
        """One of our genomes whose NCBI name carries the strain the study names, for a genome NCBI
        has since renamed ("Penicillium chrysogenum Wisconsin 54-1255" is P. rubens Wisconsin
        54-1255). The strain (the words after the binomial, or their tail) must have >= 4
        characters including a digit and appear, word for word, in exactly one of our genomes.
        A strain that is only a number ("1813" of CIRM-BRFM 1813) counts only within the same
        genus: collection numbers repeat across collections (Coemansia sp. RSA 1813)."""
        tokens = lambda text: re.findall(r"[a-z0-9]+", _plain(text))
        words = tokens(re.sub(r"\s*\([^)]*\)", "", label))
        genus = words[0] if words else ""
        extra = [w for w in words[2:] if w not in ("var", "f", "sp", "str", "strain", "subsp")]
        if not hasattr(self, "_our_names"):
            self._our_names = [(tokens(self.sci[t]), s) for s, ids in self.ours.items() for t in ids if t in self.sci]
        for i in range(len(extra)):
            strain = extra[i:]
            text = "".join(strain)
            if len(text) < 4 or not re.search(r"\d", text):
                continue
            same_genus = text.isdigit()
            hits = {s for name, s in self._our_names
                    if (not same_genus or name[0] == genus)
                    and any(name[j:j + len(strain)] == strain for j in range(len(name) - len(strain) + 1))}
            if hits:
                return hits.pop() if len(hits) == 1 else None
        return None

    @staticmethod
    def _stem(epithet):
        return re.sub(r"(ius|iae|ae|us|um|is|a|i|e)$", "", _plain(epithet).replace("-", "").rstrip("."))

    def _spelling_variant(self, label):
        """One of our species that a name NCBI does not know spells differently: same genus, and
        the epithet equal once Latin endings, hyphens and one typo are set aside ("Aspergillus
        nomius" = A. nomiae, "Volvox cartieri" = V. carteri), or cut at a hyphen ("Pyrenophora
        tritici" = P. tritici-repentis). Only a single candidate counts."""
        words = re.sub(r"\b(cf|aff|afn)\.\s*", "", label).split()
        if len(words) < 2:
            return None
        genus, epithet = _plain(words[0]), words[1]
        if _plain(epithet).rstrip(".") == "sp":
            return None  # "sp." names no particular species
        options = {}
        for name, s in self.our_species_by_name.items():
            parts = name.split()
            if len(parts) >= 2 and parts[0] == genus:
                options.setdefault(parts[1], set()).add(s)
        stems = {self._stem(e): e for e in options}
        close = get_close_matches(self._stem(epithet), list(stems), n=2, cutoff=0.85)
        hits = {s for c in close for s in options[stems[c]]}
        if not hits:
            hits = {s for e, ss in options.items() if e.startswith(_plain(epithet).rstrip(".") + "-") for s in ss}
        if not hits and not options:
            # a misspelt genus with the epithet exact ("Phyomyces blakesleeanus" = Phycomyces,
            # "Thalassiossira pseudonana" = Thalassiosira)
            genera = {name.split()[0] for name in self.our_species_by_name}
            for g in get_close_matches(genus, list(genera), n=2, cutoff=0.85):
                hits |= {s for name, s in self.our_species_by_name.items()
                         if name.split()[0] == g and len(name.split()) >= 2 and name.split()[1] == _plain(epithet)}
        return hits.pop() if len(hits) == 1 else None


# ---------------------------------------------------------------- our side

def read_stage(path):
    """accession -> set of families, from a wide table (one column per family), or from a folder
    of per-family alignments (<family>_genes_*.fasta)."""
    members = defaultdict(set)
    if os.path.isdir(path):
        by_family = defaultdict(set)
        for name in os.listdir(path):
            if not name.endswith(".fasta") or "_genes_" not in name:
                continue
            family = name.split("_genes_")[0]
            with open(os.path.join(path, name)) as handle:
                for line in handle:
                    if line.startswith(">"):
                        acc = line[1:].split()[0]
                        by_family[family].add(acc.split("|")[1] if "|" in acc else acc)
        for broad, specific in NESTED_MODEL_EXCLUSIONS.items():
            for family in specific:
                by_family[broad] -= by_family[family]
        for family, accs in by_family.items():
            for acc in accs:
                members[acc].add(family)
        return members
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            for family, acc in row.items():
                if acc:
                    members[acc].add(family)
    return members



def read_tree_members(wanted):
    """family -> (accessions in wanted that are leaves of its search tree, number of leaves)."""
    members = {}
    for family, name in TREE_FILES.items():
        search = name.split("_hmm_E1000")[0].split("_db_FAMSA")[0]
        text = open(os.path.join(SEARCH_DIR, search, name), encoding="utf-8", errors="ignore").read()
        leaves = re.findall(r"(?:sp|tr)\|([A-Z0-9]+)\|", text)
        members[family] = ({a for a in leaves if a in wanted}, len(leaves))
    return members

# ---------------------------------------------------------------- matching to our data

MIN_CONTAINED = 150  # aa: one sequence inside the other (GPCRs are mostly over 200 aa)
MIN_PIECE = 30      # aa: gap-free pieces of an alignment or contig sequence
ACCESSION = re.compile(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}")
HEADER = re.compile(r"(?:sp|tr)\|(\w+)\|(\w+)\s.*?\bOX=(\d+)")
GENE = re.compile(r"\bGN=(\S+)")


def _md5(seq):
    return hashlib.md5(seq.encode()).hexdigest()


def gap_free(seq):
    return re.sub(r"[^A-Z]", "", seq.upper())


def _fasta_gz(path):
    """(accession, taxID, sequence) from a gzipped FASTA with UniProt headers."""
    acc, taxid, chunks = None, None, []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if acc:
                    yield acc, taxid, "".join(chunks)
                m = re.match(r">(?:sp|tr)\|(\w+)\|.*?\bOX=(\d+)", line)
                acc, taxid, chunks = (m.group(1), m.group(2), []) if m else (None, None, [])
            elif acc:
                chunks.append(line.strip().upper())
    if acc:
        yield acc, taxid, "".join(chunks)


def our_sequences(species_wanted, md5_wanted, matcher, accessions_wanted=()):
    """Our full-length sequences (zenodo_data): species -> {sequence: accessions} for the species
    in species_wanted (every species if None), MD5 -> accessions for the sequences in md5_wanted
    (receptors whose study names no species), accession -> sequence for accessions_wanted, and
    accession -> taxID for every sequence kept."""
    paths = sorted(p for pattern in SEQUENCE_FILES for p in glob.glob(pattern)
                   if not os.path.basename(p).startswith(CONTROLS))
    if not paths:
        raise SystemExit("No sequences found in the archive; set GPCR_ZENODO_DIR.")
    by_species, by_md5, taxid_of, by_acc = defaultdict(lambda: defaultdict(set)), defaultdict(set), {}, {}
    species_of, seen = {}, set()
    for path in paths:
        for acc, taxid, seq in _fasta_gz(path):
            if acc in seen:
                continue
            seen.add(acc)
            if taxid not in species_of:
                species_of[taxid] = matcher.species_of(taxid)
            keep = species_of[taxid] is not None and (species_wanted is None or species_of[taxid] in species_wanted)
            if keep:
                by_species[species_of[taxid]][seq].add(acc)
            if acc in accessions_wanted:
                by_acc[acc] = seq
            if md5_wanted and _md5(seq) in md5_wanted:
                by_md5[_md5(seq)].add(acc)
                keep = True
            if keep:
                taxid_of[acc] = taxid
    print("our sequences: %d read, %d kept from %d species" % (len(seen), len(taxid_of), len(by_species)), flush=True)
    return by_species, by_md5, taxid_of, by_acc


def match_by_sequence(seq, ours):
    """(how, accessions) for one receptor against ours = {sequence: accessions} of its species."""
    full = gap_free(seq)
    if full in ours:
        return "exact sequence", set(ours[full])
    hit = set()
    if len(full) >= MIN_CONTAINED:
        for s, accs in ours.items():
            if full in s or (len(s) >= MIN_CONTAINED and s in full):
                hit |= accs
        if hit:
            return "sequence contained", hit
    # a contig translation or alignment row: gaps and stops can stand for missing exons
    pieces = [p for p in re.split(r"[^A-WYZ]+", seq.upper()) if len(p) >= MIN_PIECE]
    if len(pieces) > 1:
        for s, accs in ours.items():
            if all(p in s for p in pieces):
                hit |= accs
        if hit:
            return "sequence pieces contained", hit
    return "", set()


MIN_SIMILAR = 0.90  # identical residues of the alignment over the length of the longer sequence
MIN_PROXY_SHARE = 0.5  # receptors (with a sequence) matched in a sister genome, to compare a genome there
MIN_PROXY_RECEPTORS = 2
SIMILAR_CANDIDATES = 20  # candidates (most shared 5-mers, at least 30% of them) aligned per receptor
KMER = 5


_ALIGNER = None


def _aligner():
    """Global alignment, BLOSUM62, gap open 11, extend 1, end gaps free (built once)."""
    global _ALIGNER
    if _ALIGNER is None:
        from Bio import Align
        from Bio.Align import substitution_matrices
        _ALIGNER = Align.PairwiseAligner()
        _ALIGNER.substitution_matrix = substitution_matrices.load("BLOSUM62")
        _ALIGNER.open_gap_score, _ALIGNER.extend_gap_score = -11, -1
        _ALIGNER.target_end_gap_score = _ALIGNER.query_end_gap_score = 0
    return _ALIGNER


def identity_over_longer(q, s):
    """Identical residues of the best global alignment over the length of the longer sequence, so
    a fragment never counts as near-identical to a full-length protein, in either direction."""
    a = _aligner().align(q, s)[0]
    same = sum(1 for i, j in zip(*a.indices) if i >= 0 and j >= 0 and q[i] == s[j])
    return same / max(len(q), len(s))


class SimilarityIndex:
    """Near-identical matches within one species, for a receptor whose sequence is another gene
    model of one of our proteins (an edited, older or newer version): candidates by shared 5-mers,
    then the identical residues of a global alignment over the longer of the two sequences. A
    sequence much shorter or longer than ours (a fragment, a partial gene) cannot reach the
    threshold and is left to the MMseqs2 search."""

    def __init__(self, ours):
        self.seqs = list(ours.items())
        self.index = defaultdict(set)
        for i, (s, _) in enumerate(self.seqs):
            for j in range(len(s) - KMER + 1):
                self.index[s[j:j + KMER]].add(i)

    def match(self, seq):
        q = gap_free(seq)
        kmers = {q[j:j + KMER] for j in range(len(q) - KMER + 1)}
        if not kmers:
            return "", set()
        shared = Counter(i for k in kmers for i in self.index.get(k, ()))
        best, accs = 0.0, set()
        # ties (paralogs sharing as many 5-mers) are broken by position, so every run picks the same candidates
        for i, n in sorted(shared.items(), key=lambda kv: (-kv[1], kv[0]))[:SIMILAR_CANDIDATES]:
            s, a = self.seqs[i]
            if n < 0.3 * min(len(kmers), len(s)):
                break
            if min(len(q), len(s)) < MIN_SIMILAR * max(len(q), len(s)):
                continue  # the length difference alone keeps identity over the longer below the threshold
            score = identity_over_longer(q, s)
            if score > best + 1e-9:
                best, accs = score, set(a)
            elif abs(score - best) <= 1e-9:
                accs |= a
        if best >= MIN_SIMILAR:
            return "similar sequence (%.2f)" % best, accs
        return "", set()


# Identifiers our headers do not carry (dictyBase, Ensembl, NCBI gene, AGI or ORF names) are
# translated to UniProt accessions with the identifier table of the earlier UniProt-based run
# (id_cache.json); a translated accession counts only if it is in our search headers, in the
# species the study names. Set to False to use no UniProt-derived information. Kept on (decided
# 2026-10-03): only the gene-ID -> accession link is external; the protein itself must be in our
# data. When an ID gives several accessions, the one we have is what counts.
USE_XREF_TABLE = True
XREF_TABLE = os.path.join(SCRIPT_DIR, "literature", "id_cache.json")


def id_variants(token):
    """An identifier as the study writes it and the forms our headers may carry: the accession of
    a UniProt entry name (C0RZ86_PARBP -> C0RZ86), a locus tag without its transcript or protein
    suffix (TEQG_00345T0 -> TEQG_00345, AGAP004829-PA -> AGAP004829)."""
    t = token.strip().upper()
    forms = [t]
    m = re.fullmatch(r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})_[A-Z0-9]{3,5}", t)
    if m:
        forms.append(m.group(1))
    for pattern in (r"(.+_\d+)T\d+", r"(.+)-P[A-Z]", r"(.+)\.\d+"):
        m = re.fullmatch(pattern, t)
        if m:
            forms.append(m.group(1))
    return list(dict.fromkeys(forms))


def scan_searches(tokens, accessions):
    """One pass over our search outputs, whose headers come from the database we searched.
    Returns identifier -> {(accession, taxID)} for the identifiers in tokens (an accession, entry
    name or GN= gene name, upper case), accession -> labels of the searches that report it (for
    the accessions given and those identified), accession -> taxID, and the species code of
    entry names (NEUCR in Q1K6I3_NEUCR) -> taxIDs."""
    found, hits, taxid_of, codes = defaultdict(set), defaultdict(set), {}, defaultdict(set)
    for label, name in SEARCHES:
        with open_text(os.path.join(SEARCH_DIR, name)) as handle:
            for line in handle:
                m = HEADER.search(line)
                if not m:
                    continue
                acc, entry, taxid = m.groups()
                if "_" in entry:
                    codes[entry.rsplit("_", 1)[1].upper()].add(taxid)
                gene = GENE.search(line, m.end())
                keys = [acc, entry.upper()] + ([gene.group(1).upper()] if gene else [])
                named = [k for k in keys if k in tokens]
                for k in named:
                    found[k].add((acc, taxid))
                if named or acc in accessions:
                    hits[acc].add(label)
                    taxid_of[acc] = taxid
        print("searched: %s" % label, flush=True)
    return found, hits, taxid_of, codes


def explain(r, fams, acc, ours, stage_fams, trees, searched):
    """The step at which a listed receptor leaves our pipeline, and a one-line account of every
    step it passes: genome, match to our data, searches, its family's search and tree, clade,
    7TM filter, topology filter, overlap resolution and final family. fams is None for a study
    family with no counterpart, ["*"] for "not a GPCR"; then any family of ours counts."""
    if not ours:
        if r["named"] or acc:
            return "genome not in our set", "Genome not in our set."
        if r["genome"]:
            return "genome not in our set (name not in the NCBI taxonomy)", (
                "The species name is in neither the NCBI taxonomy nor our genome set, and the sequence matches "
                "no protein of ours exactly.")
        return "no species given; not found in our data", (
            "The study gives no species, and neither the sequence (exact) nor an identifier matches our data.")
    if not acc:
        return "not found in our data", ("Genome in our set; neither the sequence nor an identifier matches a protein in "
                                         "our clades, our search trees or our search hits.")
    parts = ["genome in our set", "matched to %s (%s)" % (acc, r["match"])]
    if not searched:
        return "no search hit", "; ".join(parts + ["hit by none of our searches"]) + "."
    target = None if fams in (None, ["*"]) else set(fams)
    own = [s for s in searched if target is None or set(SEARCH_FAMILIES[s]) & target]
    others = [s for s in searched if s not in own]
    stopped = None
    parts.append("hit by %d of our searches" % len(searched))
    if target is not None:
        if own:
            parts.append("including its family's search (%s)" % ", ".join(own))
        else:
            parts.append("but not by its family's search")
            stopped = "not hit by its family's search"
        if own:
            with_tree = [f for f in sorted(target) if f in TREE_FILES and f in {x for s in own for x in SEARCH_FAMILIES[s]}]
            leaf = [f for f in with_tree if acc in trees[f][0]]
            if leaf:
                parts.append("a leaf of the %s tree" % "/".join(leaf))
            elif with_tree and not target & SAMPLED_TREES:
                parts.append("not among the %s sequences of the %s tree" % (
                    "/".join(str(trees[f][1]) for f in with_tree), "/".join(with_tree)))
                stopped = stopped or "outside its family's tree (sequence limit)"
            elif target & SAMPLED_TREES:
                parts.append("%s membership from the complete list (sampled tree)" % "/".join(sorted(target & SAMPLED_TREES)))
    same = (lambda f: bool(set(f) - {"SLT_TRNS"})) if target is None else (lambda f: bool(set(f) & target))
    steps = [("clade", "in the %s clade", "not in its family's clade" if target else "in none of our clades"),
             ("7TM", "passed the 7TM filter", "removed by the 7TM filter"),
             ("topology", "passed the topology filter (0.8)", "removed by the topology filter (0.8)"),
             ("final", "kept after overlap resolution: final %s", "removed at overlap resolution")]
    for stage, yes, no in steps:
        if same(stage_fams[stage]):
            fam = ",".join(sorted(f for f in stage_fams[stage] if (target is None and f != "SLT_TRNS") or
                                  (target is not None and f in target)))
            parts.append(yes % fam if "%s" in yes else yes)
        else:
            parts.append(no)
            stopped = stopped or no
            break
    if target is not None:
        elsewhere = sorted(set(stage_fams["final"]) - target - {"SLT_TRNS"})
        if elsewhere and not same(stage_fams["final"]):
            parts.append("assigned to our %s family instead" % ",".join(elsewhere))
            stopped = "assigned to another family of ours (%s)" % ",".join(elsewhere)
    if others:
        parts.append("also hit by: %s" % ", ".join(others))
    if stopped is None:
        stopped = ("final in its family" if target is not None else
                   "final in our %s" % ",".join(sorted(set(stage_fams["final"]) - {"SLT_TRNS"})))
    return stopped, "; ".join(parts) + "."


def read_counts():
    counts = {}
    with open(COUNTS) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            counts[row["taxID"]] = {k[:-6]: int(v) for k, v in row.items() if k.endswith("_count")}
    return counts


def _write(name, rows, fields=None):
    with open(os.path.join(OUT_DIR, name), "w", newline="", encoding="utf-8") as handle:
        w = csv.DictWriter(handle, fieldnames=fields or list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


# ---------------------------------------------------------------- tables

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    records = read_all()
    seqs = [r for r in records if r["kind"] == "sequence"]
    claims = [r for r in records if r["kind"] == "claim" and r["their_family"] in FAMILY_MAP]
    unmapped_claims = Counter((r["study"], r["their_family"]) for r in records
                              if r["kind"] == "claim" and r["their_family"] not in FAMILY_MAP)
    print("%d records: %d claims compared, %d receptors traced" % (len(records), len(claims), len(seqs)), flush=True)

    matcher = GenomeMatcher({r["genome"] for r in records if r["genome"]}, set())
    species_cache = {}

    def species_of(taxid):
        if taxid not in species_cache:
            species_cache[taxid] = matcher.species_of(taxid)
        return species_cache[taxid]

    # Our sequences of every species are read, so a receptor whose genome is learnt from its
    # accession or from its study's other receptors can still be matched by sequence.
    for r in seqs:
        r["named"] = matcher.match(r["genome"])[0] if r["genome"] else None
        r["accs"], r["match"] = set(), ""
        r["genome_how"] = "manual decision (literature/genome_name_decisions.tsv)" if r.get("genome_decision") else ""
    md5_wanted = {_md5(gap_free(r["sequence"])) for r in seqs if r["sequence"] and not r["named"]}
    acc_ids = {v for r in seqs for t in r["ids"] for v in id_variants(t) if ACCESSION.fullmatch(v)}
    by_species, by_md5, taxid_of, by_acc = our_sequences(None, md5_wanted, matcher, acc_ids)
    species_wanted = set(by_species)

    xref = {}
    if USE_XREF_TABLE and os.path.exists(XREF_TABLE):
        xref = {k.upper(): {a for a, _ in v} for k, v in json.load(open(XREF_TABLE)).items()}
    tokens = {v for r in seqs for t in r["ids"] for v in id_variants(t)}
    tokens |= {a for r in seqs for t in r["ids"] for a in xref.get(t.upper(), ())}
    found, hits, more, codes = scan_searches(tokens, set(taxid_of))
    taxid_of.update(more)

    # 1. by accession, first: a UniProt accession the study gives (alone or inside an entry name)
    # that is in our search headers is that protein, and its genome is the receptor's genome,
    # whatever name the study gives the species (names change, accessions stay).
    renamed = Counter()
    for r in seqs:
        for v in (v for t in r["ids"] for v in id_variants(t) if ACCESSION.fullmatch(v)):
            pairs = {p for p in found.get(v, ()) if p[0] == v}
            if not pairs:
                continue
            acc, taxid = pairs.pop()
            species = species_of(taxid)
            r["accs"], r["match"] = {acc}, "accession %s" % acc
            if r["sequence"] and acc in by_acc:
                q = gap_free(r["sequence"])
                if not (q in by_acc[acc] or by_acc[acc] in q):
                    r["match"] += " (sequence differs)"
            if species != r["named"]:
                renamed[(r["study"], r["genome"], matcher.sci.get(species, taxid))] += 1
                r["named"], r["genome_how"] = species, "accession"
            break
    for (study, genome, name), n in sorted(renamed.items()):
        print("genome from accession: %s | %r -> %s (%d)" % (study, genome, name, n))

    # 2. by sequence, in the species the study names (or anywhere, exact only, if it names none)
    def by_sequence(r):
        if r["accs"] or not r["sequence"]:
            return
        if r["named"] in species_wanted:
            r["match"], r["accs"] = match_by_sequence(r["sequence"], by_species[r["named"]])
        elif not r["named"] and by_md5.get(_md5(gap_free(r["sequence"]))):
            r["match"], r["accs"] = "exact sequence", set(by_md5[_md5(gap_free(r["sequence"]))])
    for r in seqs:
        by_sequence(r)

    # 3. a genome name the taxonomy does not resolve to one of ours takes the genome of the study's
    # own receptors under that name, when every matched one is in the same species of ours
    # ("Endogone sp FLAS 59071" = Endogone sp. FLAS-F59071). Receptors listed with no species and
    # still unmatched take the species of their entry-name code in our headers (Q8NIR1_NEUCR).
    label_species = defaultdict(set)
    for r in seqs:
        if r["genome"] and r["accs"] and not matcher.ours.get(matcher.match(r["genome"])[0] or ""):
            label_species[(r["study"], r["genome"])] |= {species_of(taxid_of.get(a, "")) for a in r["accs"]}
    genome_alias = {k: v.pop() for k, v in label_species.items() if len(v) == 1 and None not in v}
    # taxIDs of the UniProt accessions (and UniParc IDs) listed with no species, from the UniProt
    # entry or, for retired entries, its last version in UniSave (literature/uniprot_accession_taxa.json)
    acc_taxa_path = os.path.join(SCRIPT_DIR, "literature", "uniprot_accession_taxa.json")
    acc_taxa = json.load(open(acc_taxa_path)) if os.path.exists(acc_taxa_path) else {}
    for (study, genome), s in sorted(genome_alias.items()):
        print("genome from the study's receptors: %s | %r -> %s" % (study, genome, matcher.sci.get(s)))
    for r in seqs:
        alias = genome_alias.get((r["study"], r["genome"]))
        if alias and not r["accs"]:
            r["named"], r["genome_how"] = alias, "the study's receptors under this name"
            by_sequence(r)
        if not r["genome"] and not r["accs"]:
            for t in r["ids"] + [r["seq_id"]]:
                hit = next((acc_taxa[v] for v in id_variants(t) + [re.split(r"[\s_]", t.strip())[0]]
                            if acc_taxa.get(v, {}).get("taxid")), None)
                if hit and species_of(hit["taxid"]):
                    r["named"], r["genome"] = species_of(hit["taxid"]), hit["organism"]
                    r["genome_how"] = "UniProt record of the study's accession"
                    by_sequence(r)
                    break
        if not r["genome"] and not r["accs"]:
            for t in r["ids"]:
                m = re.fullmatch(r"\w+_([A-Z0-9]{3,5})", t.strip().upper())
                species = {species_of(x) for x in codes.get(m.group(1), ())} if m else set()
                if len(species) == 1 and None not in species:
                    r["named"], r["genome_how"] = species.pop(), "entry-name species code"
                    by_sequence(r)
                    break

    # 4. by identifier, for the rest; an identifier counts only in the species of the receptor,
    # or, when it has none, only if every protein carrying it is in one species
    def identify(r, lookup, label):
        for token in r["ids"]:
            pairs = set().union(*[found.get(v, set()) for v in lookup(token)])
            if r["named"]:
                pairs = {p for p in pairs if species_of(p[1]) == r["named"]}
            elif len({species_of(t) for _, t in pairs}) > 1:
                pairs = set()
            if pairs:
                r["accs"] |= {a for a, _ in pairs}
                r["match"] = r["match"] or label % token

    unmatched = [r for r in seqs if not r["accs"]]
    in_ours = lambda r: not (r["named"] and not matcher.ours.get(r["named"]))
    for r in unmatched:
        if in_ours(r):
            identify(r, id_variants, "identifier %s")
    # 5. a near-identical sequence (another gene model of the same protein), in the named species
    similarity = {}
    for r in unmatched:
        if not r["accs"] and r["sequence"] and r["named"] in species_wanted:
            if r["named"] not in similarity:
                similarity[r["named"]] = SimilarityIndex(by_species[r["named"]])
            r["match"], r["accs"] = similarity[r["named"]].match(r["sequence"])
    # 6. identifiers our headers do not carry, translated by the earlier run's identifier table
    for r in unmatched:
        if not r["accs"] and in_ours(r):
            identify(r, lambda t: xref.get(t.upper(), ()), "identifier %s (cross-reference table)")

    # 7. receptors not placed in their own genome: an identical protein of ours in any other species,
    # or else a near-identical one (MIN_SIMILAR) in another species of the same genus, with the
    # nearest taxon the two species share. Reported only (other_species_* columns): it changes no
    # comparison, since a protein in another species says nothing about the study's genome.
    def ancestors(t):
        path = []
        while t in matcher.parent and len(path) < 80:
            path.append(t)
            if matcher.parent[t] == t:
                break
            t = matcher.parent[t]
        return path

    def genus_of(t):
        return next((a for a in ancestors(t) if matcher.rank.get(a) == "genus"), None) if t else None

    md5_all = defaultdict(set)
    for sp, d in by_species.items():
        for s, accs in d.items():
            md5_all[_md5(s)] |= {(a, sp) for a in accs}
    our_by_genus = defaultdict(set)
    for sp in by_species:
        our_by_genus[genus_of(sp)].add(sp)
    for r in seqs:
        r["other"] = None
        if r["accs"] or not r["sequence"]:
            continue
        identical = {(a, sp) for a, sp in md5_all.get(_md5(gap_free(r["sequence"])), ()) if sp != r["named"]}
        if identical:
            r["other"] = ("identical sequence", identical)
            continue
        best = (0.0, None, set())
        for sp in sorted(our_by_genus.get(genus_of(r["named"]), set()) - {r["named"]} if genus_of(r["named"]) else ()):
            if sp not in similarity:
                similarity[sp] = SimilarityIndex(by_species[sp])
            how, accs = similarity[sp].match(r["sequence"])
            if accs and float(how.split("(")[1].rstrip(")")) > best[0]:
                best = (float(how.split("(")[1].rstrip(")")), how, {(a, sp) for a in accs})
        if best[1]:
            r["other"] = (best[1], best[2])
    shared_taxon, wanted = {}, set()
    for r in seqs:
        if r["other"] and r["named"]:
            own = ancestors(r["named"])
            for _, sp in r["other"][1]:
                lca = next((a for a in ancestors(sp) if a in own), None)
                shared_taxon[(r["named"], sp)] = lca
                wanted.add(lca)
    lca_name = {}
    with open_text(os.path.join(NCBI_DIR, "names.dmp")) as handle:
        for line in handle:
            f = line.split("	|	")
            if f[0] in wanted and f[3].startswith("scientific name"):
                lca_name[f[0]] = f[1]

    def relation(named, sp):
        lca = shared_taxon.get((named, sp))
        if not lca:
            return "study gives no species"
        return "%s %s" % (matcher.rank.get(lca, "?"), lca_name.get(lca, lca))

    # 7b. a genome not in our set is compared in a sister genome of ours: the same genus, holding an
    # identical or near-identical protein for at least MIN_PROXY_SHARE of the study's receptors with
    # a sequence (and at least MIN_PROXY_RECEPTORS of them). Its claims are compared there and its
    # receptors followed through that genome's proteins. Every pairing is written, with its evidence,
    # to literature_genome_proxies.tsv, and labelled "sister genome" wherever it is used.
    our_name = lambda sp: matcher.sci.get(sp, sp)
    label_receptors = defaultdict(list)
    for r in seqs:
        if r["genome"] and not r["accs"] and not matcher.ours.get(r["named"] or ""):
            label_receptors[(r["study"], r["genome"])].append(r)
    proxy, proxy_rows = {}, []
    for (study, label), rs in sorted(label_receptors.items()):
        with_seq = [r for r in rs if r["sequence"]]
        sister_hits = defaultdict(list)
        for r in with_seq:
            for sp in {sp for _, sp in (r["other"][1] if r["other"] else ())}:
                sister_hits[sp].append(r)
        if not sister_hits:
            continue
        sp, hit = max(sorted(sister_hits.items()), key=lambda kv: len(kv[1]))
        named = with_seq[0]["named"]
        same_genus = (genus_of(named) == genus_of(sp)) if named else \
            _plain(label.split()[0]) == _plain(our_name(sp).split()[0])
        if not same_genus or len(hit) < MIN_PROXY_RECEPTORS or len(hit) < MIN_PROXY_SHARE * len(with_seq):
            continue
        proxy[(study, label)] = sp
        identical = sum(r["other"][0] == "identical sequence" for r in hit)
        proxy_rows.append({"study": study, "genome_as_reported": label, "our_genome": our_name(sp),
                           "receptors_listed": len(rs), "receptors_with_sequence": len(with_seq),
                           "identical": identical, "near_identical": len(hit) - identical,
                           "shared_taxon": relation(named, sp) if named else "genus " + label.split()[0]})
        for r in rs:
            r["proxy"], r["genome_how"] = sp, "sister genome (same genus; literature_genome_proxies.tsv)"
            if r in hit:
                r["accs"] = {a for a, s in r["other"][1] if s == sp}
                r["match"] = "sister genome %s: %s" % (our_name(sp), r["other"][0])
    if proxy_rows:
        _write("literature_genome_proxies.tsv", proxy_rows)
    print("sister genomes: %d study genomes compared in a genome of ours of the same genus" % len(proxy_rows))

    # 8. our side
    stages = [(name, read_stage(path)) for name, path in STAGES]
    counts = read_counts()
    trees = read_tree_members({a for r in seqs for a in r["accs"]})

    def reach(acc, target):
        final = stages[-1][1].get(acc, set())
        return bool(target and final & target), bool(final), sum(bool(m.get(acc)) for _, m in stages)

    seq_rows, groups = [], defaultdict(lambda: {"claims": [], "seqs": []})
    destinations = defaultdict(Counter)
    for r in claims:
        species, ours, how, name = matcher.match(r["genome"])
        alias = genome_alias.get((r["study"], r["genome"]))
        if not ours and alias:
            species, ours, how, name = alias, matcher.ours[alias], "the study's receptors under this name", matcher.sci.get(alias)
        if not ours and (r["study"], r["genome"]) in proxy:
            species = proxy[(r["study"], r["genome"])]
            ours, how, name = matcher.ours[species], "sister genome (same genus; literature_genome_proxies.tsv)", our_name(species)
        if ours and r.get("genome_decision"):
            how = "manual decision (literature/genome_name_decisions.tsv)"
        g = groups[(r["study"], r["their_family"], species or "?" + r["genome"])]
        g["claims"].append(r)
        g.setdefault("match", (species, ours, how, name or r["genome"]))
    for r in seqs:
        fams = FAMILY_MAP.get(r["their_family"])
        target = None if fams in (None, ["*"]) else set(fams)
        acc = max(sorted(r["accs"]), key=lambda a: reach(a, target)) if r["accs"] else None
        taxid = taxid_of.get(acc, "") if acc else ""
        species = r.get("proxy") or r["named"] or (species_of(taxid) if taxid else None)
        ours = matcher.ours.get(species, []) if species else []
        name = matcher.sci.get(species, r["genome"]) if species else r["genome"]
        how = (r["genome_how"] or ("name" if r["named"] else "taxID of the matched protein")) if species else "unresolved"
        stage_fams = {s: sorted(m.get(acc, set())) for s, m in stages} if acc else {}
        searched = sorted(hits.get(acc, ())) if acc else None
        stopped, explanation = explain(r, fams, acc, ours, stage_fams, trees, searched)
        r.update(acc=acc, species=species, stage_fams=stage_fams, searched=searched, ours=bool(ours))
        seq_rows.append({
            "study": r["study"], "their_family": r["their_family"],
            "our_families": "any" if fams == ["*"] else (",".join(fams) if fams else "none (no counterpart)"),
            "id": r["seq_id"], "genome": name or r["genome"], "genome_how": how,
            "sequence_source": r["sequence_source"], "our_taxid": taxid, "our_accession": acc or "",
            "other_matching_accessions": ",".join(sorted(r["accs"] - {acc})),
            "match": r["match"] or ("not found" if ours else ""),
            "our_search_hits": "" if searched is None else (",".join(searched) or "-"),
            "searches_not_assigned": "" if searched is None else (",".join(
                s for s in searched if not set(SEARCH_FAMILIES[s]) & set(stage_fams["clade"])) or "-"),
            **{"our_%s" % s: ",".join(stage_fams.get(s, [])) or ("-" if stage_fams else "") for s, _ in stages},
            "stopped_at": stopped, "explanation": explanation,
            "other_species_match": "; ".join("%s (%s)" % (a, matcher.sci.get(sp, sp)) for a, sp in sorted(r["other"][1])) if r["other"] else "",
            "other_species_how": r["other"][0] if r["other"] else "",
            "other_species_shared_taxon": "; ".join(sorted({relation(r["named"], sp) for _, sp in r["other"][1]})) if r["other"] else "",
        })
        if fams is None or r["their_family"] in INFERRED_MAP:
            d = destinations[(r["study"], r["their_family"])]
            d["listed"] += 1
            if not ours:
                d["genome not in our set" if (r["named"] or r["genome"]) else "no species given"] += 1
            elif not acc:
                d["not found in our data"] += 1
            elif set(stage_fams["final"]) - {"SLT_TRNS"}:
                for f in set(stage_fams["final"]) - {"SLT_TRNS"}:
                    d["final " + f] += 1
            elif any(stage_fams[s] for s, _ in stages):
                d["in a clade, removed by filters"] += 1
            else:
                d["in none of our clades"] += 1
        if fams is not None:
            g = groups[(r["study"], r["their_family"], species or "?" + (r["genome"] or "unknown"))]
            g["seqs"].append(r)
            g.setdefault("match", (species, ours, how, name or r["genome"]))

    # 9. one row per study, family and genome
    rows = []
    for (study, family, _), g in sorted(groups.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[1]["match"][3] or "")):
        species, ours, how, name = g["match"]
        mapped, cl, members = FAMILY_MAP[family], g["claims"], g["seqs"]
        present = [c for c in cl if c["status"] == "present"]
        if cl:
            status = "present" if present or members else "absent"
            counted = [c["count"] for c in (present or cl) if c["count"] != ""]
            count = max(counted) if counted else ""
            source = "; ".join(dict.fromkeys(c["source"] for c in cl))
            if len(cl) > 1:
                source += " (%d reports for this species, %d with the family)" % (len(cl), len(present))
            if members and not present and not any(c["status"] == "searched" for c in cl):
                source += "; receptors listed although the claim says absent"
            if not members and any(c["status"] == "searched" for c in cl):
                source += "; none of their receptors is in this genome"
        else:
            status = "absent" if family == NOT_A_GPCR else "present"
            count, source = len(members), "receptors listed in the supplement"
        reported = "; ".join(dict.fromkeys(c["genome"] for c in cl)) if cl else (
            members[0]["genome"] or members[0]["seq_id"])
        evidence = ("receptor sequences or IDs" if members else
                    next((c["evidence"] for c in cl if c.get("evidence")), "") or "count only")
        row = {"study": study, "their_family": family, "our_families": ",".join(mapped) if mapped != ["*"] else "any",
               "genome_as_reported": reported, "genome_matched": name if ours else "", "our_taxids": ",".join(ours),
               "match": how if ours else "", "their_status": status if family != NOT_A_GPCR else "not a GPCR",
               "their_count": count, "their_source": source, "evidence": evidence}
        traced = [m for m in members if m["stage_fams"]]
        if mapped == ["*"]:
            same = lambda fams: bool(set(fams) - {"SLT_TRNS"})  # SLT_TRNS is a background, not a GPCR family
        else:
            same = lambda fams: bool(set(fams) & set(mapped))
        listed = lambda value: value if members else ""
        row["receptors_listed"] = listed(len(members))
        row["receptors_in_our_genomes"] = listed(sum(m["ours"] for m in members))
        row["receptors_found_in_our_data"] = listed(len(traced))
        row["in_our_searches"] = listed(sum(bool(m["searched"]) for m in traced))
        matching = [s for s, fams in SEARCH_FAMILIES.items() if mapped == ["*"] or set(fams) & set(mapped)]
        hit = [m for m in traced if set(m["searched"]) & set(matching)]
        outside = [m for m in hit if not same(m["stage_fams"]["clade"])]
        held = Counter(f for m in outside for f in (m["stage_fams"]["clade"] or ["no clade"]))
        row["hit_by_our_matching_search"] = listed(len(hit))
        row["hit_but_not_in_our_clade"] = listed(len(outside))
        row["their_clade_instead"] = "; ".join("%s=%d" % kv for kv in sorted(held.items()))
        for s, _ in stages:
            row["in_our_%s" % s] = listed(sum(same(m["stage_fams"][s]) for m in traced))
        other = Counter(f for m in traced if not same(m["stage_fams"]["final"]) for f in m["stage_fams"]["final"])
        row["final_in_other_family"] = "; ".join("%s=%d" % kv for kv in sorted(other.items()))
        row["not_in_our_clades"] = listed(sum(not any(m["stage_fams"][s] for s, _ in stages) for m in traced))
        if not ours:
            row["our_status"], row["our_count"], row["comparison"] = "", "", "genome not in our set"
        else:
            if mapped == ["*"]:
                our_count = sum(same(m["stage_fams"]["final"]) for m in traced)
            else:
                our_count = sum(counts.get(t, {}).get(f, 0) for t in ours for f in mapped)
            row["our_status"], row["our_count"] = ("present" if our_count else "absent"), our_count
            if family == NOT_A_GPCR:
                row["comparison"] = "they: not a GPCR; we assign a family" if our_count else "agree: not in our families"
            elif status == "present":
                if our_count:
                    row["comparison"] = "agree: present"
                elif not members:
                    row["comparison"] = "they report it, we do not: no receptor to trace"
                elif not traced:
                    row["comparison"] = "they report it, we do not: their receptors not found in our data"
                elif other:
                    row["comparison"] = "they report it, we do not: their receptors are in another family of ours"
                elif any(m["stage_fams"]["clade"] and not m["stage_fams"]["final"] for m in traced):
                    row["comparison"] = "they report it, we do not: their receptors removed by our filters"
                else:
                    row["comparison"] = "they report it, we do not: their receptors in no clade of ours"
            else:
                row["comparison"] = "agree: absent" if not our_count else "they report it absent, we find it"
        rows.append(row)

    _write("literature_comparison.tsv", rows)
    _write("literature_sequences.tsv", seq_rows)

    # Genome names still not matched to one of ours, with the closest of our species by name.
    # For review only: the candidates are not used as matches.
    our_names = {_plain(matcher.sci[s]): matcher.sci[s] for s in matcher.ours if s in matcher.sci}
    unresolved = defaultdict(lambda: {"studies": set(), "rows": 0})
    for row in rows:
        if row["comparison"] == "genome not in our set":
            u = unresolved[row["genome_as_reported"]]
            u["studies"].add(row["study"])
            u["rows"] += 1
    review = []
    for label, u in sorted(unresolved.items()):
        species = matcher.match(label)[0]
        words = re.sub(r"\s*\([^)]*\)", "", label).replace("_", " ").split()
        binomial = _plain(" ".join(words[:2]))
        close = get_close_matches(binomial, list(our_names), n=3, cutoff=0)  # the three closest, however far
        genus = [our_names[n] for n in our_names if words and n.split()[0] == _plain(words[0])]
        review.append({
            "genome_as_reported": label, "studies": "; ".join(sorted(u["studies"])), "rows": u["rows"],
            "ncbi": ("species taxID %s, not in our set" % species) if species else "not resolved in the NCBI taxonomy",
            "closest_our_species": "; ".join("%s (%.2f)" % (our_names[n], SequenceMatcher(None, binomial, n).ratio())
                                             for n in close),
            "our_species_same_genus": "; ".join(sorted(genus)),
        })
    _write("literature_unresolved_genomes.tsv", review)

    # 10. study families with no counterpart (and those mapped by where their receptors end)
    columns = ["genome not in our set", "no species given", "not found in our data", "in none of our clades",
               "in a clade, removed by filters"]
    dest_rows = []
    for (study, family), d in sorted(destinations.items()):
        finals = {k[6:]: n for k, n in d.items() if k.startswith("final ")}
        top = max(finals, key=finals.get) if finals else ""
        dest_rows.append({
            "study": study, "their_family": family,
            "status": "mapped to %s by where its receptors end" % ",".join(FAMILY_MAP[family]) if family in INFERRED_MAP
                      else "not compared (no counterpart)",
            "receptors_listed": d["listed"], **{c.replace(" ", "_").replace(",", ""): d[c] for c in columns},
            "final_in_our_families": "; ".join("%s=%d" % kv for kv in sorted(finals.items(), key=lambda kv: -kv[1])),
            "main_family": top, "main_family_share_of_final": "%.2f" % (finals[top] / sum(finals.values())) if top else ""})
    _write("literature_unmapped_families.tsv", dest_rows)

    summary = Counter((row["study"], row["their_family"], row["our_families"], row["comparison"]) for row in rows)
    with open(os.path.join(OUT_DIR, "literature_summary.tsv"), "w", newline="", encoding="utf-8") as handle:
        w = csv.writer(handle, delimiter="\t")
        w.writerow(["study", "their_family", "our_families", "comparison", "genomes"])
        for key, n in sorted(summary.items()):
            w.writerow(list(key) + [n])
        w.writerow([])
        w.writerow(["genome claims not compared (no counterpart among our families; their receptors, if listed, "
                    "are traced in literature_unmapped_families.tsv)", "", "", "", ""])
        for (study, family), n in sorted(unmapped_claims.items()):
            w.writerow([study, family, "", "claims", n])
    print(Counter(re.sub(r" \(.*", "", r["stopped_at"]) for r in seq_rows).most_common())


if __name__ == "__main__":
    main()
