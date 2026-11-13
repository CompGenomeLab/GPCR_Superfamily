"""Families scored and key nodes of the NCBI taxonomy used to name family origins (stage 06).

KEY_NODES lists (label, node name in the NCBI taxonomy tree) from the root outwards; an origin
falling between two key nodes is named by the most recent key node containing it.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import CONTROLS, FAMILIES as GPCR_FAMILIES, stage  # noqa: E402

FAMILIES = GPCR_FAMILIES + CONTROLS
OUT_DIR = str(stage("06_family_origins"))
TAXONOMY_DIR = str(stage("06_family_origins", "taxonomy"))

# (label, NCBI node name), from the root outwards.
KEY_NODES = [
    ("LECA", "root"),
    ("Unikonta", "clade__Unikonta"),
    ("Opisthokonta", "clade 1__Opisthokonta"),
    ("Holozoa", "clade 2__Holozoa"),
    ("Choanoflagellata", "family__Salpingoecidae"),
    ("Metazoa", "kingdom__Metazoa"),
    ("Porifera", "subclass__Heteroscleromorpha"),
    ("Placozoa", "genus__Trichoplax"),
    ("Eumetazoa", "clade 3__Eumetazoa"),
    ("Cnidaria", "phylum__Cnidaria"),
    ("Bilateria", "clade 4__Bilateria"),
    ("Protostomia", "clade 5__Protostomia"),
    ("Ecdysozoa", "clade 6__Ecdysozoa"),
    ("Nematoda", "phylum__Nematoda"),
    ("Lophotrochozoa", "clade 7__Lophotrochozoa"),
    ("Mollusca", "phylum__Mollusca"),
    ("Deuterostomia", "clade 5__Deuterostomia"),
    ("Echinodermata", "phylum__Echinodermata"),
    ("Chordata", "phylum__Chordata"),
    ("Craniata", "subphylum__Craniata"),
    ("Gnathostomata", "clade 7__Gnathostomata"),
    ("Euteleostomi", "clade 9__Euteleostomi"),
    ("Nucletmycea", "clade 2__Holomycota"),
    ("Fungi", "kingdom__Fungi"),
    ("Dikarya", "subkingdom__Dikarya"),
    ("Ascomycota", "phylum__Ascomycota"),
    ("Basidiomycota", "phylum__Basidiomycota"),
    ("Glomeromycota", "class__Glomeromycetes"),
    ("Chytridiomycota", "no rank 1__Chytridiomycota%20incertae%20sedis"),
    ("Microsporidia", "phylum__Microsporidia"),
    ("Amoebozoa", "clade 1__Amoebozoa"),
    ("Dictyostelia", "class__Eumycetozoa"),
    ("Archaeplastida", "clade__Archaeplastida"),
    ("Viridiplantae", "kingdom__Viridiplantae"),
    ("Embryophyta", "clade 1__Embryophyta"),
    ("Chlorophyta", "phylum__Chlorophyta"),
    ("Rhodophyta", "phylum__Rhodophyta"),
    ("SAR", "clade__Sar"),
    ("Stramenopiles", "clade 1__Stramenopiles"),
    ("Alveolata", "clade 1__Alveolata"),
    ("Ciliophora", "phylum__Ciliophora"),
    ("Rhizaria", "clade 1__Rhizaria"),
    ("Discoba", "clade__Discoba"),
    ("Metamonada", "clade__Metamonada"),
    ("Haptophyta", "phylum__Haptophyta"),
]