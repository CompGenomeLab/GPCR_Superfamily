"""Shared locations of the repository and the Zenodo data archive.

Every script finds its inputs and writes its outputs through this module, so no path is
hard-coded. The archive is expected beside the repository,

    project/
    |-- GPCR_Superfamily/   (this repository)
    `-- zenodo_data/        (the unpacked Zenodo archive)

and is found elsewhere through the GPCR_ZENODO_DIR environment variable. The archive has one
folder per stage, named like the stage folders of `code_and_data/`, so `stage("04_overlap_assignment")`
is where stage 04 reads and writes its data.

Scripts import it with

    import sys; from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from gpcr_paths import CODE, DATA, stage
"""
import gzip
import os
from pathlib import Path

CODE = Path(__file__).resolve().parent
REPOSITORY = CODE.parent
DATA = Path(os.environ.get("GPCR_ZENODO_DIR", REPOSITORY.parent / "zenodo_data")).resolve()

# The 19 GPCR families, in the order of the manuscript figures, under the internal names
# used in file names; DISPLAY gives the names used in the text and figures.
FAMILIES = ["classA", "Olfactory", "Vomeronasal1", "classT", "Nematode", "GPR1", "cAMP",
            "classB1", "classB2", "Mth", "classF", "FSLB", "GP143", "TM116", "GP157",
            "STE3", "classC", "Vomeronasal2", "STE2"]
CONTROLS = ["Galpha", "Gbeta", "Ggamma", "TBP", "RBP1"]
DISPLAY = {"classA": "Class A", "classB1": "Class B1", "classB2": "Class B2", "classC": "Class C",
           "classF": "Class F", "classT": "Taste-2", "FSLB": "FSL", "Vomeronasal1": "Vomeronasal 1",
           "Vomeronasal2": "Vomeronasal 2", "Galpha": "G alpha", "Gbeta": "G beta",
           "Ggamma": "G gamma", "RBP1": "RPB1"}


def open_text(path, mode="rt"):
    """Open a text file for reading, or its gzip-compressed copy (path + '.gz') when the
    archive holds only that, so the archive can be used without decompressing it."""
    path = Path(path)
    if not path.exists() and Path(str(path) + ".gz").exists():
        path = Path(str(path) + ".gz")
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", errors="replace")
    return open(path, mode, encoding="utf-8", errors="replace")


def exists(path):
    """True when the file or its gzip-compressed copy is present."""
    return Path(path).exists() or Path(str(path) + ".gz").exists()


def stage(name, *parts):
    """A path inside one stage folder of the data archive."""
    return DATA.joinpath(name, *parts)


def code(name, *parts):
    """A path inside one stage folder of this repository."""
    return CODE.joinpath(name, *parts)
