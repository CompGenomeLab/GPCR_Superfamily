"""Writes literature/parsed/<study>.tsv, what each reader in literature_readers.py takes from a study's
supplementary files, so the comparison runs without the files themselves (they are not distributed).

Needs the supplementary files (literature/SOURCES.md), in a folder laid out like literature/, given by
GPCR_LITERATURE_ORIGINALS. read_all() then reads the parsed tables; read_all(from_originals=True) the files.

Usage:
    GPCR_LITERATURE_ORIGINALS=<folder> python write_parsed_literature.py
"""
import literature_readers as L

for reader in L.READERS:
    if reader is L.read_curated:
        continue
    records = reader()
    L.write_parsed(reader, records)
    print("%-28s %6d records" % (reader.__name__[len("read_"):], len(records)))
