# 08 · Superfamily alignments

Methods: *Superfamily Alignment*.

`build_representative_fastas.py` trims every representative to its 7TM core, from the first to the
seventh helix of its single `OWIMOWIMOWIMOWI` pattern with five flanking residues on either side
(`extract_tm_region.py`), and writes the four sequence sets: the three sets of stage 07 and the core set
with the 18 literature sequences of stage 05 (493). Each set is aligned with MAFFT G-INS-i
(`--globalpair --maxiterate 1000`), once with `--ep 0.123` (main) and once with the default gap
extension penalty, and trimmed with ClipKIT `kpic-gappy` (`superfamily_mafft.sh`). The topology of every aligned
residue is projected with `../03_topology_filter/project_topology_onto_alignment.py`.

```bash
python build_representative_fastas.py
sbatch superfamily_mafft.sh        # MAFFT G-INS-i (--ep 0.123 and default) and ClipKIT kpic-gappy
```

Verified: the four trimmed input sets written from the archive are identical to those the alignments
were built from.

**Archive** (`08_superfamily_alignment/`): `fasta/<set>/` (full-length sequences, 7TM cores, core
boundaries; `representatives_and_literature.3line` with every topology) and `alignments/`
(`<set>_ginsi_ep0.123.fasta`, `<set>_ginsi.fasta`, their `_kpicgappy` trimmed versions and
`_topology_converted` projections). The alignments are also on gpcrevolution.org.
