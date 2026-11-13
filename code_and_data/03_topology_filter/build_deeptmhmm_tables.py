"""Combine DeepTMHMM predictions into one table per family, then keep the GPCR-like cores.

Each family directory holds a DeepTMHMM run, either directly or split across shard
subdirectories. Every TMRs.gff3 and predicted_topologies.3line found beneath a family
directory is treated as part of that family's run and merged on accession.

Selection is by topology rather than by helix count. The membrane state is split by the
direction the helix crosses, M for inside to outside and W for outside to inside, the
same alphabet project_topology_onto_alignment.py uses. A GPCR read in that alphabet is
OWIMOWIMOWIMOWI once repeats are collapsed, which is seven helices with the N terminus
outside and the C terminus inside. A sequence is kept when that pattern occurs exactly
once, so the seven helices forming the core are uniquely located.

Keeping on the pattern rather than on a count of seven admits the sequences whose signal
peptide was predicted as an eighth helix, which are a real and homogeneous class, and it
refuses the sequences with nine or more helices where the pattern matches in two phases
and the core cannot be pinned down. Sequences missing a terminal helix are refused as
before, but the table now records which terminus is missing, because most of them are
truncated records rather than failed predictions.

The accessions are also written as two wide family-by-accession tables, one holding every
sequence with a DeepTMHMM prediction and one holding only the kept subset.

A combined 3line file is written across all families, carrying the oriented M/W labels,
so the whole set can be projected onto an alignment from one file. That file is what the
archive holds (03_topology_filter/deeptmhmm_predictions_oriented.3line); the per-family
selection used downstream is made from it by write_core_membership.py.

Usage:
    python build_deeptmhmm_tables.py --root <DeepTMHMM output, one folder per family>
"""

import argparse
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import stage  # noqa: E402

parser = argparse.ArgumentParser(description="Build per-family DeepTMHMM tables and GPCR-like-core subsets")
parser.add_argument('--root', required=True,
                    help='DeepTMHMM output directory holding one subdirectory per family')
parser.add_argument('--output-dir', default=str(stage('03_topology_filter', 'deeptmhmm_tables')))
parser.add_argument('--families', nargs='*', default=None,
                    help='restrict to these family directory names (default: all)')
parser.add_argument('--merged-table', default=None,
                    help='path of the wide family-by-accession table of kept sequences '
                         '(default: <output-dir>/accessions_by_family_single_core.tsv)')
parser.add_argument('--prefilter-table', default=None,
                    help='path of the wide family-by-accession table of all sequences before '
                         'selection '
                         '(default: <output-dir>/accessions_by_family_with_prediction.tsv)')
parser.add_argument('--combined-3line', default=None,
                    help='path of the combined 3line file across all families, with oriented '
                         'M/W labels (default: 03_topology_filter/deeptmhmm_predictions_oriented.3line)')
parser.add_argument('--combined-3line-kept-only', action='store_true',
                    help='write only the kept sequences to the combined 3line instead of all')
args = parser.parse_args()

COLUMNS = ['accession', 'prediction', 'TM_count', 'core_count', 'topology_class',
           'core_start', 'core_stop', 'sequence', 'label_string', 'oriented_label_string',
           'collapsed']

count_line = re.compile(r'^# (\S+) Number of predicted TMRs: (\d+)\s*$')

SIDES = ('I', 'O')
# A GPCR once repeats are collapsed: N terminus outside, seven helices alternating
# W (outside to inside) and M (inside to outside), C terminus inside.
GPCR_PATTERN = 'OWIMOWIMOWIMOWI'
# Where TM1 and TM7 sit inside that pattern, used to map a match back to residues.
CORE_FIRST_HELIX = 1
CORE_LAST_HELIX = 13


def accession_of(protein_id):
    """tr|A0A8H5TWS9|A0A8H5TWS9_FUSHE -> A0A8H5TWS9; anything else passes through."""
    fields = protein_id.split('|')
    return fields[1] if len(fields) >= 3 else protein_id


def orient_helices(topology):
    """Split the membrane state M into two, by the direction the helix crosses.

    Verbatim in behaviour with project_topology_onto_alignment.py: 'IMO' becomes M and
    'OMI' becomes W. A run at either terminus is decided from whichever side it has, and
    a re-entrant helix flanked by the same state is left as M rather than guessed at.
    """
    out = list(topology)
    for match in re.finditer(r'M+', topology):
        start, end = match.start(), match.end()
        before = topology[start - 1] if start > 0 else ''
        after = topology[end] if end < len(topology) else ''
        if before in SIDES and after in SIDES:
            if before == 'O' and after == 'I':
                direction = 'W'
            elif before == 'I' and after == 'O':
                direction = 'M'
            else:
                continue
        elif before in SIDES:
            direction = 'W' if before == 'O' else 'M'
        elif after in SIDES:
            direction = 'W' if after == 'I' else 'M'
        else:
            continue
        out[start:end] = direction * (end - start)
    return ''.join(out)


def topology_runs(labels):
    """Return [(label, start, end)] for every maximal run, end exclusive, 0-based."""
    return [(match.group()[0], match.start(), match.end())
            for match in re.finditer(r'(.)\1*', labels)]


def find_cores(oriented):
    """Return (collapsed string, [(tm1_start, tm7_end)]) for every GPCR-like core found.

    Occurrences are counted overlapping, because a sequence with more helices than one
    core can match the pattern in two phases and neither phase is the right answer.
    """
    runs = topology_runs(oriented)
    collapsed = ''.join(label for label, _, _ in runs)
    cores = []
    index = collapsed.find(GPCR_PATTERN)
    while index != -1:
        cores.append((runs[index + CORE_FIRST_HELIX][1], runs[index + CORE_LAST_HELIX][2]))
        index = collapsed.find(GPCR_PATTERN, index + 1)
    return collapsed, cores


def classify(collapsed, n_cores):
    """Name the topology, so the rejected sequences can be told apart in the stats."""
    if n_cores == 1:
        return 'GPCR_like_core'
    if n_cores > 1:
        return 'multi_core'
    if collapsed.endswith('IMOWIMOWIMOWI') and collapsed.lstrip('S') == 'IMOWIMOWIMOWI':
        return 'missing_TM1'
    if collapsed.lstrip('S') == 'OWIMOWIMOWIMO':
        return 'missing_TM7'
    if 'IMOWIMOWIMOWIMO' in collapsed:
        return 'inverted'
    return 'other'


def parse_gff3(path, counts, duplicates):
    """Read 'Number of predicted TMRs' comment lines into counts[accession]."""
    with open(path) as infile:
        for line in infile:
            match = count_line.match(line)
            if not match:
                continue
            accession = accession_of(match.group(1))
            n_tm = int(match.group(2))
            if accession in counts:
                duplicates.append((accession, 'gff3', str(path)))
                continue
            counts[accession] = n_tm


def parse_3line(path, records, duplicates):
    """Read header/sequence/label triples into records[accession]."""
    with open(path) as infile:
        lines = [line.rstrip('\n') for line in infile if line.strip()]
    if len(lines) % 3 != 0:
        sys.exit(f'{path}: expected a multiple of three non-blank lines, found {len(lines)}')
    for index in range(0, len(lines), 3):
        header, sequence, labels = lines[index:index + 3]
        if not header.startswith('>'):
            sys.exit(f'{path}: expected a header at line {index + 1}, found {header[:40]!r}')
        if len(sequence) != len(labels):
            sys.exit(f'{path}: sequence and label lengths differ for {header[:60]!r}')
        parts = header[1:].rsplit('|', 1)
        protein_id = parts[0].strip() if len(parts) == 2 else header[1:].strip()
        prediction = parts[1].strip() if len(parts) == 2 else ''
        accession = accession_of(protein_id)
        if accession in records:
            duplicates.append((accession, '3line', str(path)))
            continue
        records[accession] = (protein_id, prediction, sequence, labels)


def write_merged_table(path, family_to_accessions):
    """Write one column per family, accessions as rows, short columns padded blank.

    The wide layout every stage reads: each column name is a family.
    """
    families = {family: accessions for family, accessions in family_to_accessions.items()
                if accessions}
    if not families:
        return None
    height = max(len(accessions) for accessions in families.values())
    names = sorted(families, key=str.casefold)
    with open(path, 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t', lineterminator='\n')
        writer.writerow(names)
        for index in range(height):
            writer.writerow([families[name][index] if index < len(families[name]) else ''
                             for name in names])
    return path


def write_tsv(path, rows):
    with open(path, 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t', lineterminator='\n')
        writer.writerow(COLUMNS)
        writer.writerows(rows)


root = Path(args.root)
if not root.is_dir():
    sys.exit(f'{root} is not a directory')

output_dir = Path(args.output_dir)
output_dir.mkdir(parents=True, exist_ok=True)

families = sorted(path.name for path in root.iterdir() if path.is_dir())
if args.families:
    missing = [name for name in args.families if name not in families]
    if missing:
        sys.exit(f'no such family directory under {root}: {", ".join(missing)}')
    families = [name for name in families if name in args.families]

combined_path = Path(args.combined_3line) if args.combined_3line else (
    stage('03_topology_filter', 'deeptmhmm_predictions_oriented.3line'))
combined_path.parent.mkdir(parents=True, exist_ok=True)
combined = open(combined_path, 'w', newline='')
n_combined = 0

stats_rows = []
kept_by_family = {}
all_by_family = {}
for family in families:
    family_dir = root / family
    gff3_paths = sorted(family_dir.rglob('TMRs.gff3'))
    line3_paths = sorted(family_dir.rglob('predicted_topologies.3line'))
    if not gff3_paths or not line3_paths:
        print(f'{family}: no DeepTMHMM output found, skipping', file=sys.stderr)
        continue

    counts, records, duplicates = {}, {}, []
    for path in gff3_paths:
        parse_gff3(path, counts, duplicates)
    for path in line3_paths:
        parse_3line(path, records, duplicates)

    shared = sorted(set(counts) & set(records))
    gff3_only = sorted(set(counts) - set(records))
    line3_only = sorted(set(records) - set(counts))

    rows = []
    classes = {}
    for accession in shared:
        protein_id, prediction, sequence, labels = records[accession]
        oriented = orient_helices(labels)
        collapsed, cores = find_cores(oriented)
        topology_class = classify(collapsed, len(cores))
        classes[topology_class] = classes.get(topology_class, 0) + 1
        # 1-based inclusive, the span of the seven core helices themselves.
        core_start = cores[0][0] + 1 if len(cores) == 1 else ''
        core_stop = cores[0][1] if len(cores) == 1 else ''
        rows.append([accession, prediction, counts[accession], len(cores), topology_class,
                     core_start, core_stop, sequence, labels, oriented, collapsed])
        if not args.combined_3line_kept_only or len(cores) == 1:
            combined.write(f'>{protein_id} | {prediction}\n{sequence}\n{oriented}\n')
            n_combined += 1

    kept = [row for row in rows if row[3] == 1]

    write_tsv(output_dir / f'{family}.tsv', rows)
    write_tsv(output_dir / f'{family}_singlecore.tsv', kept)
    with open(output_dir / f'{family}_singlecore_accessions.txt', 'w', newline='') as outfile:
        for row in kept:
            outfile.write(f'{row[0]}\n')

    kept_by_family[family] = [row[0] for row in kept]
    all_by_family[family] = [row[0] for row in rows]

    # The helix count is never used to select, only reported, so that the composition
    # of the kept set stays visible without a count being what built it.
    kept_counts = {}
    for row in kept:
        kept_counts[row[2]] = kept_counts.get(row[2], 0) + 1

    stats_rows.append({
        'family': family,
        'n_shards': len(gff3_paths),
        'n_total': len(rows),
        'n_kept': len(kept),
        'n_removed': len(rows) - len(kept),
        'percent_kept': f'{100 * len(kept) / len(rows):.2f}' if rows else 'NA',
        'n_kept_7TM': kept_counts.get(7, 0),
        'n_kept_8TM': kept_counts.get(8, 0),
        'n_kept_9plus': sum(n for tm, n in kept_counts.items() if tm >= 9),
        'n_missing_TM1': classes.get('missing_TM1', 0),
        'n_missing_TM7': classes.get('missing_TM7', 0),
        'n_multi_core': classes.get('multi_core', 0),
        'n_inverted': classes.get('inverted', 0),
        'n_other': classes.get('other', 0),
        'n_gff3_only': len(gff3_only),
        'n_3line_only': len(line3_only),
        'n_duplicate_accessions': len(duplicates),
    })

    print(f'{family}: {len(rows)} sequences from {len(gff3_paths)} file set(s), '
          f'{len(kept)} kept with one GPCR-like core, {len(rows) - len(kept)} removed')
    for accession, source, path in duplicates:
        print(f'  duplicate accession {accession} in {source} {path}, first occurrence kept',
              file=sys.stderr)
    if gff3_only:
        print(f'  {len(gff3_only)} accession(s) in gff3 but not 3line, dropped: '
              f'{", ".join(gff3_only[:5])}', file=sys.stderr)
    if line3_only:
        print(f'  {len(line3_only)} accession(s) in 3line but not gff3, dropped: '
              f'{", ".join(line3_only[:5])}', file=sys.stderr)

combined.close()

if not stats_rows:
    sys.exit('no families produced a table')

stats_fields = list(stats_rows[0])
totals = {field: sum(row[field] for row in stats_rows) for field in
          ['n_shards', 'n_total', 'n_kept', 'n_removed', 'n_kept_7TM', 'n_kept_8TM',
           'n_kept_9plus', 'n_missing_TM1', 'n_missing_TM7', 'n_multi_core', 'n_inverted',
           'n_other', 'n_gff3_only', 'n_3line_only', 'n_duplicate_accessions']}
totals['family'] = 'ALL'
totals['percent_kept'] = (f'{100 * totals["n_kept"] / totals["n_total"]:.2f}'
                          if totals['n_total'] else 'NA')

with open(output_dir / 'filtering_stats.tsv', 'w', newline='') as outfile:
    writer = csv.DictWriter(outfile, fieldnames=stats_fields, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(stats_rows)
    writer.writerow(totals)

prefilter_path = Path(args.prefilter_table) if args.prefilter_table else (
    output_dir / 'accessions_by_family_with_prediction.tsv')
prefilter_path.parent.mkdir(parents=True, exist_ok=True)
write_merged_table(prefilter_path, all_by_family)

merged_path = Path(args.merged_table) if args.merged_table else (
    output_dir / 'accessions_by_family_single_core.tsv')
merged_path.parent.mkdir(parents=True, exist_ok=True)
write_merged_table(merged_path, kept_by_family)

print(f'\nTotal: {totals["n_total"]} sequences across {len(stats_rows)} families, '
      f'{totals["n_kept"]} kept with one GPCR-like core, {totals["n_removed"]} removed '
      f'({totals["percent_kept"]}% kept)')
print(f'  kept:    {totals["n_kept_7TM"]} with 7 helices, {totals["n_kept_8TM"]} with 8, '
      f'{totals["n_kept_9plus"]} with 9 or more')
print(f'  removed: {totals["n_missing_TM1"]} missing TM1, {totals["n_missing_TM7"]} missing TM7, '
      f'{totals["n_multi_core"]} multi-core, {totals["n_inverted"]} inverted, '
      f'{totals["n_other"]} other')
print(f'Wrote tables and filtering_stats.tsv to {output_dir}')
print(f'Wrote prefilter family-by-accession table: {prefilter_path}')
print(f'Wrote single-core family-by-accession table: {merged_path}')
print(f'Wrote combined oriented 3line ({n_combined} sequences): {combined_path}')
