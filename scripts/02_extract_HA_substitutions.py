#!/usr/bin/env python3
"""Count substitutions in an explicitly curated Tunisian subclade K set.

Python 3 standard library only. Inputs: Nextclade TSV and a text file containing
the retained K isolate accessions (one per line or comma-separated). This script
does not infer which passage record to retain: curate that list beforehand.
Region labels are preserved; no conversion to HA0 coordinates is performed.
"""
import argparse
import collections
import csv
import hashlib
import json
import re
from pathlib import Path

AA = re.compile(r'^(SigPep|HA1|HA2):([A-Z])(\d+)([A-Z])$')
REGIONS = {'SigPep': 0, 'HA1': 1, 'HA2': 2}

def parse_mutation(value):
    match = AA.fullmatch(value)
    if not match:
        raise ValueError(f'Unexpected amino-acid substitution: {value}')
    return match[1], int(match[3])

def run(nextclade, accessions, outdir):
    raw_ids = re.split(r'[\s,;]+', accessions.read_text().strip())
    if not raw_ids or any(not re.fullmatch(r'EPI_ISL_\d+', x) for x in raw_ids):
        raise ValueError('Accession file must contain EPI_ISL identifiers only.')
    if len(raw_ids) != len(set(raw_ids)):
        raise ValueError('Duplicate accession in curated list.')
    wanted = set(raw_ids)
    with nextclade.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'seqName', 'clade', 'aaSubstitutions', 'unknownAaRanges',
                    'aaDeletions', 'aaInsertions', 'qc.overallStatus', 'errors'}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f'Missing Nextclade columns: {sorted(missing)}')
        selected = []
        seen = set()
        for row in reader:
            matches = set(re.findall(r'EPI_ISL_\d+', row['seqName'])) & wanted
            if not matches:
                continue
            if len(matches) != 1:
                raise ValueError('Multiple curated identifiers in one sequence header.')
            accession = matches.pop()
            if accession in seen:
                raise ValueError(f'Duplicate Nextclade record: {accession}')
            seen.add(accession)
            if 'A/Tunisia/' not in row['seqName'] or row['clade'] != 'K':
                raise ValueError(f'Curated record is not Tunisian subclade K: {accession}')
            if row['qc.overallStatus'] != 'good' or row['errors']:
                raise ValueError(f'Nextclade quality requires review: {accession}')
            if any(row[field] for field in ['unknownAaRanges', 'aaDeletions', 'aaInsertions']):
                raise ValueError(f'Unknown residues or amino-acid indels require review: {accession}')
            mutations = set(filter(None, row['aaSubstitutions'].split(',')))
            for mutation in mutations:
                parse_mutation(mutation)
            selected.append((accession, row['seqName'], mutations))
    if wanted != seen:
        raise ValueError(f'Curated accessions absent from TSV: {sorted(wanted-seen)}')
    n = len(selected)
    if n < 2:
        raise ValueError('At least two sequences are required.')
    counts = collections.Counter(m for _, _, mutations in selected for m in mutations)
    mutations = sorted(counts, key=lambda m: (REGIONS[parse_mutation(m)[0]], parse_mutation(m)[1], m))
    outdir.mkdir(parents=True, exist_ok=True)
    def write_csv(name, header, records):
        with (outdir/name).open('w', newline='') as handle:
            writer = csv.writer(handle)
            writer.writerow(header)
            writer.writerows(records)
    header = ['Mutation_full', 'Mutation', 'Region', 'Position', 'Count',
              'Denominator', 'Frequency_percent', 'Main_substitution']
    records = [[m, m.split(':')[1], *parse_mutation(m), counts[m], n,
                round(100*counts[m]/n, 6), 'Yes' if counts[m] >= n-1 else 'No'] for m in mutations]
    write_csv('HA_Tunisia_mutation_counts.csv', header, records)
    write_csv('HA_Tunisia_main_mutations.csv', header, [r for r in records if r[-1] == 'Yes'])
    write_csv('HA_Tunisia_mutation_heatmap_matrix.csv', ['Accession', 'Sequence']+mutations,
              [[a, name]+[int(m in muts) for m in mutations] for a, name, muts in sorted(selected)])
    summary = {'sequences': n, 'distinct_substitutions': len(counts), 'main_threshold_count': n-1,
               'main_substitutions': sum(c >= n-1 for c in counts.values()),
               'reference': 'The reference and dataset version must be documented from the original Nextclade dataset.',
               'input_sha256': hashlib.sha256(nextclade.read_bytes()).hexdigest(),
               'curated_accessions_sha256': hashlib.sha256(accessions.read_bytes()).hexdigest()}
    (outdir/'extraction_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nextclade', type=Path, required=True)
    parser.add_argument('--accessions', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, default=Path('results'))
    args = parser.parse_args()
    run(args.nextclade, args.accessions, args.outdir)
