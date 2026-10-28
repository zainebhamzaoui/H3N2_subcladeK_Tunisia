#!/usr/bin/env python3
"""Annotate HA1 site locations using Table 3 of Lees et al. (2010).

DOI: 10.1093/bioinformatics/btq160. Previously identified residues and proposed
additional assignments are reported separately. These are site locations, not
experimental evidence of altered antigenicity. Requires only Python 3.
Input: HA_Tunisia_mutation_counts.csv from corrected script 02.
"""
import argparse
import csv
import math
import re
from pathlib import Path

def residues(text):
    result = set()
    for token in text.split(','):
        bounds = [int(n) for n in token.strip().split('-')]
        result.update(range(bounds[0], bounds[-1]+1))
    return result

PREVIOUS = {k: residues(v) for k, v in {
    'A': '122,124,126,130-133,135,137,138,140,142-146,150,152,168',
    'B': '128,129,155-160,163-165,186-190,192-194,196-198',
    'C': '44-48,50,51,53,54,273,275,276,278-280,294,297,299,300,304,305,307-312',
    'D': '96,102,103,117,121,167,170-177,179,182,201,203,207-209,212-219,226-230,238,240,242,244,246-248',
    'E': '57,59,62,63,67,75,78,80-83,86-88,91,92,94,109,260-262,265',
}.items()}
ADDITIONAL = {k: residues(v) for k, v in {
    'A': '71-72,98,127,141,148-149,151,255',
    'B': '161-162,199',
    'C': '41-43,49,52,55,271-272,274,282,284-285,287-293,295-296,298,301-303,313-314',
    'D': '95,97,99-101,104-105,107,118,120,166,169,178,180,183-184,200,204-206,210-211,220-225,231-236,239,243,245,257-258',
    'E': '56,58,60,64-65,68-70,73-74,76-77,79,84-85,89-90,93,110-115,119,259,263-264,267-270',
}.items()}

def annotate(region, position):
    if region != 'HA1':
        return '', 'Not applicable to HA1 A-E scheme'
    for scheme, category in [(PREVIOUS, 'Previously identified'), (ADDITIONAL, 'Proposed additional assignment')]:
        hits = [site for site, positions in scheme.items() if position in positions]
        if hits:
            if len(hits) != 1:
                raise ValueError('Overlapping site definitions.')
            return hits[0], category
    return '', 'Unassigned in selected scheme'

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=Path('results/HA_Tunisia_mutation_counts.csv'))
    parser.add_argument('--output', type=Path, default=Path('results/HA_Tunisia_antigenic_annotation.csv'))
    args = parser.parse_args()
    with args.input.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames
        required = {'Mutation_full', 'Region', 'Position', 'Count', 'Denominator', 'Frequency_percent'}
        if not required.issubset(fields or []):
            raise ValueError('Input must be the count table from corrected script 02.')
        rows = list(reader)
    seen = set()
    denominators = set()
    for row in rows:
        label = row['Mutation_full']
        match = re.fullmatch(r'(SigPep|HA1|HA2):[A-Z](\d+)[A-Z]', label)
        if not match or label in seen:
            raise ValueError(f'Invalid or duplicate mutation: {label}')
        seen.add(label)
        position, n, count = int(row['Position']), int(row['Denominator']), int(row['Count'])
        if row['Region'] != match[1] or position != int(match[2]):
            raise ValueError(f'Inconsistent region/position for {label}')
        if n < 2 or not 1 <= count <= n or not math.isclose(float(row['Frequency_percent']),100*count/n,abs_tol=1e-5):
            raise ValueError(f'Inconsistent frequency for {label}')
        denominators.add(n)
        row['Main_substitution'] = 'Yes' if count >= n-1 else 'No'
        row['Antigenic_site'], row['Site_assignment_category'] = annotate(row['Region'], position)
        row['Annotation_source_DOI'] = '10.1093/bioinformatics/btq160'
    if not rows or len(denominators) != 1:
        raise ValueError('A nonempty table with a common denominator is required.')
    extra = ['Main_substitution', 'Antigenic_site', 'Site_assignment_category', 'Annotation_source_DOI']
    fields = list(fields)+[f for f in extra if f not in fields]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f'Annotated {len(rows)} substitutions; {sum(r["Main_substitution"]=="Yes" for r in rows)} main substitutions.')
    print(args.output)

if __name__ == '__main__':
    main()
