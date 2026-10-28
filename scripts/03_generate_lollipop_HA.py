#!/usr/bin/env python3
"""Plot HA substitutions from the corrected script 02 count table.

Requires matplotlib. By default, shows all observed substitutions with a dashed
threshold at (N-1)/N. Use --main-only to show only main substitutions. Region
labels and positions are preserved; no antigenic effects are inferred.
"""
import argparse
import csv
import math
import re
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ORDER = {'SigPep': 0, 'HA1': 1, 'HA2': 2}
COLORS = {'SigPep': '#777777', 'HA1': '#23689B', 'HA2': '#B35D24'}

def load_counts(path):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle)
        required = {'Mutation_full', 'Region', 'Position', 'Count', 'Denominator', 'Frequency_percent'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError('Use HA_Tunisia_mutation_counts.csv from corrected script 02.')
        rows = list(reader)
    if not rows:
        raise ValueError('The input table is empty.')
    seen, denominators = set(), set()
    for row in rows:
        label = row['Mutation_full']
        match = re.fullmatch(r'(SigPep|HA1|HA2):([A-Z])(\d+)([A-Z])', label)
        if not match or label in seen:
            raise ValueError(f'Invalid or duplicate mutation: {label}')
        seen.add(label)
        row['Count'] = int(row['Count'])
        row['Denominator'] = int(row['Denominator'])
        row['Position'] = int(row['Position'])
        n, count = row['Denominator'], row['Count']
        if n < 2 or not 1 <= count <= n:
            raise ValueError(f'Invalid count or denominator: {label}')
        if row['Region'] != match[1] or row['Position'] != int(match[3]):
            raise ValueError(f'Region or position disagrees with mutation label: {label}')
        frequency = 100 * count / n
        if not math.isclose(float(row['Frequency_percent']), frequency, abs_tol=1e-5):
            raise ValueError(f'Frequency disagrees with count: {label}')
        row['frequency'] = frequency
        denominators.add(n)
    if len(denominators) != 1:
        raise ValueError('This plot requires a common denominator for all positions.')
    return sorted(rows, key=lambda r: (ORDER[r['Region']], r['Position'], r['Mutation_full'])), denominators.pop()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=Path('results/HA_Tunisia_mutation_counts.csv'))
    parser.add_argument('--out-prefix', type=Path, default=Path('results/Figure_HA_lollipop'))
    parser.add_argument('--main-only', action='store_true')
    args = parser.parse_args()
    rows, n = load_counts(args.input)
    main_count = sum(r['Count'] >= n-1 for r in rows)
    if args.main_only:
        rows = [r for r in rows if r['Count'] >= n-1]
    if not rows:
        raise ValueError('No substitutions meet the selected threshold.')
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'svg.fonttype': 'none'})
    fig, ax = plt.subplots(figsize=(8, max(5, .245*len(rows)+1.4)))
    for i, row in enumerate(rows):
        value = row['frequency']
        color = COLORS[row['Region']]
        ax.hlines(i, 0, value, color=color, linewidth=1.4)
        ax.plot(value, i, 'o', color=color, markersize=5)
        ax.text(value+1.5, i, f"{row['Count']}/{n}", va='center', fontsize=8)
    ax.set_yticks(range(len(rows)), [r['Mutation_full'] for r in rows])
    ax.invert_yaxis()
    ax.set_xlim(0, 119)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel('Frequency among Tunisian subclade K sequences (%)')
    ax.axvline(100*(n-1)/n, color='#555555', linestyle='--', linewidth=.8,
               label=f'Main substitution threshold ({n-1}/{n})')
    for side in ['top', 'right']:
        ax.spines[side].set_visible(False)
    ax.legend(loc='lower right', frameon=False, fontsize=8)
    fig.tight_layout()
    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    for suffix in ['png', 'pdf', 'svg']:
        destination = Path(str(args.out_prefix)+'.'+suffix)
        fig.savefig(destination, dpi=300, facecolor='white')
        print(destination)
    plt.close(fig)
    print(f'N={n}; main substitutions={main_count}; displayed={len(rows)}')

if __name__ == '__main__':
    main()
