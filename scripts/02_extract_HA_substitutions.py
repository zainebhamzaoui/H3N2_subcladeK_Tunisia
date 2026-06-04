#!/usr/bin/env python3

"""
Extract HA amino-acid substitutions from Nextclade output.

Input:
    results/nextclade_HA_Tunisia_plus_global.tsv

Outputs:
    results/HA_Tunisia_mutation_counts.csv
    results/HA_Tunisia_mutation_heatmap_matrix.csv
    results/HA_Tunisia_main_mutations_min4.csv

Notes:
    Raw GISAID FASTA files are not redistributed in this repository.
"""

import os
import re
import pandas as pd


NEXTCLADE_FILE = "results/nextclade_HA_Tunisia_plus_global.tsv"
OUT_COUNTS = "results/HA_Tunisia_mutation_counts.csv"
OUT_MATRIX = "results/HA_Tunisia_mutation_heatmap_matrix.csv"
OUT_MAIN = "results/HA_Tunisia_main_mutations_min4.csv"

os.makedirs("results", exist_ok=True)


def find_column(df, candidates):
    for cand in candidates:
        for col in df.columns:
            if col.lower() == cand.lower():
                return col
    for cand in candidates:
        for col in df.columns:
            if cand.lower() in col.lower():
                return col
    return None


def parse_ha_mutations(value):
    if pd.isna(value):
        return []

    text = str(value).strip()
    if text in ["", "nan", "None"]:
        return []

    parts = re.split(r"[,; ]+", text)
    mutations = []

    for part in parts:
        part = part.strip()
        if part.startswith(("SigPep:", "HA1:", "HA2:")):
            mutations.append(part)

    return mutations


def clean_mutation(mut):
    return (
        str(mut)
        .replace("SigPep:", "")
        .replace("HA1:", "")
        .replace("HA2:", "")
    )


def mutation_region(mut):
    if str(mut).startswith("SigPep:"):
        return "SP"
    if str(mut).startswith("HA1:"):
        return "HA1"
    if str(mut).startswith("HA2:"):
        return "HA2"
    return "Unknown"


def mutation_position(mut):
    match = re.search(r"(\d+)", str(mut))
    if match:
        return int(match.group(1))
    return 9999


def short_sequence_name(seq):
    seq = str(seq)

    parts = seq.split("|")
    for part in parts:
        if "A/Tunisia/" in part:
            seq = part
            break

    seq = seq.replace("HA|", "")
    seq = seq.split("|")[0]
    return seq


def main():
    if not os.path.exists(NEXTCLADE_FILE):
        raise FileNotFoundError(
            f"Input file not found: {NEXTCLADE_FILE}\n"
            "Please run 01_run_nextclade_HA.sh first or place the Nextclade TSV in results/."
        )

    df = pd.read_csv(NEXTCLADE_FILE, sep="\t", low_memory=False)

    seq_col = find_column(df, ["seqName", "sequence", "strain"])
    aa_col = find_column(
        df,
        [
            "aaSubstitutions",
            "substitutions",
            "privateAaMutations.substitutions",
        ],
    )

    if seq_col is None:
        raise ValueError("Could not identify sequence-name column in Nextclade TSV.")

    if aa_col is None:
        raise ValueError("Could not identify amino-acid substitution column in Nextclade TSV.")

    tun = df[df[seq_col].astype(str).str.contains("A/Tunisia/", regex=False)].copy()

    if len(tun) == 0:
        raise ValueError("No Tunisian sequences detected in Nextclade output.")

    tun["Sequence_short"] = tun[seq_col].apply(short_sequence_name)
    tun["HA_mutations_full"] = tun[aa_col].apply(parse_ha_mutations)

    all_mutations = sorted(
        set(m for muts in tun["HA_mutations_full"] for m in muts),
        key=lambda x: (mutation_region(x), mutation_position(x), x),
    )

    # Mutation count table
    rows = []
    for mut in all_mutations:
        count = tun["HA_mutations_full"].apply(lambda x: mut in x).sum()
        rows.append(
            {
                "Mutation_full": mut,
                "Mutation": clean_mutation(mut),
                "Region": mutation_region(mut),
                "Position": mutation_position(mut),
                "Count": int(count),
                "Frequency_percent": 100 * count / len(tun),
            }
        )

    counts = pd.DataFrame(rows)
    counts = counts.sort_values(["Count", "Region", "Position"], ascending=[False, True, True])
    counts.to_csv(OUT_COUNTS, index=False)

    # Presence/absence matrix
    matrix = pd.DataFrame(index=tun["Sequence_short"])
    for mut in all_mutations:
        matrix[clean_mutation(mut)] = tun["HA_mutations_full"].apply(lambda x: int(mut in x)).values

    matrix.to_csv(OUT_MATRIX)

    # Main mutations, threshold count >= 4
    main = counts[counts["Count"] >= 4].copy()
    main.to_csv(OUT_MAIN, index=False)

    print("Tunisian sequences:", len(tun))
    print("Total HA substitutions:", len(counts))
    print("Main HA substitutions with count >= 4:", len(main))
    print("Created:")
    print(OUT_COUNTS)
    print(OUT_MATRIX)
    print(OUT_MAIN)


if __name__ == "__main__":
    main()
