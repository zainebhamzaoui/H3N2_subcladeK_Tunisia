#!/usr/bin/env python3

"""
Generate supplementary tables for the Tunisian A/H3N2 variant K HA analysis.

Outputs:
    supplementary_tables/Supplementary_Table_S1_complete_HA_substitutions.csv
    supplementary_tables/Supplementary_Table_S2_closest_international_HA_neighbors.csv
    supplementary_tables/Supplementary_Table_S3_Tunisian_HA_sequence_metadata.csv
    supplementary_tables/Supplementary_Tables_S1_S2_S3.xlsx
"""

import os
import re
import pandas as pd
from pathlib import Path


OUTDIR = Path("supplementary_tables")
OUTDIR.mkdir(parents=True, exist_ok=True)

S1_OUT = OUTDIR / "Supplementary_Table_S1_complete_HA_substitutions.csv"
S2_OUT = OUTDIR / "Supplementary_Table_S2_closest_international_HA_neighbors.csv"
S3_OUT = OUTDIR / "Supplementary_Table_S3_Tunisian_HA_sequence_metadata.csv"
XLSX_OUT = OUTDIR / "Supplementary_Tables_S1_S2_S3.xlsx"

S1_INPUT = "results/HA_Tunisia_mutation_counts.csv"
S2_INPUT = "results/HA_Tunisia_closest_international_neighbors.csv"
S3_INPUT = "results/nextclade_HA_Tunisia_plus_global.tsv"


def get_position(mut):
    match = re.search(r"(\d+)", str(mut))
    return int(match.group(1)) if match else 9999


def short_sequence_name(seq):
    seq = str(seq)
    parts = seq.split("|")

    for part in parts:
        if "A/Tunisia/" in part:
            return part

    return seq.split()[0]


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


def make_s1():
    if not os.path.exists(S1_INPUT):
        raise FileNotFoundError(S1_INPUT)

    df = pd.read_csv(S1_INPUT)

    required = ["Mutation", "Position", "Count", "Frequency_percent", "Region"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns for S1: {missing}")

    # Add epitope column if absent
    if "Epitope" not in df.columns:
        df["Epitope"] = "Not assigned"

    s1 = df[["Mutation", "Position", "Count", "Frequency_percent", "Region", "Epitope"]].copy()
    s1 = s1.rename(columns={"Region": "HA_region"})
    s1 = s1.sort_values(["Position", "Mutation"])

    s1.to_csv(S1_OUT, index=False)
    return s1


def make_s2():
    if not os.path.exists(S2_INPUT):
        raise FileNotFoundError(S2_INPUT)

    df = pd.read_csv(S2_INPUT)

    expected = [
        "Tunisian_sequence",
        "Closest_international_HA_neighbor",
        "Neighbor_country",
        "Neighbor_region",
        "Neighbor_collection_date",
        "Nucleotide_distance",
        "Compared_positions",
    ]

    existing = [c for c in expected if c in df.columns]
    s2 = df[existing].copy()

    s2 = s2.sort_values(
        ["Nucleotide_distance", "Neighbor_country", "Tunisian_sequence"],
        na_position="last"
    )

    s2.to_csv(S2_OUT, index=False)
    return s2


def make_s3():
    if not os.path.exists(S3_INPUT):
        raise FileNotFoundError(S3_INPUT)

    df = pd.read_csv(S3_INPUT, sep="\t", low_memory=False)

    seq_col = find_column(df, ["seqName", "sequence", "strain"])
    clade_col = find_column(df, ["clade"])
    subclade_col = find_column(df, ["proposedSubclade", "subclade"])

    if seq_col is None:
        raise ValueError("Could not identify sequence-name column for S3.")

    tun = df[df[seq_col].astype(str).str.contains("A/Tunisia/", regex=False)].copy()

    rows = []
    for _, r in tun.iterrows():
        seq = str(r[seq_col])
        rows.append({
            "Sequence_name": short_sequence_name(seq),
            "Country": "Tunisia",
            "Segment": "HA",
            "Clade": r[clade_col] if clade_col else "",
            "Subclade": r[subclade_col] if subclade_col else "",
        })

    s3 = pd.DataFrame(rows).drop_duplicates()
    s3 = s3.sort_values("Sequence_name")

    s3.to_csv(S3_OUT, index=False)
    return s3


def main():
    s1 = make_s1()
    s2 = make_s2()
    s3 = make_s3()

    with pd.ExcelWriter(XLSX_OUT, engine="openpyxl") as writer:
        s1.to_excel(writer, sheet_name="Table S1 HA substitutions", index=False)
        s2.to_excel(writer, sheet_name="Table S2 closest neighbors", index=False)
        s3.to_excel(writer, sheet_name="Table S3 metadata", index=False)

    print("Created supplementary tables:")
    print(S1_OUT)
    print(S2_OUT)
    print(S3_OUT)
    print(XLSX_OUT)
    print()
    print("Rows:")
    print("S1:", len(s1))
    print("S2:", len(s2))
    print("S3:", len(s3))


if __name__ == "__main__":
    main()
