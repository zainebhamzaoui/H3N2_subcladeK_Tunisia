#!/usr/bin/env python3

"""
Identify closest international HA neighbors for Tunisian A/H3N2 sequences.

Input:
    alignment/HA_Tunisia_plus_global_2025Dec_2026Jan.aligned.fasta

Outputs:
    results/HA_Tunisia_closest_international_neighbors.csv
    results/HA_Tunisia_closest_neighbor_country_summary.csv
    results/HA_Tunisia_closest_neighbor_region_summary.csv
    results/HA_Tunisia_closest_neighbor_distance_summary.csv

Notes:
    Raw GISAID sequences are not redistributed in this repository.
    The aligned FASTA should contain Tunisian sequences and contemporaneous global HA sequences.
"""

import os
import re
from collections import Counter
import pandas as pd
from Bio import SeqIO


ALIGNMENT_FILE = "alignment/HA_Tunisia_plus_global_2025Dec_2026Jan.aligned.fasta"

OUT_NEIGHBORS = "results/HA_Tunisia_closest_international_neighbors.csv"
OUT_COUNTRY = "results/HA_Tunisia_closest_neighbor_country_summary.csv"
OUT_REGION = "results/HA_Tunisia_closest_neighbor_region_summary.csv"
OUT_DISTANCE = "results/HA_Tunisia_closest_neighbor_distance_summary.csv"

os.makedirs("results", exist_ok=True)


def is_tunisia(header):
    return "A/Tunisia/" in str(header)


def extract_strain(header):
    text = str(header)

    if "|" in text:
        parts = text.split("|")
        for part in parts:
            if part.startswith("A/"):
                return part

    match = re.search(r"(A/[A-Za-z0-9_\- ]+/[^| ]+/\d{4})", text)
    if match:
        return match.group(1)

    return text.split()[0]


def short_tunisia_name(header):
    strain = extract_strain(header)
    strain = strain.replace("A/Tunisia/", "A/Tunisia/")
    return strain


def extract_country(strain):
    s = str(strain)

    # Strain format usually: A/Country-or-region/ID/year
    parts = s.split("/")
    if len(parts) >= 2:
        location = parts[1].strip()
    else:
        location = s

    loc_low = location.lower()

    mapping = {
        "france": "France",
        "russia": "Russia",
        "saint petersburg": "Russia",
        "saint peterburg": "Russia",
        "st petersburg": "Russia",
        "china": "China",
        "jilin": "China",
        "khabarovsk": "Russia",
        "tianjin": "China",
        "canada": "Canada",
        "nova scotia": "Canada",
        "ontario": "Canada",
        "united states": "United States",
        "usa": "United States",
        "delaware": "United States",
        "netherlands": "Netherlands",
        "sweden": "Sweden",
        "morocco": "Morocco",
        "poland": "Poland",
        "lebanon": "Lebanon",
        "norway": "Norway",
        "united kingdom": "United Kingdom",
        "england": "United Kingdom",
        "australia": "Australia",
        "perth": "Australia",
    }

    for key, value in mapping.items():
        if key in loc_low:
            return value

    # fallback: use location string
    return location


def country_to_region(country):
    mapping = {
        "Tunisia": "North Africa",
        "Morocco": "North Africa",
        "France": "Europe",
        "Netherlands": "Europe",
        "Sweden": "Europe",
        "Poland": "Europe",
        "Norway": "Europe",
        "United Kingdom": "Europe",
        "Russia": "Asia/Middle East",
        "China": "Asia/Middle East",
        "Lebanon": "Asia/Middle East",
        "Canada": "North America",
        "United States": "North America",
        "Australia": "Oceania",
    }

    return mapping.get(country, "Other")


def extract_date(strain):
    # Try to extract year or full date-like part
    text = str(strain)

    m = re.search(r"(\d{4}-\d{2}-\d{2})", text)
    if m:
        return m.group(1)

    # fallback: final /2025 or /2026
    parts = text.split("/")
    if len(parts) >= 4:
        return parts[-1]

    return ""


def nucleotide_distance(seq1, seq2):
    s1 = str(seq1).upper()
    s2 = str(seq2).upper()

    distance = 0
    compared = 0

    for a, b in zip(s1, s2):
        if a in ["N", "-", "?"] or b in ["N", "-", "?"]:
            continue
        compared += 1
        if a != b:
            distance += 1

    return distance, compared


def main():
    if not os.path.exists(ALIGNMENT_FILE):
        raise FileNotFoundError(
            f"Input alignment not found: {ALIGNMENT_FILE}\n"
            "Please create the HA alignment before running this script."
        )

    records = list(SeqIO.parse(ALIGNMENT_FILE, "fasta"))

    if len(records) == 0:
        raise ValueError("No sequences found in alignment.")

    tunisian = [r for r in records if is_tunisia(r.description)]
    international = [r for r in records if not is_tunisia(r.description)]

    if len(tunisian) == 0:
        raise ValueError("No Tunisian sequences found in alignment.")

    if len(international) == 0:
        raise ValueError("No international sequences found in alignment.")

    rows = []

    for t in tunisian:
        best = None

        for g in international:
            dist, compared = nucleotide_distance(t.seq, g.seq)

            if compared == 0:
                continue

            if best is None or dist < best["Nucleotide_distance"]:
                strain = extract_strain(g.description)
                country = extract_country(strain)
                region = country_to_region(country)

                best = {
                    "Tunisian_sequence": short_tunisia_name(t.description),
                    "Closest_international_HA_neighbor": strain,
                    "Neighbor_country": country,
                    "Neighbor_region": region,
                    "Neighbor_collection_date": extract_date(strain),
                    "Nucleotide_distance": dist,
                    "Compared_positions": compared,
                }

        if best is not None:
            rows.append(best)

    neighbors = pd.DataFrame(rows)

    neighbors = neighbors.sort_values(
        ["Nucleotide_distance", "Neighbor_country", "Tunisian_sequence"]
    )

    neighbors.to_csv(OUT_NEIGHBORS, index=False)

    country_summary = (
        neighbors["Neighbor_country"]
        .value_counts()
        .reset_index()
    )
    country_summary.columns = ["Country", "Number_of_Tunisian_sequences"]
    country_summary.to_csv(OUT_COUNTRY, index=False)

    region_summary = (
        neighbors["Neighbor_region"]
        .value_counts()
        .reset_index()
    )
    region_summary.columns = ["Region", "Number_of_Tunisian_sequences"]
    region_summary.to_csv(OUT_REGION, index=False)

    distance_summary = (
        neighbors["Nucleotide_distance"]
        .value_counts()
        .sort_index()
        .reset_index()
    )
    distance_summary.columns = ["Nucleotide_distance", "Number_of_Tunisian_sequences"]
    distance_summary.to_csv(OUT_DISTANCE, index=False)

    print("Tunisian sequences:", len(tunisian))
    print("International sequences:", len(international))
    print("Closest-neighbor rows:", len(neighbors))
    print("Created:")
    print(OUT_NEIGHBORS)
    print(OUT_COUNTRY)
    print(OUT_REGION)
    print(OUT_DISTANCE)


if __name__ == "__main__":
    main()
