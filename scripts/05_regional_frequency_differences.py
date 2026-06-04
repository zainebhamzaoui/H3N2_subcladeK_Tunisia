#!/usr/bin/env python3

"""
Calculate regional HA substitution frequency differences relative to Tunisia
and generate the regional-difference heatmap.

Inputs:
    results/nextclade_HA_Tunisia_plus_global.tsv
    results/HA_Tunisia_main_mutations_min4.csv

Outputs:
    results/HA_global_region_mutation_frequencies.csv
    results/HA_global_region_mutation_delta_vs_Tunisia.csv
    results/Figure_HA_regional_frequency_differences.png
    results/Figure_HA_regional_frequency_differences.pdf

Notes:
    Values correspond to:
    regional mutation frequency minus Tunisian mutation frequency.
"""

import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


NEXTCLADE_FILE = "results/nextclade_HA_Tunisia_plus_global.tsv"
MAIN_MUT_FILE = "results/HA_Tunisia_main_mutations_min4.csv"

OUT_FREQ = "results/HA_global_region_mutation_frequencies.csv"
OUT_DELTA = "results/HA_global_region_mutation_delta_vs_Tunisia.csv"
OUT_PNG = "results/Figure_HA_regional_frequency_differences.png"
OUT_PDF = "results/Figure_HA_regional_frequency_differences.pdf"

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
            mutations.append(clean_mutation(part))

    return mutations


def clean_mutation(mut):
    return (
        str(mut)
        .replace("SigPep:", "")
        .replace("HA1:", "")
        .replace("HA2:", "")
    )


def extract_country_from_seqname(seq):
    text = str(seq)

    parts = text.split("|")
    strain = None

    for p in parts:
        if p.startswith("A/"):
            strain = p
            break

    if strain is None:
        strain = text.split()[0]

    pieces = strain.split("/")
    if len(pieces) >= 2:
        return pieces[1].strip()

    return "Unknown"


def country_to_region(country):
    c = str(country).strip().lower()

    if "tunisia" in c:
        return "Tunisia"

    if "morocco" in c:
        return "Morocco"

    other_africa = [
        "kenya", "ghana", "cameroon", "south africa", "senegal",
        "egypt", "nigeria", "ethiopia", "algeria", "libya"
    ]
    if any(x in c for x in other_africa):
        return "Other Africa"

    europe = [
        "france", "netherlands", "sweden", "poland", "norway",
        "united kingdom", "england", "germany", "spain", "italy",
        "belgium", "denmark", "iceland", "ireland", "finland",
        "austria", "switzerland", "portugal", "czech", "slovenia"
    ]
    if any(x in c for x in europe):
        return "Europe"

    asia_middle_east = [
        "china", "russia", "saint petersburg", "saint peterburg",
        "jilin", "khabarovsk", "tianjin", "lebanon", "oman",
        "bahrain", "kazakhstan", "qatar", "saudi", "japan",
        "korea", "india", "thailand", "vietnam", "israel",
        "turkey", "iran", "iraq", "jordan", "uae", "emirates"
    ]
    if any(x in c for x in asia_middle_east):
        return "Asia/Middle East"

    north_america = [
        "canada", "united states", "usa", "delaware",
        "nova scotia", "ontario", "toronto"
    ]
    if any(x in c for x in north_america):
        return "North America"

    oceania = ["australia", "new zealand", "perth"]
    if any(x in c for x in oceania):
        return "Oceania"

    return "Other"


def main():
    if not os.path.exists(NEXTCLADE_FILE):
        raise FileNotFoundError(NEXTCLADE_FILE)

    if not os.path.exists(MAIN_MUT_FILE):
        raise FileNotFoundError(MAIN_MUT_FILE)

    df = pd.read_csv(NEXTCLADE_FILE, sep="\t", low_memory=False)
    main_mut = pd.read_csv(MAIN_MUT_FILE)

    seq_col = find_column(df, ["seqName", "sequence", "strain"])
    aa_col = find_column(
        df,
        ["aaSubstitutions", "substitutions", "privateAaMutations.substitutions"]
    )

    if seq_col is None:
        raise ValueError("Could not identify sequence-name column.")

    if aa_col is None:
        raise ValueError("Could not identify amino-acid substitution column.")

    selected_mutations = main_mut["Mutation"].astype(str).tolist()

    df["Country"] = df[seq_col].apply(extract_country_from_seqname)
    df["Region_Group"] = df["Country"].apply(country_to_region)
    df["HA_mutations"] = df[aa_col].apply(parse_ha_mutations)

    region_order = [
        "Tunisia",
        "Morocco",
        "Other Africa",
        "Europe",
        "Asia/Middle East",
        "North America",
        "Oceania",
        "Other",
    ]

    region_order = [r for r in region_order if r in df["Region_Group"].unique()]

    # Frequencies
    freq_table = pd.DataFrame(index=selected_mutations)

    region_counts = {}

    for region in region_order:
        sub = df[df["Region_Group"] == region].copy()
        n = len(sub)
        region_counts[region] = n

        values = []
        for mut in selected_mutations:
            if n == 0:
                values.append(np.nan)
            else:
                count = sub["HA_mutations"].apply(lambda x: mut in x).sum()
                values.append(100 * count / n)

        freq_table[region] = values

    freq_table.to_csv(OUT_FREQ)

    if "Tunisia" not in freq_table.columns:
        raise ValueError("Tunisia column not found. Cannot calculate delta.")

    delta = freq_table.subtract(freq_table["Tunisia"], axis=0)

    # Remove Tunisia column for final figure, but keep it in output table
    delta.to_csv(OUT_DELTA)

    plot_regions = [r for r in region_order if r != "Tunisia"]
    plot_delta = delta[plot_regions].copy()

    # Better labels with n
    xlabels = [f"{r}\n(n={region_counts.get(r, 0)})" for r in plot_regions]

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 12,
        "axes.labelsize": 10,
        "xtick.labelsize": 8,
        "ytick.labelsize": 9,
    })

    fig_h = max(5.4, 0.28 * len(selected_mutations) + 2.0)
    fig_w = 9.6

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=300)

    im = ax.imshow(
        plot_delta.values,
        cmap="RdBu_r",
        vmin=-10,
        vmax=10,
        aspect="auto"
    )

    ax.set_xticks(np.arange(len(plot_regions)))
    ax.set_xticklabels(xlabels, rotation=45, ha="right")

    ax.set_yticks(np.arange(len(selected_mutations)))
    ax.set_yticklabels(selected_mutations)

    ax.set_title(
        "Regional frequency differences of Tunisian HA substitutions",
        fontweight="bold",
        pad=14
    )

    ax.set_xlabel("Geographic region")
    ax.set_ylabel("Main Tunisian HA substitutions")

    # grid
    ax.set_xticks(np.arange(-0.5, len(plot_regions), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(selected_mutations), 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)

    # numbers
    for i in range(plot_delta.shape[0]):
        for j in range(plot_delta.shape[1]):
            val = plot_delta.iloc[i, j]
            if pd.isna(val):
                txt = "NA"
            elif round(val) > 0:
                txt = f"+{int(round(val))}"
            else:
                txt = f"{int(round(val))}"

            ax.text(j, i, txt, ha="center", va="center", fontsize=8)

    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label("Frequency difference vs Tunisia (%)")

    plt.tight_layout()

    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    fig.savefig(OUT_PDF, bbox_inches="tight")
    plt.close(fig)

    print("Created:")
    print(OUT_FREQ)
    print(OUT_DELTA)
    print(OUT_PNG)
    print(OUT_PDF)


if __name__ == "__main__":
    main()
