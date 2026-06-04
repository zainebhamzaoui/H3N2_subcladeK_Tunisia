#!/usr/bin/env python3

"""
Generate lollipop plot of main HA amino-acid substitutions.

Input:
    results/HA_Tunisia_main_mutations_min4.csv

Outputs:
    results/Figure_HA_lollipop_main_substitutions.png
    results/Figure_HA_lollipop_main_substitutions.pdf
"""

import os
import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

INPUT = "results/HA_Tunisia_main_mutations_min4.csv"

OUT_PNG = "results/Figure_HA_lollipop_main_substitutions.png"
OUT_PDF = "results/Figure_HA_lollipop_main_substitutions.pdf"

os.makedirs("results", exist_ok=True)


def get_position(mut):
    m = re.search(r"(\d+)", str(mut))
    return int(m.group(1)) if m else 9999


def plot_position(row):
    region = str(row["Region"])
    pos = int(row["Position"])

    if region == "SP":
        return pos
    if region == "HA1":
        return 16 + pos
    if region == "HA2":
        return 344 + pos

    return pos


def assign_epitope(mut, region):
    pos = get_position(mut)

    if region != "HA1":
        return "Other"

    epitope_A = set(list(range(121, 146)) + list(range(153, 159)) + list(range(171, 174)))
    epitope_B = set(list(range(155, 161)) + list(range(186, 199)))
    epitope_C = set(list(range(44, 58)) + list(range(276, 281)))
    epitope_D = set(list(range(96, 103)) + list(range(201, 225)))
    epitope_E = set(list(range(62, 83)) + list(range(260, 270)))

    if pos in epitope_A:
        return "Epitope A"
    if pos in epitope_B:
        return "Epitope B"
    if pos in epitope_C:
        return "Epitope C"
    if pos in epitope_D:
        return "Epitope D"
    if pos in epitope_E:
        return "Epitope E"

    return "Other"


def main():
    if not os.path.exists(INPUT):
        raise FileNotFoundError(
            f"Input file not found: {INPUT}. "
            "Please run 02_extract_HA_substitutions.py first."
        )

    df = pd.read_csv(INPUT)

    required = ["Mutation", "Region", "Position", "Count"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in input file: {missing}")

    df["Plot_position"] = df.apply(plot_position, axis=1)
    df["Epitope"] = df.apply(lambda r: assign_epitope(r["Mutation"], r["Region"]), axis=1)

    df = df.sort_values("Plot_position")

    epitope_colors = {
        "Epitope A": "#1f77b4",
        "Epitope B": "#ff7f0e",
        "Epitope C": "#2ca02c",
        "Epitope D": "#d62728",
        "Epitope E": "#9467bd",
        "Other": "#bdbdbd",
    }

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 12,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
    })

    fig, ax = plt.subplots(figsize=(12.5, 4.8), dpi=300)

    # Background regions
    ax.axvspan(0, 16, color="#eaf4fb", alpha=0.8)
    ax.axvspan(16, 344, color="#fdeaea", alpha=0.55)
    ax.axvspan(344, 560, color="#eaf6e9", alpha=0.8)

    # Region labels
    ymax = max(df["Count"]) + 6
    ax.text(8, ymax - 1.0, "SP", color="#1f77b4", ha="center", va="top", fontsize=12, fontweight="bold")
    ax.text(180, ymax - 1.0, "HA1", color="#d62728", ha="center", va="top", fontsize=12, fontweight="bold")
    ax.text(455, ymax - 1.0, "HA2", color="#2ca02c", ha="center", va="top", fontsize=12, fontweight="bold")

    # Boundaries
    ax.axvline(16, color="0.45", linestyle="--", linewidth=0.8)
    ax.axvline(344, color="0.45", linestyle="--", linewidth=0.8)

    # Lollipop
    for _, r in df.iterrows():
        x = r["Plot_position"]
        y = r["Count"]
        color = epitope_colors.get(r["Epitope"], "#bdbdbd")

        ax.vlines(x, 0, y, color="#d62728" if r["Region"] == "HA1" else color, linewidth=1.4, alpha=0.9)
        ax.scatter(x, y, s=80, color=color, edgecolor="black", linewidth=0.5, zorder=3)

        label = str(r["Mutation"])
        ax.text(x, y + 0.7, label, ha="center", va="bottom", rotation=90, fontsize=7.5)

    ax.set_title(
        "Main HA amino-acid substitutions in Tunisian A/H3N2 variant K sequences",
        fontweight="bold",
        pad=12
    )
    ax.set_xlabel("Position along HA protein")
    ax.set_ylabel("Number of Tunisian sequences")

    ax.set_xlim(0, 560)
    ax.set_ylim(0, ymax)

    ax.grid(axis="y", linestyle="--", alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    legend_handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=epitope_colors[e],
               markeredgecolor="black", markersize=7, label=e.replace("Epitope ", "Epitope "))
        for e in ["Epitope A", "Epitope B", "Epitope C", "Epitope D", "Epitope E", "Other"]
    ]

    ax.legend(
        handles=legend_handles,
        ncol=6,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.18)
    )

    plt.tight_layout()

    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    fig.savefig(OUT_PDF, bbox_inches="tight")
    plt.close(fig)

    print("Created:")
    print(OUT_PNG)
    print(OUT_PDF)


if __name__ == "__main__":
    main()
