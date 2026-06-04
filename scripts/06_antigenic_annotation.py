#!/usr/bin/env python3

"""
Generate antigenic annotation figure for main HA substitutions.

Inputs:
    results/HA_Tunisia_main_mutations_min4.csv
    results/HA_global_region_mutation_delta_vs_Tunisia.csv

Outputs:
    results/HA_antigenic_annotation_main_substitutions.csv
    results/Figure_HA_antigenic_annotation_main_substitutions.png
    results/Figure_HA_antigenic_annotation_main_substitutions.pdf
"""

import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches


MAIN_MUT_FILE = "results/HA_Tunisia_main_mutations_min4.csv"
DELTA_FILE = "results/HA_global_region_mutation_delta_vs_Tunisia.csv"

OUT_CSV = "results/HA_antigenic_annotation_main_substitutions.csv"
OUT_PNG = "results/Figure_HA_antigenic_annotation_main_substitutions.png"
OUT_PDF = "results/Figure_HA_antigenic_annotation_main_substitutions.pdf"

os.makedirs("results", exist_ok=True)


def get_position(mut):
    match = re.search(r"(\d+)", str(mut))
    if match:
        return int(match.group(1))
    return 9999


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


def classify_regional_difference(value):
    if pd.isna(value):
        return "Not assessed"
    if value <= 2:
        return "Highly conserved"
    if value <= 5:
        return "Minor variation"
    return "Regional variation"


def main():
    if not os.path.exists(MAIN_MUT_FILE):
        raise FileNotFoundError(MAIN_MUT_FILE)

    if not os.path.exists(DELTA_FILE):
        raise FileNotFoundError(DELTA_FILE)

    df = pd.read_csv(MAIN_MUT_FILE)
    delta = pd.read_csv(DELTA_FILE, index_col=0)

    required = ["Mutation", "Region", "Position", "Count", "Frequency_percent"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in main mutation file: {missing}")

    df["Epitope"] = df.apply(lambda r: assign_epitope(r["Mutation"], r["Region"]), axis=1)

    max_diffs = []
    for mut in df["Mutation"]:
        if mut in delta.index:
            vals = pd.to_numeric(delta.loc[mut], errors="coerce").dropna()

            if "Tunisia" in vals.index:
                vals = vals.drop("Tunisia")

            if len(vals) > 0:
                max_diffs.append(abs(vals).max())
            else:
                max_diffs.append(pd.NA)
        else:
            max_diffs.append(pd.NA)

    df["Max_abs_regional_difference"] = max_diffs
    df["International_pattern"] = df["Max_abs_regional_difference"].apply(classify_regional_difference)

    df = df.sort_values(["Count", "Position"], ascending=[False, True]).copy()

    df.to_csv(OUT_CSV, index=False)

    epitope_colors = {
        "Epitope A": "#1f77b4",
        "Epitope B": "#ff7f0e",
        "Epitope C": "#2ca02c",
        "Epitope D": "#d62728",
        "Epitope E": "#9467bd",
        "Other": "#bdbdbd",
    }

    pattern_colors = {
        "Highly conserved": "#2166ac",
        "Minor variation": "#fddbc7",
        "Regional variation": "#b2182b",
        "Not assessed": "#eeeeee",
    }

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8.5,
        "axes.titlesize": 11.5,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
    })

    n = len(df)

    fig_w = 8.6
    fig_h = max(5.0, 0.32 * n + 1.5)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=300)

    ax.set_xlim(0, 10)
    ax.set_ylim(-0.5, n - 0.5)
    ax.invert_yaxis()

    x_mut = 0.25
    x_bar0 = 2.25
    x_bar1 = 5.10
    x_freq_text = 5.25
    x_epitope = 6.60
    x_delta = 8.10

    ax.text(x_mut, -0.95, "Mutation", ha="left", va="bottom", fontsize=9.2, fontweight="bold")
    ax.text((x_bar0 + x_bar1) / 2, -0.95, "Frequency in Tunisia", ha="center", va="bottom", fontsize=9.2, fontweight="bold")
    ax.text(x_epitope, -0.95, "Epitope", ha="center", va="bottom", fontsize=9.2, fontweight="bold")
    ax.text(x_delta, -0.95, "Max regional\ndifference", ha="center", va="bottom", fontsize=9.2, fontweight="bold")

    for i in range(n):
        if i % 2 == 0:
            ax.add_patch(
                patches.Rectangle(
                    (0, i - 0.5),
                    10,
                    1,
                    facecolor="#f7f7f7",
                    edgecolor="none",
                    zorder=0
                )
            )

    for i, (_, r) in enumerate(df.iterrows()):
        mut = str(r["Mutation"])
        count = int(round(r["Count"]))
        freq = float(r["Frequency_percent"])
        epitope = str(r["Epitope"])
        delta_value = float(r["Max_abs_regional_difference"])
        pattern = str(r["International_pattern"])

        ax.text(x_mut, i, mut, ha="left", va="center", fontsize=8.8)

        total_width = x_bar1 - x_bar0
        bar_width = total_width * freq / 100

        ax.add_patch(
            patches.Rectangle(
                (x_bar0, i - 0.17),
                total_width,
                0.34,
                facecolor="white",
                edgecolor="0.80",
                linewidth=0.35
            )
        )

        ax.add_patch(
            patches.Rectangle(
                (x_bar0, i - 0.17),
                bar_width,
                0.34,
                facecolor="#4c78a8",
                edgecolor="none"
            )
        )

        ax.text(
            x_freq_text,
            i,
            f"{count}/39 ({freq:.0f}%)",
            ha="left",
            va="center",
            fontsize=8.0
        )

        epi_color = epitope_colors.get(epitope, "#bdbdbd")
        epi_label = epitope.replace("Epitope ", "")

        ax.scatter(
            x_epitope,
            i,
            s=135,
            color=epi_color,
            edgecolor="black",
            linewidth=0.35,
            zorder=3
        )

        ax.text(
            x_epitope,
            i,
            epi_label,
            ha="center",
            va="center",
            fontsize=6.7,
            color="black",
            zorder=4
        )

        delta_color = pattern_colors.get(pattern, "#eeeeee")
        size = 75 + 18 * delta_value

        ax.scatter(
            x_delta,
            i,
            s=size,
            color=delta_color,
            edgecolor="black",
            linewidth=0.35,
            zorder=3
        )

        ax.text(
            x_delta,
            i,
            f"{delta_value:.0f}%",
            ha="center",
            va="center",
            fontsize=6.8,
            color="white" if pattern == "Regional variation" else "black",
            zorder=4
        )

    for x in [1.95, 5.95, 7.35]:
        ax.axvline(x, color="0.84", linewidth=0.65)

    ax.set_title(
        "Antigenic annotation of main HA substitutions in Tunisian A/H3N2 variant K sequences",
        fontweight="bold",
        pad=18
    )

    ax.set_xticks([])
    ax.set_yticks([])

    for spine in ax.spines.values():
        spine.set_visible(False)

    legend_y = n + 0.22
    legend_y2 = n + 0.72

    ax.text(0.25, legend_y, "Epitope:", ha="left", va="center", fontsize=8, fontweight="bold")

    legend_epitopes = ["Epitope A", "Epitope B", "Epitope C", "Epitope D", "Epitope E", "Other"]

    for j, epi in enumerate(legend_epitopes):
        xx = 1.05 + j * 0.82
        ax.scatter(xx, legend_y, s=72, color=epitope_colors[epi], edgecolor="black", linewidth=0.3)
        ax.text(xx + 0.16, legend_y, epi.replace("Epitope ", ""), ha="left", va="center", fontsize=7.2)

    ax.text(0.25, legend_y2, "Max regional difference:", ha="left", va="center", fontsize=8, fontweight="bold")

    legend_items = [
        ("Highly conserved", "#2166ac"),
        ("Minor variation", "#fddbc7"),
        ("Regional variation", "#b2182b")
    ]

    for j, (lab, col) in enumerate(legend_items):
        xx = 2.15 + j * 1.78
        ax.scatter(xx, legend_y2, s=72, color=col, edgecolor="black", linewidth=0.3)
        ax.text(xx + 0.16, legend_y2, lab, ha="left", va="center", fontsize=7.2)

    ax.set_ylim(n + 1.0, -1.25)

    plt.tight_layout()

    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    fig.savefig(OUT_PDF, bbox_inches="tight")
    plt.close(fig)

    print("Created:")
    print(OUT_CSV)
    print(OUT_PNG)
    print(OUT_PDF)


if __name__ == "__main__":
    main()
