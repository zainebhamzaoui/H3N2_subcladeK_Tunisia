# Tunisian influenza A(H3N2) subclade K HA analysis

This repository accompanies the genetic characterization of Tunisian influenza A(H3N2) viruses during the 2025/2026 season, with a focus on subclade K.

The repository is being updated to reflect the revised manuscript. Existing scripts and supplementary tables correspond to an earlier analysis and should not be used to reproduce the revised results until this update is complete.

The earlier analysis included:

- HA sequence processing and Nextclade analysis
- extraction of HA amino-acid substitutions
- generation of HA substitution lollipop plots
- identification of closest international HA genetic neighbors
- regional comparison of HA substitution frequencies
- antigenic annotation of main HA substitutions
- generation of supplementary tables

## Repository structure

```text
github_H3N2_variantK_Tunisia/
├── README.md
├── data_example/
│   └── README_data.md
├── environment/
│   └── requirements.txt
├── results_example/
├── scripts/
│   ├── 01_run_nextclade_HA.sh
│   ├── 02_extract_HA_substitutions.py
│   ├── 03_generate_lollipop_HA.py
│   ├── 04_closest_HA_neighbors.py
│   ├── 05_regional_frequency_differences.py
│   ├── 06_antigenic_annotation.py
│   └── 07_generate_supplementary_tables.py
└── supplementary_tables/
    ├── Supplementary_Table_S1_complete_HA_substitutions.csv
    ├── Supplementary_Table_S2_closest_international_HA_neighbors.csv
    ├── Supplementary_Table_S3_Tunisian_HA_sequence_metadata.csv
    └── Supplementary_Tables_S1_S2_S3.xlsx
