#!/usr/bin/env bash

# Run Nextclade on combined Tunisian and global A/H3N2 HA sequences.
# Raw GISAID FASTA files are not redistributed in this repository.
# Users should prepare the input FASTA according to the README instructions.

set -euo pipefail

INPUT_FASTA="data/HA_Tunisia_plus_global_2025Dec_2026Jan.fasta"
DATASET="nextclade_h3n2_ha_dataset"
OUTDIR="results"

mkdir -p "${OUTDIR}"

nextclade run \
  --input-dataset "${DATASET}" \
  --output-tsv "${OUTDIR}/nextclade_HA_Tunisia_plus_global.tsv" \
  --output-tree "${OUTDIR}/nextclade_HA_Tunisia_plus_global.tree.json" \
  "${INPUT_FASTA}"
