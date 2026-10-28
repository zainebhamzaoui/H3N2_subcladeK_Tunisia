# Input data and revised HA substitution workflow

## Sequence dataset

The phylogenetic sequence dataset is identified by GISAID EPI_SET_261027se:
https://doi.org/10.55876/gis8.261027se

Only HA was analyzed. Raw GISAID sequences are not supplied in this repository.
Retrieve sequences through GISAID with appropriate access credentials and
acknowledge their contributors. The EPI_SET identifies the phylogenetic dataset;
it does not identify the Nextclade reference dataset.

## Curated analysis subset

`data_example/Tunisian_K_37_accessions.txt` lists the 37 retained Tunisian
subclade K isolate accessions. The phylogeny contains 39 Tunisian viruses, but
the Tunisian J.2.2 and J.2.4.2 viruses are excluded from K substitution frequencies.
The accession list selects one curated record per strain. It does not remove
identical HA sequences from distinct strains.

## Required input

Place the Nextclade TSV export containing the curated sequences at:
`results/nextclade_HA.tsv`.

Sequence headers must contain both the strain name and EPI_ISL accession.
The extraction script checks the curated accessions, K assignments, quality
status and region-specific substitution labels. It stops if unknown amino-acid
positions or amino-acid indels require review. It does not run Nextclade itself.
The original Nextclade software version, reference dataset name and dataset
version still need to be documented. Do not assume the dataset currently
available online is identical to the one used for the existing results.

## Run from the repository root

Python 3 is required. The plotting script also requires matplotlib.

```bash
python3 scripts/02_extract_HA_substitutions.py \
  --nextclade results/nextclade_HA.tsv \
  --accessions data_example/Tunisian_K_37_accessions.txt \
  --outdir results

python3 scripts/03_generate_lollipop_HA.py

python3 scripts/06_antigenic_annotation.py

python3 scripts/07_generate_supplementary_tables.py
```

The lollipop displays all observed substitutions and marks the main-substitution
threshold. Main substitutions are present in all sequences or all but one
sequence (at least 36 of 37). Counts refer to differences from the Nextclade
reference, not necessarily to changes specific to K or differences from a vaccine.

## Expected results for the checked export

37 sequences; 33 distinct amino-acid substitutions; 17 main substitutions
(15 in 37/37 sequences and two in 36/37 sequences).

Script 02 outputs the complete count table, main-substitution table,
presence/absence matrix and an extraction summary with input SHA256 hashes.
Script 03 outputs PNG, PDF and editable SVG plots.
Script 06 outputs site annotations following Table 3 of Lees et al. (2010),
https://doi.org/10.1093/bioinformatics/btq160, distinguishing previously
identified residues from proposed additional assignments. The A-E scheme applies
to HA1, not to signal peptide or HA2. Location alone does not demonstrate altered
antigenicity.
Script 07 exports the substitution table (provisionally Supplementary Table S4)
and the main-substitution subset. Surveillance and sequence metadata tables are
prepared separately from their respective sources.

## Scope of this update

Closest-neighbor and regional-frequency analyses are excluded from the revised
manuscript. These commands document the revised substitution workflow only.
Phylogenetic inference, ancestral-state processing and final tree figure
preparation require separate documentation before the repository can be described
as reproducing all manuscript analyses. The existing script 01 and unversioned
environment requirements have not yet been fully updated.
