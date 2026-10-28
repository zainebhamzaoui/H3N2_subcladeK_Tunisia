# HA phylogeny and ancestral-state reconstruction

## Input and scope

The final phylogenetic analysis used 67 HA coding sequences: 39 Tunisian viruses,
26 international comparators and two vaccine reference viruses. The aligned
input was HA_67_aligned.fasta, with 1,701 nucleotide positions.
Sequence source: GISAID EPI_SET_261027se,
https://doi.org/10.55876/gis8.261027se . Only HA was analyzed.
Raw sequences and alignments are not supplied here.

MAFFT v7.526 was used for alignment. The exact alignment command and sequence
preparation steps are not documented in this file and still need to be recovered.
The commands below were copied from the original IQ-TREE log files.

## Software

IQ-TREE v3.1.2, executable iqtree3.
Run commands in a working directory containing the required inputs.
These are records of the original executions, not commands to rerun merely
for updating the repository. Preserve original outputs before any new execution.

## Maximum-likelihood inference

Source: H3N2_67.log.

```bash
iqtree3 -s HA_67_aligned.fasta -st DNA -m MFP -B 1000 -bnni -alrt 1000 -T AUTO -seed 20261026 --prefix H3N2_67
```

ModelFinder selected HKY+F+I under the Bayesian information criterion, as
reported in H3N2_67.iqtree. Branch support used 1,000 ultrafast bootstrap
replicates with bootstrap-tree NNI refinement and 1,000 SH-aLRT replicates.
The inference output is an unrooted tree.

## Ancestral nucleotide states

Source: H3N2_67_ASR.log.

```bash
iqtree3 -s HA_67_aligned.fasta -te H3N2_67.treefile -m HKY+F+I -asr -keep-ident -T 2 --prefix H3N2_67_ASR
```

This analysis fixes the previously inferred topology and estimates ancestral
nucleotide states under HKY+F+I. Identical sequences are retained. The original
command contains no explicit random seed; do not imply that the seed used for
the first analysis was also specified for this command.

The state output is H3N2_67_ASR.state. The associated report and tree are
H3N2_67_ASR.iqtree and H3N2_67_ASR.treefile.

## Rooting and branch annotation in the manuscript figure

The displayed tree was rooted on the branch separating A/Croatia/10136RV/2023
and A/District of Columbia/27/2023 from the remaining viruses. This display
rooting is separate from the unrooted IQ-TREE inference and is not specified
by either command above.

Ancestral nucleotide states were translated in the HA coding frame and selected
changes were mapped between adjacent reconstructed internal nodes using
separate HA1 and HA2 coordinates. The processing and figure-generation code
for these steps still needs to be added. These two commands alone do not
reproduce the annotated publication figure.

The ancestral-state report flags numerous near-zero internal branches.
Fine-scale branching and mapped ancestral changes require cautious
interpretation. Neither the tree nor branch annotations alone establish direct
transmission, geographic origin or the number of introductions.
