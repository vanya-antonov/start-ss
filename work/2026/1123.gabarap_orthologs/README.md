# GABARAP ortholog transcript analysis

## Introduction

This analysis continues the GABARAP ortholog workflow using the NCBI
ortholog table and Ensembl release 116 canonical transcript tables prepared
in the previous analysis.

## Data

Run these commands from the repository root to copy the NCBI ortholog table
and link the canonical transcript tables from the 1118 analysis.

```sh
mkdir -p work/2026/1123.gabarap_orthologs
cd work/2026/1123.gabarap_orthologs
cp ../1118.gabarap_orthologs/01_GABARAP_ncbi_orthologs.tsv .
ls ../1118.gabarap_orthologs/*.canonical_transcripts.tsv | xargs -tI{} ln -s {}
cd -
```

## Commands

Activate the project analysis environment from the repository root:

```sh
conda activate start-ss-scripts
```
