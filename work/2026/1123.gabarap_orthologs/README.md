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
cp ../1118.gabarap_orthologs/01_GABARAP_ncbi_orthologs.tsv 00_GABARAP_orthologs.tsv
ls ../1118.gabarap_orthologs/*.canonical_transcripts.tsv | xargs -tI{} ln -s {}
cd -
```

## Commands

Activate the project analysis environment from the repository root:

```sh
conda activate start-ss-scripts
```

### Step 1. Add Ensembl canonical transcript annotations

Match each ortholog gene symbol to the corresponding species' local canonical
transcript table. The output includes Ensembl gene, transcript, and protein
IDs and transcript-relative CDS coordinates.

```sh
python work/2026/1123.gabarap_orthologs/01_add_canonical_transcripts.py GABARAP
```

The resulting file is `01_GABARAP_orthologs.tsv`.
