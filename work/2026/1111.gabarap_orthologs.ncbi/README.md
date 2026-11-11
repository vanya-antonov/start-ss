# NCBI ortholog and transcript workflow

## Introduction

Retrieves NCBI orthologs for a human gene and selects each gene's longest
annotated protein-coding transcript.

## Downloads

From this experiment directory, download the gzip-compressed NCBI mapping file:

```sh
mkdir -p data/large_files/ncbi

wget -O data/large_files/ncbi/gene2ensembl.gz \
    https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz
```

This relatively large file maps NCBI Gene IDs to Ensembl Gene IDs and includes
transcript and protein accession mappings. It is stored locally for reuse in
NCBI-to-Ensembl ID conversion.

## Workflow

Run the commands from the repository root.

### Step 1: Retrieve orthologs

```sh
python work/2026/1111.gabarap_orthologs.ncbi/01_download_orthologs_from_ncbi.py GABARAP
```

Writes `01_GABARAP_ncbi_orthologs.tsv` beside the script.

### Step 2: Add longest-protein transcript

```sh
python work/2026/1111.gabarap_orthologs.ncbi/02_add_transcript_with_longest_protein.py work/2026/1111.gabarap_orthologs.ncbi/01_GABARAP_ncbi_orthologs.tsv
```

Reads the Step 1 table and writes a `02_`-prefixed TSV beside the input.
