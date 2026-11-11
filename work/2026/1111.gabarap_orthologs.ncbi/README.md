# NCBI ortholog and transcript workflow

Retrieves NCBI orthologs for a human gene and selects each gene's longest
annotated protein-coding transcript. Run commands from the repository root.

## Step 0: Download NCBI Gene–Ensembl mapping

From the experiment directory, run:

```sh
mkdir -p data/large_files/ncbi

wget -O data/large_files/ncbi/gene2ensembl.gz \
    https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz
```

This gzip-compressed file maps NCBI Gene IDs to Ensembl Gene IDs and includes
transcript and protein accession mappings. It is relatively large and is kept
locally for reuse in NCBI-to-Ensembl ID conversion.

## Step 1: Retrieve orthologs

```sh
python work/2026/1111.gabarap_orthologs.ncbi/01_download_orthologs_from_ncbi.py GABARAP
```

The command writes `01_GABARAP_ncbi_orthologs.tsv` beside the script.

## Step 2: Add longest-protein transcript

```sh
python work/2026/1111.gabarap_orthologs.ncbi/02_add_transcript_with_longest_protein.py work/2026/1111.gabarap_orthologs.ncbi/01_GABARAP_ncbi_orthologs.tsv
```
