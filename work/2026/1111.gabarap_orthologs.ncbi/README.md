# NCBI ortholog and transcript workflow

Retrieves NCBI orthologs for a human gene and selects each gene's longest
annotated protein-coding transcript. Run commands from the repository root.

## Step 1: Retrieve orthologs

```sh
python work/2026/1111.gabarap_orthologs.ncbi/01_download_orthologs_from_ncbi.py GABARAP
```

The command writes `01_GABARAP_ncbi_orthologs.tsv` beside the script.

## Step 2: Add longest-protein transcript

```sh
python work/2026/1111.gabarap_orthologs.ncbi/02_add_transcript_with_longest_protein.py work/2026/1111.gabarap_orthologs.ncbi/01_GABARAP_ncbi_orthologs.tsv
```
