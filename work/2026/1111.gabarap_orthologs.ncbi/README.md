# NCBI ortholog and transcript workflow

Retrieves NCBI orthologs for a human gene and selects each gene's longest
annotated protein-coding transcript. Run commands from the repository root.

## Step 1: Retrieve orthologs

```sh
python work/2026/1111.gabarap_orthologs.ncbi/01_download_orthologs_from_ncbi.py GABARAP
python work/2026/1111.gabarap_orthologs.ncbi/01_download_orthologs_from_ncbi.py GABARAP work/2026/1111.gabarap_orthologs.ncbi/results/01_GABARAP_ncbi_orthologs.tsv
```

The first command writes the default TSV beside the script; the second sets an
explicit output path.

## Step 2: Add longest-protein transcript

```sh
python work/2026/1111.gabarap_orthologs.ncbi/02_add_transcript_with_longest_protein.py work/2026/1111.gabarap_orthologs.ncbi/results/01_GABARAP_ncbi_orthologs.tsv
python work/2026/1111.gabarap_orthologs.ncbi/02_add_transcript_with_longest_protein.py work/2026/1111.gabarap_orthologs.ncbi/results/01_GABARAP_ncbi_orthologs.tsv work/2026/1111.gabarap_orthologs.ncbi/results/GABARAP_selected_transcripts.tsv
```

The first command automatically changes a numbered input prefix to `02_`; the
second writes to the specified output path.
