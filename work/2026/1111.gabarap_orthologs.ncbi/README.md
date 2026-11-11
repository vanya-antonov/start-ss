# NCBI ortholog and transcript workflow

Retrieves NCBI orthologs for a human gene, then selects each gene's longest
annotated protein-coding transcript. Run commands from the repository root.

## Step 1: Retrieve orthologs

```sh
python work/2026/1111.gabarap_orthologs.ncbi/01_download_orthologs_from_ncbi.py GABARAP
python work/2026/1111.gabarap_orthologs.ncbi/01_download_orthologs_from_ncbi.py GABARAP work/2026/1111.gabarap_orthologs.ncbi/results/01_GABARAP_ncbi_orthologs.tsv
```

The first command writes the default TSV beside the script; the optional second
argument in the other command sets an explicit output path.

## Step 2: Add longest-protein transcript

```sh
python work/2026/1111.gabarap_orthologs.ncbi/02_add_transcript_with_longest_protein.py work/2026/1111.gabarap_orthologs.ncbi/results/01_GABARAP_ncbi_orthologs.tsv
python work/2026/1111.gabarap_orthologs.ncbi/02_add_transcript_with_longest_protein.py work/2026/1111.gabarap_orthologs.ncbi/results/01_GABARAP_ncbi_orthologs.tsv work/2026/1111.gabarap_orthologs.ncbi/results/GABARAP_selected_transcripts.tsv
```

The default output changes a numbered filename prefix to `02_`; an optional
second argument overrides the output path. It uses NCBI gene-ID product reports
cached in `data/`, selecting the greatest annotated protein length, preferring
`NM_` accessions in ties, then sorting tied accessions alphabetically.

## Output columns

**Step 1:**

- `human_gene_symbol`, `human_gene_id`: human reference symbol and NCBI GeneID.
- `species`, `taxon_id`: target species and NCBI taxonomy ID.
- `ortholog_gene_id`, `ortholog_gene_symbol`: target NCBI GeneID and symbol.
- `ortholog_description`, `ortholog_type`: NCBI description and reported type.
- `refseq_accession`: available RefSeq accession; empty if unavailable.

**Step 2 appends:**

- `transcript_accession`, `transcript_length`: selected transcript and its length in nucleotides.
- `protein_accession`, `protein_length`: annotated protein and its length in amino acids.
- `transcript_selection_method`: selection label (`longest_protein`); new fields are empty when no suitable transcript is found.
