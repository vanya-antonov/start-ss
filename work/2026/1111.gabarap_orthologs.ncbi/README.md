# NCBI ortholog and transcript workflow

Retrieves NCBI orthologs for a human gene, then selects each gene's longest
annotated protein-coding transcript. Run commands from the repository root.

## Step 1: Retrieve orthologs

```sh
python work/2026/1111.gabarap_orthologs.ncbi/01_download_orthologs_from_ncbi.py GABARAP
python work/2026/1111.gabarap_orthologs.ncbi/01_download_orthologs_from_ncbi.py GABARAP work/2026/1111.gabarap_orthologs.ncbi/results/01_GABARAP_ncbi_orthologs.tsv
```

The first command writes the default TSV beside the script; the optional second
argument in the other command sets an explicit output path. `human_gene_symbol`
and `human_gene_id` identify the human reference; `species` and `taxon_id`
identify the target; `ortholog_gene_id`, `ortholog_gene_symbol`,
`ortholog_description`, and `ortholog_type` describe the NCBI ortholog;
`refseq_accession` gives an available RefSeq accession.

## Step 2: Add longest-protein transcript

```sh
python work/2026/1111.gabarap_orthologs.ncbi/02_add_transcript_with_longest_protein.py work/2026/1111.gabarap_orthologs.ncbi/results/01_GABARAP_ncbi_orthologs.tsv
python work/2026/1111.gabarap_orthologs.ncbi/02_add_transcript_with_longest_protein.py work/2026/1111.gabarap_orthologs.ncbi/results/01_GABARAP_ncbi_orthologs.tsv work/2026/1111.gabarap_orthologs.ncbi/results/GABARAP_selected_transcripts.tsv
```

The default output changes a numbered filename prefix to `02_`; an optional
second argument overrides the output path. It uses NCBI gene-ID product reports
cached in `data/`, selecting the greatest annotated protein length, preferring
`NM_` accessions in ties, then sorting tied accessions alphabetically.
`transcript_accession` identifies the selected transcript;
`transcript_length` is in nucleotides; `protein_accession` identifies its
annotated protein; `protein_length` is in amino acids; and
`transcript_selection_method` is `longest_protein`. Transcript-related fields
are empty when no suitable transcript is found.
