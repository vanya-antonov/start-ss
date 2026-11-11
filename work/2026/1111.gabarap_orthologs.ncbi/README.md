# NCBI ortholog and transcript workflow

Retrieves NCBI orthologs for a human gene and selects each gene's longest
annotated protein-coding transcript.

## Commands

Run from the repository root. Step 1 writes its default TSV beside the script;
the optional second argument sets an explicit output path.

```sh
python work/2026/1111.gabarap_orthologs.ncbi/01_download_orthologs_from_ncbi.py GABARAP
python work/2026/1111.gabarap_orthologs.ncbi/01_download_orthologs_from_ncbi.py GABARAP work/2026/1111.gabarap_orthologs.ncbi/results/01_GABARAP_ncbi_orthologs.tsv
```

Step 2 automatically changes a numbered input prefix to `02_`; an optional
output path overrides that name.

```sh
python work/2026/1111.gabarap_orthologs.ncbi/02_add_transcript_with_longest_protein.py work/2026/1111.gabarap_orthologs.ncbi/results/01_GABARAP_ncbi_orthologs.tsv
python work/2026/1111.gabarap_orthologs.ncbi/02_add_transcript_with_longest_protein.py work/2026/1111.gabarap_orthologs.ncbi/results/01_GABARAP_ncbi_orthologs.tsv work/2026/1111.gabarap_orthologs.ncbi/results/GABARAP_selected_transcripts.tsv
```

Step 2 uses NCBI Datasets gene-ID product reports, cached under `data/`. It
selects the longest annotated protein, prefers `NM_` accessions in length ties,
then orders tied transcript accessions alphabetically. Missing transcript
fields are left empty.

## Output columns

Step 1's `01_...tsv` contains the human gene symbol and GeneID, target species
and taxon ID, ortholog GeneID, symbol, description, orthology type, and an
available RefSeq accession. One row is written per ortholog; absent species
use `-` for ortholog fields.

Step 2 preserves those columns and appends:

- `transcript_accession`: selected NCBI transcript accession.
- `transcript_length`: transcript length in nucleotides.
- `protein_accession`: annotated protein accession.
- `protein_length`: annotated protein length in amino acids.
- `transcript_selection_method`: selection label (`longest_protein`).
