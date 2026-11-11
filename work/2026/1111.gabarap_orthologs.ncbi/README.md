# NCBI ortholog and transcript workflow

This experiment retrieves NCBI orthologs for a human gene, then selects the
protein-coding transcript with the longest annotated protein for each gene.
Raw API responses are cached in `data/` for reproducibility.

## Step 1: Retrieve orthologs

From the repository root:

```sh
python work/2026/1111.gabarap_orthologs.ncbi/01_download_orthologs_from_ncbi.py GABARAP
```

Writes `GABARAP_ncbi_orthologs.tsv` beside the script. Replace `GABARAP` with
another human gene symbol as needed.

```sh
python work/2026/1111.gabarap_orthologs.ncbi/01_download_orthologs_from_ncbi.py GABARAP work/2026/1111.gabarap_orthologs.ncbi/results/01_GABARAP_ncbi_orthologs.tsv
```

The optional second argument sets the output TSV path. The table includes all
12 configured target species, preserves multiple orthologs, and uses `-` in
ortholog fields when no ortholog is found.

## Step 2: Add longest-protein transcript

```sh
python work/2026/1111.gabarap_orthologs.ncbi/02_add_transcript_with_longest_protein.py \
    work/2026/1111.gabarap_orthologs.ncbi/results/01_GABARAP_ncbi_orthologs.tsv
```

The output defaults to a sibling filename with its numeric prefix changed to
`02_` (or `_with_transcripts` appended when there is no numeric prefix). An
explicit output path overrides that default:

```sh
python work/2026/1111.gabarap_orthologs.ncbi/02_add_transcript_with_longest_protein.py \
    work/2026/1111.gabarap_orthologs.ncbi/results/01_GABARAP_ncbi_orthologs.tsv \
    work/2026/1111.gabarap_orthologs.ncbi/results/GABARAP_selected_transcripts.tsv
```

The script uses `taxon_id` and `ortholog_gene_id` (never gene symbols) with
`GET /gene/id/{gene_ids}/product_report`. It requests the NCBI Datasets gene
product report, equivalent to the documented CLI form
`datasets summary gene gene-id <GeneID> --report product --as-json-lines`.
Responses are cached individually in `data/`; repeated taxon/GeneID pairs do
not trigger repeated requests. Requests use a 60-second timeout, retries, and
the default NCBI Datasets limit of five requests per second.

Only transcripts typed `PROTEIN_CODING` or `PROTEIN_CODING_MODEL` with an
annotated protein accession and protein length are eligible. The transcript
with the greatest annotated protein length is selected. Ties prefer `NM_`
transcripts over `XM_` transcripts, then use transcript accession order
alphabetically. Transcript length is reported in nucleotides and protein
length in amino acids. No protein length is inferred from sequence length.

The output retains every input row and column, adding `transcript_accession`,
`transcript_length`, `protein_accession`, `protein_length`, and
`transcript_selection_method` (`longest_protein`). The row count is preserved.
Rows without a suitable transcript retain empty new fields and are listed in
the console summary. API failures are reported and counted; failed rows remain
in the output with empty transcript fields.

## Limitations

Selection reflects NCBI's current gene product annotation and available
RefSeq records. A longest annotated protein is not necessarily the biologically
dominant isoform. `PROTEIN_CODING_MODEL` records may be computationally
predicted. Cached JSON is reused on later runs; remove the relevant cache file
from `data/` to request a fresh response. API content and annotations can
change over time.
