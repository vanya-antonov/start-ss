# NCBI ortholog and Ensembl ID workflow

Retrieves NCBI orthologs for a human gene and maps target NCBI Gene IDs to
Ensembl Gene IDs using NCBI's local mapping file.

## Downloads

From the repository root, download the gzip-compressed mapping file to
`data/large_files/ncbi/gene2ensembl.gz`. It includes NCBI Gene-to-Ensembl Gene
mappings and transcript/protein accession mappings; this large file is kept
locally for reuse.

```sh
mkdir -p data/large_files/ncbi

wget -O data/large_files/ncbi/gene2ensembl.gz \
    https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz
```

## Workflow

Run the commands from the repository root.

### Step 1: Retrieve orthologs

```sh
python work/2026/1111.gabarap_orthologs.ncbi/01_download_orthologs_from_ncbi.py GABARAP
```

Writes `01_GABARAP_ncbi_orthologs.tsv` beside the script.

### Step 2: Add Ensembl gene IDs

```sh
python work/2026/1111.gabarap_orthologs.ncbi/02_add_ensembl_gene_ids.py GABARAP
```

The gene name selects `01_GABARAP_ncbi_orthologs.tsv` beside the script. The
script maps rows by `(taxon_id, ortholog_gene_id)` and writes
`02_GABARAP_ncbi_orthologs.tsv`; `ensembl_gene_id` contains distinct matches
separated by semicolons, or is empty when no mapping exists.
