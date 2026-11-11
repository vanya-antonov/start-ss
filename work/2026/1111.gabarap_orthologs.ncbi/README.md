# NCBI ortholog workflow

Resolves a human gene symbol to its NCBI Gene ID, retrieves NCBI's available
orthologs, and writes rows for configured nonhuman species. Raw API responses
are saved in `data/` for reproducibility.

Run from the repository root:

```sh
python work/2026/1111.gabarap_orthologs.ncbi/01_download_orthologs_from_ncbi.py GABARAP
```

Writes `GABARAP_ncbi_orthologs.tsv` beside the script. Replace `GABARAP` with
another human gene symbol as needed.

```sh
python work/2026/1111.gabarap_orthologs.ncbi/01_download_orthologs_from_ncbi.py GABARAP work/2026/1111.gabarap_orthologs.ncbi/01_GABARAP_ncbi_orthologs.tsv
```

The optional second argument specifies the output TSV path; relative paths are
resolved from the current directory. The table contains one row per returned
ortholog matched by NCBI taxonomy ID.
