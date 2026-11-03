# Human transcript workflow

Use the `start-ss-scripts` environment described in the [project README](../README.md).
Run the following commands from the repository root.

## Step 1. Download human MANE Select transcripts

```bash
python scripts/human/01_download_mane.py
```

The script downloads the latest NCBI MANE summary and writes the filtered
transcript table to `data/human/01_mane_select_transcripts.tsv`.

## Step 2. Download MANE transcript GenBank records

After creating the MANE Select table, set an NCBI contact email (recommended)
and optionally an API key, then run:

```bash
export NCBI_EMAIL="you@example.org"
# Optional: export NCBI_API_KEY="your-key"
python scripts/human/02_download_genbank.py
```

The script reads `data/human/01_mane_select_transcripts.tsv`, downloads versioned
RefSeq records in batches, and writes successful records to
`data/human/02_mane_transcripts.gb.gz`. Failed or missing accessions are listed in
`data/human/genbank_download_failures.tsv`.
