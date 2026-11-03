# start-ss
RNA secondary structures around the start codon.

## Environments

The Django project and bioinformatics scripts use separate Conda environments.
Create and activate the Django environment from the repository root with:

```bash
mamba env create -f django-env.yml
mamba activate start-ss
```

It uses Python 3.14 and Django 6.1.2. To apply changes to `django-env.yml`:

```bash
mamba env update -f django-env.yml --prune
```

Create and activate the shared scripts environment with:

```bash
mamba env create -f scripts-env.yml
mamba activate start-ss-scripts
```

It uses Python 3.11 with pandas, Biopython, requests, and ViennaRNA. Update it
after changing `scripts-env.yml` with:

```bash
mamba env update -f scripts-env.yml --prune
```

## Download human MANE Select transcripts

With `start-ss-scripts` active, run this from the repository root:

```bash
python scripts/human/01_download_mane.py
```

The script downloads the latest NCBI MANE summary and writes the filtered
transcript table to `data/human/mane_select_transcripts.tsv`.
```
