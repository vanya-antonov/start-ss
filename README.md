# start-ss
RNA secondary structures around the start codon.

## Development environment

Create the Conda environment from the project root:

```bash
conda env create -f environment.yml
conda activate start-ss
```

The environment uses Python 3.15 and Django 6.1.2. To apply later changes to
`environment.yml`, update the existing environment with:

```bash
conda env update -f environment.yml --prune
```
