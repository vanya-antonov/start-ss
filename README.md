# start-ss
RNA secondary structures around the start codon.

## Environments

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

## Data Log
### 2026
* [1103.human](docs/2026/1103.human/README.md) - instructions to
download MANE Select transcripts and their GenBank records
* [1109.mammals](docs/2026/1109.mammals/README.md) -- cross-species workflow

