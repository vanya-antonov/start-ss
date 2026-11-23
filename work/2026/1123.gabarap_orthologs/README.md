# GABARAP ortholog transcript analysis

## Introduction

This analysis continues the GABARAP ortholog workflow using the NCBI
ortholog table and Ensembl release 116 canonical transcript tables prepared
in the previous analysis.

## Data

Run these commands from the repository root to copy the NCBI ortholog table
and link the canonical transcript tables from the 1118 analysis.

```sh
mkdir -p work/2026/1123.gabarap_orthologs
cd work/2026/1123.gabarap_orthologs
cp ../1118.gabarap_orthologs/01_GABARAP_ncbi_orthologs.tsv 00_GABARAP_orthologs.tsv
ls ../1118.gabarap_orthologs/*.canonical_transcripts.tsv | xargs -tI{} ln -s {}
cd -
```

## Commands

Activate the project analysis environment from the repository root:

```sh
conda activate start-ss-scripts
```

### Step 1. Add Ensembl canonical transcript annotations

Match each ortholog gene symbol to the corresponding species' local canonical
transcript table. The output includes Ensembl gene, transcript, and protein
IDs and transcript-relative CDS coordinates.

```sh
python work/2026/1123.gabarap_orthologs/01_add_canonical_transcripts.py GABARAP
```

The resulting file is `01_GABARAP_orthologs.tsv`.

### Step 2. Extract start-codon windows

Use the canonical transcript IDs and CDS coordinates to retrieve each transcript
from its Ensembl release 116 cDNA FASTA. Extract 23 nt upstream, the annotated
3-nt start codon, and 24 nt downstream; invalid or incomplete windows are
listed in the errors TSV. The cDNA FASTAs are read from
`data/large_files/ensembl/cdna/`.

```sh
python work/2026/1123.gabarap_orthologs/02_extract_start_codon_sequences.py GABARAP
```

This writes `02_GABARAP_start_codon_windows.fasta` and
`02_GABARAP_start_codon_windows.tsv` for successful windows, plus
`02_GABARAP_extraction_errors.tsv` for failed rows.

### Step 3. Calculate minimum free energies

Fold each successful 50-nt sequence as RNA using ViennaRNA's default `RNA.fold`
parameters, without constraints or custom model settings. The script records
the predicted structure, MFE, and GC percentage; failed folds are written to an
error report.

```sh
python work/2026/1123.gabarap_orthologs/03_calculate_mfe.py GABARAP
```

This writes `03_GABARAP_rnafold_results.tsv`,
`03_GABARAP_start_codon_structure.tsv`, and, if any folds fail,
`03_GABARAP_rnafold_errors.tsv`.
