# GABARAP ortholog transcript analysis

## Introduction

This analysis continues the GABARAP ortholog workflow using the NCBI
ortholog table and Ensembl release 116 canonical transcript tables prepared
in the previous analysis.

## Data

Run these commands from the repository root to copy the NCBI ortholog table.
Canonical transcript tables are read directly from
`data/large_files/ensembl/gtf/` and must be generated first; see
`data/large_files/README.md`.

```sh
mkdir -p work/2026/1123.gabarap_orthologs
cd work/2026/1123.gabarap_orthologs
cp ../1118.gabarap_orthologs/01_GABARAP_ncbi_orthologs.tsv 00_GABARAP_orthologs.tsv
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

### Step 2. Retrieve full transcripts and extract start-codon windows

Use the canonical transcript IDs to retrieve each full cDNA sequence from its
Ensembl release 116 FASTA. The FASTA output includes the full transcript even
when CDS annotations do not permit window extraction, and preserves ambiguous
bases. Separately, extract 23 nt upstream, the annotated 3-nt start codon, and
24 nt downstream for the results TSV; invalid or incomplete windows are listed
in the errors TSV. In the full-transcript FASTA, the annotated CDS is uppercase
and the 5′ and 3′ UTRs are lowercase. A transcript with missing or invalid CDS
coordinates is kept entirely lowercase. The cDNA FASTAs are read from
`data/large_files/ensembl/cdna/`.

```sh
python work/2026/1123.gabarap_orthologs/02_extract_start_codon_sequences.py GABARAP
```

This writes full-length sequences to `02_GABARAP_transcripts.fasta`,
`02_GABARAP_orthologs.tsv` with the 50-nt window in its `sequence` column (or
`-` when extraction fails), and `02_GABARAP_extraction_errors.tsv` for failed
window extractions. The TSV retains the annotated ortholog columns and adds the
extracted sequence and GC percentage.

### Step 3. Calculate minimum free energies

Fold each successful 50-nt sequence from the step 2 results TSV as RNA using
ViennaRNA's default `RNA.fold` parameters, without constraints or custom model
settings. The script records the predicted structure, MFE, and GC percentage;
failed folds are written to an error report. An optional `--fasta` argument can
provide a legacy 50-nt window FASTA as fallback for rows missing a TSV sequence;
the full-transcript FASTA from step 2 is not a folding input.

```sh
python work/2026/1123.gabarap_orthologs/03_calculate_mfe.py GABARAP
```

This writes `03_GABARAP_rnafold_results.tsv`,
`03_GABARAP_start_codon_structure.tsv`, and, if any folds fail,
`03_GABARAP_rnafold_errors.tsv`.

### Step 4. Estimate structure conservation

Compare each valid nonhuman dot-bracket structure with the human transcript.
All valid structures must have the same length. Set the maximum number of
character mismatches allowed for a structure to count as conserved.

```sh
python work/2026/1123.gabarap_orthologs/04_estimate_structure_conservation.py GABARAP
```

The summary reports the number of valid ortholog structures, conserved
orthologs, and the conservation percentage. It is written to
`04_GABARAP_structure_conservation.tsv`.
