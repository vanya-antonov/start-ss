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

## Step 3. Extract sequences around translation start codons

With `start-ss-scripts` active, run:

```bash
python scripts/human/03_extract_sequences.py
```

The script reads `data/human/02_mane_transcripts.gb.gz` and extracts validated
50-nt transcript-oriented windows (23 nt upstream, the annotated ATG, and 24 nt
downstream). Successful windows are written to
`data/human/03_start_codon_windows.fasta` and
`data/human/03_start_codon_windows.tsv`. Transcripts that cannot be processed
are listed in `data/human/03_extraction_errors.tsv`.

## Step 4. ViennaRNA Validation
Check the ViennaRNA command-line installation:

```bash
RNAfold --version
```

### Validate the GABARAP MFE calculation

Convert the reference DNA sequence to RNA by replacing `T` with `U`. The
expected minimum free energy (MFE) is `-14.80 kcal/mol`. With the environment
active, run this Python API check:

```bash
python - <<'PY'
import RNA

sequence_dna = "TGTCCAGCCCGGTTCCCGGGAGGATGAAGTTCGTGTACAAAGAAGAGCAT"
sequence_rna = sequence_dna.replace("T", "U")
expected_mfe = -14.80
tolerance = 0.01

structure, mfe = RNA.fold(sequence_rna)
print(f"Structure: {structure}")
print(f"Minimum free energy: {mfe:.2f} kcal/mol")
if abs(mfe - expected_mfe) <= tolerance:
    print(f"MFE matches expected value ({expected_mfe:.2f} kcal/mol).")
else:
    print(f"MFE does not match expected value ({expected_mfe:.2f} kcal/mol).")
PY
```

## Step 5. Calculate RNA secondary structure and MFE

With the `start-ss-scripts` environment active, run from the repository root:

```bash
python scripts/human/05_calculate_mfe.py
```

The default is one worker process. To enable multiprocessing, set the worker
count with `--num_cpu`; for example:

```bash
python scripts/human/05_calculate_mfe.py --num_cpu 4
```

The script reads metadata from `data/human/03_start_codon_windows.tsv` and
sequences from `data/human/03_start_codon_windows.fasta`. If an input TSV has a
`sequence` column, it uses that; otherwise it joins sequences by transcript ID
from the FASTA. This keeps nucleotide sequences out of TSV files intended for
GitHub.

The outputs are:

- `data/human/05_rnafold_results.tsv`: `gene_name`, `transcript`, `structure`,
  and `MFE` (kcal/mol, rounded to two decimals).
- `results/human/05_mane_start_codon_structure.tsv`: exactly `gene_name`,
  `transcript`, `GC` (rounded to one decimal), and `free_energy` (kcal/mol,
  rounded to two decimals).
- `data/human/05_rnafold_errors.tsv`: transcript, gene name, error type, and
  description for sequences that could not be folded.

The environment pins ViennaRNA 2.7.2 (`viennarna=2.7.2`). The script logs the
Python binding version and uses `RNA.fold()` with default RNAfold settings:
RNA input, default temperature and model details, and no constraints or custom
parameters. It does not change settings to force the reference result.

### Independently validate GABARAP

The reference transcript is `NM_007278.2`; its 50-nt DNA window has expected GC
content 54% and expected MFE `-14.80 kcal/mol`. The following independent API
check converts T to U, verifies the sequence length, and compares the MFE with a
0.01 kcal/mol tolerance:

```bash
python - <<'PY'
import RNA

sequence_dna = "TGTCCAGCCCGGTTCCCGGGAGGATGAAGTTCGTGTACAAAGAAGAGCAT"
assert len(sequence_dna) == 50
sequence_rna = sequence_dna.replace("T", "U")
structure, mfe = RNA.fold(sequence_rna)
expected_mfe = -14.80
tolerance = 0.01
print(f"Structure: {structure}")
print(f"MFE: {mfe:.2f} kcal/mol")
print(f"Matches expected MFE: {abs(mfe - expected_mfe) <= tolerance}")
PY
```

