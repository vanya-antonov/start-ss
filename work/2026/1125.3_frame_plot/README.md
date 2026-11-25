# Three-frame transcript visualization

## Introduction

Create a standalone interactive Plotly page showing all FASTA transcripts in
three translated reading frames. The full transcript is shown; frame 0 uses the
reading phase anchored at each transcript's annotated CDS start when available,
and frames +1 and +2 shift that phase by one and two nucleotides. Frame colors
are red, green, and blue. ATG start codons are marked with white
dashes, stop codons (TAA, TAG, TGA) with black dashes, and the annotated CDS
interval is shaded on frame 0. FASTA IDs are expected as
`species_name|transcript_id`.

The frame and start/stop marker conventions follow the
[RiboCrypt tutorial](https://bioconductor.statistik.uni-dortmund.de/packages/3.23/bioc/vignettes/RiboCrypt/inst/doc/RiboCrypt_app_tutorial.html).

## Data

The inputs below are copied from the GABARAP ortholog experiment. The transcript
FASTA contains Ensembl release 116 canonical cDNA sequences with CDS in uppercase
and UTRs in lowercase; the TSV contains ortholog IDs and transcript-relative CDS
coordinates.

Run from the repository root:

```sh
mkdir -p work/2026/1125.3_frame_plot
cp work/2026/1123.gabarap_orthologs/02_GABARAP_transcripts.fasta work/2026/1125.3_frame_plot/00_GABARAP_transcripts.fasta
cp work/2026/1123.gabarap_orthologs/04_GABARAP_orthologs.tsv work/2026/1125.3_frame_plot/00_GABARAP_orthologs.tsv
```

## Commands

Activate the project scripts environment from the repository root:

```sh
conda activate start-ss-scripts
```

### Step 1. Plot three transcript reading frames

Translate every input transcript in three frames and write a self-contained
interactive HTML page. The default annotation TSV is the copied `00_` file.

```sh
python work/2026/1125.3_frame_plot/01_plot_transcript_frames.py work/2026/1125.3_frame_plot/00_GABARAP_transcripts.fasta
```

The page is written to
`work/2026/1125.3_frame_plot/01_GABARAP_transcript_frames.html`. Use `--output`
to choose another path or `--annotations` to provide a different CDS TSV.

Plotly is included in `scripts-env.yml` from the existing `conda-forge` channel
to generate the interactive, self-contained HTML visualization.
