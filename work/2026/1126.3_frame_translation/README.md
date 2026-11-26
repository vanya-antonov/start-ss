# Data

```sh
mkdir -p work/2026/1126.3_frame_translation
cp work/2026/1123.gabarap_orthologs/02_GABARAP_orthologs.tsv work/2026/1126.3_frame_translation/00_GABARAP_orthologs.tsv
cp work/2026/1123.gabarap_orthologs/02_GABARAP_transcripts.fasta work/2026/1126.3_frame_translation/00_GABARAP_transcripts.fasta
```

## Commands

Activate the project scripts environment from the repository root:

```sh
conda activate start-ss-scripts
```

### Step 1. Plot three-frame translations

Generate an interactive HTML plot with the annotated CDS and the three transcript-based reading frames:

```sh
python work/2026/1126.3_frame_translation/01_plot_3_frames.py work/2026/1126.3_frame_translation/00_GABARAP_transcripts.fasta
```

The page is written to `work/2026/1126.3_frame_translation/01_GABARAP_3_frames.html`.
