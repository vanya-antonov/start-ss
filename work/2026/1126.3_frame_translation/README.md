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

Generate an interactive HTML plot with the annotated CDS and three transcript-based reading frames. CDS coordinates are 1-based and inclusive:

```sh
python work/2026/1126.3_frame_translation/01_plot_3_frames.py \
  --fasta work/2026/1126.3_frame_translation/00_GABARAP_transcripts.fasta \
  --cds_start 120 \
  --cds_end 470 \
  --output work/2026/1126.3_frame_translation/01_GABARAP_3_frames.html
```

The FASTA must contain one transcript record.
