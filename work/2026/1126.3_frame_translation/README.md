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

Generate vertically stacked interactive plots for every transcript in the FASTA. CDS coordinates are looked up from the canonical transcript TSVs.

```sh
python work/2026/1126.3_frame_translation/01_plot_3_frames.py \
  --fasta work/2026/1126.3_frame_translation/00_GABARAP_transcripts.fasta \
  --output work/2026/1126.3_frame_translation/GABARAP_all_species.html
```

### Step 2. Find the first two CDS candidates

Search every transcript for its first two ATG codons and their first in-frame stop codons. Optionally translate complete CDS 2 candidates to protein FASTA.

```sh
python work/2026/1126.3_frame_translation/02_find_cds.py \
  --fasta work/2026/1126.3_frame_translation/00_GABARAP_transcripts.fasta \
  --output work/2026/1126.3_frame_translation/02_GABARAP_cds.tsv \
  --cds_2_faa work/2026/1126.3_frame_translation/02_GABARAP_cds_2.faa
```
