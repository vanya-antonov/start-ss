# Large reference files

Download Ensembl annotation and sequence files to this directory. Keep the
files compressed; large reference files are excluded from Git and must be
downloaded separately on each machine.

## Ensembl release 116 GTF annotations

The GTF files provide gene and canonical transcript annotations, including
CDS features used to calculate transcript-relative coordinates. Run from the
repository root:

```sh
mkdir -p data/large_files/ensembl/gtf
cd data/large_files/ensembl/gtf
wget https://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz
wget https://ftp.ensembl.org/pub/release-116/gtf/mus_musculus/Mus_musculus.GRCm39.116.gtf.gz
wget https://ftp.ensembl.org/pub/release-116/gtf/bos_taurus/Bos_taurus.ARS-UCD2.0.116.gtf.gz
wget https://ftp.ensembl.org/pub/release-116/gtf/ovis_aries/Ovis_aries.ARS-UI_Ramb_v3.0.116.gtf.gz
wget https://ftp.ensembl.org/pub/release-116/gtf/ailuropoda_melanoleuca/Ailuropoda_melanoleuca.ASM200744v2.116.gtf.gz
wget https://ftp.ensembl.org/pub/release-116/gtf/loxodonta_africana/Loxodonta_africana.loxAfr3.116.gtf.gz
wget https://ftp.ensembl.org/pub/release-116/gtf/canis_lupus_familiaris/Canis_lupus_familiaris.ROS_Cfam_1.0.116.gtf.gz
wget https://ftp.ensembl.org/pub/release-116/gtf/equus_caballus/Equus_caballus.EquCab3.0.116.gtf.gz
wget https://ftp.ensembl.org/pub/release-116/gtf/oryctolagus_cuniculus/Oryctolagus_cuniculus.OryCun2.0.116.gtf.gz
wget https://ftp.ensembl.org/pub/release-116/gtf/callithrix_jacchus/Callithrix_jacchus.mCalJac1.pat.X.116.gtf.gz
cd -
```

These ten species are the ones with GTF annotations used by this workflow.

### Step 1. Extract canonical transcripts from GTF files

Convert each GTF into a TSV of canonical protein-coding transcripts. The
tables include versioned gene, transcript, and protein IDs, transcript
lengths, and CDS coordinates relative to the spliced transcript. Run from the
repository root:

```sh
python work/2026/1118.gabarap_orthologs/00_get_canonical_transcripts_from_gtf.py data/large_files/ensembl/gtf/Homo_sapiens.GRCh38.116.gtf.gz
python work/2026/1118.gabarap_orthologs/00_get_canonical_transcripts_from_gtf.py data/large_files/ensembl/gtf/Mus_musculus.GRCm39.116.gtf.gz
python work/2026/1118.gabarap_orthologs/00_get_canonical_transcripts_from_gtf.py data/large_files/ensembl/gtf/Bos_taurus.ARS-UCD2.0.116.gtf.gz
python work/2026/1118.gabarap_orthologs/00_get_canonical_transcripts_from_gtf.py data/large_files/ensembl/gtf/Ovis_aries.ARS-UI_Ramb_v3.0.116.gtf.gz
python work/2026/1118.gabarap_orthologs/00_get_canonical_transcripts_from_gtf.py data/large_files/ensembl/gtf/Ailuropoda_melanoleuca.ASM200744v2.116.gtf.gz
python work/2026/1118.gabarap_orthologs/00_get_canonical_transcripts_from_gtf.py data/large_files/ensembl/gtf/Loxodonta_africana.loxAfr3.116.gtf.gz
python work/2026/1118.gabarap_orthologs/00_get_canonical_transcripts_from_gtf.py data/large_files/ensembl/gtf/Canis_lupus_familiaris.ROS_Cfam_1.0.116.gtf.gz
python work/2026/1118.gabarap_orthologs/00_get_canonical_transcripts_from_gtf.py data/large_files/ensembl/gtf/Equus_caballus.EquCab3.0.116.gtf.gz
python work/2026/1118.gabarap_orthologs/00_get_canonical_transcripts_from_gtf.py data/large_files/ensembl/gtf/Oryctolagus_cuniculus.OryCun2.0.116.gtf.gz
python work/2026/1118.gabarap_orthologs/00_get_canonical_transcripts_from_gtf.py data/large_files/ensembl/gtf/Callithrix_jacchus.mCalJac1.pat.X.116.gtf.gz
```

The output TSVs are written beside their input GTF files in
`data/large_files/ensembl/gtf/`. Later workflow steps can use these
preprocessed tables instead of parsing the GTF files again.

## Ensembl release 116 cDNA FASTA sequences

The cDNA FASTA files contain spliced transcript sequences in transcript
orientation. Step 2 of the 1123 workflow retrieves the selected canonical
transcript sequences from these files. Run from the repository root:

```sh
mkdir -p data/large_files/ensembl/cdna
cd data/large_files/ensembl/cdna
wget https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
wget https://ftp.ensembl.org/pub/release-116/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz
wget https://ftp.ensembl.org/pub/release-116/fasta/bos_taurus/cdna/Bos_taurus.ARS-UCD2.0.cdna.all.fa.gz
wget https://ftp.ensembl.org/pub/release-116/fasta/ovis_aries/cdna/Ovis_aries.ARS-UI_Ramb_v3.0.cdna.all.fa.gz
wget https://ftp.ensembl.org/pub/release-116/fasta/ailuropoda_melanoleuca/cdna/Ailuropoda_melanoleuca.ASM200744v2.cdna.all.fa.gz
wget https://ftp.ensembl.org/pub/release-116/fasta/loxodonta_africana/cdna/Loxodonta_africana.loxAfr3.cdna.all.fa.gz
wget https://ftp.ensembl.org/pub/release-116/fasta/canis_lupus_familiaris/cdna/Canis_lupus_familiaris.ROS_Cfam_1.0.cdna.all.fa.gz
wget https://ftp.ensembl.org/pub/release-116/fasta/equus_caballus/cdna/Equus_caballus.EquCab3.0.cdna.all.fa.gz
wget https://ftp.ensembl.org/pub/release-116/fasta/oryctolagus_cuniculus/cdna/Oryctolagus_cuniculus.OryCun2.0.cdna.all.fa.gz
wget https://ftp.ensembl.org/pub/release-116/fasta/callithrix_jacchus/cdna/Callithrix_jacchus.mCalJac1.pat.X.cdna.all.fa.gz
cd -
```
