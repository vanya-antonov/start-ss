# NCBI orthologs and Ensembl canonical transcripts

NCBI will be used to identify orthologous genes, and Ensembl will provide
canonical transcript annotations, CDS coordinates, and transcript sequences.
The eventual goal is to extract 50-nt transcript regions around annotated
translation start codons: 23 nt upstream, the 3-nt start codon, and 24 nt
downstream. All Ensembl annotations and sequences will come from release 116
for consistency.

The initial analysis will include these ten species:

| Common name | Scientific name |
| --- | --- |
| Human | *Homo sapiens* |
| Mouse | *Mus musculus* |
| Cow | *Bos taurus* |
| Sheep | *Ovis aries* |
| Panda | *Ailuropoda melanoleuca* |
| Elephant | *Loxodonta africana* |
| Dog | *Canis lupus familiaris* |
| Horse | *Equus caballus* |
| Rabbit | *Oryctolagus cuniculus* |
| Marmoset | *Callithrix jacchus* |

## Downloads

Run the commands from the repository root. Keep all downloaded reference files
compressed and store them under `data/large_files/`.

### 1. Ensembl GTF annotations

The GTF files will be used to extract canonical transcript IDs and CDS
start/end coordinates in transcript coordinates.

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

### 2. Ensembl cDNA FASTA sequences

The cDNA FASTA files contain spliced transcript sequences in transcript
orientation. Canonical transcript sequences will be retrieved from these files
using transcript IDs obtained from the GTF annotations.

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

The initial workflow covers ten species because Ensembl release 116 GTF
annotations are not available in the current workflow for killer whale,
manatee, and aardvark.

## Workflow

### Step 1. Prepare NCBI ortholog table

Reuse the NCBI ortholog table generated in the previous analysis (`work/2026/1111.gabarap_orthologs.ncbi/01_GABARAP_ncbi_orthologs.tsv`), retaining the nine species included in the current Ensembl workflow and adding species keys from `constants.py`.

The resulting file is `01_GABARAP_ncbi_orthologs.tsv`.

### Step 2. Add Ensembl canonical transcript IDs

Match NCBI ortholog gene symbols against `gene_name` attributes in the Ensembl GTF annotations.

```sh
python work/2026/1118.gabarap_orthologs/02_add_canonical_transcripts.py
```

The resulting file is `02_GABARAP_ncbi_orthologs.tsv`.
