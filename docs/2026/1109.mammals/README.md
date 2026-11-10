# Mammalian analysis

This document describes the reproducible cross-species analysis of RNA
secondary structures around annotated translation start codons.

## Species and assemblies

The species are listed in analysis order. Assemblies are UCSC assembly names.

| No. | Species ID | Common name | Scientific name | UCSC assembly |
| --- | --- | --- | --- | --- |
| 1 | human | Human | *Homo sapiens* | hg38 |
| 2 | mouse | Mouse | *Mus musculus* | mm10 |
| 3 | whale | Killer whale | *Orcinus orca* | orcOrc1 |
| 4 | cow | Cow | *Bos taurus* | bosTau8 |
| 5 | sheep | Sheep | *Ovis aries* | oviAri3 |
| 6 | panda | Panda | *Ailuropoda melanoleuca* | ailMel1 |
| 7 | manatee | Manatee | *Trichechus manatus latirostris* | triMan1 |
| 8 | elephant | Elephant | *Loxodonta africana* | loxAfr3 |
| 9 | dog | Dog | *Canis lupus familiaris* | canFam3 |
| 10 | aardvark | Aardvark | *Orycteropus afer afer* | oryAfe1 |
| 11 | horse | Horse | *Equus caballus* | equCab2 |
| 12 | rabbit | Rabbit | *Oryctolagus cuniculus* | oryCun2 |
| 13 | marmoset | Marmoset | *Callithrix jacchus* | calJac3 |

## Workflow

### Step 0. Download Ensembl GTF annotations

Run from the repository root:

```bash
wget -P data/mammals https://ftp.ensembl.org/pub/current/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz
wget -P data/mammals https://ftp.ensembl.org/pub/current/gtf/mus_musculus/Mus_musculus.GRCm39.116.gtf.gz
wget -P data/mammals https://ftp.ensembl.org/pub/current/gtf/bos_taurus/Bos_taurus.ARS-UCD2.0.116.gtf.gz
wget -P data/mammals https://ftp.ensembl.org/pub/current/gtf/ovis_aries/Ovis_aries.ARS-UI_Ramb_v3.0.116.gtf.gz
wget -P data/mammals https://ftp.ensembl.org/pub/current/gtf/ailuropoda_melanoleuca/Ailuropoda_melanoleuca.ASM200744v2.116.gtf.gz
wget -P data/mammals https://ftp.ensembl.org/pub/current/gtf/loxodonta_africana/Loxodonta_africana.loxAfr3.116.gtf.gz
wget -P data/mammals https://ftp.ensembl.org/pub/current/gtf/canis_lupus_familiaris/Canis_lupus_familiaris.ROS_Cfam_1.0.116.gtf.gz
wget -P data/mammals https://ftp.ensembl.org/pub/current/gtf/equus_caballus/Equus_caballus.EquCab3.0.116.gtf.gz
wget -P data/mammals https://ftp.ensembl.org/pub/current/gtf/oryctolagus_cuniculus/Oryctolagus_cuniculus.OryCun2.0.116.gtf.gz
wget -P data/mammals https://ftp.ensembl.org/pub/current/gtf/callithrix_jacchus/Callithrix_jacchus.mCalJac1.pat.X.116.gtf.gz
```

Ensembl 116 does not provide GTFs here for killer whale, manatee, or aardvark.

### Step 1. Extract canonical transcripts and CDS coordinates

Run the parser once for each downloaded annotation:

```bash
time python scripts/mammals/01_get_canonical_transcripts_from_gtf.py --gtf data/mammals/Homo_sapiens.GRCh38.116.gtf.gz -o data/mammals/01_human --exclude_genes_without_name
time python scripts/mammals/01_get_canonical_transcripts_from_gtf.py --gtf data/mammals/Mus_musculus.GRCm39.116.gtf.gz -o data/mammals/01_mouse
time python scripts/mammals/01_get_canonical_transcripts_from_gtf.py --gtf data/mammals/Bos_taurus.ARS-UCD2.0.116.gtf.gz -o data/mammals/01_cow
time python scripts/mammals/01_get_canonical_transcripts_from_gtf.py --gtf data/mammals/Ovis_aries.ARS-UI_Ramb_v3.0.116.gtf.gz -o data/mammals/01_sheep
time python scripts/mammals/01_get_canonical_transcripts_from_gtf.py --gtf data/mammals/Ailuropoda_melanoleuca.ASM200744v2.116.gtf.gz -o data/mammals/01_panda
time python scripts/mammals/01_get_canonical_transcripts_from_gtf.py --gtf data/mammals/Loxodonta_africana.loxAfr3.116.gtf.gz -o data/mammals/01_elephant
time python scripts/mammals/01_get_canonical_transcripts_from_gtf.py --gtf data/mammals/Canis_lupus_familiaris.ROS_Cfam_1.0.116.gtf.gz -o data/mammals/01_dog
time python scripts/mammals/01_get_canonical_transcripts_from_gtf.py --gtf data/mammals/Equus_caballus.EquCab3.0.116.gtf.gz -o data/mammals/01_horse
time python scripts/mammals/01_get_canonical_transcripts_from_gtf.py --gtf data/mammals/Oryctolagus_cuniculus.OryCun2.0.116.gtf.gz -o data/mammals/01_rabbit
time python scripts/mammals/01_get_canonical_transcripts_from_gtf.py --gtf data/mammals/Callithrix_jacchus.mCalJac1.pat.X.116.gtf.gz -o data/mammals/01_marmoset
```

Each `-o` prefix produces `<prefix>.canonical_transcripts.tsv` and
`<prefix>.errors.txt`; for example, the human command writes
`data/mammals/01_human.canonical_transcripts.tsv` and
`data/mammals/01_human.errors.txt`. Without `-o`, filenames use the GTF prefix.
The table includes versioned IDs,
transcript and CDS lengths, and CDS coordinates relative to the spliced
transcript.

### Step 2. Identify human orthologs

```bash
python scripts/mammals/02_identify_orthologous_genes.py
```

Reads `gene_id` and `gene_name` from `01_Homo_sapiens.canonical_transcripts.tsv`;
writes one row per human gene to `03_human_orthologs.tsv`. It queries Ensembl
Compara for the nine other species with Ensembl 116 GTFs, caches responses in
`data/mammals/.cache/`, and retains multiple orthologs in semicolon-separated
fields. Use `--input`, `--output`, or `--cache-dir` to override paths.

The analysis identifies orthologous genes across these species, retrieves
annotated transcripts, and extracts transcript-oriented 50-nt windows
around each annotated start codon: 23 nt upstream, the ATG start codon, and
24 nt downstream. The resulting windows will be used to analyze conservation
of predicted RNA secondary structures.

The analysis records orthology criteria, transcript and isoform selection,
annotation and sequence sources, missing-data handling, and structure
prediction software, version, and settings. Predicted structure is a model;
conservation or association in these data does not establish a causal effect
on translation initiation.

Mammalian scripts and inputs are stored under `scripts/mammals/` and
`data/mammals/`. Derived results are stored under `results/mammals/`.
