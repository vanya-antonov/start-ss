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

### Download canonical transcript annotations

The workflow uses the Ensembl REST API at `/info/species` to match each
configured common and scientific name to an Ensembl production name, and
`/info/data` to determine the current release. Download the matching GTF
annotation manually before running the script. The script checks for that file
under `data/mammals/`, parses it in streaming mode, and never downloads,
deletes, or modifies the GTF. UCSC assembly IDs are not treated as Ensembl
production names. Only transcripts tagged `Ensembl_canonical` with both gene
and transcript biotypes `protein_coding` are included. CDS length is the sum
of the transcript's CDS feature lengths. Genome FASTA and GenBank files are not
downloaded.

### Step 1. Download the Ensembl GTF file

The Ensembl 116 FTP release contains GTF files for 10 of the 13 planned
species. Download the available files with these commands (run from the
repository root):

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

These commands use Ensembl's `current` FTP path and the release 116 filenames.
`wget -P` retains each filename in `data/mammals/`, where the script expects
it. The release 116 GTF directory does not contain files for killer whale,
manatee, or aardvark, so those three species cannot be downloaded from these
Ensembl paths. The downloaded `*.gtf.gz` files are ignored by Git.

### Step 2. Parse canonical transcript annotations

Pass a species ID from `SPECIES_ORDER` to process one species, or pass
`__ALL__` after all required species GTF files are available locally. Run from
the repository root or another working directory with the `start-ss-scripts`
environment active:

```bash
python scripts/mammals/01_download_canonical_transcripts.py __ALL__
```

For example, process human annotations only with
`python scripts/mammals/01_download_canonical_transcripts.py human`.

It writes one tab-separated file per species to
`data/mammals/01_canonical_transcripts_{species_id}.tsv`, using the species IDs
in `SPECIES_ORDER`. Each file contains gene ID/name and biotype, transcript ID
and biotype, protein ID, and CDS length. REST, FTP, missing-species, malformed
GTF, and missing-canonical-annotation failures are logged and cause a non-zero
exit status. A species output is written only after its GTF has been completely
downloaded and validated.

### Ortholog and structure analysis

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
